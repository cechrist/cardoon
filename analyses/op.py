"""
:mod:`op` -- Operating Point Analysis
-------------------------------------

.. module:: op
.. moduleauthor:: Carlos Christoffersen

This is an attempt to port/adapt the existing analyses in pycircuit
(https://github.com/henjo/pycircuit). For now most of the code is cut
and pasted from there.

************** This is experimental/incomplete ****************                

"""

from __future__ import print_function
import numpy as np
from globalVars import glVar
from paramset import ParamSet
from analysis import AnalysisError

class Analysis(ParamSet):
    """
    DC Operating Point Calculation

    Calculates the DC operating point of a circuit. This is also
    useful for many other analyses such as DC, AC, Transient, HB,
    etc. Charge sources are ignored in this calculation.
    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "OP"

    # Define parameters as follows
    paramDict = dict(
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)


    def run(self, circuit):
        """
        Calculates the operating point by solving nodal equations

        The state of all devices is determined by the values of the
        voltages at the controlling ports.
        """
        # 
        # 3. Call modified pycircuit's functions to solve equations
        #
        # 4. if not very difficult use pycircuit's output routines

        pass
        

    def solve(self):
        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        convergence_helpers = [self._simple, self._homotopy_gmin, 
                               self._homotopy_source, 
                               None]

        x0 = np.zeros(self.cir.n) # Would be good with a better initial guess

        for algorithm in convergence_helpers:
            if algorithm == None:
                raise last_e
            else:
                if algorithm.__doc__:
                    logging.info('Trying ' + algorithm.__doc__)
                try:
                    x = algorithm(x0)
                except (NoConvergenceError, SingularMatrix), last_e:
                    logging.warning('Problems encoutered: ' + str(last_e))
                else:
                    break

        self.result = CircuitResultDC(self.cir, x)

        return self.result

    def _simple(self, x0):
        """Simple Newton's method"""
        def func(x):
            return self.cir.i(x) + self.cir.u(0,analysis='dc'), self.cir.G(x)

        return self._newton(func, x0)

    def _homotopy_gmin(self, x0):
        """Newton's method with gmin stepping"""
        x = x0
        for gmin in (1, 1e-1, 1e-2, 0):
            n_nodes = len(self.cir.nodes)
            Ggmin = np.zeros((self.cir.n, self.cir.n))
            Ggmin[0:n_nodes, 0:n_nodes] = gmin * np.eye(n_nodes)

            def func(x):
                return self.cir.i(x) + self.cir.u(0,analysis='dc'), \
                       self.cir.G(x) + Ggmin

            x, x0 = self._newton(func, x0), x

        return x

    def _homotopy_source(self, x0):
        """Newton's method with source stepping"""
        x = x0
        for lambda_ in (0, 1e-2, 1e-1, 1):
            def func(x):
                f = self.cir.i(x) + lambda_ * self.cir.u(0,analysis='dc')
                dFdx = self.cir.G(x)
                return f, dFdx            
            x, x0 = self._newton(func, x0), x

        return x

    def _newton(self, func, x0):
        ones_nodes = np.ones(len(self.cir.nodes))
        ones_branches = np.ones(len(self.cir.branches))

        abstol = np.concatenate((self.par.iabstol * ones_nodes,
                                 self.par.vabstol * ones_branches))
        xtol = np.concatenate((self.par.vabstol * ones_nodes,
                                 self.par.iabstol * ones_branches))

        (x0, abstol, xtol) = remove_row_col((x0, abstol, xtol), self.irefnode, self.toolkit)

        try:
            result = fsolve(refnode_removed(func, self.irefnode,self.toolkit), 
                            x0, 
                            full_output = True, 
                            reltol = self.par.reltol,
                            abstol = abstol, xtol=xtol,
                            maxiter = self.par.maxiter,
                            toolkit = self.toolkit)
        except np.linearsolverError(), e:
            raise SingularMatrix(e.message)

        x, infodict, ier, mesg = result

        if ier != 1:
            raise NoConvergenceError(mesg)

        # Insert reference node voltage
        return np.concatenate((x[:self.irefnode], np.array([0.0]), x[self.irefnode:]))





def fsolve(f, x0, args=(), full_output=False, maxiter=200,
           xtol=1e-6, reltol=1e-4, abstol=1e-12, toolkit='Numeric'):
    """Solve a multidimensional non-linear equation with Newton-Raphson's method

    In each iteration the linear system

    M{J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained x_{n+1}
    
    """
    
    converged = False
    ier = 2
    for i in xrange(maxiter):
        F, J = f(x0, *args) # TODO: Make sure J is never 0, e.g. by gmin (stepping)
        xdiff = toolkit.linearsolver(J, -F)# TODO: Limit xdiff to improve convergence

        x = x0 + xdiff

        if toolkit.alltrue(abs(xdiff) < reltol * toolkit.maximum(x, x0) + xtol):
            ier = 1
            mesg = "Success"
            break
        if toolkit.alltrue(abs(F) < reltol * max(F) + abstol):
            ier = 1
            mesg = "Success"
            break
            
        x0 = x

    if ier == 2:
        mesg = "No convergence. xerror = "+str(xdiff)
    
    infodict = {}
    if full_output:
        return x, infodict, ier, mesg
    else:
        return x
