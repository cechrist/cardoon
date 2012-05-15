"""
:mod:`tran` -- Transient Analysis
---------------------------------

.. module:: tran
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

from paramset import ParamSet
from analysis import AnalysisError, ipython_drop
from integration import BEuler, Trapezoidal
from globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodalSP
import nodal

class Analysis(ParamSet):
    """
    Transient Analysis
    ------------------

    Solves nodal equations starting from ``t=0`` to ``tstop`` with a
    fixed time step (at least for now) equal to ``tstep``. Two
    integration methods are supported: Backwards Euler (``im = BE``)
    and trapezoidal (``im=trap``). Support for frequency-defined
    elements and time delays is not yet included.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The type of matrix used in
    this analysis is controlled by the ``sparse`` option. Global
    options are documented in :doc:`global_vars`. 

    One plot window is generated for each ``.plot`` statement. Use
    ``tran`` request type for this analysis. By default, only results
    for nodes listed in ``.plot`` statements are saved. To save all
    nodal variables set ``saveall`` to 1.

    Transient analysis formulation is documented in :doc:`analysis`,
    and internal classes and functions used in this analysis are
    documented in :doc:`analyses_classes`.

    Example::

        .analysis tran tstop=1ms tstep=.01ms im=BE

        .plot tran vin vout

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "tran"

    # Define parameters as follows
    paramDict = dict(
        tstop = ('Simulation stop time', 's', float, 1e-3),
        tstep = ('Time step size', 's', float, 1e-5),
        im = ('Integration method', '', str, 'BE'),
        verbose = ('Show iterations for each point', '', bool, False),
        saveall = ('Save all nodal voltages', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )

    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calculates transient analysis by solving nodal equations
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('                 Transient analysis')
        print('******************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        if glVar.sparse:
            nd = nodalSP
        else:
            nd = nodal
            print('Using dense matrices\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Select integration method
        if self.im == 'BE':
            imo = BEuler()
        elif self.im == 'trap':
            imo = Trapezoidal()
        else:
            raise AnalysisError(
                'Unknown integration method: {0}'.format(self.im))

        # Create nodal objects and solve for initial state
        nd.make_nodal_circuit(circuit)
        dc = nd.DCNodal(circuit)
        tran = nd.TransientNodal(circuit, imo)
        x = dc.get_guess()
        # Use sources including transient values for t == 0
        sV = tran.get_source(0.)
        # solve DC equations
        try: 
            print('Calculating DC operating point ... ', end='')
            sys.stdout.flush()
            (x, res, iterations) = solve(x, sV, dc.convergence_helpers)
            print('Succeded.\n')
        except NoConvergenceError as ce:
            print('Failed.\n')
            print(ce)
            return
        dc.save_OP(x)
        tran.set_IC(self.tstep)
        # Release memory in dc object?
        del(dc)

        # Create time vector
        timeVec = np.arange(start=0., stop = self.tstop, step = self.tstep, 
                            dtype=float)
        nsamples = len(timeVec)
        circuit.tran_timevec = timeVec

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.tran_v = np.empty(nsamples)
                term.tran_v[0] = x[term.nD_namRC]                
            circuit.nD_ref.tran_v = np.zeros(nsamples)
        else:
            # Only save requested nodes
            for outreq in circuit.outReqList:
                if outreq.type == 'tran':
                    for termname in outreq.varlist:
                        term = circuit.termDict[termname]
                        term.tran_v = np.empty(nsamples)
                        term.tran_v[0] = x[term.nD_namRC]

        # Save initial values
        xOld = x
        tIter = 0
        tRes = 0.
        dots = 50
        print('System dimension: {0}'.format(circuit.nD_dimension))
        print('Number of samples: {0}'.format(nsamples))
        if self.verbose:
            print('-------------------------------------------------')
            print(' Step    | Time (s)     | Iter.    | Residual    ')
            print('-------------------------------------------------')
        else:
            print('Printing one dot every {0} samples:'.format(dots))
            print('.', end='')
            sys.stdout.flush()

        for i in xrange(1, nsamples):
            qVec = tran.update_q(xOld)
            imo.accept(qVec)
            sV = imo.f_n1()
            sV += tran.get_source(timeVec[i])
            # solve equations: use previous time-step solution as an
            # initial guess
            try: 
                (x, res, iterations) = solve(xOld, sV, 
                                             tran.convergence_helpers)
            except NoConvergenceError as ce:
                print(ce)
                return
            # Save results
            xOld[:] = x
            if self.saveall:
                for term in circuit.nD_termList:
                    term.tran_v[i] = x[term.nD_namRC]                
            else:
                # Only save requested nodes
                for outreq in circuit.outReqList:
                    if outreq.type == 'tran':
                        for termname in outreq.varlist:
                            term = circuit.termDict[termname]
                            term.tran_v[i] = x[term.nD_namRC]
            # Keep some info about iterations
            tIter += iterations
            tRes += res
            if self.verbose:
                print('{0:8} | {1:12} | {2:8} | {3:12}'.format(
                        i, timeVec[i], iterations, res))
            elif not i%dots:
                print('.', end='')
                sys.stdout.flush()

        # Calculate average residual and iterations
        avei = float(tIter) / nsamples
        aver = tRes / nsamples
        print('\nAverage iterations: {0}'.format(avei))
        print('Average residual: {0}\n'.format(aver))

        # Process output requests.  In the future this may be moved
        # to a common module that processes output requests such as
        # plot, print, save, etc.
        flag = False
        for outreq in circuit.outReqList:
            if outreq.type == 'tran':
                flag = True
                plt.figure()
                plt.grid(True)
                plt.xlabel('Time [s]')
                for termname in outreq.varlist:
                    term = circuit.termDict[termname]
                    plt.plot(timeVec, term.tran_v, 
                             label = 'Term: {0} [{1}]'.format(
                            term.nodeName, term.unit)) 
                if len(outreq.varlist) == 1:
                    plt.ylabel(
                        'Term: {0} [{1}]'.format(term.nodeName, term.unit))
                else:
                    plt.legend()
        if flag:
            plt.show()

        def getvec(termname):
            return circuit.termDict[termname].tran_v

        if self.shell:
            ipython_drop("""
Available commands:
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
    plt.* to access pyplot commands (plt.plot(x,y), plt.show(), etc.)
""", globals(), locals())

