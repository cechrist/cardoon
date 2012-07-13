"""
:mod:`tran` -- Transient Analysis
---------------------------------

.. module:: tran
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import sys
import numpy as np

from paramset import ParamSet
import analysis 
from integration import BEuler, Trapezoidal
from globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodal
import nodalSP

class Transient(ParamSet):
    """
    Transient Analysis
    ------------------

    Solves nodal equations starting from ``t=0`` to ``tstop`` with a
    fixed time step equal to ``tstep``. Two integration methods are
    supported: Backwards Euler (``im = BE``) and trapezoidal
    (``im=trap``). Support for frequency-defined elements and time
    delays is not yet included.

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
        im = ('Integration method', '', str, 'trap'),
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
            raise analysis.AnalysisError(
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

        # Get terminals to plot/save from circuit. 
        termSet = circuit.get_requested_terms('tran')

        # Special treatment for ground terminal
        termSet1 = set(termSet)
        if circuit.nD_ref in termSet1:
            termSet1.remove(circuit.nD_ref)
            circuit.nD_ref.tran_v = np.zeros(nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.tran_v = np.empty(nsamples)
                term.tran_v[0] = x[term.nD_namRC]                
            circuit.nD_ref.tran_v = np.zeros(nsamples)
        else:
            # Only save requested nodes
            for term in termSet1:
                term.tran_v = np.empty(nsamples)
                term.tran_v[0] = x[term.nD_namRC]

        # Save initial values
        xOld = x
        tIter = 0
        tRes = 0.
        dots = 50
        print('System dimension: {0}'.format(circuit.nD_dimension))
        print('Number of samples: {0}'.format(nsamples))
        print('Integration method: {0}'.format(self.im))
        if self.verbose:
            print('-------------------------------------------------')
            print(' Step    | Time (s)     | Iter.    | Residual    ')
            print('-------------------------------------------------')
        else:
            print('Printing one dot every {0} samples:'.format(dots))
            sys.stdout.flush()

        for i in xrange(1, nsamples):
            tran.accept(xOld)
            sV = tran.get_rhs(timeVec[i])
            # solve equations: use previous time-step solution as an
            # initial guess
            if i > 1 and glVar.sparse:
                # Re-use factorized Jacobian: This saves the time to
                # evaluate the function and Jacobian plus the time for
                # factorization. Only sparse implementation stores
                # factorized Jacobian
                xOld -= tran.get_chord_deltax(sV)
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
                for term in termSet1:
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

        # Process output requests.  
        analysis.process_requests(circuit, 'tran', 
                                  timeVec, 'Time [s]', 'tran_v')

        def getvec(termname):
            return circuit.termDict[termname].tran_v

        if self.shell:
            analysis.ipython_drop("""
Available commands:
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())



aClass = Transient

