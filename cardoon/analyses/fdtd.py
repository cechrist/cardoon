"""
:mod:`fdtd` -- Finite differences time domain steady-state analysis
-------------------------------------------------------------------

.. module:: fdtd
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

"""

from __future__ import print_function
import os.path
import numpy as np
import scipy.sparse as sp
import sys

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from analysis import AnalysisError, ipython_drop
from fsolve import solve, NoConvergenceError
import nodalSP as nd
import nodalFDTD
import analysis 

# Valid request types
reqTypes = ['fdtd']

class FDTD(ParamSet):
    r"""
    FDTD Steady-State Analysis
    --------------------------

    This is an experimental analysis. It solves the steady-state
    problem using a finite-diference time-domain approach (FDTD).

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The type of matrix used in
    this analysis is controlled by the ``sparse`` option. Global
    options are documented in :doc:`global_vars`. 

    Example::

        .analysis fdtd T = 1us nsamples = 100  

        .plot fdtd in out

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "fdtd"

    # Define parameters as follows
    paramDict = dict(
        T = ('Fundamental period (overrides frequency)', 's', float, 0.),
        f = ('Fundamental frequency', 'Hz', float, 0.),
        nsamples = ('Number of samples in period (must be > 5)', '', int, 100),
        dcguess = ('Use DC operating point as initial guess', '', bool, False),
	saveall = ('Save all nodal voltages', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calculates steady state using FDTD method
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('                  FDTD analysis')
        print('******************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        # Check netlist parameters
        if self.T <= 0.:
            if self.f <=0.:
                raise analysis.AnalysisError(
                    'Invalid period/frequency: {0}/{1}'.format(self.T, self.f))
            else:
                self.T = 1./self.f
        if self.nsamples < 4:
            raise analysis.AnalysisError(
                'Number of samples is too low: {0}'.format(self.nsamples))

        if not glVar.sparse:
            print('Warning: dense matrices not supported in this analysis, using sparse matrices.\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodal and FDTD objects
        nd.make_nodal_circuit(circuit)
	fdtd = nodalFDTD.FDTDNodal(circuit, self.nsamples, self.T)

        if (self.dcguess):
            dc = nd.DCNodal(circuit)
            xDC = dc.get_guess()
            sVDC = dc.get_source()
            # solve DC equations
            try: 
                print('Calculating DC operating point ... ', end='')
                sys.stdout.flush()
                (xDC, res, iterations) = solve(xDC, sVDC, 
                                               dc.convergence_helpers)
                print('Succeded.\n')
            except NoConvergenceError as ce:
                print('Failed.\n')
                print(ce)
                return
            dc.save_OP(xDC)
	    # Use DC solution as initial guess
            x = np.tile(xDC, self.nsamples)
        else:
            # Start from zero instead (should be an option)
            x = np.zeros(fdtd.dim)

        # Solve equations
        print('System size: ', fdtd.dim)
	sV = fdtd.get_source()
        # solve equations 
        try: 
            print('Solving system ... ', end='')
            sys.stdout.flush()
            (x, res, iterations) = solve(x, sV, fdtd.convergence_helpers)
            print('Succeded.\n')
        except NoConvergenceError as ce:
            print('Failed.\n')
            print(ce)
            return

        # re-order circuit variables (samples for one variable in one
        # block)

        p, pt = nodalFDTD.calculate_P_PT(circuit.nD_dimension, self.nsamples)
        # Re-shape so that xHat[1] returns all samples for nodal variable 1
        xHat = np.reshape(x[p], (circuit.nD_dimension, self.nsamples))

        xMatrix = xHat*0
        for i in range(fdtd.xMatrix.shape[0]):
            xMatrix = np.vstack((xMatrix,np.reshape(fdtd.xMatrix[i,:][p], (circuit.nD_dimension, self.nsamples))[1,:]))

        # fdtd.timeVec contains the time vector
        timeVec = fdtd.timeVec
	
        # Save results from fdtd analysis ********************************

        # Get terminals to plot/save from circuit. 
        termSet = circuit.get_requested_terms('fdtd')

        # Special treatment for ground terminal
        termSet1 = set(termSet)
        if circuit.nD_ref in termSet1:
            termSet1.remove(circuit.nD_ref)
            circuit.nD_ref.fdtd_v = np.zeros(self.nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.fdtd_v = xHat[term.nD_namRC]
        else:
            # Only save requested nodes
            for term in termSet1:
                term.fdtd_v = xHat[term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        import matplotlib.pylab as plt
        plt.ion()
        plt.figure(3)
        for i in range(xMatrix.shape[0]):
            plt.plot(timeVec, xMatrix[i,:])
        #import pdb; pdb.set_trace()
        np.savetxt('fdtdoutputs.txt', xMatrix, delimiter=',')

        # Process output requests.  
        analysis.process_requests(circuit, 'fdtd', 
                                  timeVec, 'Time [s]', 'fdtd_v')

        def getvec(termname):
            return circuit.find_term(termname).fdtd_v

        if self.shell:
            import matplotlib.pylab as plt
            plt.ion()
            analysis.ipython_drop("""
Available commands:
    matplotlib.pylab imported as 'plt'
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())

aClass = FDTD



