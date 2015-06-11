"""
:mod:`op` -- Finite differences time domain steady-state analysis
-----------------------------------------------------------------

.. module:: op
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

"""

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import sys

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from analysis import AnalysisError, ipython_drop
from fsolve import solve, NoConvergenceError
import nodalSP as nd
import nodal_FDTD
import analysis 

class FDTD(ParamSet):
    r"""
    FDTD Steady-State Analysis
    --------------------------

    Write documentation here ...

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The type of matrix used in
    this analysis is controlled by the ``sparse`` option. Global
    options are documented in :doc:`global_vars`. 

    Example::

        .analysis fdtd T = 1us nsamples = 100  

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "fdtd"

    # Define parameters as follows
    paramDict = dict(
        T = ('Fundamental period', 's', float, 0.),
        nsamples = ('Number of samples in period (must be > 5)', '', int, 100),
	saveall = ('Save all nodal voltages', '', bool, False)
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
            raise analysis.AnalysisError(
                'Invalid period: {0}'.format(self.T))
        if self.nsamples < 4:
            raise analysis.AnalysisError(
                'Number of samples is too low: {0}'.format(self.nsamples))

        if not glVar.sparse:
            print('Warning: dense matrices not supported in this analysis, using sparse matrices.\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodal object
        nd.make_nodal_circuit(circuit)
#	dc = nd.DCNodal(circuit)
	fdtd = nodal_FDTD.FDTDNodal(circuit, self.nsamples, self.T)

	# For now just start with zeros as initial guess
	x = np.zeros(fdtd.dim)
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
        p, pt = nodal_FDTD.calculate_P_PT(circuit.nD_dimension, self.nsamples)
        xHat = x[p]

        # fdtd.timeVec contains the time vector
        timeVec = fdtd.timeVec
	
#        # Get terminals to plot/save from circuit. 
#        termSet = circuit.get_requested_terms('fdtd')
#
#        # Special treatment for ground terminal
#        termSet1 = set(termSet)
#        if circuit.nD_ref in termSet1:
#            termSet1.remove(circuit.nD_ref)
#            circuit.nD_ref.tran_v = np.zeros(self.nsamples)
#
#        # Allocate vectors for results
#        if self.saveall:
#            for term in circuit.nD_termList:
#                term.tran_v = np.empty(self.nsamples)
#                term.tran_v[0] = x[term.nD_namRC]                
#            circuit.nD_ref.tran_v = np.zeros(self.nsamples)
#        else:
#            # Only save requested nodes
#            for term in termSet1:
#                term.tran_v = np.empty(self.nsamples)
#                term.tran_v[0] = x[term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        # Save results from fdtd analysis
        # to be implemented ...

        import matplotlib.pylab as plt
        plt.ion()
        # Plot first nodal variable
        plt.plot(timeVec, xHat[:self.nsamples])

        analysis.ipython_drop("\nfdtd shell:", globals(), locals())

aClass = FDTD



