"""
:mod:`ssw` -- Steady-state analysis based on wavelets
-----------------------------------------------------

.. module:: ssw
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

"""

from __future__ import print_function
import numpy as np
import sys

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodalSP as nd
import nodalWavelet 
import analysis 

# Valid request types
reqTypes = ['ssw']

class SSW(ParamSet):
    r"""
    SSW Steady-State Analysis using Wavelets
    ----------------------------------------

    Solve the steady-state problem in wavelet domain. This analysis
    not implemented for efficiency but just as a proof of
    concept. Frequency-defined elements are supported.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The ``sparse`` option is
    ignored in this analysis (sparse matrices always used). Global
    options are documented in :doc:`global_vars`.

    Example::

        .analysis ssw T = 1us nsamples = 100  

        .plot ssw in out

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "ssw"

    # Define parameters as follows
    paramDict = dict(
        T = ('Fundamental period (overrides frequency)', 's', float, 0.),
        f = ('Fundamental frequency', 'Hz', float, 0.),
        nsamples = ('Number of samples in period (always even, power of 2 for multilevel)', '', int, 64),
        wavelet = ('Wavelet family (none=FDTD)', '', str, 'db4'),
        multilevel = ('Use multilevel transform', '', bool, False),
        deriv = ('Derivative type: d2, d4, Fourier', '', str, 'd2'),
#        dcguess = ('Use DC operating point as initial guess', '', bool, False),
        step = ('Directly try conservative convergence helpers', '', bool, False),
        ssfactor = ('Initial source stepping factor', '', float, 0.5),
	saveall = ('Save all nodal voltages', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calculates steady state using SSW method
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('           Steady-State Wavelet analysis')
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
        if self.nsamples < 16:
            raise analysis.AnalysisError(
                'Number of samples is too low: {0}'.format(self.nsamples))
        if np.mod(self.nsamples,2):
            raise analysis.AnalysisError(
                'Number of samples must be even: {0}'.format(self.nsamples))
        if self.multilevel and (2**int(np.log2(self.nsamples))
                                != self.nsamples):
            raise analysis.AnalysisError(
                'multilevel set but nsamples not a power of two: {0}'.format(self.nsamples))

        if not glVar.sparse:
            print('Warning: dense matrices not supported, using sparse matrices.\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodalwav object
        nd.make_nodal_circuit(circuit)
	nodalwav = nodalWavelet.WaveletNodal(circuit, self.nsamples, self.T,
                                             self.wavelet, self.multilevel,
                                             self.deriv, self.ssfactor)

        # Retrieve wavelet transform matrix
        W = nodalwav.W
        Wi = nodalwav.Wi
        # The code below disabled because it may not work well
#        if (self.dcguess):
#            dc = nd.DCNodal(circuit)
#            xDC = dc.get_guess()
#            sVDC = dc.get_source()
#            # solve DC equations
#            try: 
#                print('Calculating DC operating point ... ', end='')
#                sys.stdout.flush()
#                (xDC, res, iterations) = solve(xDC, sVDC, 
#                                               dc.convergence_helpers)
#                print('Succeded.\n')
#            except NoConvergenceError as ce:
#                print('Failed.\n')
#                print(ce)
#                return
#            dc.save_OP(xDC)
#	    # Use DC solution as initial guess
#            u = np.atleast_2d(W.dot(np.ones(self.nsamples))).T
#            x = np.dot(u, np.atleast_2d(xDC)).ravel()
#        else:
        # Start from zero instead (should be an option)
        x = np.zeros(nodalwav.dim)

        # Solve equations
        print('System size: ', nodalwav.dim)
	sV = nodalwav.get_source()
        # solve equations 
        try: 
            print('Solving system ... ', end='')
            sys.stdout.flush()
            if self.step:
                # Skip direct solution and try conservative
                # convergence helpers directly
                (x, res, iterations) = solve(x, sV,
                                             nodalwav.convergence_helpers[2:])
            else:
                (x, res, iterations) = solve(x, sV,
                                             nodalwav.convergence_helpers)
            print('Succeded.\n')
        except NoConvergenceError as ce:
            print('Failed.\n')
            print(ce)
            return

        print('\nFinal coefficients density:',
              100. * np.sum(np.abs(x)>1e-2)/len(x), '%')
        print('Final Jacobian density: ',
              100. * nodalwav.Jac.nnz / nodalwav.dim**2,'%')
        
        # re-shape circuit variables and convert to time domain for display
        xHat = Wi.dot(x.reshape((circuit.nD_dimension, self.nsamples)).T)

        # nodalwav.timeVec contains the time vector
        timeVec = nodalwav.timeVec
	
        # Save results from analysis ********************************

        # Get terminals to plot/save from circuit. 
        termSet = circuit.get_requested_terms('ssw')

        # Special treatment for ground terminal
        termSet1 = set(termSet)
        if circuit.nD_ref in termSet1:
            termSet1.remove(circuit.nD_ref)
            circuit.nD_ref.ssw_v = np.zeros(self.nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.ssw_v = xHat[:,term.nD_namRC]
        else:
            # Only save requested nodes
            for term in termSet1:
                term.ssw_v = xHat[:,term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        # Process output requests.  
        analysis.process_requests(circuit, 'ssw', 
                                  timeVec, 'Time [s]', 'ssw_v')

        def getvec(termname):
            return circuit.find_term(termname).ssw_v

        if self.shell:
            import matplotlib.pylab as plt
            plt.ion()
            analysis.ipython_drop("""
Available commands:
    matplotlib.pylab imported as 'plt'
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())

aClass = SSW


