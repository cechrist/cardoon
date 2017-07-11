"""
:mod:`sssw` -- Steady-state Sparse Wavelet Analysis
-----------------------------------------------------

.. module:: sssw
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

"""

from __future__ import print_function
import numpy as np
import sys

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodalSP as nd
import nodalSparseWavelet 
import analysis 

# Valid request types
reqTypes = ['sssw', 'ssswcoeff']

class SSSW(ParamSet):
    r"""
    SSSW Steady-State Sparse Wavelet Analysis
    ----------------------------------------

    Solve the steady-state problem in wavelet domain. Analysis efficiency
    is increased by sensing the positions of the significant amplitude wavelet
    coefficients and creating a reduced size Jacobian that will solve for only
    those coefficients.  This results in a sparse solution that can be solved
    more efficiently than the standard steady state wavelet analysis.
    This analysis is not implemented for efficiency but just as a proof of
    concept. Frequency-defined elements are supported.
    
    There are several methods included which are used for ongoing research.
    These methods can be accessed with the 'method' option:
            Classical - Just do the standard ssw analysis but with the least
            squares solver.  This is just for speed/coefficient density
            comparison and testing.
            Threshold - Threshold off (Set to zero) low amplitude coefficients
            of the solution before returning the result.  This is for testing
            if returning a solution that doesn't have all of its coefficients
            will result in a solution that is similar/the same as the standard
            ssw method.
            RecalcRect - Determines where the significant amplitude non-zero
            coefficients are by first calculating the solution and then forms
            a tall rectangular Jacobian baised on the locations.  The reduced
            tall Jacobian is then solved and this new result is returned.  This
            method is used for determining if the simulation will converge on
            a similar/the same result as the standard ssw method if a perfect
            estimate of the coefficient locations can be determined per
            iteration.
            SenseRect (Default) - Senses the locations of the significant
            amp Waveletlitude non-zero coefficients of the solution each iteration,
            forms a reduced tall Jacobian baised on these locations, solves the
            reduced system, and returns the solution of the reduced system each
            iteration.  This method is for testing the current research idea.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The ``sparse`` option is
    ignored in this analysis (sparse matrices always used). Global
    options are documented in :doc:`global_vars`.

    Example::

        .analysis ssw T = 1us nsamples = 100 method = 'SenseRect' 

        .plot ssw in out

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "sssw"

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
        shell = ('Drop to ipython shell after calculation', '', bool, False), 
        method = ('Classical, Threshold, RecalcRect, SenseRect', '', str, 'SenseRect')
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calculates steady state using SSSW method
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('           Steady-State Sparse Wavelet analysis')
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
	nodalwav = nodalSparseWavelet.SparseWaveletNodal(circuit, self.nsamples, self.T,
                                             self.wavelet, self.multilevel,
                                             self.deriv, self.ssfactor, self.method)

        # Retrieve wavelet transform matrix
        W = nodalwav.W
        Wi = nodalwav.Wi

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
                                             nodalwav.convergence_helpers[1:])
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
        xTilde = x.reshape((circuit.nD_dimension, self.nsamples)).T
        xHat = Wi.dot(xTilde)

        # nodalwav.timeVec contains the time vector
        timeVec = nodalwav.timeVec
	
        # Save results from analysis ********************************

        # Get terminals to plot/save from circuit. 
        termSetw = circuit.get_requested_terms('ssswcoeff')
        termSet = circuit.get_requested_terms('sssw')

        # Special treatment for ground terminal
        termSet1w = set(termSetw)
        if circuit.nD_ref in termSet1w:
            termSet1w.remove(circuit.nD_ref)
            circuit.nD_ref.ssswcoeff_v = np.zeros(self.nsamples)
        termSet1 = set(termSet)
        if circuit.nD_ref in termSet1:
            termSet1.remove(circuit.nD_ref)
            circuit.nD_ref.sssw_v = np.zeros(self.nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.ssswcoeff_v = xTilde[:,term.nD_namRC]
                term.sssw_v = xHat[:,term.nD_namRC]
        else:
            # Only save requested nodes
            for term in termSet1w:
                term.ssswcoeff_v = xTilde[:,term.nD_namRC]
            for term in termSet1:
                term.sssw_v = xHat[:,term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        # Process output requests.  
        cnVec = np.arange(0, self.nsamples)
        analysis.process_requests(circuit, 'ssswcoeff', 
                                  cnVec, 'Coefficient number', 'ssswcoeff_v',
                                  style = '-o')

        analysis.process_requests(circuit, 'sssw', 
                                  timeVec, 'Time [s]', 'sssw_v')
        
        def getvec(termname):
            return circuit.find_term(termname).sssw_v

        if self.shell:
            import matplotlib.pylab as plt
            plt.ion()
            analysis.ipython_drop("""
Available commands:
    matplotlib.pylab imported as 'plt'
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())

aClass = SSSW



