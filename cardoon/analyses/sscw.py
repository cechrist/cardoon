"""
:mod:`sscw` -- Steady-state analysis based on compressed wavelet coefficients
-----------------------------------------------------------------------------

.. module:: sscw
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

"""

from __future__ import print_function
import numpy as np
import sys

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodalSP as nd
import nodalCompressed
import analysis 

# Valid request types
reqTypes = ['sscw', 'sscwcoeff']

class SSCW(ParamSet):
    r"""
    SSCW Steady-State Analysis using Compressed Wavelet Coefficients
    ----------------------------------------------------------------

    Solve the steady-state problem in compressed wavelet domain using
    compressed sampling theory. This analysis is experimental and
    expected to change. Frequency-defined elements are supported.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The ``sparse`` option is
    ignored in this analysis (sparse matrices always used). Global
    options are documented in :doc:`global_vars`.

    Example::

        .analysis sscw T = 1us nsamples = 128 ncoeff = 32

        .plot sscw in out

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "sscw"

    # Define parameters as follows
    paramDict = dict(
        T = ('Fundamental period (overrides frequency)', 's', float, 0.),
        f = ('Fundamental frequency', 'Hz', float, 0.),
        nsamples = ('Number of samples in period (power of 2)', '', int, 128),
        wavelet = ('Wavelet family', '', str, 'db4'),
        deriv = ('Derivative type: d2, d4, Fourier', '', str, 'd2'),
        ncoeff = ('Number of compressed coefficients', '', int, None),        
        step = ('Directly try conservative convergence helpers', '', bool,
                False),
        ssfactor = ('Initial source stepping factor', '', float, 0.5),
	saveall = ('Save all nodal voltages', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calculates steady state using SSCW method
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('       Steady-State Compressed Wavelet analysis')
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
        if 2**int(np.log2(self.nsamples)) != self.nsamples:
            raise analysis.AnalysisError(
                'nsamples not a power of two: {0}'.format(self.nsamples))
        if not self.ncoeff:
            self.ncoeff = self.nsamples / 2
        elif self.ncoeff > self.nsamples / 2:
             raise analysis.AnalysisError(
                'ncoeff too large for given nsamples: {0}'.format(self.ncoeff))
        
        if not glVar.sparse:
            print('Warning: dense matrices not supported, using sparse matrices.\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodalwav object
        nd.make_nodal_circuit(circuit)
	nodalcomp = nodalCompressed.CompressedNodal(circuit, self.nsamples,
                                                    self.T, self.wavelet,
                                                    self.deriv, self.ncoeff,
                                                    self.ssfactor)

        # Start from zero 
        x = np.zeros(nodalcomp.dim)

        # Solve equations
        print('System size: ', nodalcomp.dim)
	sV = nodalcomp.get_source()
        # solve equations 
        try: 
            print('Solving system ... ', end='')
            sys.stdout.flush()
            if self.step:
                # Skip direct solution and try conservative
                # convergence helpers directly
                (x, res, iterations) = solve(x, sV,
                                             nodalcomp.convergence_helpers[1:])
            else:
                (x, res, iterations) = solve(x, sV,
                                             nodalcomp.convergence_helpers)
            print('Succeded.\n')
        except NoConvergenceError as ce:
            print('Failed.\n')
            print(ce)
            return

        print('\nFinal coefficients density:',
              100. * np.sum(np.abs(x)>1e-2)/len(x), '%')
        print('Final Jacobian density: ',
              100. * nodalcomp.Jac.nnz / nodalcomp.dim**2,'%')
        
        Wi = nodalcomp.Wi
        # re-shape circuit variables and convert to time domain for display
        xTilde = x.reshape((circuit.nD_dimension, self.nsamples)).T
        xHat = Wi.dot(xTilde)

        # nodalcomp.timeVec contains the time vector
        timeVec = nodalcomp.timeVec
	
        # Save results from analysis ********************************

        # Get terminals to plot/save from circuit. 
        termSetw = circuit.get_requested_terms('sscwcoeff')
        termSet = circuit.get_requested_terms('sscw')

        #import pdb; pdb.set_trace()
        # Special treatment for ground terminal
        if circuit.nD_ref in termSetw:
            termSetw.remove(circuit.nD_ref)
            circuit.nD_ref.sscwcoeff_v = np.zeros(self.nsamples)
        if circuit.nD_ref in termSet:
            termSet.remove(circuit.nD_ref)
            circuit.nD_ref.sscw_v = np.zeros(self.nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.sscwcoeff_v = xTilde[:,term.nD_namRC]
                term.sscw_v = xHat[:,term.nD_namRC]
        else:
            # Only save requested nodes
            for term in termSetw:
                term.sscwcoeff_v = xTilde[:,term.nD_namRC]
            for term in termSet:
                term.sscw_v = xHat[:,term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        # Process output requests.  
        cnVec = np.arange(0, self.nsamples)
        analysis.process_requests(circuit, 'sscwcoeff', 
                                  cnVec, 'Coefficient number', 'sscwcoeff_v',
                                  style = '-o')
        analysis.process_requests(circuit, 'sscw', 
                                  timeVec, 'Time [s]', 'sscw_v')

        def getvec(termname):
            return circuit.find_term(termname).sscw_v

        if self.shell:
            import matplotlib.pylab as plt
            plt.ion()
            analysis.ipython_drop("""
Available commands:
    matplotlib.pylab imported as 'plt'
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())

aClass = SSCW



