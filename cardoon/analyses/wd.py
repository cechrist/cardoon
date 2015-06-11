"""
:mod:`wd` -- Wavelet Domain steady-state analysis
-------------------------------------------------------------------

.. module:: wd
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
import nodalWD
import analysis
import CSFunctions as CSF

# Valid request types
reqTypes = ['wd']

class WD(ParamSet):
    r"""
    WD Steady-State Analysis
    --------------------------

    This is an experimental analysis. It solves the steady-state
    problem using a Wavelet Transform approach.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``. The type of matrix used in
    this analysis is controlled by the ``sparse`` option. Global
    options are documented in :doc:`global_vars`. 

    Example::

        .analysis wd T = 1us nsamples = 100  

        .plot wd in out

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "wd"

    # Define parameters as follows
    paramDict = dict(
        T = ('Fundamental period (overrides frequency)', 's', float, 0.),
        f = ('Fundamental frequency', 'Hz', float, 0.),
        m = ('Desired Number of Wavelet Domain Coefficients', '', float, 0.1),
        nsamples = ('Number of samples in period (must be a power of 2 and > 5)', '', int, 100),
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
        print('*******************************************************')
        print('                Wavelet Domain analysis                ')
        print('*******************************************************')
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
        elif (np.log(self.nsamples)/np.log(2)).is_integer() == False:
            raise analysis.AnalysisError(
                'Number of samples must be a power of 2')
            return 0
        if not glVar.sparse:
            print('Warning: dense matrices not supported in this analysis, using sparse matrices.\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodal and WD objects
        nd.make_nodal_circuit(circuit)
	wd = nodalWD.WDNodal(circuit, self.nsamples, self.T, self.m)
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
            x = np.zeros(wd.dim)

        # Solve equations
        print('System size: ', wd.dim)
	sV = wd.get_source()
        # solve equations 
        try: 
            print('Solving system ... ', end='')
            sys.stdout.flush()
            
            #xcomp = np.empty(cs.csdim)
            #for r in range(circuit.nD_dimension):
            #    xcomp[r*cs.m:r*cs.m+cs.m] = np.dot(cs.PSIWM, x[r*self.nsamples:r*self.nsamples+self.nsamples])            
            #x = np.copy(xcomp)
            #del xcomp
            x = np.dot(wd.IkWM.todense(),x).A1
            #import pdb; pdb.set_trace()
            print(x.shape)
            (x, res, iterations) = solve(x, sV, wd.convergence_helpers)
            print('Succeded.\n')
        except NoConvergenceError as ce:
            print('Failed.\n')
            print(ce)
            return

        # re-order circuit variables (samples for one variable in one
        # block)
        p, pt = nodalWD.calculate_P_PT(circuit.nD_dimension, self.nsamples)
        
        
        #x = cs.sVecuncomp
        #xcomp = np.empty(cs.csdim)        


        #x = np.dot(cs.p.todense(), x)
        #x = x.A1
        #for res in range(circuit.nD_dimension):
        #    xcomp[res*cs.m:res*cs.m+cs.m] = np.dot(cs.PSIWM, x[res*self.nsamples:res*self.nsamples+self.nsamples])
        #for r in range(circuit.nD_dimension):
        #    x[r*self.nsamples:r*self.nsamples+self.nsamples] = np.dot(cs.PSIWMinv, xcomp[r*cs.m:r*cs.m+cs.m])
        #x = np.dot(cs.pt.todense(), x)
        #x = x.A1
        
        #x =   np.reshape(x[pt], (circuit.nD_dimension* self.nsamples))        
        
        #import CSFunctions as csf
        #x, SupportMatrix = csf.CSRecover (x, circuit.nD_dimension, cs.PSIWM, True, SearchSpeed = 50)
        x = wd.IkWMinv*x        
        #import pdb; pdb.set_trace()
        # Re-shape so that xHat[1] returns all samples for nodal variable 1
        xHat = np.reshape(x[p], (circuit.nD_dimension, self.nsamples))
        #xHat = np.reshape(x[p], (circuit.nD_dimension, self.nsamples))
        #xHat = np.reshape(sV, (circuit.nD_dimension, self.nsamples))
        xMatrix = xHat*0
        for i in range(wd.xMatrix.shape[0]):
            xMatrix = np.vstack((xMatrix,np.reshape(wd.xMatrix[i,:][p], (circuit.nD_dimension, self.nsamples))[1,:]))

        #xHat[1] = np.dot(cs.PSIWMinv, np.dot(cs.PSIWM, xHat[1]))

        
        # cs.timeVec contains the time vector
        timeVec = wd.timeVec
	
        # Save results from cs analysis ********************************

        # Get terminals to plot/save from circuit. 
        termSet = circuit.get_requested_terms('wd')

        # Special treatment for ground terminal
        termSet1 = set(termSet)
        if circuit.nD_ref in termSet1:
            termSet1.remove(circuit.nD_ref)
            circuit.nD_ref.cs_v = np.zeros(self.nsamples)

        # Allocate vectors for results
        if self.saveall:
            for term in circuit.nD_termList:
                term.cs_v = xHat[term.nD_namRC]
        else:
            # Only save requested nodes
            for term in termSet1:
                term.cs_v = xHat[term.nD_namRC]

        print('Number of iterations = ', iterations)
        print('Residual = ', res)
        import matplotlib.pylab as plt
        plt.ion()
        plt.figure(1)
        for i in range(xMatrix.shape[0]):
            plt.plot(timeVec, xMatrix[i,:])
#        import pdb; pdb.set_trace()
        np.savetxt('wdoutputs.txt', xMatrix, delimiter=',')
        # Process output requests.  
        analysis.process_requests(circuit, 'wd', 
                                  timeVec, 'Time [s]', 'cs_v')

        def getvec(termname):
            return circuit.find_term(termname).cs_v

        if self.shell:
            import matplotlib.pylab as plt
            plt.ion()
            analysis.ipython_drop("""
Available commands:
    matplotlib.pylab imported as 'wd'
    timeVec: time vector
    getvec(<terminal>) to retrieve results (if result saved)
""", globals(), locals())

aClass = WD



