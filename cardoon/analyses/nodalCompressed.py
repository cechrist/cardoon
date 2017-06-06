"""
:mod:`nodalCompressed` -- Steady-state Nodal Analysis Using Compressed Sensing
------------------------------------------------------------------------------

.. module:: nodalCompressed
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

This module contains specific classes/functions for steady state nodal
analysis using wavelet coefficients.

"""

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import scipy.sparse.linalg 
from analysis import AnalysisError
from cardoon.globalVars import glVar
from fsolve import fsolve_Newton, NoConvergenceError
from nodalWavelet import WaveletNodal
import csUtil
import pywt

def get_wv(nsamples, levels = None, wtype='db4'):
    # Create weight vector (wv) to give preference to low-subspace coefficients
    wt = pywt.Wavelet('db4')
    maxlev = pywt.dwt_max_level(nsamples, wt)
    if not levels:
        levels = maxlev
    wv = np.zeros(nsamples)
    end = nsamples/2**levels
    base = 2.
    wv[0:end] = 1.
    for j in range(levels):
        beg = end
        end = 2*beg
        wv[beg:end] = base**j
    return wv

#----------------------------------------------------------------------

class CompressedNodal(WaveletNodal):
    """
    Steady-state using compressed wavelets coefficients

    This class only sets up equations. Equations are solved elsewhere.

    For now only the highest wavelet subspace level is compressed.

    Matrices and vectors (G, C, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit()) and number of samples. 
    """

    def __init__(self, ckt, nsamples, T, wavelet, deriv, ncc, sstep):
        """
        Arguments:

        ckt: circuit instance
        nsamples: number of samples in period
        T: period
        wavelet: wavelet type, string
        multilevel: it True, use multilevel Wavelet transform
        deriv: derivative type, string
        ncc: number of compressed coefficients
        sstep: source stepping factor

        """
        # Init base class: Use multilevel wavelet transform to
        # guarantee that the second half of the coefficients
        # constitute the highest subspace level.
        # Set Jacobian matrix format to 'bsr'
        super(CompressedNodal, self).__init__(ckt, nsamples, T, wavelet,
                                              True, deriv, sstep, matformat='bsr')

        self.convergence_helpers = [self.solve_compressed, 
#                                    self.solve_homotopy_source, 
                                    None]

        self.wcw = get_wv(nsamples, wtype=wavelet)
        
        # Number of wavelet coefficients to be compressed
        self.nwcoeff = nsamples / 2
        # Number of compressed coefficients 
        self.m = ncc
        # Create compression matrix
        self.Fi = csUtil.get_C(self.m, self.nwcoeff)
        # To be allocated in factor_and_solve()
        self.data = None
        self.__count = 0

# IMPORTANT: a better idea for the reduced system is to just create
# another wavelet object: problem is that we need to be careful with
# attributes written to ckt instances. A better approach would be to
# allow updating nsamples in NodalWavelet.

    def solve_compressed(self, x0, sV):
        # This is the main flow for compressed analysis
        #
        print("Solving low-res system ...")
        (x, res, iter1, success) = self.solve_reduced(x0, sV)
        print("Iterations: {0}, Residual: {1}".format(iter1, res)) 
        # 2. Refine solution and find support
        print("Finding support with CS ...")
        (x, res, iter2, success) = self.solve_cs(x, sV)
        print("Iterations: {0}, Residual: {1}".format(iter2, res)) 
#        # 3. Further refine solution using fixed support
#        print("Refining solution ...")
#        (x, res, iter3, success) = self.solve_cs(x, sV)
#        print("Iterations: {0}, Residual: {1}".format(iter2, res)) 
        iterations = iter1 + iter2 #+ iter3
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError(
                'No convergence. iter = {0} res = {1}'.format(iterations, res))

    def solve_reduced(self, x0, sV):
        # Approximately solve a low-resolution problem

        # Relax tolerances
        glVar.abstol *= 1e2
        glVar.reltol *= 1e2

        xExpA = np.empty((self.ckt.nD_dimension, self.nsamples))
        xExp = xExpA.ravel()
        def get_deltax(x):
            xExpA[:,:self.nwcoeff] = x.reshape((self.ckt.nD_dimension,
                                                self.nwcoeff))
            (iVec, Jac) = self.get_i_Jac(xExp) 
            return self.factor_and_solve_reduced(sV - iVec, Jac)
        def f_eval(x):
            xExpA[:,:self.nwcoeff] = x.reshape((self.ckt.nD_dimension,
                                                self.nwcoeff))
            iVec = self.get_i(xExp)
            errFunc = (iVec - sV).reshape((self.ckt.nD_dimension,
                                           self.nsamples))[:,
                                                0:self.nwcoeff].ravel()
            return errFunc
        # Reduce size of initial guess vector
        x0red = x0.reshape((self.ckt.nD_dimension,
                            self.nsamples))[:, 0:self.nwcoeff].ravel()
        (x, res, iterations, success) = \
            fsolve_Newton(x0red, get_deltax, f_eval)
        # Expand solution vector
        xExpA[:,:self.nwcoeff] = x.reshape((self.ckt.nD_dimension,
                                            self.nwcoeff))
        # Restore tolerances
        glVar.abstol *= 1e-2
        glVar.reltol *= 1e-2

        if success:
            return (xExp, res, iterations, success)
        else:
            return (xExp, res, iterations, True)
#            raise NoConvergenceError(
#                'No convergence. iter = {0} res = {1}'.format(iterations, res))
        

    def factor_and_solve_reduced(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc using half resolution

        errFunc and Jac are uncompressed 

        Returns the opposite of the Newton correction, but that
        is what fsolve() is expecting (-deltax)
        """
        # Create reduced system
        # No need to allocate memory, just use a view of existing Jac
        data = Jac.data[:, :self.nwcoeff, :self.nwcoeff]
        sysShape = (self.dim/2, self.dim/2)
        reducedJac = sp.bsr_matrix((data,
                                    Jac.indices, Jac.indptr),
                                   shape = sysShape).tocsc()
        # One variable per row
        reducedRhs = errFunc.reshape((self.ckt.nD_dimension,
                                         self.nsamples))[:,
                                                         0:self.nwcoeff].ravel()
        factorized = sp.linalg.factorized(reducedJac)
        return factorized(reducedRhs)


    def solve_cs(self, x0, sV):
        # Approximately solve using cs
        # Relax tolerances
        looseFactor = 1.
        #glVar.abstol = .005
        #glVar.reltol = 1e-3
        #glVar.maxiter = 10

        blocksize = self.nwcoeff + self.m
        def get_deltax(x):
            #import pdb; pdb.set_trace()
            x1 = x.copy()
            xRef = x.copy()
            delta_xc = x.copy()
            # Array view of input vector
            xA = x1.reshape((self.ckt.nD_dimension, blocksize))
            xARef = xRef.reshape((self.ckt.nD_dimension, blocksize))
            xA2 = delta_xc.reshape((self.ckt.nD_dimension, blocksize))
            # Compressed-coefficient-only view
            xcC = xA[:,self.nwcoeff:]
            xcCo = xcC.copy()
            xcCRef = xARef[:,self.nwcoeff:]
            xcC2 = xA2[:,self.nwcoeff:]
            # Array view of previously allocated vector
            deltaxA = self.deltaxVec.reshape((self.ckt.nD_dimension,
                                              self.nsamples))
            deltaxA[:, 0:self.nwcoeff] = xA[:, 0:self.nwcoeff]
            # Decompress coefficients
            xHigh = csUtil.recover(self.Fi, xcC.T)
            # Copy into deltax
            deltaxA[:, self.nwcoeff:] = xHigh.T
            (iVec, Jac) = self.get_i_Jac(self.deltaxVec)
            b1 = sV - iVec
            # Solve decompressed system
            factorized = sp.linalg.factorized(Jac.tocsc())
            self.deltaxVec[:] = factorized(b1)
            # Compress solution
            xARef[:, 0:self.nwcoeff] = deltaxA[:, 0:self.nwcoeff]
            xcCRef[:,:] = np.dot(self.Fi, deltaxA[:, self.nwcoeff:].T).T

            b2 = Jac.dot(self.deltaxVec)
            b3 = b1
            x1[:] = 0.
            for i in range(2):
                delta_xc[:] = self.factor_and_solve_compressed(Jac, b3,
                                                               xcC + xcCo,
                                                               blocksize)
                err = max(abs(delta_xc))
                x1 += delta_xc
                # Decompress coefficients
                xHigh = csUtil.recover(self.Fi, xcC2.T)
                # Copy into deltax
                deltaxA[:, 0:self.nwcoeff] = xA2[:, :self.nwcoeff]
                deltaxA[:, self.nwcoeff:] = xHigh.T
                # Calculate error
                b3 -= Jac.dot(self.deltaxVec)
                nb3 = max(abs(b3))
                print('delta_xc:', err, 'b3: ', nb3)
                if nb3 < .01:
                    break;

            #x1[:] = self.factor_and_solve_compressed(Jac, b1, xcC, blocksize)
            print(max(abs(xRef)) , max(abs(x1 - xRef)))
            return x1

        def f_eval(x):
            # Array view of input vector
            xA = x.reshape((self.ckt.nD_dimension, blocksize))
            # Compressed-coefficient-only view
            xcC = xA[:,self.nwcoeff:]
            # Array view of previously allocated vector
            deltaxA = self.deltaxVec.reshape((self.ckt.nD_dimension,
                                              self.nsamples))
            deltaxA[:, 0:self.nwcoeff] = xA[:, 0:self.nwcoeff]
            # Decompress coefficients
            xHigh = csUtil.recover(self.Fi, xcC.T)
            # Copy into deltax
            deltaxA[:, self.nwcoeff:] = xHigh.T
            iVec = self.get_i(self.deltaxVec) 
            return iVec - sV

        x0A = x0.reshape((self.ckt.nD_dimension, self.nsamples))
        x0C = x0A[:,:blocksize].ravel()
        (x, res, iterations, success) = fsolve_Newton(x0C, get_deltax, f_eval)
        xA = x.reshape((self.ckt.nD_dimension, blocksize))
        xcC = xA[:,self.nwcoeff:]
        # Array view of previously allocated vector
        deltaxA = self.deltaxVec.reshape((self.ckt.nD_dimension, self.nsamples))
        deltaxA[:, 0:self.nwcoeff] = xA[:, 0:self.nwcoeff]
        # Decompress coefficients
        xHigh = csUtil.recover(self.Fi, xcC.T)
        # Copy into deltax
        deltaxA[:, self.nwcoeff:] = xHigh.T        
        # Restore tolerances
        glVar.abstol /= looseFactor
        glVar.reltol /= looseFactor

        if success:
            return (self.deltaxVec, res, iterations, success)
        else:
            return (self.deltaxVec, res, iterations, True)
#            raise NoConvergenceError(
#                'No convergence. iter = {0} res = {1}'.format(iterations, res))


    def factor_and_solve_compressed(self, Jac, b1, xcC, blocksize):
        """
        Solves a compressed system

        b1: rhs vector
        xcC: array view of current compressed coefficients (one per row)
        """
        # Allocate memory for compressed Jacobians only if needed
        if self.data == None:
            # Allocate and save for next iteration
            #self.data = np.empty((Jac.data.shape[0], blocksize, blocksize))
            self.data = np.zeros((Jac.data.shape[0], self.nsamples, blocksize))
        # 2. Find decompression Jacobian for each variable
        nvars = self.ckt.nD_dimension
        # Number of samples to be compressed
        ncsamples = self.nsamples - self.nwcoeff
        Jc = np.empty((nvars, ncsamples, self.m))
        varSupport = np.empty((nvars, ncsamples), dtype=bool)
        for i in range(0, nvars):
            Jc[i] = csUtil.get_Jc(self.Fi, xcC[i])
            r1 = abs(Jc[i]).max(axis=1)
            varSupport[i] = (r1 != 0)

        # 3. Find required support for each block row
        rowSupport = np.zeros((nvars, ncsamples), dtype=bool)
        j = 0
        for i,k in enumerate(Jac.indices):
            if i == Jac.indptr[j+1]:
                j += 1
            rowSupport[j] += varSupport[k]
           
        # 4. Setup linear system
        # Right-hand-side vector
        reducedb1 = np.empty(nvars * blocksize)
        b1A = b1.reshape((nvars, self.nsamples))
        rb1A = reducedb1.reshape((nvars, blocksize))
        j = -1
        block = np.zeros((self.nsamples, blocksize))
        # Left compression matrix for this row
        Fij = np.zeros((self.m, ncsamples))
        for i,k in enumerate(Jac.indices):
            if i == Jac.indptr[j+1]:
                # Go to next row
                j += 1
                Fij.fill(0.)
                Fij[:,rowSupport[j]] = self.Fi[:,rowSupport[j]]
                # Fill entry in reducedb1
                rb1A[j,:self.nwcoeff] = b1A[j,:self.nwcoeff]
                rb1A[j,self.nwcoeff:] = np.dot(Fij, b1A[j,self.nwcoeff:])
            # Multiply in blocks: first post-mult
            # Identity part
            block[:,:self.nwcoeff] = Jac.data[i,:,:self.nwcoeff]
            # Decompression Jacobian
            block[:,self.nwcoeff:] = np.dot(Jac.data[i,:,self.nwcoeff:],
                                            Jc[k])
            self.data[i,:,:] = block
            # Now pre-mult by compression matrix
            #self.data[i,:self.nwcoeff,:] = block[:self.nwcoeff,:]
            #self.data[i,self.nwcoeff:,:] = np.dot(Fij,
            #                                      block[self.nwcoeff:,:])
            
        sysShape = (nvars*self.nsamples, nvars*blocksize)
        reducedJac = sp.bsr_matrix((self.data,
                                    Jac.indices, Jac.indptr),
                                   shape = sysShape).tocsc()
        #reducedb1 = reducedJac.T.dot(b1)
        #reducedJac = reducedJac.T.dot(reducedJac).tocsc()

        result = scipy.sparse.linalg.lsmr(reducedJac, b1,
                                          atol=1e-8, btol=1e-8, maxiter=2000)
#        print('code:',result[1],'iter:',result[2])
#        #result = la.lstsq(Jac.todense(), b1)
        return result[0]

####        return np.dot(la.pinv(reducedJac.todense()), b1)
#        try:
#            factorized = sp.linalg.factorized(reducedJac)
#            x = factorized(reducedb1)
#        except:
#            result = scipy.sparse.linalg.lsmr(reducedJac, reducedb1,
#                                          atol=1e-7, btol=1e-7, maxiter=1000)
#            print('code:',result[1],'iter:',result[2])
#            x = result[1]
#        return x



        
#----------------------------------------------------------------------
#  
#  
#      
#          self.deltaxVec.fill(0.)
#          deltaxA = self.deltaxVec.reshape((self.ckt.nD_dimension,
#                                               self.nsamples))
#          deltaxA[:, 0:self.nwcoeff] = x.reshape((self.ckt.nD_dimension,
#                                                     self.nwcoeff))
#  
#          # Calculate residual for high subspace levels
#          b1 = Jac.dot(self.deltaxVec) - errFunc
#          self.__count += 1
#          if max(abs(b1)) > 1e-6 and self.__count > 50:            
#              blocksize = self.nwcoeff + self.m
#              xc = np.zeros(self.ckt.nD_dimension * blocksize)
#              # Array view (row-oriented)
#              xcA = xc.reshape(self.ckt.nD_dimension, blocksize)
#              # Copy initial guess
#              xcA[:, :self.nwcoeff] = deltaxA[:, 0:self.nwcoeff] 
#              # Compressed-coefficient-only view
#              xcC = xcA[:,self.nwcoeff:]
#              #import pdb; pdb.set_trace()
#              # Use fixed number of iterations for now
#              factor = 1.
#              for i in range(10):
#                  delta_xc = self.solve_compressed(Jac, factor * b1,
#                                                   xcC, blocksize)
#                  xc -= delta_xc / factor
#                  # Decompress coefficients
#                  xHigh = csUtil.recover(self.Fi, xcC.T)
#                  # Copy into deltax
#                  deltaxA[:, 0:self.nwcoeff] = xcA[:, :self.nwcoeff]
#                  deltaxA[:, self.nwcoeff:] = xHigh.T
#                  # Calculate error
#                  b1 = Jac.dot(self.deltaxVec) - errFunc
#                  err1 = la.norm(b1)
#                  err2 = la.norm(delta_xc) / factor
#                  if err2 < 10.:
#                      factor = min(10. / err2, 1e3)
#                  print(err1, err2, factor)
#                  if i>1 and (err1 < .05 or err2 < 1e-2):
#                      break
#  
#          print('-------------------------')
#          return self.deltaxVec
#  
  
  



        

#        fullSupport = np.zeros(ncsamples, dtype=bool)
#        for sup in varSupport:
#            fullSupport += sup
#        Fij[:,fullSupport] = self.Fi[:,fullSupport]
