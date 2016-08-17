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
from analysis import AnalysisError
from cardoon.globalVars import glVar
from nodalWavelet import WaveletNodal
import csUtil



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

        # Number of wavelet coefficients to be compressed
        self.nwcoeff = nsamples / 2
        # Number of compressed coefficients 
        self.m = ncc
        # Create compression matrix
        self.Fi = csUtil.get_C(self.m, self.nwcoeff)
        # To be allocated in factor_and_solve()
        self.data = None
        

    def factor_and_solve(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc

        errFunc and Jac are uncompressed 

        So it returns the opposite of the Newton correction, but that
        is what fsolve() is expecting (-deltax)

        The idea is to trick the nonlinear solver to think that we 
        are solving the uncompressed system
        """
        # Steps:
        #
        # 1. Solve reduced system (preconditioner), substract result
        # to main system.

        #import pdb; pdb.set_trace()
        # Create reduced system
        # No need to allocate memory, just use a view of existing Jac
        data = Jac.data[:, :self.nwcoeff, :self.nwcoeff]
        sysShape = (self.dim/2, self.dim/2)
        reducedJac = sp.bsr_matrix((data,
                                    Jac.indices, Jac.indptr),
                                   shape = sysShape).tocsc()
        # One variable per row
        reducedRhs = np.empty(self.ckt.nD_dimension * self.nwcoeff)
        reducedRhs[:] = errFunc.reshape((self.ckt.nD_dimension,
                                         self.nsamples))[:,
                                                         0:self.nwcoeff].ravel()
        factorized = sp.linalg.factorized(reducedJac)
        x = factorized(reducedRhs)

        self.deltaxVec.fill(0.)
        deltaxview = self.deltaxVec.reshape((self.ckt.nD_dimension,
                                             self.nsamples))
        deltaxview[:, 0:self.nwcoeff] = x.reshape((self.ckt.nD_dimension,
                                                   self.nwcoeff))

        # Calculate residual for high subspace levels
        b1 = Jac.dot(self.deltaxVec) - errFunc

        blocksize = self.nwcoeff + self.m
        # Allocate memory for compressed Jacobians only if needed
        if not self.data:
            # Allocate and save for next iteration
            self.data = np.empty((Jac.data.shape[0], blocksize, blocksize))

        xc = np.zeros(self.ckt.nD_dimension * blocksize)
        # Array view (row-oriented)
        xcA = xc.reshape(self.ckt.nD_dimension, blocksize)
        # Compressed-coefficient-only view
        xcC = xcA[:,self.nwcoeff:]
        # Use fixed number of iterations for now
        for i in range(3):
            delta_xc = self.solve_compressed(Jac, b1, xcC, blocksize)
            xc -= delta_xc
            # Decompress coefficients
            xHigh = csUtil.recover(self.Fi, xcC.T)
            # Copy into deltax
            deltaxview[:, 0:self.nwcoeff] = xcA[:, :self.nwcoeff]
            deltaxview[:, self.nwcoeff:] = xHigh.T
            # Calculate error
            b1 = Jac.dot(self.deltaxVec) - errFunc

        print(la.norm(b1))
        return self.deltaxVec


    def solve_compressed(self, Jac, b1, xcC, blocksize):
        """
        Solves a compressed system

        b1: rhs vector
        xcC: array view of current compressed coefficients (one per row)
        """
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
        block = np.empty((self.nsamples, blocksize))
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
            # Now pre-mult by compression matrix
            self.data[i,:self.nwcoeff,:] = block[:self.nwcoeff,:]
            self.data[i,self.nwcoeff:,:] = np.dot(Fij, block[self.nwcoeff:,:])
            
        sysShape = (nvars*blocksize, nvars*blocksize)
        reducedJac = sp.bsr_matrix((self.data,
                                    Jac.indices, Jac.indptr),
                                   shape = sysShape).tocsc()

        factorized = sp.linalg.factorized(reducedJac)
        x = factorized(reducedb1)

        return x




        

