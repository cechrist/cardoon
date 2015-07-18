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
import scipy.linalg as linalg
from nodal import set_xin, set_i, _NLFunction
import nodalWavelet as nw
import nodalSP as nd
from analysis import AnalysisError
from cardoon.globalVars import glVar
import spams

#----------------------------------------------------------------------

class CompressedNodal(nd._NLFunctionSP):
    """
    Steady-state equations using Compressed Wavelet Coefficients

    This class only sets up equations. Equations are solved elsewhere.

    Matrices and vectors (G, C, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit()) and number of samples. 

    Depends on NodalWavelet
    """

    def __init__(self, ckt, nsamples, ncc, T, wavelet, deriv, sstep):
        """
        Arguments:

        ckt: circuit instance
        nsamples: number of samples in period
        ncc: number of compressed coefficients
        T: period
        wavelet: wavelet type, string
        deriv: derivative type, string
        sstep: source stepping factor

        """
        # Init base class
        super(CompressedNodal, self).__init__()
        # Source stepping factor
        self._sstep = sstep
        # Override convergence helpers (generic gmin does not work with this)
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_source,
                                    self.solve_homotopy_source2, 
                                    None]

        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref
        # Save circuit
        self.ckt = ckt
        self.nsamples = nsamples
        self.ncc = ncc
        self.dim = self.ckt.nD_dimension*self.ncc
        # T is the period
        self.T = T
        # Allocate matrices/vectors
        # G and C (done in refresh())

        # Vectors are allocated as matrices, each row is the set of
        # samples for one variable. Example: iVec[2] is the set of
        # samples for the second current in iVec.

        # iVec: total compressed current
        self.iVecA = np.zeros((self.ckt.nD_dimension, self.ncc))
        self.iVec = self.iVecA.ravel()
        # nonlinear current only
        self.inlArray = np.zeros((self.nsamples, self.ckt.nD_dimension))
        # Nonlinear charge only
        self.qArray = np.zeros((self.nsamples, self.ckt.nD_dimension))
        # Decompressed Source vector 
        self.sVecA = np.zeros((self.ckt.nD_dimension, self.nsamples))

        # Delta_x is required as workspace by the nonlinear solver
        self.deltaxVec = np.zeros(self.ckt.nD_dimension * self.ncc)

        # Create time vector
        self.timeVec = np.linspace(0., T, self.nsamples, endpoint=False)

        # TODO: include time delays
        #self.tdVecList = []

        # Create derivative matrix (need in dense format)
        if deriv=='d2':
            D = nw.deriv_2(self.nsamples, self.T).toarray()
        elif deriv=='d4':
            D = nw.deriv_4(self.nsamples, self.T)
        elif deriv=='Fourier':
            D = nw.deriv_fourier(self.nsamples, self.T)
        else:
            raise AnalysisError(
                'Invalid deriv value: {0}'.format(deriv))
            
        # Create compression matrix
        Dc = np.asfortranarray(np.random.normal(scale=1.,
                                                size=(self.ncc, self.nsamples)))
        # Norm of columns of Dc must be 1
        self.Dc = np.asfortranarray(Dc / np.tile(np.sqrt((Dc*Dc).sum(axis=0)),
                                                 (Dc.shape[0],1)))
        # for testing purposes ...
        #self.Dc /= self.ncc
        #self.Dc += np.asfortranarray(np.eye(self.ncc, self.nsamples))
        
        # Wavelet matrices. It may be possible to avoid its creation
        # (using the fast wavelet transform) but for the moment is it
        # not clear that we would save much with that.
        if wavelet == 'none':
            raise AnalysisError(
                'Compressed analysis requires a Wavelet transform')
        else:
            # Use multilevel transforms (W dense, Wi sparse)
            W = nw.wavelet_f(self.nsamples, wavelet)
            self.Wi = sp.csr_matrix(nw.wavelet_i(self.nsamples, wavelet))

        self.DcW = np.dot(self.Dc, W)
        self.DcWD = np.dot(self.Dc, np.dot(W, D))

        # The following used in get_Jac_i only
        self.mtmp = np.empty((self.nsamples, self.ncc))
                             
        # Signal reconstruction parameters
        # w: wavelet weight coefficients
        maxlev = nw.wavelet_levels(self.nsamples, wavelet)
        wv = np.zeros(self.nsamples)
        end = self.nsamples/2**maxlev
        base = 2.
        wv[0:end] = base**(-maxlev+1)
        for j in range(maxlev):
            beg = end
            end = 2*beg
            wv[beg:end] = base**(-maxlev+j+1)
        nsig = self.ckt.nD_dimension
        wArray = np.asfortranarray(np.repeat(wv,nsig).reshape((self.nsamples,
                                                               nsig)))
        self.lassoParam = {
            # This is a waste, but lassoWeighted expects it
            'W' : wArray,
            'lambda1' : 1e-10,
            # -1 means use as many cores as possible
            'numThreads' : -1,
            'mode' : 1
            }

        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        # set this flag to request adaptive recovery matrix 
        self.adaptiveDcRec = False
        # Create decompression matrix list (all pointing to same matrix)
        self.WiDcRec = self.Wi.dot(linalg.pinv2(self.Dc))
        self.WiDcRecList = self.ckt.nD_dimension * [self.WiDcRec]

        # Big matrices with rectangular blocks (current compression
        # only) named *Hat. Full compression: *HatC

        # GTriplet stores Jacobian matrix for a single time sample
        # Format is (data, row, col)
        GTriplet = ([], [], [])
        # Insert linear contributions: Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                nd.set_quad(GTriplet, *vccs)
        self.G = sp.coo_matrix((GTriplet[0], GTriplet[1:]),
                               (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                               dtype = float).tocsr()
        # Create matrix for linear part of circuit
        eyeNsamples = sp.eye(self.nsamples, self.nsamples, format='csr')
        # Note if the pseudo-inverse is used for recovery, blocks are
        # diagonal, but for flexibility, let's assume decompression
        # matrix is different
        GHat = sp.kron(self.G, self.DcW, 'bsr')

        # CTriplet stores Jacobian matrix for a single time sample
        CTriplet = ([], [], [])
        # Add C matrix
        for elem in self.ckt.nD_elemList:
            for vcqs in elem.nD_linVCQS:
                nd.set_quad(CTriplet, *vcqs)
        self.C = sp.coo_matrix((CTriplet[0], CTriplet[1:]),
                               (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                               dtype = float).tocsr()
        # Rectangular dense matrix
        CHat = sp.kron(self.C, self.DcWD, 'bsr')
        
        # Frequency-defined elements: not the optimum way to create
        # the matrix but should do for now. The underlying assumption
        # is that matrix structure is the same for all frequencies.
        ycontrib = False
        if self.ckt.nD_freqDefinedElem:
            ycontrib = True
            # Frequency vector (nsamples is even)
            nfreq = self.nsamples // 2 + 2
            fvec = np.arange(0, nfreq) / self.T            
            YTriplet = ([], [], [])
            # Special case for DC values
            for elem in self.ckt.nD_freqDefinedElem:
                nd.set_Jac_triplet(YTriplet, elem.nD_fpos, elem.nD_fneg, 
                                   elem.nD_fpos, elem.nD_fneg,
                                   elem.get_G_matrix())
            Y0 = sp.coo_matrix((YTriplet[0], YTriplet[1:]),
                               (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                               dtype = complex).tocsr()
            ydata = np.empty((Y0.nnz, nfreq), dtype=complex)
            ydata[:,0] = Y0.data
            # Create additional matrices for other frequencies
            n = 1
            for f in fvec[1:]:
                YTriplet = ([], [], [])
                for elem in self.ckt.nD_freqDefinedElem:
                    nd.set_Jac_triplet(YTriplet, elem.nD_fpos, elem.nD_fneg, 
                                       elem.nD_fpos, elem.nD_fneg,
                                       elem.get_Y_matrix(f))
                Yn = sp.coo_matrix((YTriplet[0], YTriplet[1:]),
                                  (self.ckt.nD_dimension,
                                   self.ckt.nD_dimension), 
                                   dtype = complex).tocsr()
                # Save data vector only
                ydata[:,n] = Yn.data
                n += 1
            # import pdb; pdb.set_trace()
            # Create convolution matrices
            # Inverse-transform Y: for now use numpy. Later switch to pyfftw
            Yblocks = np.empty((Y0.nnz, self.ncc, self.nsamples))
            Ytmp1 = np.empty((self.nsamples, self.nsamples))
            for Yv, Yb in zip(ydata, Yblocks):
                yv = np.fft.irfft(Yv, self.nsamples)
                # Reverse order of elements in yv
                yi = yv[::-1]
                # Create convolution matrix
                Ytmp1.fill(0.)
                for i in range(0, self.nsamples):
                    Ytmp1[i] = np.roll(yi, i+1)
                # Add Compression and Wavelet Transform
                np.dot(self.DcW, Ytmp1, out=Yb)

            YHat = sp.bsr_matrix((Yblocks, Y0.indices, Y0.indptr),
                                 shape=(self.dim,
                                        self.ckt.nD_dimension * self.nsamples))

        # Linear response matrix
        if ycontrib:
            self.MHat = GHat + CHat + YHat
        else:
            self.MHat = GHat + CHat
        # Compressed linear matrix (default recovery)
        self.MHatC = self.MHat.dot(sp.block_diag(self.WiDcRecList)).tobsr(
            blocksize=(self.ncc, self.ncc))

            
        # Create nonlinear structure for one block
        JiTriplet = ([], [], []) 
        JqTriplet = ([], [], []) 

        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # The values that we insert do not matter, we are just
            # interested in the structure
            outJac = np.empty((len(elem.csOutPorts), elem.nD_nxin),
                              dtype = float)
            qJac = np.empty((len(elem.qsOutPorts), elem.nD_nxin),
                            dtype = float)
            nd.set_Jac_triplet(JiTriplet, elem.nD_cpos, elem.nD_cneg, 
                               elem.nD_vpos, elem.nD_vneg, outJac)
            nd.set_Jac_triplet(JqTriplet, elem.nD_qpos, elem.nD_qneg, 
                               elem.nD_vpos, elem.nD_vneg, qJac)

        # Create CSR matrices used to find the nonzero structure for
        # the corresponding BSR *-hat.
        self.Ji = sp.coo_matrix((JiTriplet[0], JiTriplet[1:]),
                                (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                                dtype = float).tocsr()
        self.Jq = sp.coo_matrix((JqTriplet[0], JqTriplet[1:]),
                                (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                                dtype = float).tocsr()

        # Create vectors and dictionary to store nonlinear response.
        # Access to [i,k] vector: data[dl[i][j]]: this is achieved
        # using a list for dictionaries (dl), the list index is the
        # row and the dictionary keys are the column numbers.
        self.Ji.nw_dl = []
        for i in range(0, self.Ji.shape[0]):
            d=dict()
            baseIdx = self.Ji.indptr[i]
            endIdx = self.Ji.indptr[i+1]
            for j,k in enumerate(self.Ji.indices[baseIdx:endIdx]):
                d[k] = baseIdx + j
            self.Ji.nw_dl.append(d)
        # Reserve space for entries in time domain (diagonal matrices)
        self.Ji.nw_dataDiag = np.zeros((self.Ji.nnz, self.nsamples))
        # Use dense matrices for compressed blocks
        self.Ji.nw_dataBlock = np.empty((self.Ji.nnz, self.ncc, self.ncc))
        
        # Same for Jq
        self.Jq.nw_dl = []
        for i in range(0, self.Jq.shape[0]):
            d=dict()
            baseIdx = self.Jq.indptr[i]
            endIdx = self.Jq.indptr[i+1]
            for j,k in enumerate(self.Jq.indices[baseIdx:endIdx]):
                d[k] = baseIdx + j
            self.Jq.nw_dl.append(d)
        # Reserve space for entries in time domain (diagonal matrices)
        self.Jq.nw_dataDiag = np.zeros((self.Jq.nnz, self.nsamples))
        # Use dense matrices for compressed blocks
        self.Jq.nw_dataBlock = np.empty((self.Jq.nnz, self.ncc, self.ncc))

        # Create help arrays in elements
        for elem in self.ckt.nD_nlinElem:
            pos, neg = nw.get_Jac_idx(self.Ji, elem,
                                      elem.nD_cpos, elem.nD_cneg, 
                                      elem.nD_vpos, elem.nD_vneg)
            elem.nw_posJi = pos
            elem.nw_negJi = neg
            pos, neg = nw.get_Jac_idx(self.Jq, elem,
                                      elem.nD_qpos, elem.nD_qneg, 
                                      elem.nD_vpos, elem.nD_vneg)
            elem.nw_posJq = pos
            elem.nw_negJq = neg


    def get_source(self):
        """
        Get the source vector considering DC and TD source components
        """
        # Erase vector first. sVecA is the array view of the same thing
        self.sVecA.fill(0.)
        for elem in self.ckt.nD_sourceDCElem:
            for j, ctime in enumerate(self.timeVec):
                # first get the destination row/columns 
                outTerm = elem.nD_sourceOut
                current = elem.get_DCsource()
                # This may not need optimization because we usually do not
                # have too many independent sources
                # import pdb; pdb.set_trace()
                if outTerm[0] >= 0:
                    self.sVecA[outTerm[0], j] -= current
                if outTerm[1] >= 0:
                    self.sVecA[outTerm[1], j] += current
        for elem in self.ckt.nD_sourceTDElem:
            for j, ctime in enumerate(self.timeVec):            
                # first get the destination row/columns 
                outTerm = elem.nD_sourceOut
                current = elem.get_TDsource(ctime)
                # This may not need optimization because we usually do not
                # have too many independent sources
                # import pdb; pdb.set_trace()
                if outTerm[0] >= 0:
                    self.sVecA[outTerm[0], j] -= current
                if outTerm[1] >= 0:
                    self.sVecA[outTerm[1], j] += current
        # Convert to Wavelet Domain
        s = np.dot(self.DcW, self.sVecA.T)
        return s.T.flatten()
        

    def get_rhs(self):
        """
        Returns system rhs vector: s' 

        s' = source current - currents generated by other circuit components

        """
        # TODO: add convolution
        return self.get_source() - self.get_i()

    def get_i(self, xVec):
        """
        Calculate total current

        returns currents generated by non-sources

        xVec: input vector of nodal voltages (compressed). 

        iVec: output vector of currents (compressed)
        """
        #import pdb; pdb.set_trace()
        # Reconstruct signals
        xHat = self.recover(xVec, updateMHatC = False)
        # Linear contribution
        self.iVec[:] = self.MHat.dot(xHat.T.ravel())
        
        # Nonlinear contribution
        # Erase arrays
        self.inlArray.fill(0.)
        self.qArray.fill(0.)
        for elem in self.ckt.nD_nlinElem:
            for j in range(self.nsamples):
                # first have to retrieve port voltages from xVec
                xin = np.zeros(elem.nD_nxin)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xHat[j,:])
                # We can not yet treat delays
                assert elem.nDelays == 0
                outV = elem.eval(xin)
                # Update iVec, qVec
                set_i(self.inlArray[j], elem.nD_cpos, elem.nD_cneg, outV)
                set_i(self.qArray[j], elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
        # iVec = M x + DcW i(x) + DcW D q(x)
        self.iVecA += np.dot(self.DcW, self.inlArray).T
        self.iVecA += np.dot(self.DcWD, self.qArray).T
        return self.iVec

        
    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)

        xVec: input vector of compressed nodal voltages. 

        iVec: output vector of compressed currents

        Jac: system Jacobian
        """
        # Recover signals and matrix
        xHat = self.recover(xVec, updateMHatC = True)
        # Linear contribution
        self.iVec[:] = self.MHat.dot(xHat.T.ravel())
        # Alternative Linear contribution (not sure if better):
        # self.iVec[:] = self.MHatC.dot(xVec)

        # Nonlinear contribution
        # Erase arrays
        self.inlArray.fill(0.)
        self.qArray.fill(0.)
        self.Ji.nw_dataDiag.fill(0.)
        self.Jq.nw_dataDiag.fill(0.)
        for elem in self.ckt.nD_nlinElem:
            xin = np.zeros(elem.nD_nxin)
            for j in range(self.nsamples):
                # first have to retrieve port voltages from xVec
                xin.fill(0.)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xHat[j,:])
                # We can not yet treat delays
                assert elem.nDelays == 0
                (outV, outJac) = elem.eval_and_deriv(xin)
                # Update iVec, qVec
                set_i(self.inlArray[j], elem.nD_cpos, elem.nD_cneg, outV)
                set_i(self.qArray[j], elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
                qJac = outJac[len(elem.csOutPorts):,:]
                # Set Jacobian entries for all samples at once
                nw.set_Jac(self.Ji.nw_dataDiag,
                           elem.nw_posJi, elem.nw_negJi, outJac, j)
                nw.set_Jac(self.Jq.nw_dataDiag,
                           elem.nw_posJq, elem.nw_negJq, qJac, j)
        # iVec = M x + DcW i(x) + DcW D q(x)
        self.iVecA += np.dot(self.DcW, self.inlArray).T
        self.iVecA += np.dot(self.DcWD, self.qArray).T
        # Form system Jacobian
        #import pdb; pdb.set_trace()
        for i in range(self.Ji.nnz):
            # Retrieve recovery matrix (different for each variable)
            WiDcRec = self.WiDcRecList[self.Ji.indices[i]]
            # Use element-wise multiplication for diagonal matrix
            np.multiply(WiDcRec.T, self.Ji.nw_dataDiag[i], out = self.mtmp.T)
            np.dot(self.DcW, self.mtmp, out=self.Ji.nw_dataBlock[i])
        JiHatC = sp.bsr_matrix((self.Ji.nw_dataBlock,
                                self.Ji.indices, self.Ji.indptr),
                               shape=(self.dim, self.dim))
        for i in range(self.Jq.nnz):
            # Retrieve recovery matrix (different for each variable)
            WiDcRec = self.WiDcRecList[self.Jq.indices[i]]
            # Use element-wise multiplication for diagonal matrix
            np.multiply(WiDcRec.T, self.Jq.nw_dataDiag[i], out = self.mtmp.T)
            np.dot(self.DcWD, self.mtmp, out=self.Jq.nw_dataBlock[i])
        JqHatC = sp.bsr_matrix((self.Jq.nw_dataBlock,
                                self.Jq.indices, self.Jq.indptr),
                               shape=(self.dim, self.dim))
        # Jac: system matrix. In the future we can allocate memory
        # just once and operate directly on blocks.
        self.Jac = self.MHatC + JiHatC + JqHatC
        if glVar.verbose:
            print('Jacobian density= ',
                  100. * self.Jac.nnz / self.dim**2,'%')
        return (self.iVec, self.Jac)


    def recover(self, xVec, updateMHatC = False):
        """
        Recover samples from compressed coefficients

        Returns signal samples

        xVec: input vector of compressed nodal voltages. 

        if updateMHatC is True, update self.MHatC according to reconstruction

        """
        # Reconstruct signals
        xArray = xVec.reshape((self.ckt.nD_dimension, self.ncc)).T

        # At this point we have to decide what kind of
        # recovery matrix to use. Initially could use fixed PI, but at
        # some poing we have to switch. For now use a flag to indicate
        # to update the compressive matrices or not. Alternative Idea:
        # use an internal flag controlled by the nonlinear solver to
        # set the policy: 1. general PI; 2. adaptive custom PI;
        # 3. fixed custom PI.
        if self.adaptiveDcRec:
            xTilde = spams.lassoWeighted(xArray, D=self.Dc, **self.lassoParam)
            # Calculate samples of nodal voltages
            xHat = self.Wi.dot(xTilde).toarray()
            if glVar.verbose:
                print('x_Tilde density=', 100. * xTilde.nnz/self.dim, '%')
            if updateMHatC:
                D1 = np.empty((self.ncc, self.nsamples))
                # Update decompression matrices
                # Use nonzero structure from each column of xTilde 
                indptr = xTilde.indptr
                colidx = xTilde.indices
                for j in range(xTilde.shape[1]):
                    D1.fill(0.)
                    nzidx = colidx[indptr[j]:indptr[j+1]]
                    # Check if any elements are considered at all
                    if len(nzidx) > 0:
                        D1[:, nzidx] = self.Dc[:, nzidx]
                        # Try pinv and pinv2. Here we have to point the
                        # elements of the list to new matrices since
                        # originally we only allocate one
                        self.WiDcRecList[j] = self.Wi.dot(linalg.pinv2(D1))
                    else:
                        # Use default 
                        self.WiDcRecList[j] = self.WiDcRec
                self.MHatC = self.MHat.dot(sp.block_diag(
                    self.WiDcRecList)).tobsr(blocksize=(self.ncc, self.ncc))
        else:
            xHatv = sp.block_diag(self.WiDcRecList, format='bsr').dot(xVec)
            # Calculate samples of nodal voltages
            xHat = xHatv.reshape((self.ckt.nD_dimension, self.nsamples)).T
                    
        return xHat

        
