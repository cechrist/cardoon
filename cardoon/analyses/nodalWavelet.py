"""
:mod:`nodalWavelet` -- Steady-state Nodal Analysis Using Wavelets
-----------------------------------------------------------------

.. module:: nodalSS
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

This module contains specific classes/functions for steady state nodal
analysis using wavelet coefficients.

Adds attributes to elements starting with "nw_"

"""

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import pywt
from nodal import set_xin, set_i, _NLFunction
import nodalSP as nd
from analysis import AnalysisError
from cardoon.globalVars import glVar


#----------------------------------------------------------------------
# Jacobian handling

def get_Jac_idx(M, elem, posRows, negRows, posCols, negCols):
    """
    Set index matrix in elem for quick Jacobian filling

    M: sparse big Jacobian matrix (with dictionary list: .nw_dl)
    elem: element instance
    """
    # Lists with element Jacobian indices
    posContrib = np.empty((2, len(posRows)*len(posCols)
                           + len(negRows)*len(negCols)), dtype=int)
    negContrib = np.empty((2, len(posRows)*len(negCols)
                           + len(negRows)*len(posCols)), dtype=int)
    ncols = elem.nD_nxin
    posCount = 0
    negCount = 0
    # i,j large Jacobian indices
    # i1, j1 element Jacobian indices
    for i1, i in posRows:
        for j1, j in posCols:
            posContrib[:, posCount] = [i1*ncols+j1, M.nw_dl[i][j]]
            posCount += 1
        for j1, j in negCols:
            negContrib[:, negCount] = [i1*ncols+j1, M.nw_dl[i][j]]
            negCount += 1
    for i1, i in negRows:
        for j1, j in negCols:
            posContrib[:, posCount] = [i1*ncols+j1, M.nw_dl[i][j]]
            posCount += 1
        for j1, j in posCols:
            negContrib[:, negCount] = [i1*ncols+j1, M.nw_dl[i][j]]
            negCount += 1
    # Return arrays to be saved
    return posContrib, negContrib

def set_Jac(dataVec, posContrib, negContrib, Jac, j):
    """
    Sets big Jacobian elements from Jac into dataVec

    posContrib is one of the .nw_posJ? attributes in element 
    (same for negContrib) 

    Jac is element Jacobian
    
    j is sample number
    """
    # I have no choice but to go element by element as numpy has
    # trouble assigning 2 values to the same destination.
    for k,l in posContrib.T:
        dataVec[l, j] += Jac.flat[k]
    for k,l in negContrib.T:
        dataVec[l, j] -= Jac.flat[k]

    # This stuff below does not work:
    #    dataVec[posContrib[1], j] += Jac.flat[posContrib[0]]
    #    dataVec[negContrib[1], j] -= Jac.flat[negContrib[0]]

#----------------------------------------------------------------------
# Derivatives

def deriv_1(nsamples, T):
    """
    Returns D matrix using first order derivatives (BE formula)
    
    nsamples: number of samples
    T: period
    """
    data = np.ones((3, nsamples))
    data *= nsamples / T
    data[[0,2],:] *= -1.
    offsets = np.array([-1, 0, nsamples-1])
    return sp.dia_matrix((data,offsets), shape=(nsamples, nsamples))

def deriv_2(nsamples, T):
    """
    Returns D matrix using second order derivatives (2-point formula)
    
    nsamples: number of samples
    T: period
    """
    data = np.ones((4, nsamples))
    data *= 0.5 * nsamples / T
    data[[1,3],:] *= -1.
    offsets = np.array([-nsamples+1, -1, 1, nsamples-1])
    return sp.dia_matrix((data,offsets), shape=(nsamples, nsamples))

def deriv_alpha(nsamples, T, alpha):
    """
    Returns D matrix using second order derivatives as follows:
    
    alpha * FE + (1-alpha) * BE , 

    where FE and BE stand for forward and backward Euler.

    alpha = 0 is BE formula
    alpha = 0.5 is 2-point formula

    nsamples: number of samples
    T: period
    """
    data = np.ones((5, nsamples))
    data *= nsamples / T
    data[2,:] *= (1. - 2. * alpha)
    data[[1,4],:] *= (-1. + alpha)
    data[[0,3],:] *= alpha
    offsets = np.array([-nsamples+1, -1, 0, 1, nsamples-1])
    return sp.dia_matrix((data,offsets), shape=(nsamples, nsamples))

def deriv_4(nsamples, T):
    """
    Returns D matrix using 4th order derivatives (5-point formula)
    
    nsamples: number of samples
    T: period
    """
    # Create a dense matrix first (easier)
    row = np.zeros(nsamples)
    row[1:3] = [8.,-1.]
    row[-2:] = [1.,-8.]
    row *= nsamples / T / 12.
    Dd = np.zeros((nsamples, nsamples))
    for i in range(0, nsamples):
        Dd[i] = np.roll(row, i)
    # create diagonal sparse matrix
    return Dd


def deriv_fourier(nsamples, T):
    """
    Returns D matrix using Fourier derivatives (*jw in freq. domain)
    
    (this matrix is dense)
    nsamples: number of samples
    T: period
    """
    # Assume number of samples is even
    nfreq = nsamples // 2 + 2
    w0 = 2. * np.pi / T
    Y = 1j * w0 * np.arange(0,nfreq)
    # Inverse-transform Y: for now use numpy. Later switch to pyfftw
    y = np.fft.irfft(Y, nsamples)
    # Reverse order of elements in y
    yi = y[::-1]
    # Create convolution matrix
    D = np.zeros((nsamples, nsamples))
    for i in range(0, nsamples):
        D[i] = np.roll(yi, i+1)
    return D

#----------------------------------------------------------------------
# Wavelet functions

def wavelet_levels(nsamples, wtype='db4'):
    """
    Wrapper for the pywt function
    """
    w = pywt.Wavelet(wtype)
    maxlev = pywt.dwt_max_level(nsamples, w)
    return maxlev

# Multilevel transforms: require nsamples = 2^n
def wavelet_f(nsamples, wtype='db4'):
    """
    Returns wavelet transform matrix in dense format

    nsamples must be a power of two
    """
    w = pywt.Wavelet(wtype)
    levels = pywt.dwt_max_level(nsamples, w)
    row = np.zeros(nsamples)
    W = np.empty((nsamples, nsamples))
    for i in range(nsamples):
        row[i] = 1.
        coeffs = pywt.wavedec(row, wtype, level=levels, mode='per')
        row[i] = 0.
        W[:,i] = np.concatenate(coeffs)
    # Return sparse matrix as most entries should be zero
    return W

# Soveiko-style transforms
def wavelet_fs(nsamples, wtype='db4'):
    """
    Returns wavelet transform matrix in dense format
    """
    row = np.zeros(nsamples)
    W = np.empty((nsamples, nsamples))
    for i in range(nsamples):
        row[i] = 1.
        (cA, cD) = pywt.dwt(row, wtype, 'per')
        row[i] = 0.
        # Arrange coefficients for banded-diagonal structure. This
        # will change if we use multi-resolution.
        W[::2,i] = cA
        W[1::2,i] = cD
    # Return sparse matrix as most entries should be zero
    return W


#----------------------------------------------------------------------

class WaveletNodal(nd._NLFunctionSP):
    """
    Builds steady-state system of equations using Wavelets

    This class only sets up equations. Equations are solved elsewhere.

    Matrices and vectors (G, C, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit()) and number of samples. 
    """

    def __init__(self, ckt, nsamples, T, wavelet, multilevel, deriv, alpha,
                 sstep, matformat='csr'):
        """
        Arguments:

        ckt: circuit instance
        nsamples: number of samples in period
        T: period
        wavelet: wavelet type, string
        multilevel: it True, use multilevel Wavelet transform
        deriv: derivative type, string
        sstep: source stepping factor
        matformat: Jacobian matrix format ('bsr' for compressed analysis)

        """
        # Init base class
        super(WaveletNodal, self).__init__()
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
        # T is the period
        self.T = T
        # Sparse matrix format for Jacobian
        self.matformat = matformat
        # Save input parameters for refresh()
        self._deriv = deriv
        self._alpha = alpha
        self._wavelet = wavelet
        self._multilevel = multilevel
        # Allocate matrices/vectors (done in refresh())
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        self.dim = self.ckt.nD_dimension*self.nsamples

        # Vectors are allocated as matrices, each row is the set of
        # samples for one variable. Example: iVec[2] is the set of
        # samples for the second current in iVec.

        # iVec: total current
        self.iVecA = np.zeros((self.ckt.nD_dimension, self.nsamples))
        self.iVec = self.iVecA.ravel()
        
        # nonlinear current only
        self.inlArray = np.zeros((self.nsamples, self.ckt.nD_dimension))
        # Nonlinear charge only
        self.qArray = np.zeros((self.nsamples, self.ckt.nD_dimension))
        # Source vector 
        self.sVecA = np.zeros((self.ckt.nD_dimension, self.nsamples))

        # Delta_x is required as workspace by the nonlinear solver
        self.deltaxVec = np.zeros(self.ckt.nD_dimension * self.nsamples)

        # Create time vector
        self.timeVec = np.linspace(0., self.T, self.nsamples, endpoint=False)

        # TODO: include time delays
        #self.tdVecList = []

        deriv = self._deriv
        # Create derivative matrix
        if deriv=='alpha':
            D = deriv_alpha(self.nsamples, self.T, self._alpha)
        elif deriv=='d2':
            D = deriv_2(self.nsamples, self.T)
        elif deriv=='d4':
            D = sp.dia_matrix(deriv_4(self.nsamples, self.T))
        elif deriv=='Fourier':
            D = deriv_fourier(self.nsamples, self.T)
        else:
            raise AnalysisError(
                'Invalid deriv value: {0}'.format(deriv))
        #import pdb; pdb.set_trace()
            
        # Wavelet matrices. It may be possible to avoid its creation
        # (using the fast wavelet transform) but for the moment is it
        # not clear that we would save much with that.
        wavelet = self._wavelet
        if wavelet == 'none':
            # Use identity matrices for wavelet transform: not the
            # most efficient way for FDTD but useful for testing
            self.W = sp.eye(self.nsamples, self.nsamples, format='csr')
            self.Wi = self.W
        else:
            if self._multilevel:
                # Use multilevel transforms
                self.W = sp.csr_matrix(wavelet_f(self.nsamples, wavelet))
                # Inverse transform is just the transpose
                self.Wi = self.W.T
            else:
                # Use Soveiko-style matrices (single level, banded diagonal)
                self.W = sp.csr_matrix(wavelet_fs(self.nsamples, wavelet))
                # Inverse transform is just the transpose
                self.Wi = self.W.T
        self.WD = self.W.dot(D)

        # The following used in get_Jac_i only
        self.WiT = self.Wi.todense().getA().T
        self.mtmp = np.empty(self.WiT.shape)        

        # Big matrices named *Hat

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
                               dtype = float)
        # Create matrix for linear part of circuit
        eyeNsamples = sp.eye(self.nsamples, self.nsamples, format='csr')
        GHat = sp.kron(self.G, np.eye(self.nsamples))

        # CTriplet stores Jacobian matrix for a single time sample
        CTriplet = ([], [], [])
        # Add C matrix
        for elem in self.ckt.nD_elemList:
            for vcqs in elem.nD_linVCQS:
                nd.set_quad(CTriplet, *vcqs)
        self.C = sp.coo_matrix((CTriplet[0], CTriplet[1:]),
                               (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                               dtype = float)
        # If we use Fourier derivatives, WD is dense and needs special
        # treatment
        WDWi = self.W.dot(self.WD.T).T
        CHat = sp.kron(self.C, WDWi)
        
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
            # Create convolution matrices
            # Inverse-transform Y: for now use numpy. Later switch to pyfftw
            Yblocks = np.empty((Y0.nnz, self.nsamples, self.nsamples))
            Ytmp1 = np.empty((self.nsamples, self.nsamples))
            Ytmp2 = np.empty((self.nsamples, self.nsamples))
            for Yv, Yb in zip(ydata, Yblocks):
                yv = np.fft.irfft(Yv, self.nsamples)
                # Reverse order of elements in yv
                yi = yv[::-1]
                # Create convolution matrix
                Ytmp1.fill(0.)
                for i in range(0, self.nsamples):
                    Ytmp1[i] = np.roll(yi, i+1)
                # Add Wavelet Transform
                Ytmp2 = self.W.dot(Ytmp1.T).T
                Yb[:,:] = self.W.dot(Ytmp2)

            YHat = sp.bsr_matrix((Yblocks, Y0.indices, Y0.indptr),
                                 shape=(self.dim, self.dim))

        # Linear response matrix
        if ycontrib:
            MHat = GHat + CHat + YHat
        else:
            MHat = GHat + CHat

        if self.matformat == 'bsr':
            self.MHat = sp.bsr_matrix(MHat, blocksize=(self.nsamples,
                                                       self.nsamples))
        elif self.matformat == 'csr':
            self.MHat = MHat
        else:
            raise AnalysisError(
                'Unsupported format: {0}'.format(self.matformat))
            
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
        # Use dense matrices just as a proof of concept. 
        self.Ji.nw_dataBlock = np.empty((self.Ji.nnz,
                                         self.nsamples, self.nsamples))
        self.JiHat = sp.bsr_matrix((self.Ji.nw_dataBlock,
                                    self.Ji.indices, self.Ji.indptr),
                                   shape=self.MHat.shape)
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
        # Use dense matrices just as a proof of concept. 
        self.Jq.nw_dataBlock = np.empty((self.Jq.nnz,
                                         self.nsamples, self.nsamples))
        self.JqHat = sp.bsr_matrix((self.Jq.nw_dataBlock,
                                    self.Jq.indices, self.Jq.indptr),
                                   shape=self.MHat.shape)

        # Create help arrays in elements
        for elem in self.ckt.nD_nlinElem:
            pos, neg = get_Jac_idx(self.Ji, elem,
                                   elem.nD_cpos, elem.nD_cneg, 
                                   elem.nD_vpos, elem.nD_vneg)
            elem.nw_posJi = pos
            elem.nw_negJi = neg
            pos, neg = get_Jac_idx(self.Jq, elem,
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
        s = self.W.dot(self.sVecA.T)
        return s.T.flatten()
        

    def get_rhs(self, xVec):
        """
        Returns system rhs vector: s' 

        s' = source current - currents generated by other circuit components

        """
        return self.get_source() - self.get_i(xVec)


    def get_i(self, xVec):
        """
        Calculate total current

        returns currents generated by non-sources

        xVec: input vector of nodal voltages (for all samples). 

        iVec: output vector of currents (for all samples)
        """
        # Convert to time domain
        xArray = self.Wi.dot(xVec.reshape((self.ckt.nD_dimension,
                                           self.nsamples)).T)
        # Linear contribution
        self.iVec[:] = self.MHat.dot(xVec)
        # Nonlinear contribution
        # Erase arrays
        self.inlArray.fill(0.)
        self.qArray.fill(0.)
        for elem in self.ckt.nD_nlinElem:
            for j in range(self.nsamples):
                # first have to retrieve port voltages from xVec
                xin = np.zeros(elem.nD_nxin)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xArray[j,:])
                # We can not yet treat delays
                assert elem.nDelays == 0
                outV = elem.eval(xin)
                # Update iVec, qVec
                set_i(self.inlArray[j], elem.nD_cpos, elem.nD_cneg, outV)
                set_i(self.qArray[j], elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
        # iVec = M x + W i(x) + W D q(x)
        self.iVecA += self.W.dot(self.inlArray).T
        self.iVecA += self.WD.dot(self.qArray).T
        return self.iVec

        
    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        if glVar.verbose:
            print('x_hat density=',
                  100. * np.sum(np.abs(xVec)>1e-2)/len(xVec), '%')
        # Convert to time domain
        xArray = self.Wi.dot(xVec.reshape((self.ckt.nD_dimension,
                                           self.nsamples)).T)
        #import pdb; pdb.set_trace()
        # Linear contribution
        self.iVec[:] = self.MHat.dot(xVec)
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
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xArray[j,:])
                # We can not yet treat delays
                assert elem.nDelays == 0
                (outV, outJac) = elem.eval_and_deriv(xin)
                # Update iVec, qVec
                set_i(self.inlArray[j], elem.nD_cpos, elem.nD_cneg, outV)
                set_i(self.qArray[j], elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
                qJac = outJac[len(elem.csOutPorts):,:]
                # Set Jacobian entries for all samples at once
                set_Jac(self.Ji.nw_dataDiag,
                        elem.nw_posJi, elem.nw_negJi, outJac, j)
                set_Jac(self.Jq.nw_dataDiag,
                        elem.nw_posJq, elem.nw_negJq, qJac, j)
        # iVec = M x + W i(x) + W D q(x)
        self.iVecA += self.W.dot(self.inlArray).T
        self.iVecA += self.WD.dot(self.qArray).T
        # Form system Jacobian
        #import pdb; pdb.set_trace()
        offsets = np.zeros(1, dtype=int)
        for i in range(self.Ji.nnz):
            #  # Experimental: try sparse matrix block calculation
            #  d = sp.dia_matrix((self.Ji.nw_dataDiag[i],offsets),
            #                    shape=(self.nsamples, self.nsamples))
            #  self.Ji.nw_dataBlock[i] = self.W.dot(d.dot(self.Wi)).todense()

            # Use element-wise multiplication for diagonal matrix
            np.multiply(self.WiT, self.Ji.nw_dataDiag[i], out = self.mtmp)
            self.Ji.nw_dataBlock[i] = self.W.dot(self.mtmp.T)
            # Optional: threshold matrix to eliminate near-zero entries
            # pywt.thresholding.hard(self.Ji.nw_dataBlock[i], 1e-13)
        for i in range(self.Jq.nnz):
            # Use element-wise multiplication for diagonal matrix
            np.multiply(self.WiT, self.Jq.nw_dataDiag[i], out = self.mtmp)
            self.Jq.nw_dataBlock[i] = self.WD.dot(self.mtmp.T)
        # Jac: system matrix. In the future we can allocate memory
        # just once and operate directly on blocks.
        self.Jac = self.MHat + self.JiHat + self.JqHat
        if glVar.verbose:
            print('Jacobian density= ',
                  100. * self.Jac.nnz / self.dim**2,'%')
        return (self.iVec, self.Jac)


