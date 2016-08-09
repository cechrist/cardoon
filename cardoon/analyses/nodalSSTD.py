"""
:mod:`nodaSSTD` -- Steady-state Nodal Analysis Time Domain
----------------------------------------------------------

.. module:: nodalSSTD
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

This module contains specific classes/functions for steady state nodal
analysis in time domain

Adds attributes to elements starting with "nsstd_"

"""

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg 
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

    M: sparse big Jacobian matrix (with dictionary list: .nsstd_dl)
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
            posContrib[:, posCount] = [i1*ncols+j1, M.nsstd_dl[i][j]]
            posCount += 1
        for j1, j in negCols:
            negContrib[:, negCount] = [i1*ncols+j1, M.nsstd_dl[i][j]]
            negCount += 1
    for i1, i in negRows:
        for j1, j in negCols:
            posContrib[:, posCount] = [i1*ncols+j1, M.nsstd_dl[i][j]]
            posCount += 1
        for j1, j in posCols:
            negContrib[:, negCount] = [i1*ncols+j1, M.nsstd_dl[i][j]]
            negCount += 1
    # Return arrays to be saved
    return posContrib, negContrib

def set_Jac(dataVec, posContrib, negContrib, Jac, j):
    """
    Sets big Jacobian elements from Jac into dataVec

    posContrib is one of the .nsstd_posJ? attributes in element 
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
    # Return dense array to be consistent with the other D matrices
    return sp.dia_matrix((data,offsets),
                         shape=(nsamples, nsamples)).todense().getA()

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

class _NLFunctionIter(_NLFunction):
    """
    Nonlinear function interface class for large sparse matrices
    Also works for non-square systems

    This is an abstract class to be used as a base class by other
    nodal-based classes such as DCNodal. Only methods that differ from
    nodal._NLFunction are defined here.

    """

    def __init__(self):
        # Override convergence helpers (generic gmin does not work with this)
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_source,
                                    self.solve_homotopy_source2, 
                                    None]
        # Used in source stepping helpers
        self._sstep = 0.5

    def factor_and_solve(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc

        So it returns the opposite of the Newton correction, but that
        is what fsolve() is expecting (-deltax)
        """
        # import pdb; pdb.set_trace()
#        result = scipy.sparse.linalg.lsmr(Jac, errFunc,
#                                          atol=1e-6, btol=1e-6)#, iter_lim=30)
#        print('code:',result[1],'iter:',result[2])
        result = scipy.sparse.linalg.gmres(Jac, errFunc, maxiter=20)
        print('code:',result[1])
        self.deltaxVec[:] = result[0]

        return self.deltaxVec


#--------------------------------------------------------------------------
    
class SSTDNodal(nd._NLFunctionSP):
    """
    Builds steady-state system of equations in time domain

    This class only sets up equations. Equations are solved elsewhere.

    Matrices and vectors (G, C, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit()) and number of samples. 
    """

    def __init__(self, ckt, nsamples, T, deriv, sstep):
        """
        Arguments:

        ckt: circuit instance
        nsamples: number of samples in period
        T: period
        deriv: derivative type, string
        sstep: source stepping factor

        """
        # Init base class
        super(SSTDNodal, self).__init__()
        # Override convergence helpers (generic gmin does not work with this)
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_source,
                                    self.solve_homotopy_source2, 
                                    None]
        # Source stepping factor
        self._sstep = sstep

        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref
        # Save circuit
        self.ckt = ckt
        self.nsamples = nsamples
        self.dim = self.ckt.nD_dimension*self.nsamples
        # T is the period
        self.T = T
        # Allocate matrices/vectors
        # G and C (done in refresh())

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
        self.timeVec = np.linspace(0., T, self.nsamples, endpoint=False)

        # TODO: include time delays
        #self.tdVecList = []

        # Create derivative matrix
        if deriv=='d2':
            D = deriv_2(self.nsamples, self.T)
        elif deriv=='d4':
            D = deriv_4(self.nsamples, self.T)
        elif deriv=='Fourier':
            D = deriv_fourier(self.nsamples, self.T)
        else:
            raise AnalysisError(
                'Invalid deriv value: {0}'.format(deriv))
            
        self.D = D
        self.refresh()


    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
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
                               dtype = float).tocsr()
        # Create matrix for linear part of circuit
        eyeNsamples = sp.eye(self.nsamples, self.nsamples, format='csr')
        GHat = sp.kron(self.G, np.eye(self.nsamples), 'csr')

        # CTriplet stores Jacobian matrix for a single time sample
        CTriplet = ([], [], [])
        # Add C matrix
        for elem in self.ckt.nD_elemList:
            for vcqs in elem.nD_linVCQS:
                nd.set_quad(CTriplet, *vcqs)
        self.C = sp.coo_matrix((CTriplet[0], CTriplet[1:]),
                               (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                               dtype = float).tocsr()

        CHat = sp.kron(self.C, self.D, 'csr')
        
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
                Yb[:,:] = Ytmp1

            YHat = sp.bsr_matrix((Yblocks, Y0.indices, Y0.indptr),
                                 shape=(self.dim, self.dim)).tocsr()

        # Linear response matrix
        if ycontrib:
            self.MHat = GHat + CHat + YHat
        else:
            self.MHat = GHat + CHat
            
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
        self.Ji.nsstd_dl = []
        for i in range(0, self.Ji.shape[0]):
            d=dict()
            baseIdx = self.Ji.indptr[i]
            endIdx = self.Ji.indptr[i+1]
            for j,k in enumerate(self.Ji.indices[baseIdx:endIdx]):
                d[k] = baseIdx + j
            self.Ji.nsstd_dl.append(d)
        # Reserve space for entries in time domain (diagonal matrices)
        self.Ji.nsstd_dataDiag = np.zeros((self.Ji.nnz, self.nsamples))
        # Use dense matrices for now
        self.Ji.nsstd_dataBlock = np.zeros((self.Ji.nnz,
                                            self.nsamples, self.nsamples))
        self.Ji.nsstd_diagIdx = np.diag_indices(self.nsamples)
        
        # Similar thing for Jq
        self.Jq.nsstd_dl = []
        for i in range(0, self.Jq.shape[0]):
            d=dict()
            baseIdx = self.Jq.indptr[i]
            endIdx = self.Jq.indptr[i+1]
            for j,k in enumerate(self.Jq.indices[baseIdx:endIdx]):
                d[k] = baseIdx + j
            self.Jq.nsstd_dl.append(d)
        # Reserve space for entries in time domain (diagonal matrices)
        self.Jq.nsstd_dataDiag = np.zeros((self.Jq.nnz, self.nsamples))
        # Use dense matrices just as a proof of concept. 
        self.Jq.nsstd_dataBlock = np.empty((self.Jq.nnz,
                                         self.nsamples, self.nsamples))

        # Create help arrays in elements
        for elem in self.ckt.nD_nlinElem:
            pos, neg = get_Jac_idx(self.Ji, elem,
                                   elem.nD_cpos, elem.nD_cneg, 
                                   elem.nD_vpos, elem.nD_vneg)
            elem.nsstd_posJi = pos
            elem.nsstd_negJi = neg
            pos, neg = get_Jac_idx(self.Jq, elem,
                                   elem.nD_qpos, elem.nD_qneg, 
                                   elem.nD_vpos, elem.nD_vneg)
            elem.nsstd_posJq = pos
            elem.nsstd_negJq = neg


        

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

        #import pdb; pdb.set_trace()
        return self.sVecA.flatten()
        

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

        xVec: input vector of nodal voltages (for all samples). 

        iVec: output vector of currents (for all samples)
        """
        # Convert to time domain
        xArray = xVec.reshape((self.ckt.nD_dimension, self.nsamples)).T
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
        # iVec = M x + i(x) + D q(x)
        self.iVecA += self.inlArray.T
        self.iVecA += self.D.dot(self.qArray).T
        return self.iVec

        
    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        # Convert to time domain
        xArray = xVec.reshape((self.ckt.nD_dimension, self.nsamples)).T
        #import pdb; pdb.set_trace()
        # Linear contribution
        self.iVec[:] = self.MHat.dot(xVec)
        # Nonlinear contribution
        # Erase arrays
        self.inlArray.fill(0.)
        self.qArray.fill(0.)
        self.Ji.nsstd_dataDiag.fill(0.)
        self.Jq.nsstd_dataDiag.fill(0.)
        self.Jq.nsstd_dataBlock[:] = self.D
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
                set_Jac(self.Ji.nsstd_dataDiag,
                        elem.nsstd_posJi, elem.nsstd_negJi, outJac, j)
                set_Jac(self.Jq.nsstd_dataDiag,
                        elem.nsstd_posJq, elem.nsstd_negJq, qJac, j)
        # iVec = M x + i(x) + D q(x)
        self.iVecA += self.inlArray.T
        self.iVecA += self.D.dot(self.qArray).T
        # Form system Jacobian
        #import pdb; pdb.set_trace()
        # Copy diagonal into dense matrix block
        self.Ji.nsstd_dataBlock[:,
                                self.Ji.nsstd_diagIdx[0],
                                self.Ji.nsstd_diagIdx[1]] = \
                                                self.Ji.nsstd_dataDiag
        JiHat = sp.bsr_matrix((self.Ji.nsstd_dataBlock,
                               self.Ji.indices, self.Ji.indptr),
                              shape=self.MHat.shape).tocsr()
        
        for i in range(self.Jq.nnz):
            self.Jq.nsstd_dataBlock[i] *= self.Jq.nsstd_dataDiag[i]
                   
        JqHat = sp.bsr_matrix((self.Jq.nsstd_dataBlock,
                               self.Jq.indices, self.Jq.indptr),
                              shape=self.MHat.shape).tocsr()
        # Jac: system matrix. In the future we can allocate memory
        # just once and operate directly on blocks.
        self.Jac = self.MHat + JiHat + JqHat
        #self.Jac.prune()
        if glVar.verbose:
            print('Jacobian density= ',
                  100. * self.Jac.nnz / self.dim**2,'%')
        return (self.iVec, self.Jac)


