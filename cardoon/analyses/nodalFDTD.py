"""
:mod:`nodal` -- Nodal Analysis using FDTD 
-----------------------------------------

.. module:: nodalFDTD
.. moduleauthor:: Carlos Christoffersen, Kris Fedick

This module contains specific classes/functions for FDTD nodal
analysis. 

"""

from __future__ import print_function
import os.path
from warnings import warn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg 
from fsolve import fsolve_Newton, NoConvergenceError
from integration import BEuler
from nodal import set_xin, set_i, restore_RCnumbers, delay_interp, _NLFunction
from cardoon.globalVars import glVar
import nodalSP as nd

#----------------------------------------------------------------------
# Auxiliary functions

def calculate_P_PT(numofnodes, numofsamples):
    """
    Calculates both P and the transpose in compressed format:

    P = [j0, j1, ... , j(n*m)] ; j: row index of P.T
    PT = [i0, i1, ... , i(n*m)] ; i: col index of P.T

    n = numofnodes
    m = numofsamples

    returns (P, PT)
    """
    P = np.zeros(numofnodes*numofsamples, dtype=int)
    PT = np.zeros(numofnodes*numofsamples, dtype=int)
    j = 0
    # sn: sample number
    for sn in range(numofsamples):
        i = sn
        # nn: node number
        for nn in range(numofnodes):
            P[i] = j
            PT[j] = i
            j += 1
            i += numofsamples
        
    return (P, PT)
            

def set_Jac(JacX, Jac, mpidx, mnidx, jacpidx, jacnidx):
    """
    Set current contributions of element Jacobian into JacX 

    JacX: JacI or JacQ (block Jacobians for each time sample)
    Jac: element Jacobian
    mpidx, mnidx, jacpidx, jacnidx: relative indexes (stored in element)

    JacX._base_idx: base index to keep track of current position in
    coo matrix (different for for JacI or JacQ). This attribute must
    exist and must be initially set to zero.

    """
#    import pdb; pdb.set_trace()
    JacX.data[mpidx + JacX._base_idx] = Jac.flat[jacpidx]
    JacX._base_idx += len(mpidx)
    JacX.data[mnidx + JacX._base_idx] = -Jac.flat[jacnidx]
    JacX._base_idx += len(mnidx)


#----------------------------------------------------------------------

class FDTDNodal(nd._NLFunctionSP):
    """
    Builds FDTD system of equations

    This class only sets up FDTD equations. Equations are solved
    elsewhere. 

    Matrices and vectors (G, C, JacI, JacQ, s, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit()) and number of samples. 
    """

    def __init__(self, ckt, nsamples, T):
        """
        Arguments:

        ckt: circuit instance
        nsamples: number of samples in period
        T: period

        """
        # Init base class
        super(FDTDNodal, self).__init__()
        # Override convergence helpers (generic gmin does not work with this)
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_source, 
                                    None]

        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref
        # Save circuit
        self.ckt = ckt
        self.nsamples = nsamples
        self.dim = self.ckt.nD_dimension*self.nsamples
        self.T = T
        # Allocate matrices/vectors
        # G and C (done in refresh())

        # iVec = G x + i(x)   total current
        self.iVec = np.zeros(self.dim)

        # Total charge: C x + q(x)
        self.qVec = np.zeros(self.dim)
        # Source vector 
        self.sVec = np.zeros(self.dim)
        # Delta_x includes samples for all nodes in circuit
        self.deltaxVec = np.zeros(self.ckt.nD_dimension * self.nsamples)

        # Create time vector
        self.timeVec = np.linspace(0., T, self.nsamples, endpoint=False)

        # TODO: include time delays
        # List of time-delay vectors: they are stored here to allow
        # sharing of one circuit amongst several TransientNodal
        # objects, to use in hierarchical simulation.
        #self.tdVecList = []
        self.xMatrix = self.sVec
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        # Names with 'B' refer to the 'big' matrices

        # GTriplet stores Jacobian matrix for a single time sample
        # Format is (data, row, col)
        GTriplet = ([], [], [])
        # Insert linear contributons: Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                nd.set_quad(GTriplet, *vccs)
#        # Frequency-defined elements not included for now
#        for elem in self.ckt.nD_freqDefinedElem:
#            set_Jac_triplet(JacTriplet, elem.nD_fpos, elem.nD_fneg, 
#                            elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        
        # CTriplet stores Jacobian matrix for a single time sample
        CTriplet = ([], [], [])
        # Add C matrix
        for elem in self.ckt.nD_elemList:
            for vcqs in elem.nD_linVCQS:
                nd.set_quad(CTriplet, *vcqs)
        
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


        # Allocate big Jacobian matrices and big 'D'
        (P, PT) = calculate_P_PT(self.ckt.nD_dimension, self.nsamples)
        # Create D (using centred derivatives)
        data = np.ones((4, self.nsamples))
        data *= 0.5 * self.nsamples / self.T
        data[[1,3],:] *= -1.
        offsets = np.array([-self.nsamples+1, -1, 1, self.nsamples-1])
        D = sp.dia_matrix((data,offsets), shape=(self.nsamples, self.nsamples))
        Id = sp.identity(self.ckt.nD_dimension, format='dia')
        BDcoo = sp.kron(Id, D).tocoo()
#        import pdb; pdb.set_trace()
        # Pre-multiply by P.T
        vtmp = np.array(BDcoo.row, copy=True)
        for j,val in enumerate(vtmp):
            BDcoo.row[j] = P[val]
        # Post-multiply by P
        vtmp = np.array(BDcoo.col, copy=True)
        for i,val in enumerate(vtmp):
            BDcoo.col[i] = P[val]
        # Convert to compressed row
        self.BD = BDcoo.tocsc()

	# Create overall system Jacobian matrices: 'big' JI and JQ
        # Step 1: Create data and index arrays for big matrices
        gnnz = len(GTriplet[0])
        nnzBG = self.nsamples * (gnnz + len(JiTriplet[0]))
        dataBG = np.empty(shape = nnzBG, dtype = float)
        ijBG = np.empty(shape = (2, nnzBG), dtype = int)

        qnnz = len(CTriplet[0])
        nnzBC = self.nsamples * (qnnz + len(JqTriplet[0]))
        dataBC = np.empty(shape = nnzBC, dtype = float)
        ijBC = np.empty(shape = (2, nnzBC), dtype = int)
	
        # Step 2: Add linear contributions
        dataG = np.array(GTriplet[0], dtype=float)
        ijG = np.array(GTriplet[1:], dtype=int)
        for j in range(self.nsamples):
            base = j * gnnz
            dataBG[base : base+gnnz] = dataG
            ijBG[:, base : base+gnnz] = ijG
            ijBG[:, base : base+gnnz] += self.ckt.nD_dimension * j
        dataQ = np.array(CTriplet[0], dtype=float)
        ijQ = np.array(CTriplet[1:], dtype=int)
        for j in range(self.nsamples):
            base = j * qnnz
            dataBC[base : base+qnnz] = dataQ
            ijBC[:, base : base+qnnz] = ijQ
            ijBC[:, base : base+qnnz] += self.ckt.nD_dimension * j

        # Step 3: Add only index values for nonlinear contributions
	# Update the conductance/current jacobian
        baseG = self.nsamples * gnnz
        ji_nnz = len(JiTriplet[0])
        ijJi = np.array(JiTriplet[1:], dtype=int)
        for j in range(self.nsamples):
            ijBG[:, baseG : baseG + ji_nnz] = ijJi
            ijBG[:, baseG : baseG + ji_nnz] += self.ckt.nD_dimension * j
            baseG += ji_nnz

        baseC = self.nsamples * qnnz
        jq_nnz = len(JqTriplet[0])
        ijJq = np.array(JqTriplet[1:], dtype=int)
        for j in range(self.nsamples):
            ijBC[:, baseC : baseC + jq_nnz] = ijJq
            ijBC[:, baseC : baseC + jq_nnz] += self.ckt.nD_dimension * j 
            baseC += jq_nnz
        
        # Step 4: create coo matrices
        # 'big' G matrix
        baseG = self.nsamples * gnnz
        baseC = self.nsamples * qnnz
        self.BGcoo = sp.coo_matrix((dataBG, ijBG),
                                  shape = (self.dim, self.dim), 
                                  dtype = float)
        # This is the starting index for nonlinear contributions
        self.BGcoo._base_idx = baseG
        self.baseG = baseG
        # 'big' C matrix
        self.BCcoo = sp.coo_matrix((dataBC, ijBC),
                                   (self.dim, self.dim), 
                                   dtype = float)
        # This is the starting index for nonlinear contributions
        self.BCcoo._base_idx = baseC
        self.baseC = baseC

        # Save G and C in csr format for fast multiplication
        Gcoo = sp.coo_matrix((GTriplet[0], (GTriplet[1], GTriplet[2])),
                             (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                             dtype = float)
                             
        # Convert to compressed-row (csr) format for efficient
        # matrix-vector multiplication
        self.G = Gcoo.tocsr()
        Ccoo = sp.coo_matrix((CTriplet[0], (CTriplet[1], CTriplet[2])),
                             (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                             dtype = float)
        # Convert to compressed-row (csr) format for efficient
        # matrix-vector multiplication
        self.C = Ccoo.tocsr()


    def get_source(self):
        """
        Get the source vector considering DC and TD source components
        """
        # Erase vector first. 
        self.sVec.fill(0.)
        for j, ctime in enumerate(self.timeVec):
            offset = j * self.ckt.nD_dimension
            for elem in self.ckt.nD_sourceDCElem:
                # first get the destination row/columns 
                outTerm = elem.nD_sourceOut
                current = elem.get_DCsource()
                # This may not need optimization because we usually do not
                # have too many independent sources
                # import pdb; pdb.set_trace()
                if outTerm[0] >= 0:
                    self.sVec[outTerm[0]+offset] -= current
                if outTerm[1] >= 0:
                    self.sVec[outTerm[1]+offset] += current
            for elem in self.ckt.nD_sourceTDElem:
                # first get the destination row/columns 
                outTerm = elem.nD_sourceOut
                current = elem.get_TDsource(ctime)
                # This may not need optimization because we usually do not
                # have too many independent sources
                # import pdb; pdb.set_trace()
                if outTerm[0] >= 0:
                    self.sVec[outTerm[0]+offset] -= current
                if outTerm[1] >= 0:
                    self.sVec[outTerm[1]+offset] += current
        return self.sVec
        
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
        for j in range(self.nsamples):
            base = self.ckt.nD_dimension * j
            itmp = self.iVec[base:base+self.ckt.nD_dimension]
            qtmp = self.qVec[base:base+self.ckt.nD_dimension]
            xtmp = xVec[base:base+self.ckt.nD_dimension]
            # Linear contribution
            itmp[:] = self.G * xtmp
            qtmp[:] = self.C * xtmp
            dcounter = 0
            # Nonlinear contribution
            for elem in self.ckt.nD_nlinElem:
                # first have to retrieve port voltages from xVec
                xin = np.zeros(elem.nD_nxin)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xtmp)
                # We can not yet treat delays
                assert elem.nDelays == 0
                outV = elem.eval(xin)
                # Update iVec, qVec
                set_i(itmp, elem.nD_cpos, elem.nD_cneg, outV)
                set_i(qtmp, elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
        # iVec = [G x + i(x)] + bigD [C x + q(x)] 
        self.iVec += self.BD * self.qVec
        return self.iVec
        
    def storedata(self, xVec):
        self.xMatrix = np.vstack((self.xMatrix,xVec))
        
    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        self.BGcoo._base_idx = self.baseG
        self.BCcoo._base_idx = self.baseC
        # print(max(abs(xVec)))
        self.storedata(xVec)

        for j in range(self.nsamples):

            base = self.ckt.nD_dimension * j

            itmp = self.iVec[base:base+self.ckt.nD_dimension]
            qtmp = self.qVec[base:base+self.ckt.nD_dimension]
            xtmp = xVec[base:base+self.ckt.nD_dimension]
            
            # Linear contribution
            itmp[:] = self.G * xtmp
            qtmp[:] = self.C * xtmp
            
            # Nonlinear contribution
            for elem in self.ckt.nD_nlinElem:
                # first have to retrieve port voltages from xVec
                xin = np.zeros(elem.nD_nxin)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xtmp)
                # We can not yet treat delays
                assert elem.nDelays == 0

                (outV, outJac) = elem.eval_and_deriv(xin)
                # Update iVec, qVec and block Jacobians now. 
                set_i(itmp, elem.nD_cpos, elem.nD_cneg, outV)
                set_i(qtmp, elem.nD_qpos, elem.nD_qneg, 
                      outV[len(elem.csOutPorts):])
                qJac = outJac[len(elem.csOutPorts):,:]
                set_Jac(self.BGcoo, outJac, *elem.nD_csidx)
                set_Jac(self.BCcoo, qJac, *elem.nD_qsidx)
          
        # iVec = [G x + i(x)] + bigD [C x + q(x)] 
        self.iVec += self.BD * self.qVec
        # print(self.iVec)
        BC = self.BCcoo.tocsc()
        Jac =  self.BGcoo.tocsc() + self.BD * BC
        # TODO: make sure result is in csc format
        #import pdb; pdb.set_trace()
        return (self.iVec, Jac)


    def get_guess(self):
        """
        Retrieve guesses from vPortGuess in each nonlinear device
        
        Returns a guess vector for DC operating point used for initial guess
	Format is nodal voltages for each sample.
        """
        x0 = np.zeros(self.dim)
        for elem in self.ckt.nD_nlinElem:
            try:
                # Only set positive side. This is not the only way
                # and may not work well in some cases but this is a
                # guess anyway
                for i,j in elem.nD_vpos:
                    x0[j] = elem.vPortGuess[i]
            except AttributeError:
                # if vPortGuess not given just leave things unchanged
                pass
        return x0

    def set_IC(self):
        """
        Set initial conditions (ICs)

        Retrieves ICs from DC operating point info in elements
        """
        # Nodal variables
        xVec = np.zeros(self.ckt.nD_dimension*self.nsamples)
        # Get initial nodal voltages
        for i,term in enumerate(self.ckt.nD_termList):
            xVec[i*self.nsamples] = term.nD_vOP
        # Calculate total charges
	self.refresh()
        # Initialize time delay structures
	# TODO:  Add time delay initialization

    def save_OP(self, xVec):
        """
        Save nodal voltages in terminals and set OP in elements

        The following information is saved:

          * The nodal voltage vector for the circuit (self.xop)

          * The nodal voltage in each terminal (term.nD_vOP)

          * The port voltages in each nonlinear device (elem.nD_xOP)

          * The operating point (OP) information in nonlinear devices

        """
        # Set nodal voltage of reference to zero
        self.ckt.nD_ref.nD_vOP = 0.
        # Save nodal vector
        self.xop = xVec
        for v,term in zip(xVec, self.ckt.nD_termList):
            term.nD_vOP = v
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            # Set OP in element (discard return value)
            elem.nD_xOP = xin
            elem.get_OP(xin)
