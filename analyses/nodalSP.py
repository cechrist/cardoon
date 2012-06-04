"""
:mod:`nodalSP` -- Nodal Analysis using Scipy sparse matrices
------------------------------------------------------------

.. module:: nodalSP
.. moduleauthor:: Carlos Christoffersen

This module contains basic classes/functions for nodal analysis. These
are part of the standard netlist analyses but they can also be used
independently.

This implementation uses Scipy sparse matrices. By default SuperLU is
used to factor the matrix. UMFPack could be used, but SuperLU seems to
be faster. The main matrix is built in triplet format (Scipy
coo_matrix format).

Many of the functions are directly imported from the ``nodal`` module
to avoid redundancy. Some of these could be optimized for better
performance.

"""

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg 
from fsolve import fsolve_Newton, NoConvergenceError
from integration import BEuler
import nodal
from nodal import set_xin, set_i, restore_RCnumbers, _NLFunction

# ****************** Stand-alone functions to be optimized ****************

def triplet_append(G, val, row, col):
    """
    Add val at [row, col] position into G matrix (triplet format)
    """
    G[0].append(val)
    G[1].append(row)
    G[2].append(col)


def set_quad(G, row1, col1, row2, col2, g):
    """
    Set transconductance/transcapacitance quad

    G: target matrix in (dataVec, rowVec, colVec) format (triplet format)
    """
    #import pdb; pdb.set_trace()
    if col1 >= 0:
        if row1 >= 0:
            triplet_append(G, g, row1, col1)
        if row2 >= 0:
            triplet_append(G, -g, row2, col1)
    if col2 >= 0:
        if row1 >= 0:
            triplet_append(G, -g, row1, col2)
        if row2 >= 0:
            triplet_append(G, g, row2, col2)


def set_Jac_triplet(M, posRows, negRows, posCols, negCols, Jac):
    """
    Set elements in M
    
    M: target matrix in (dataVec, rowVec, colVec) format (triplet format)

    The looping is not optimal but has the advantage that all positive
    values are inserted first. This simplify index calculation in
    create_additional_indexes()
    """
    for i1, i in posRows:
        for j1, j in posCols:
            triplet_append(M, Jac[i1,j1], i, j)
    for i1, i in negRows:
        for j1, j in negCols:
            triplet_append(M, Jac[i1,j1], i, j)
    for i1, i in posRows:
        for j1, j in negCols:
            triplet_append(M, -Jac[i1,j1], i, j)
    for i1, i in negRows:
        for j1, j in posCols:
            triplet_append(M, -Jac[i1,j1], i, j)

# ********************  Regular functions *******************************

def make_nodal_circuit(ckt, reference='gnd'):
    """
    Add attributes to Circuit/Elements/Terminals for nodal analysis

    Similar to nodal.make_nodal_circuit but in addition calls
    create_additional_indexes()
    """
    nodal.make_nodal_circuit(ckt, reference)
    # Generate sparse-matrix index mappings
    for elem in ckt.nD_nlinElem:
        create_additional_indexes(elem)
    

def process_nodal_element(elem):
    """
    Process element for nodal analysis
    """
    nodal.process_nodal_element(elem)
    if elem.isNonlinear:
        create_additional_indexes(elem)


def create_additional_indexes(elem):
    """
    Creates pre-calculated indexes for nonlinear devices

    The following indexes are calculated:

    nD_csidx = (cpidx, cnidx, jacpidx, jacnidx)
    nD_qsidx = (qpidx, qnidx, jacqpidx, jacqnidx)

    cpidx, cnidx, qpidx and qnidx are relative to 0 in this
    element. The final position in the coo matrix is unknown at this
    time.

    Important note: these indexes depend on the order used to insert
    elements in set_Jac_triplet()
    """
    lvpos = len(elem.nD_vpos)
    lvneg = len(elem.nD_vneg)
    lcpos = len(elem.nD_cpos)
    lcneg = len(elem.nD_cneg)
    lqpos = len(elem.nD_qpos)
    lqneg = len(elem.nD_qneg)
    # Create vectors
    cpidx = np.arange(0, lvpos * lcpos + lvneg * lcneg, dtype = int)
    cnidx = np.arange(0, lvpos * lcneg + lvneg * lcpos, dtype = int) 
    qpidx = np.arange(0, lvpos * lqpos + lvneg * lqneg, dtype = int)
    qnidx = np.arange(0, lvpos * lqneg + lvneg * lqpos, dtype = int) 

    ncols = len(elem.controlPorts)
    # Current source Jac indexes
    jaci = [i1 * ncols + j1 
            for i1, i in elem.nD_cpos for j1, j in elem.nD_vpos]
    jaci += [i1 * ncols + j1 
            for i1, i in elem.nD_cneg for j1, j in elem.nD_vneg]
    jacpidx = np.array(jaci, dtype=int)
    jaci = [i1 * ncols + j1 
            for i1, i in elem.nD_cpos for j1, j in elem.nD_vneg]
    jaci += [i1 * ncols + j1 
            for i1, i in elem.nD_cneg for j1, j in elem.nD_vpos]
    jacnidx = np.array(jaci, dtype=int)
    # Charge source Jac indexes
    jaci = [i1 * ncols + j1 
            for i1, i in elem.nD_qpos for j1, j in elem.nD_vpos]
    jaci += [i1 * ncols + j1 
            for i1, i in elem.nD_qneg for j1, j in elem.nD_vneg]
    jacqpidx = np.array(jaci, dtype=int)
    jaci = [i1 * ncols + j1 
            for i1, i in elem.nD_qpos for j1, j in elem.nD_vneg]
    jaci += [i1 * ncols + j1 
            for i1, i in elem.nD_qneg for j1, j in elem.nD_vpos]
    jacqnidx = np.array(jaci, dtype=int)

    elem.nD_csidx = (cpidx, cnidx, jacpidx, jacnidx)
    elem.nD_qsidx = (qpidx, qnidx, jacqpidx, jacqnidx)


def _set_Jac(self, M, Jac, mpidx, mnidx, jacpidx, jacnidx):
    """
    Set current contributions of Jac into M (in coo format)
    """
    # Apparently fancy indexing is a little faster than using
    # ``nimpy.put()`` and ``numpy.take()`` to fill the matrix.

#   import pdb; pdb.set_trace()
#    np.put(M.data, mpidx + self._mbase, np.take(Jac, jacpidx))
    M.data[mpidx + self._mbase] = Jac.flat[jacpidx]
    self._mbase += len(mpidx)
#    np.put(M.data, mnidx + self._mbase, -np.take(Jac, jacnidx))
    M.data[mnidx + self._mbase] = -Jac.flat[jacnidx]
    self._mbase += len(mnidx)


#-------------------------------------------------------------------------
# ****************************** Classes *********************************
#-------------------------------------------------------------------------

class _NLFunctionSP(_NLFunction):
    """
    Nonlinear function interface class using sparse matrices

    This is an abstract class to be used as a base class by other
    nodal-based classes such as DCNodal. Only methods that differ from
    nodal._NLFunction are defined here.

    """

    def __init__(self):
        # List here the functions that can be used to solve equations
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_gmin2, 
                                    self.solve_homotopy_source, 
                                    None]
        sp.linalg.use_solver(useUmfpack = False, assumeSortedIndices = True)
        self._factorized = None

    def _get_deltax(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc
        """
        # import pdb; pdb.set_trace()
        M = Jac.tocsc()
        #print(M.todense())
        #print('------------------------------------------------')
        # Decompose matrix
        try:
            self._factorized = sp.linalg.factorized(M)
        except RuntimeError as re:
            # In the future we could try a least-squares solution here
            print('\nProblem factoring matrix: ', re)
            exit()
        # Solve linear system
        self.deltaxVec[:] = self._factorized(errFunc)
        return self.deltaxVec

    def get_chord_deltax(self, sV, iVec=None):
        """
        Get deltax for sV, iVec using existing factored Jacobian

        Useful for the first iteration of transient analysis. If iVec
        not given the stored value is used.
        """
        if iVec == None:
            iVec = self.iVec
        return self._factorized(iVec - sV)

    def solve_homotopy_gmin2(self, x0, sV):
        """Newton's method with gmin stepping (nonlinear ports)"""
        # Adds gmin in parallel with nonlinear element (external) ports
        # Create Gmin matrix (structure is fixed)
        G1triplet = ([], [], [])
        for elem in self.ckt.nD_nlinElem:
            Gadd = np.eye(len(elem.nD_extPorts))
            # import pdb; pdb.set_trace()
            set_Jac_triplet(G1triplet, elem.nD_epos, elem.nD_eneg, 
                            elem.nD_epos, elem.nD_eneg, Gadd)
        G1coo = sp.coo_matrix((G1triplet[0], G1triplet[1:]),
                              (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                              dtype = float)
        Gones = G1coo.tocsr()
        tmpVec =  np.empty(self.ckt.nD_dimension)
        # Must perform symbolic factorization first 
        self.__doSymbolic = True

        def get_deltax(xVec):
            (iVec, Jac) = self.get_i_Jac(xVec) 
            tmpVec[:] = Gones * xVec
            iVec += self.gmin * tmpVec
            Jac1 = Jac.tocsc() + self.gmin * Gones
            assert Jac1.format == 'csc'
            return self._get_deltax(iVec - sV, Jac1)

        def f_eval(xVec):
            iVec = self.get_i(xVec)
            tmpVec[:] = Gones * xVec
            iVec += self.gmin * tmpVec
            return iVec - sV

        (x, res, iterations, success) = \
            self._homotopy(0.5, self._set_gmin, x0, get_deltax, f_eval)
        # Make sure next time matrix is re-factorized next time
        self.__doSymbolic = True
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError('gmin stepping did not converge')



#---------------------------------------------------------------------------

class DCNodal(_NLFunctionSP):
    """
    Calculates the DC part of currents and Jacobian

    Matrices and vectors (G, Jac, s) are allocated here. This is to
    centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt):
        # Init base class
        super(DCNodal, self).__init__()

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        # Allocate vectors
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.iVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)
        if hasattr(self.ckt, 'nD_namRClist'):
            # Allocate external currents vector
            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    # Use custom function defined in this module
    set_Jac = _set_Jac

    def refresh(self):
        """
        Re-generate linear matrices (and Jacobian structure)

        Used for parameter sweeps
        """
        # Jac is (G1 + di/dv) in doc. 
        JacTriplet = ([], [], [])
        # Insert linear contributons: Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(JacTriplet, *vccs)
        # Frequency-defined elements
        for elem in self.ckt.nD_freqDefinedElem:
            set_Jac_triplet(JacTriplet, elem.nD_fpos, elem.nD_fneg, 
                            elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        # G here is G1 = G + G0 in documentation 
        Gcoo = sp.coo_matrix((JacTriplet[0], JacTriplet[1:]),
                             (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                             dtype = float)
        # Convert to compressed-row (csr) format for efficient
        # matrix-vector multiplication
        self.G = Gcoo.tocsr()

        # Create nonlinear structure.  Mark the end of linear part in
        # coo matrix
        self._mbaseLin = len(JacTriplet[0])

        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # The values that we insert do not matter, we are just
            # interested in the structure
            outJac = np.empty((len(elem.csOutPorts), len(elem.controlPorts)),
                              dtype = float)
            set_Jac_triplet(JacTriplet, elem.nD_cpos, elem.nD_cneg, 
                            elem.nD_vpos, elem.nD_vneg, outJac)
        self.Jaccoo = sp.coo_matrix((JacTriplet[0], JacTriplet[1:]),
                                    (self.ckt.nD_dimension, 
                                     self.ckt.nD_dimension), 
                                    dtype = float)

            
    def set_ext_currents(self, extIvec):
        """
        Set external currents applied to subcircuit

        extIvec: vector of external currents. Length of this vector
        should be equal to the number of external connections. The sum
        of all currents must be equal to zero (KCL)

        This will fail if not a subcircuit
        """
        # This idea still needs some testing
        assert sum(extIvec[:ncurrents]) == 0
        ncurrents = len(self.ckt.nD_namRClist)
        # Must do the loop in case there are repeated connections
        for val,rcnum in zip(extIvec, self.ckt.nD_namRClist):
            self.extSVec[rcnum] = val

    def get_guess(self):
        """
        Retrieve guesses from vPortGuess in each nonlinear device
        
        Returns a guess vector
        """
        x0 = np.zeros(self.ckt.nD_dimension)
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

    def get_source(self):
        """
        Get the source vector considering only the DC source components
        """
        # Erase vector first. 
        try:
            # If subcircuit add external currents
            self.sVec[:] = self.extSVec
        except AttributeError:
            # Not a subcircuit
            self.sVec.fill(0.)
        for elem in self.ckt.nD_sourceDCElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_DCsource()
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        return self.sVec
        

    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G xVec + i(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents
        """
        # Linear contribution
        self.iVec[:] = self.G * xVec
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            outV = elem.eval(xin)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
        return self.iVec

    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)::

            iVec = G xVec + i(xVec)
            Jac = G + (di/dx)(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        # Linear contribution
        self.iVec[:] = self.G * xVec
        # Nonlinear contribution
        self._mbase = self._mbaseLin
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            self.set_Jac(self.Jaccoo, outJac, *elem.nD_csidx)

        return (self.iVec, self.Jaccoo)


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
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            # Set OP in element (discard return value)
            elem.nD_xOP = xin
            elem.get_OP(xin)


    
#----------------------------------------------------------------------

class TransientNodal(_NLFunctionSP):
    """
    Keeps track of transient analysis equations. 

    This class only sets up transient equations. Equations are solved
    elsewhere. Circuit history required for numerical integration is
    kept by an integration method class. 

    Matrices and vectors (G, C, JacI, JacQ, s, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt, im):
        """
        Arguments:

        ckt: circuit instance

        im: Integration method instance. 

        """
        # Init base class
        super(TransientNodal, self).__init__()

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        self.im = im

        # Allocate matrices/vectors
        # G, C and G'

        # iVec = G' x + i'(x)   total current
        self.iVec = np.empty(self.ckt.nD_dimension)

        # Total charge: C x + q(x)
        self.qVec = np.empty(self.ckt.nD_dimension)
        # Source vector at current time s(t) 
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)

#        if hasattr(self.ckt, 'nD_namRClist'):
#            # Allocate external currents vector
#            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    # Use custom function defined in this module
    set_Jac = _set_Jac

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        # Jac is (G1 + di/dv) in doc. 
        JacTriplet = ([], [], [])
        # Insert linear contributons: Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(JacTriplet, *vccs)
#        # Frequency-defined elements not included for now
#        for elem in self.ckt.nD_freqDefinedElem:
#            set_Jac_triplet(JacTriplet, elem.nD_fpos, elem.nD_fneg, 
#                            elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        # Mark beginning of C
        self._Cbase = len(JacTriplet[0])
        # Add C matrix
        for elem in self.ckt.nD_elemList:
            for vccs in elem.nD_linVCQS:
                set_quad(JacTriplet, *vccs)
        # Save C values
        Coo = sp.coo_matrix((JacTriplet[0][self._Cbase:], 
                             (JacTriplet[1][self._Cbase:], 
                              JacTriplet[2][self._Cbase:])),
                            (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                            dtype = float)
        self._Cdata = Coo.data
        self.C = Coo.tocsr()
        
        # Create nonlinear structure.  Mark the end of linear part in
        # coo matrix
        self._mbaseLin = len(JacTriplet[0])

        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # The values that we insert do not matter, we are just
            # interested in the structure
            outJac = np.empty((len(elem.csOutPorts), len(elem.controlPorts)),
                              dtype = float)
            qJac = np.empty((len(elem.qsOutPorts), len(elem.controlPorts)),
                            dtype = float)
            set_Jac_triplet(JacTriplet, elem.nD_cpos, elem.nD_cneg, 
                            elem.nD_vpos, elem.nD_vneg, outJac)
            set_Jac_triplet(JacTriplet, elem.nD_qpos, elem.nD_qneg, 
                            elem.nD_vpos, elem.nD_vneg, qJac)
        self.Jaccoo = sp.coo_matrix((JacTriplet[0], JacTriplet[1:]),
                                    (self.ckt.nD_dimension, 
                                     self.ckt.nD_dimension), 
                                    dtype = float)

    def set_IC(self, h):
        """
        Set initial conditions (ICs)

        h: (initial) time step

        Retrieves ICs from DC operating point info in elements
        """
        # Nodal variables
        xVec = np.zeros(self.ckt.nD_dimension)
        # Get nodal voltages for tstart
        for i,term in enumerate(self.ckt.nD_termList):
            xVec[i] = term.nD_vOP
        # Calculate total charges
        self.update_q(xVec)
        # initialize integration element
        self.im.init(h, self.qVec)
        # Generate Gp 
        self.update_Gp()


    def update_Gp(self):
        """
        Recalculate Gp from im information
        """
        # Gp = G + a0 C in documentation (Y0 not implemented yet)
        data = self.Jaccoo.data
        row = self.Jaccoo.row
        col = self.Jaccoo.col
        # Calculate a0 C and save into Jacobian
        data[self._Cbase:self._mbaseLin] = self.im.a0 * self._Cdata
        Gcoo = sp.coo_matrix((data[:self._mbaseLin], 
                              (row[:self._mbaseLin], col[:self._mbaseLin])),
                             (self.ckt.nD_dimension, self.ckt.nD_dimension), 
                             dtype = float)
        # Convert to compressed-row (csr) format for efficient
        # matrix-vector multiplication
        self.Gp = Gcoo.tocsr()
        return self.Gp


    def update_q(self, xVec):
        """ 
        Recalculate qVec for a given value of xVec
        """
        # Calculate total q vector
        # Calculate linear charges first
        self.qVec[:] = self.C * xVec
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)  
            outV = elem.eval(xin)
            set_i(self.qVec, elem.nD_qpos, elem.nD_qneg, 
                  outV[len(elem.csOutPorts):])
        return self.qVec
            

    def get_source(self, ctime):
        """
        Get the source vector considering DC and TD source components

        ctime: current time.
        """
        # Erase vector first. 
        try:
            # If subcircuit add external currents
            self.sVec[:] = self.extSVec
        except AttributeError:
            # Not a subcircuit
            self.sVec.fill(0.)
        for elem in self.ckt.nD_sourceDCElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_DCsource()
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        for elem in self.ckt.nD_sourceTDElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_TDsource(ctime)
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        return self.sVec
        

    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G' xVec + i'(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents
        """
        # Linear contribution
        self.iVec[:] = self.Gp * xVec
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            outV = elem.eval(xin)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_i(self.iVec, elem.nD_qpos, elem.nD_qneg, 
                  self.im.a0 * outV[len(elem.csOutPorts):])
        return self.iVec

    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)::

            iVec = G' xVec + i'(xVec)
            Jac = G' + (di'/dx)(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        # Linear contribution
        self.iVec[:] = self.Gp * xVec
        # Nonlinear contribution
        self._mbase = self._mbaseLin
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_i(self.iVec, elem.nD_qpos, elem.nD_qneg, 
                  self.im.a0 * outV[len(elem.csOutPorts):])
            self.set_Jac(self.Jaccoo, outJac, *elem.nD_csidx)
            qJac = self.im.a0 * outJac[len(elem.csOutPorts):,:]
            self.set_Jac(self.Jaccoo, qJac, *elem.nD_qsidx)

        return (self.iVec, self.Jaccoo)


