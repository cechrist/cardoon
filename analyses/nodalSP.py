"""
:mod:`nodalSP` -- Nodal Analysis using pysparse
-----------------------------------------------

.. module:: nodalSP
.. moduleauthor:: Carlos Christoffersen

This module contains basic classes/functions for nodal analysis. These
are part of the standard netlist analyses but they can also be used
independently.

This implementation uses pysparse. It is much more efficient than the
dense implementation but still a lot could be gained with finer
control over matrix factorization. At this time the main Jacobian is
(almost) created and factored from scratch at every iteration.

Also there is a lot of overlap with the ``nodal`` module. 

"""

import numpy as np
import pysparse
from fsolve import fsolve_Newton, NoConvergenceError
from integration import BEuler

# ****************** Stand-alone functions (perhaps to be optimized) **********

def set_quad(G, row1, col1, row2, col2, g):
    """
    Set transconductance/transcapacitance quad

    G: target matrix
    """
    if col1 >= 0:
        if row1 >= 0:
            G[row1, col1] += g
        if row2 >= 0:
            G[row2, col1] -= g
    if col2 >= 0:
        if row1 >= 0:
            G[row1, col2] -= g
        if row2 >= 0:
            G[row2, col2] += g


def set_xin(xin, posCols, negCols, xVec):
    """
    Calculate input port voltage xin
    """
    for i,j in posCols:
        xin[i] += xVec[j]
    for i,j in negCols:
        xin[i] -= xVec[j]


def set_i(iVec, posRows, negRows, current):
    """
    Set current contributions of current into iVec
    """
    for i,j in posRows:
        iVec[j] += current[i]
    for i,j in negRows:
        iVec[j] -= current[i]


def set_Jac(M, posRows, negRows, posCols, negCols, Jac):
    """
    Set current contributions of Jac into M
    """
    for i1, i in posRows:
        for j1, j in posCols:
            M[i,j] += Jac[i1,j1]
        for j1, j in negCols:
            M[i,j] -= Jac[i1,j1]
    for i1, i in negRows:
        for j1, j in negCols:
            M[i,j] += Jac[i1,j1]
        for j1, j in posCols:
            M[i,j] -= Jac[i1,j1]


# ********************** Regular functions *****************************

def make_nodal_circuit(ckt, reference='gnd'):
    """
    Add attributes to Circuit/Elements/Terminals for nodal analysis

    This functionality should be useful for any kind of nodal-based
    analysis (DC, AC, TRAN, HB, etc.)

    Takes a Circuit instance (ckt) as an argument. If the circuit
    contains the 'gnd' node, it is used as the reference. Otherwise a
    reference node must be indicated.

    New attributes are added in Circuit/Element/Terminal
    instances. All new attributes start with ``nD_``

    Works with subcircuits too (see ``nD_namRClist`` attribute)
    """
    # get ground node
    ckt.nD_ref = ckt.get_term(reference)

    # make a list of all non-reference terminals in circuit 
    ckt.nD_termList = ckt.termDict.values() + ckt.get_internal_terms()
    # remove ground node from terminal list
    ckt.nD_termList.remove(ckt.nD_ref)
    # Assign a number (0-inf) to all nodes. For reference nodes
    # assign -1 
    ckt.nD_ref.nD_namRC = -1
    # Make a list of all elements
    ckt.nD_elemList = ckt.elemDict.values()
    # Set RC number of reference terminals to -1
    for elem in ckt.nD_elemList:
        if elem.localReference:
            elem.neighbour[elem.localReference].nD_namRC = -1
    # For the future: use graph techniques to find the optimum
    # terminal order
    for i, term in enumerate(ckt.nD_termList):
        term.nD_namRC = i

    # Store internal RC numbers for later use
    for elem in ckt.nD_elemList:
        elem.nD_intRC = [term.nD_namRC for term in 
                         elem.neighbour[elem.numTerms:]]

    # Dimension is the number of unknowns to solve for
    ckt.nD_dimension = len(ckt.nD_termList)
    # Number of external terminals excluding reference
    ckt.nD_nterms = len(ckt.termDict.values()) - 1

    # Create specialized element lists
    ckt.nD_nlinElem = filter(lambda x: x.isNonlinear, ckt.nD_elemList)
    ckt.nD_freqDefinedElem = filter(lambda x: x.isFreqDefined, ckt.nD_elemList)
    ckt.nD_sourceDCElem = filter(lambda x: x.isDCSource, ckt.nD_elemList)
    ckt.nD_sourceTDElem = filter(lambda x: x.isTDSource, ckt.nD_elemList)
    ckt.nD_sourceFDElem = filter(lambda x: x.isFDSource, ckt.nD_elemList)

    # Map row/column numbers directly into VC*S descriptions
    for elem in ckt.nD_elemList:
        process_nodal_element(elem)

    # Subcircuit-connection processing
    try:
        connectTerms = ckt.get_connections()
        # List of RC numbers of external connections
        ckt.nD_namRClist = [term.nD_namRC for term in connectTerms]
    except AttributeError:
        # Not a subcircuit
        pass

def restore_RCnumbers(elem):
    """
    Restore RC numbers in internal terminals

    Assumption is number of internal terminals is the same
    """
    for term, i in zip(elem.neighbour[elem.numTerms:], elem.nD_intRC):
        term.nD_namRC = i

def process_nodal_element(elem):
    """
    Process element for nodal analysis
    """
    # Create list with RC numbers (choose one)
    # rcList = map(lambda x: x.nD_namRC, elem.neighbour)
    rcList = [x.nD_namRC for x in elem.neighbour]

    # Translate linear VCCS/VCQS format
    def convert_vcs(x):
        """
        Converts format of VC*S 

        input: [(contn1, contn2), (outn1, outn2), g]

        output: [row1, col1, row2, col2, g]
        """
        col1 = rcList[x[0][0]]
        col2 = rcList[x[0][1]]
        row1 = rcList[x[1][0]]
        row2 = rcList[x[1][1]]
        return [row1, col1, row2, col2, x[2]]
    elem.nD_linVCCS = map(convert_vcs, elem.linearVCCS)
    elem.nD_linVCQS = map(convert_vcs, elem.linearVCQS)

    # Translate positive and negative terminal numbers
    def create_list(portlist):
        """
        Converts an internal port list into 2 lists with (+-) nodes

        The format of each list is::

            (internal term number, namRC number)
        """
        tmp0 = [rcList[x1[0]] for x1 in portlist]
        tmp1 = [rcList[x1[1]] for x1 in portlist]
        return ([(i, j) for i,j in enumerate(tmp0) if j > -1],
                [(i, j) for i,j in enumerate(tmp1) if j > -1])

    # Convert nonlinear device descriptions to a format more
    # readily usable for the NA approach
    if elem.isNonlinear:
        # Control voltages
        (elem.nD_vpos, elem.nD_vneg) = create_list(elem.controlPorts)
        # Current source terminals
        (elem.nD_cpos, elem.nD_cneg) = create_list(elem.csOutPorts)
        # Charge source terminals
        (elem.nD_qpos, elem.nD_qneg) = create_list(elem.qsOutPorts)
# Disabled since not needed for now
#            # Create list with external ports for gmin calculation
#            nports = elem.numTerms - 1
#            elem.nD_extPorts = [(i, nports) for i in range(nports)]
#            (elem.nD_epos, elem.nD_eneg) = create_list(elem.nD_extPorts)

    # Convert frequency-defined elements
    if elem.isFreqDefined:
        (elem.nD_fpos, elem.nD_fneg) =  create_list(elem.fPortsDefinition)

    # Translate source output terms
    if elem.isDCSource or elem.isTDSource or elem.isFDSource:
        # first get the destination row/columns 
        n1 = rcList[elem.sourceOutput[0]]
        n2 = rcList[elem.sourceOutput[1]]
        elem.nD_sourceOut = (n1, n2)



#-------------------------------------------------------------------------
# ****************************** Classes *********************************
#-------------------------------------------------------------------------

class _NLFunction:
    """
    Nonlinear function interface class

    This is an abstract class to be use as a base class by other
    nodal-based classes such as DCNodal
    """

    def __init__(self):
        # List here the functions that can be used to solve equations
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_source, 
                                    self.solve_homotopy_gmin, 
                                    None]

    def _get_deltax(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc
        """
        umf = pysparse.umfpack.factorize(
            Jac, strategy="UMFPACK_STRATEGY_UNSYMMETRIC")
        umf.solve(errFunc, self.deltaxVec)
        return self.deltaxVec

    # The following functions used to solve equations, originally from
    # pycircuit
    def solve_simple(self, x0, sV):
        #"""Simple Newton's method"""
        # Docstring removed to avoid printing this all the time
        def get_deltax(x):
            (iVec, Jac) = self.get_i_Jac(x) 
            return self._get_deltax(iVec - sV, Jac)
    
        def f_eval(x):
            iVec = self.get_i(x) 
            return iVec - sV
    
        return fsolve_Newton(x0, get_deltax, f_eval)
    
    def solve_homotopy_source(self, x0, sV):
        """Newton's method with source stepping"""
        x = np.copy(x0)
        totIter = 0
        for lambda_ in np.linspace(start = .1, stop = 1., num = 10):
            def get_deltax(x):
                (iVec, Jac) = self.get_i_Jac(x) 
                return self._get_deltax(iVec - lambda_ * sV, Jac)
            
            def f_eval(x):
                iVec = self.get_i(x) 
                return iVec - lambda_ * sV
            (x, res, iterations) = fsolve_Newton(x, get_deltax, f_eval)
            print('lambda = {0}, res = {1}, iter = {2}'.format(lambda_, 
                                                               res, iterations))
            totIter += iterations

        return (x, res, totIter)

    def solve_homotopy_gmin(self, x0, sV):
        """Newton's method with gmin stepping"""
        x = np.copy(x0)
        totIter = 0
        idx = np.arange(self.ckt.nD_nterms)
        for gminexp in np.arange(-1., -8., -1):
            gmin = 10.**gminexp
            val = gmin * np.ones(self.ckt.nD_nterms)
            # Add gmin from ground to all external nodes. Assumes all
            # external nodes are sorted first in the vector. This will
            # not work if the terminal order is changed.
            def get_deltax(xvec):
                (iVec, Jac) = self.get_i_Jac(xvec) 
                iVec[idx] += gmin * xvec[idx]
                Jac.update_add_at(val, idx, idx)
                return self._get_deltax(iVec - sV, Jac)
            def f_eval(xvec):
                iVec = self.get_i(xvec)
                iVec[idx] += gmin * xvec[idx]
                return iVec - sV
            (x, res, iterations) = fsolve_Newton(x0, get_deltax, f_eval)
            print('gmin = {0}, res = {1}, iter = {2}'.format(
                    gmin, res, iterations))
            totIter += iterations

        # Call solve_simple with better initial guess
        (x, res, iterations) = self.solve_simple(x, sV)
        print('gmin = 0, res = {0}, iter = {1}'.format(res, iterations))
        totIter += iterations

        return (x, res, totIter)

# Old code for gmin homotopy: add conductances in parallel to
# nonlinear device external ports
#
#        for elem in self.ckt.nD_nlinElem:
#            Gadd = np.eye(len(elem.nD_extPorts))
#            set_Jac(Ggmin, elem.nD_epos, elem.nD_eneg, 
#                    elem.nD_epos, elem.nD_eneg, Gadd)


#---------------------------------------------------------------------------

class DCNodal(_NLFunction):
    """
    Calculates the DC part of currents and Jacobian

    Matrices and vectors (G, Jac, s) are allocated here. This is to
    centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt):
        # Init base class
        _NLFunction.__init__(self)

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        # Allocate matrices/vectors
        # G here is G1 = G + G0 in documentation (allocated in refresh())

        # Jac is (G1 + di/dv) in doc
        self.Jac = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.iVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)
        if hasattr(self.ckt, 'nD_namRClist'):
            # Allocate external currents vector
            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parameter sweeps
        """
        self.Gll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.Gll, *vccs)
        # Frequency-defined elements
        for elem in self.ckt.nD_freqDefinedElem:
            set_Jac(self.Gll, elem.nD_fpos, elem.nD_fneg, 
                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        # Free unused memory
        self.Gll.compress()        
        # Create G in csr form for efficient matrix-vector multiplication
        self.G = self.Gll.to_csr()

            
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
                # Only add to positive side. This is not the only way
                # and may not work well in some cases but this is a
                # guess anyway
                for i,j in elem.nD_vpos:
                    x0[j] += elem.vPortGuess[i]
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
        self.G.matvec(xVec, self.iVec)
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
        # Erase sparse matrix
        self.Jac.scale(0.)
        # Linear contribution
        self.G.matvec(xVec, self.iVec)
        self.Jac.shift(1., self.Gll)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_Jac(self.Jac, elem.nD_cpos, elem.nD_cneg, 
                    elem.nD_vpos, elem.nD_vneg, outJac)

        return (self.iVec, self.Jac)

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

class TransientNodal(_NLFunction):
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
        _NLFunction.__init__(self)

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        self.im = im

        # Allocate matrices/vectors
        # G, C and G'
        self.Gll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        self.Cll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # iVec = G' x + i'(x)   total current
        self.iVec = np.empty(self.ckt.nD_dimension)
        # System Jacobian: G' + di'/dx
        self.Jac = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # Total charge: C x + q(x)
        self.qVec = np.empty(self.ckt.nD_dimension)
        # Source vector at current time s(t) 
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)

#        if hasattr(self.ckt, 'nD_namRClist'):
#            # Allocate external currents vector
#            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        self.Gll.scale(0.)
        self.Cll.scale(0.)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.Gll, *vccs)
            for vccs in elem.nD_linVCQS:
                set_quad(self.Cll, *vccs)
#        # Frequency-defined elements not included for now
#        for elem in self.ckt.nD_freqDefinedElem:
#            set_Jac(self.G, elem.nD_fpos, elem.nD_fneg, 
#                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        self.Gll.compress()
        self.Cll.compress()
        self.C = self.Cll.to_csr()

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
        self.Gpll = self.Gll.copy()
        self.Gpll.shift(self.im.a0, self.Cll)
        self.Gp = self.Gpll.to_csr()
        return self.Gp


    def update_q(self, xVec):
        """ 
        Recalculate qVec for a given value of xVec
        """
        # Calculate total q vector
        # Calculate linear charges first
        self.C.matvec(xVec, self.qVec)
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
        self.Gp.matvec(xVec, self.iVec)
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
        # Erase sparse matrix
        self.Jac.scale(0.)
        # Linear contribution
        self.Gp.matvec(xVec, self.iVec)
        self.Jac.shift(1., self.Gpll)
        # Nonlinear contribution
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
            set_Jac(self.Jac, elem.nD_cpos, elem.nD_cneg, 
                    elem.nD_vpos, elem.nD_vneg, outJac)
            qJac = self.im.a0 * outJac[len(elem.csOutPorts):,:]
            set_Jac(self.Jac, elem.nD_qpos, elem.nD_qneg, 
                    elem.nD_vpos, elem.nD_vneg, qJac)

        return (self.iVec, self.Jac)


