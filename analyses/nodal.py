"""
:mod:`nodal` -- Nodal Approach
------------------------------

.. module:: nodal
.. moduleauthor:: Carlos Christoffersen

************** This is experimental/incomplete ****************

Some working comments to include in doc:

Tried advanced numpy indexing to avoid loops for each element of the
matrices but unfortunately G[obj] += Jac does not work when some of
the elements selected by obj are repeated (see tentative numpy
tutorial at http://www.scipy.org/NumPy_Tutorial)

So far there seem to be two ways to overcome this, to be tested as
soon as this code starts working: the first is to use a few cython
functions doing the inner loops. The second is to create a giant AD
tape with the whole circuit.

Note about storing nodal voltages in Terminals: currently voltages are
not stored by default. The main reason for this is efficiency as it is
less work to operate directly from the vector of unknowns in the
equation-solving routine.

It seems that it is possible to eliminate all the if statements in
set_i() and set_Jac() by eliminating the (-1)s from the row and column
vectors. The problem is that we would have to store the corresponding
non-zero indexes for current and Jac --> this is implemented now. May
need some fine tuning. For linear conductances we keep the if
statements as this is done only once.

"""

import numpy as np
from fsolve import fsolve, NoConvergenceError

# ****************** Stand-alone functions to be optimized ****************

def set_quad(G, row1, col1, row2, col2, g):
    """
    Set transconductance/transcapacitance quad

    G: target matrix
    """
    # good cython candidate
    #import pdb; pdb.set_trace()
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
    # good candidate for cython
    for i,j in posCols:
        xin[i] += xVec[j]
    for i,j in negCols:
        xin[i] -= xVec[j]


def set_i(iVec, posRows, negRows, current):
    """
    Set current contributions of current into iVec
    """
    # good candidate for cython
    for i,j in posRows:
        iVec[j] += current[i]
    for i,j in negRows:
        iVec[j] -= current[i]


def set_Jac(M, posRows, negRows, posCols, negCols, Jac):
    """
    Set current contributions of Jac into M
    """
    # good candidate for cython
    for i1, i in posRows:
        for j1, j in posCols:
            M[i,j] += Jac[i1,j1]
    for i1, i in negRows:
        for j1, j in negCols:
            M[i,j] += Jac[i1,j1]
    for i1, i in posRows:
        for j1, j in negCols:
            M[i,j] -= Jac[i1,j1]
    for i1, i in negRows:
        for j1, j in posCols:
            M[i,j] -= Jac[i1,j1]


# ********************** Regular functions *****************************

def make_nodal_circuit(ckt, reference='gnd'):
    """
    Add attributes to Circuit/Elements/Terminals for nodal analysis

    This functionalyty should be useful for any kind of nodal-based
    analysis (DC, AC, TRAN, HB, etc.)

    Takes a Circuit instance (ckt) as an argument. If the circuit
    contains the 'gnd' node, it is used as the reference. Otherwise a
    reference node must be indicated.

    New attributes are added in Circuit/Element/Terminal
    instances. All new attributes start with "nD_"

    In the future this should also work with subcircuits. External
    connections have to handled for that.
    """
    # get ground node
    ckt.nD_ref = ckt.get_term(reference)

    # get a list of all terminals in circuit (internal/external)
    ckt.nD_termList = ckt.termDict.values() + ckt.get_internal_terms()
    # remove ground node from terminal list
    ckt.nD_termList.remove(ckt.nD_ref)

    # Dimension is the number of unknowns to solve for
    ckt.nD_dimension = len(ckt.nD_termList)

    # For the future: use graph techniques to find the optimum
    # terminal order

    # Assign a number (0-inf) to all nodes. For the reference node
    # assign -1 
    ckt.nD_ref.nD_namRC = -1
    for i, term in enumerate(ckt.nD_termList):
        term.nD_namRC = i

    # Get a list of all elements and nonlinear devices/sources
    ckt.nD_elemList = ckt.elemDict.values()
    ckt.nD_nlinElem = filter(lambda x: x.isNonlinear, ckt.nD_elemList)
    ckt.nD_freqDefinedElem = filter(lambda x: x.isFreqDefined, ckt.nD_elemList)
    ckt.nD_sourceDCElem = filter(lambda x: x.isDCSource, ckt.nD_elemList)
    ckt.nD_sourceTDElem = filter(lambda x: x.isTDSource, ckt.nD_elemList)
    ckt.nD_sourceFDElem = filter(lambda x: x.isFDSource, ckt.nD_elemList)

    # Map row/column numbers directly into VC*S descriptions
    for elem in ckt.nD_elemList:
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

            The format of each list is: 

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

        # Convert frequency-defined elements
        if elem.isFreqDefined:
            (elem.nD_fpos, elem.nD_fneg) =  create_list(elem.fPortsDefinition)

        # Translate source output terms
        if elem.isDCSource or elem.isTDSource or elem.isFDSource:
            # first get the destination row/columns 
            n1 = rcList[elem.sourceOutput[0]]
            n2 = rcList[elem.sourceOutput[1]]
            elem.nD_sourceOut = (n1, n2)

        


# ****************** Classes ****************

#------------------------------------------------------------------------

class DCNodal:
    """
    Calculates the DC part of currents and Jacobian

    Matrices and vectors (G, Jac, s) are allocated here. This is to
    centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt):
        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        # List here the functions that can be used to solve equations
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_gmin, 
                                    self.solve_homotopy_source, 
                                    None]

        # Allocate matrices/vectors
        self.G = np.zeros((self.ckt.nD_dimension, self.ckt.nD_dimension))
        self.Jac = np.zeros((self.ckt.nD_dimension, self.ckt.nD_dimension))
        self.sVec = np.zeros(self.ckt.nD_dimension)
        self.iVec = np.zeros(self.ckt.nD_dimension)

        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.G, *vccs)
        # Frequency-defined elements
        for elem in self.ckt.nD_freqDefinedElem:
            set_Jac(self.G, elem.nD_fpos, elem.nD_fneg, 
                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
            
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
        self.sVec[:] = 0.
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
        # Erase arrays
        self.iVec[:] = 0.
        # Linear contribution
        self.iVec += np.dot(self.G, xVec)
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

        Returns: (iVec, Jac)

            iVec = G xVec + i(xVec)
            Jac = G + Ji(xVec)

        xVec: input vector of nodal voltages. 
        iVec: output vector of currents
        Jac: system Jacobian
        """
        # Erase arrays
        self.iVec[:] = 0.
        self.Jac[:] = 0.
        # Linear contribution
        self.iVec += np.dot(self.G, xVec)
        self.Jac += self.G
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
        """
        # Set nodal voltage of reference to zero
        self.ckt.nD_ref.nD_v = 0.
        for v,term in zip(xVec, self.ckt.nD_termList):
            term.nD_v = v
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            # Set OP in element (discard return value)
            elem.get_OP(xin)

    # The following functions used to solve equations, originally from
    # pycircuit
    def solve_simple(self, x0, sV):
        """Simple Newton's method"""
        def f_Jac_eval(x):
            (iVec, Jac) = self.get_i_Jac(x) 
            return (iVec - sV, Jac)
    
        def f_eval(x):
            iVec = self.get_i(x) 
            return iVec - sV
    
        return fsolve(x0, f_Jac_eval, f_eval)
    
        
    def solve_homotopy_source(self, x0, sV):
        """Newton's method with source stepping"""
        x = np.copy(x0)
        totIter = 0
        for lambda_ in np.linspace(start = .1, stop = 1., num = 10):
            def f_Jac_eval(x):
                (iVec, Jac) = self.get_i_Jac(x) 
                return (iVec - lambda_ * sV, Jac)
            
            def f_eval(x):
                iVec = self.get_i(x) 
                return iVec - lambda_ * sV
            (x, res, iterations) = fsolve(x, f_Jac_eval, f_eval)
            print('lambda = {0}, res = {1}, iter = {2}'.format(lambda_, 
                                                               res, iterations))
            totIter += iterations

        return (x, res, totIter)
    

    def solve_homotopy_gmin(self, x0, sV):
        """Newton's method with gmin stepping"""
        x = np.copy(x0)
        totIter = 0
        # Add a conductance in parallel with every node
        idx = np.arange(self.ckt.nD_dimension)
        for gmin in (1., .1, .01, 1e-3, 1e-4, 1e-5, 0.):
            # Add conductances to G matrix
            self.G[idx,idx] += gmin
            try:
                (x, res, iterations) = self.solve_simple(x, sV)
                print('gmin = {0}, res = {1}, iter = {2}'.format(gmin, 
                                                                 res, 
                                                                 iterations))
            except NoConvergenceError:
                # Restore original G matrix
                self.G[idx,idx] -= gmin
                raise

            totIter += iterations
            # Restore original G matrix
            self.G[idx,idx] -= gmin

        # Call solve_simple with better initial guess
        totIter += iterations

        return (x, res, totIter)
    
    
    

#----------------------------------------------------------------------

class TransientNodal(DCNodal):
    """
    Adds C and VCQS to DC part of currents and Jacobian

    Matrices and vectors (G, Jac, s) are allocated here

    Requires a NodalCircuit instance and possibly an integration method
    """

    def __init__(self, nodalCircuit):

        DCNodal.nD_init__(self, nodalCircuit)
        # Generate C here 
        self.C = np.zeros((self.ckt.nD_dimension, self.ckt.nD_dimension))
        for elem in self.ckt.nD_elemList:
            for vcqs in elem.nD_linVCQS:
                set_quad(self.C, *vcqs)

    def get_source(self):
        """
        Add time-domain component to DC
        """
        pass

    def get_i(self, xVec):
        """
        Add contribution of C and VCQS to DC. Requires integration method
        """
        pass

    def get_i_Jac(self, xVec):
        """
        Add contribution of C and VCQS to DC. Requires integration method
        """
        pass
