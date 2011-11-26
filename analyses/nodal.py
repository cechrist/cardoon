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


def set_quad(G, row1, col1, row2, col2, g):
    """
    Set transconductance/transcapacitance quad

    G: target matrix
    """
    # good cython candidate
    if col1 > 0:
        if row1 > 0:
            G[row1, col1] += g
        if row2 > 0:
            G[row2, col1] -= g
    if col2 > 0:
        if row1 > 0:
            G[row1, col2] -= g
        if row2 > 0:
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


class NodalCircuit:
    """
    Provides methods for nodal analysis of circuit

    Takes a circuit as an argument. If the circuit contains the 'gnd'
    node, it is used as the reference. Otherwise a reference node must
    be indicated.

    For now only DC is implemented.

    In the future this should also work with subcircuits. External
    connections have to handled for that.

    """
    def __init__(self, circuit, reference='gnd'):
        # Save circuit
        self.cir = circuit
        # get ground node
        self.ref = circuit.get_term(reference)

        # get a list of all terminals in circuit (internal/external)
        self.termList = circuit.termDict.values() + circuit.get_internal_terms()
        # remove ground node from terminal list
        self.termList.remove(self.ref)

        # Dimension is the number of unknowns to solve for
        self.dimension = len(self.termList)

        # For the future: use graph techniques to find the optimum
        # terminal order

        # Assign a number (0-inf) to all nodes. For the reference node
        # assign -1 
        self.ref.__namRC = -1
        for i, term in enumerate(self.termList):
            term.__namRC = i

        # Get a list of all elements and nonlinear devices/sources
        self.elemList = circuit.elemDict.values()
        self.nlinElements = filter(lambda x: x.isNonlinear, self.elemList)
        self.sourceDCElements = filter(lambda x: x.isDCSource, self.elemList)

        # Map row/column numbers directly into VC*S descriptions
        for elem in self.elemList:
            # Create list with RC numbers (choose one)
            # rcList = map(lambda x: x.__namRC, elem.neighbour)
            rcList = [x.__namRC for x in elem.neighbour]

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
            elem.__linVCCS = map(convert_vcs, elem.linearVCCS)
            elem.__linVCQS = map(convert_vcs, elem.linearVCQS)

            # Convert nonlinear device descriptions to a format more
            # readily usable for the NA approach
            if elem.isNonlinear:
                # Translate positive and negative terminal numbers
                def create_list(portlist):
                    tmp0 = [rcList[x1[0]] for x1 in portlist]
                    tmp1 = [rcList[x1[1]] for x1 in portlist]
                    return ([(i, j) for i,j in enumerate(tmp0) if j > -1],
                            [(i, j) for i,j in enumerate(tmp1) if j > -1])
                # Control voltages
                (elem.__vpos, elem.__vneg) = create_list(elem.controlPorts)
                # Current source terminals
                (elem.__cpos, elem.__cneg) = create_list(elem.csOutPorts)
                # Charge source terminals
                (elem.__qpos, elem.__qneg) = create_list(elem.qsOutPorts)

            # Translate source output terms
            if elem.isDCSource:
                # first get the destination row/columns 
                n1 = rcList[elem.sourceOutput[0]]
                n2 = rcList[elem.sourceOutput[1]]
                elem.__sourceOut = (n1, n2)

        # Generate G matrix (never changes)
        self.G = np.zeros((self.dimension, self.dimension))
        # Generate C here too
        self.C = np.zeros((self.dimension, self.dimension))
        for elem in self.elemList:
            for vccs in elem.__linVCCS:
                set_quad(self.G, *vccs)
            for vcqs in elem.__linVCQS:
                set_quad(self.C, *vcqs)



    def get_guess(self):
        """
        Retrieve guesses from vPortGuess in each nonlinear device
        
        Returns a guess vector
        """
        x0 = np.zeros(self.dimension)
        for elem in self.nlinElements:
            try:
                # Only add to positive side. This is not the only way
                # and may not work well in some cases but this is a
                # guess anyway
                for i,j in elem.__vpos:
                    x0[j] += elem.vPortGuess[i]
            except AttributeError:
                # if vPortGuess not given just leave things unchanged
                pass
        return x0


    def get_DC_source(self, sVec):
        """
        Get the source vector considering only the DC source components

        sVec is the source destination source vector. It is passed as
        an argument to avoid having to create a new vector from
        scratch each time this function is called.
        """
        # Erase vector first. 
        sVec[:] = 0.
        for elem in self.sourceDCElements:
            # first get the destination row/columns 
            outTerm = elem.__sourceOut
            current = elem.get_DCsource()
            if outTerm[0] > 0:
                sVec[outTerm[0]] += current
            if outTerm[1] > 0:
                sVec[outTerm[1]] -= current
        

    def get_DC_i(self, xVec, iVec):
        """
        Calculate total current

        iVec = G xVec + i(xVec)

        xVec: input vector of nodal voltages. 
        iVec: output vector of currents

        Output vector is passed as an argument to avoid having to
        create one from scratch each time this function is called.
        """
        # Erase vector
        iVec[:] = 0
        # Linear contribution
        iVec += np.dot(self.G, xVec)
        Jac += self.G
        # Nonlinear contribution
        for elem in self.nlinElements:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.__vpos, elem.__vneg, xVec)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(iVec, elem.__cpos, elem.__cneg, outV)


    def get_DC_i_Jac(self, xVec, iVec, Jac):
        """
        Calculate total current and Jacobian

        iVec = G xVec + i(xVec)
        Jac = G + Ji(xVec)

        xVec: input vector of nodal voltages. 
        iVec: output vector of currents
        Jac: system Jacobian

        Output vectors/matrix are passed as arguments to avoid having
        to create new ones from scratch each time this function is
        called.
        """
        # Erase vectors
        iVec[:] = 0
        Jac[:] = 0
        # Linear contribution
        iVec += np.dot(self.G, xVec)
        Jac += self.G
        # Nonlinear contribution
        for elem in self.nlinElements:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.__vpos, elem.__vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(iVec, elem.__cpos, elem.__cneg, outV)
            set_Jac(Jac, elem.__cpos, elem.__cneg, 
                    elem.__vpos, elem.__vneg, outJac)



    def save_OP(self, xVec):
        """
        Save nodal voltages in terminals and set OP in elements
        """
        # Set nodal voltage of reference to zero
        self.ref.__v = 0.
        for v,term in zip(xVec, self.termList):
            term.__v = v
        for elem in self.nlinElements:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.__vpos, elem.__vneg, xVec)
            # Set OP in element (discard return value)
            elem.get_OP(xin)

