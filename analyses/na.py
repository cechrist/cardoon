"""
:mod:`na` -- Nodal Approach
---------------------------

.. module:: na
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np

def setQuad(G, idx, g):
    """
    Set VCCS quad in G NAM

    idx are the G indexes in the following format::

        [(on1, cn1), (on2, cn2), (on2, cn1), (on1, cn2)]

    """
    # This is a good candidate function to compile with cython?
    G[idx[0]] = g
    G[idx[1]] = g
    G[idx[2]] = - g
    G[idx[3]] = - g


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
        # get a list of all terminals in circuit (internal/external)
        termList = circuit.termDict.values() + circuit.get_internal_terms()

        # get ground node
        self.ref = circuit.get_term(reference)
        # Set nodal voltage of reference to zero
        self.ref.v = 0.
        # remove ground node from terminal list
        termList.remove(self.ref)

        # Assign a number (0-inf) to all nodes except gnd=-1. For the
        # future: use graph techniques to find the optimum terminal
        # order
        self.ref.__namRC = -1
        for i, term in enumerate(termList):
            term.__namRC = i

        # Dimension is the number of unknowns to solve for
        self.dimension = len(termList)

        # Get a list of all elements and nonlinear devices/sources
        self.elemList = circuit.elemDict.values()
        self.nlinElements = filter(lambda x: x.isNonlinear, elemList)
        self.sourceDCElements = filter(lambda x: x.isDCSource, elemList)

        # Map row/column numbers directly into VC*S descriptions
        for elem in elemList:
            # Create list with RC numbers
            rcList = map(lambda x: x.__namRC, elem.neighbour)

            # Create quad pairs for linearVC*S
            def convert_vcs(x):
                """
                Converts format of VC*S 

                input: [(cn1, cn2), (on1, on2), g]
                output: [(on1, cn1), (on2, cn2), (on2, cn1), (on1, cn2)]
                """
                cn1 = rcList[x[0][0]]
                cn2 = rcList[x[0][1]]
                on1 = rcList[x[1][0]]
                on2 = rcList[x[1][1]]
                # May omit unnecesary pairs here?
                return [(on1, cn1), (on2, cn2), (on2, cn1), (on1, cn2)]

            elem.__linVCCidx = map(convert_vcs, elem.linearVCCS)
            elem.__linVCQidx = map(convert_vcs, elem.linearVCQS)

            # Convert nonlinear device descriptions
            if elem.isNonlinear:
                # Rather than doing this, here we should (somehow)
                # precalculate the indexes of all Jacobian entries
                def convert_port(x):
                    n1 = rcList[x[0]]
                    n2 = rcList[x[1]]
                    return (n1, n2)
                elem.__controlPorts = map(convert_port, elem.controlPorts)
                elem.__csOutPorts = map(convert_port, elem.csOutPorts)
                elem.__qsOutPorts = map(convert_port, elem.qsOutPorts)
                # add code for time-delayed sources

            # Translate source output terms
            if elem.isDCSource:
                # first get the destination row/columns 
                n1 = rcList[elem.sourceOutput[0]]
                n2 = rcList[elem.sourceOutput[1]]
                elem.__namSourceOut = (n1, n2)

        # Generate G matrix (never changes)
        self.G = np.zeros((self.dimension, self.dimension))
        for elem in elemList:
            # check for VCCS in element
            for i, vccs in enumerate(elem.linearVCCS):
                # Must decide which is worse: always checking for gnd
                # or keeping the gnd row/column and deleting it later.
                #
                set_quad(self.G, elem.__linVCCidx[i], vccs[2])

        # Should generate C here too
        self.C = np.zeros((self.dimension, self.dimension))
        for elem in elemList:
            # check for VCCS in element
            for i, vcqs in enumerate(elem.linearVCQS):
                # Must decide which is worse: always checking for gnd
                # or keeping the gnd row/column and deleting it later.
                #
                set_quad(self.C, elem.__linVCQidx[i], vcqs[2])

        
        
    def get_DC_source(self, sVec, time):
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
            outTerm = elem.__namSourceOut
            current = elem.get_DCsource()
            if outTerm1 > 0:
                sVec[outTerm[0]] += current
            if outTerm2 > 0:
                sVec[outTerm[1]] -= current

    def get_i_G_DC(self, xVec, iVec, Jac):
        """
        Calculate total current and Jacobian

        iVec = G xVec + i(xVec)
        Jac = G + Ji(xVec)
        """
        pass

