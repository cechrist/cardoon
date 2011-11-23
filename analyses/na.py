"""
:mod:`na` -- Nodal Approach
---------------------------

.. module:: na
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np

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
        self.ref.rcnumber = -1
        for i, term in enumerate(termList):
            term.rcnumber = i

        # Dimension is the number of unknowns to solve for
        self.dimension = len(termList)

        # Get a list of all elements and nonlinear devices/sources
        self.elemList = circuit.elemDict.values()
        self.nlinElements = filter(lambda x: x.isNonlinear, elemList)
        self.sourceDCElements = filter(lambda x: x.isDCSource, elemList)
        
        
    def get_DC_source(self, sVec, time):
        """
        Get the source vector considering only the DC source components

        sVec is the source destination source vector. It is passed as
        an argument to avoid having to create a new vector from
        scratch each time this function is called.
        """
        # Erase vector first. Not sure if this is more efficient
        sVec *= 0.
        for elem in self.sourceDCElements:
            # first get the destination row/columns 
            outTerm1 = elem.neighbour[elem.sourceOutput[0]]
            outTerm2 = elem.neighbour[elem.sourceOutput[1]]
            current = elem.get_DCsource()
            sVec[outTerm1] += current
            sVec[outTerm2] -= current

    def get_i_G(self, iVec):
        pass

