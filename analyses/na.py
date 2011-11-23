"""
:mod:`na` -- Nodal Approach
---------------------------

.. module:: na
.. moduleauthor:: Carlos Christoffersen
"""

class NodalCircuit:
    """
    Provides methods for nodal analysis of circuit

    Takes a circuit as an argument. If the circuit contains the 'gnd'
    node, it is used as the reference. Otherwise a reference node must
    be indicated.

    In the future this should also work with subcircuits. External
    connections have to handled for that.

    """

    def __init__(self, circuit, reference='gnd'):
        # 
        self.cir = circuit
        self.ref = circuit.get_term(reference)
        # Set nodal voltage to zero
        sef.ref.v = 0.
