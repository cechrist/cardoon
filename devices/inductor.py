"""
:mod:`inductor` -- Linear inductor
-----------------------------------

.. module:: inductor
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir

class Device(cir.Element):
    """
    Inductor
    --------

    Connection diagram::

                 __  __  __  _ 
        0       /  \/  \/  \/ \          1
          o----+   /\  /\  /\  +-------o    External view
                  (_/ (_/ (_/  

    Netlist example::

        ind:l1 1 0 l=3uH


    Internal Topology
    +++++++++++++++++

    Internal implementation uses a gyrator (adds il internal node)::

                                        il/gyr    Term: il
        0  o---------+            +----------------+
                     | gyr V(il)  |                |
          +         /|\          /^\               |
        Vin        ( | )        ( | ) gyr Vin    ----- gyr^2 * L
          -         \V/          \|/             -----
                     |            |                |
        1  o---------+            +----------------+
                                          |
                                         --- tref 
                                          V


    """
    # Device category
    category = "Basic components"

    # devtype is the 'model' name
    devType = "ind"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    paramDict = dict(
        l = ('Inductance', 'H', float, 0.)
        )

    def __init__(self, instanceName):
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)


    def process_params(self):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # remove any existing internal connections
        self.clean_internal_terms()

        if not self.l:
            raise cir.CircuitError(self.instanceName 
                                   + ': Inductance can not be zero')
        # Connect internal terminal
        til = self.add_internal_term('il', '{0} A'.format(glVar.gyr)) 
        tref = self.add_reference_term() 
        # Setup gyrator
        # Access to global variables is through the glVar 
        self.linearVCCS = [((0,1), (tref, til), glVar.gyr), 
                           ((til, tref), (0,1), glVar.gyr)]
        cap = self.l * glVar.gyr**2
        self.linearVCQS = [((til, tref), (til, tref), cap)]

        # Adjust according to temperature (not needed so far)
        # self.set_temp_vars(self.temp)

#    def set_temp_vars(self, temp):
#        """
#        Calculate temperature-dependent variables for temp given in C
#        """
#        # Absolute temperature (note temp is in deg. C)
#        # T = const.T0 + temp
#        deltaT = temp - self.tnom
#        self.g /= (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)


