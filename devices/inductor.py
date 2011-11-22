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
    Linear inductor::

                 __  __  __  _ 
        0       /  \/  \/  \/ \          1
          o----+   /\  /\  /\  +-------o    External view
                  (_/ (_/ (_/  

    Netlist example::

        ind:l1 1 0 l=3uH


    Internal Topology
    +++++++++++++++++

    Internal implementation uses a gyrator (adds one internal node)::

                                          2
        0  o---------+            +----------------+
                     | gyr V21    |                |
          +         /|\          /^\               |
        Vin        | | |        | | | gyr Vin    ----- gyr^2 * L
          -         \V/          \|/             -----
                     |            |                |
        1  o---------+------------+----------------+


    """

    # devtype is the 'model' name
    devType = "ind"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    # isNonlinear = True
    # needsDelays = True
    # isFreqDefined = True
    # isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Nonlinear device attributes
    # csOutPorts = ((0, 2), )
    # csContPorts = ((0, 3), (1, 3), (2, 3))
    # qsOutPorts = ( )
    # qsContPorts = ( )
    # csDelayedContPorts = ( )

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
            raise cir.CircuitError(self.nodeName 
                                   + ': Inductance can not be zero')
        # Connect internal terminal
        self.add_internal_terms(1)
        # Set unit for internal term
        self.neighbour[2].unit = '* {0} A'.format(1./glVar.gyr)
        # Setup gyrator
        # Access to global variables is through the glVar 
        self.linearVCCS = [((0,1), (1,2), glVar.gyr), 
                           ((2,1), (0,1), glVar.gyr)]
        cap = self.l * glVar.gyr * glVar.gyr
        self.linearVCQS = [((2, 1), (2, 1), cap)]

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

# Here you can add additional functions and classes that only are
# visible withing this module.

