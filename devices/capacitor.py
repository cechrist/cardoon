"""
:mod:`capacitor` -- Linear capacitor
------------------------------------

.. module:: capacitor
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from globalVars import glVar
import circuit as cir

class Device(cir.Element):
    """
    Linear Capacitor
    ----------------

    Connection diagram::

                   || C
      0 o----------||---------o 1
                   ||

    Netlist example::

        cap:c1 1 2 c=10uF

    """

    # devtype is the 'model' name
    devType = "cap"

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
        c = ('Capacitance', 'F', float, 0.)
        )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)


    def process_params(self):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # Access to global variables is through the glVar 
        if not self.c:
            raise cir.CircuitError(self.nodeName 
                                   + ': Capacitance can not be zero')

        # Adjust according to temperature (not needed so far)
        # self.set_temp_vars(self.temp)
        self.linearVCQS = [((0, 1), (0, 1), self.c)]

#    def set_temp_vars(self, temp):
#        """
#        Calculate temperature-dependent variables for temp given in C
#        """
#        # Absolute temperature (note temp is in deg. C)
#        # T = const.T0 + temp
#        deltaT = temp - self.tnom
#        self.g /= (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)


