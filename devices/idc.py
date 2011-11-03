"""
DC current source

-------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir

class Device(cir.Element):
    """
    DC current source. Includes temperature dependence
                ______ 
               /      \ idc
    0 o-------+  --->  +---------o 1
               \______/
    """

    # devtype is the 'model' name
    devType = "idc"

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
    isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Nonlinear device attributes
    # csOutPorts = ((0, 2), )
    # csContPorts = ((0, 3), (1, 3), (2, 3))
    # qsOutPorts = ( )
    # qsContPorts = ( )
    # csDelayedContPorts = ( )

    # Independent source attribute: output port
    sourceOutput = (0, 1)

    paramDict = dict(
        cir.Element.tempItem,
        idc = ('DC current', 'A', float, 0.),
        tnom = ('Nominal temperature', 'C', float, 27.),
        tc1 = ('Current temperature coefficient 1', '1/C', float, 0.),
        tc2 = ('Current temperature coefficient 2', '1/C^2', float, 0.)
        )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)


    def process_params(self, circuit):
        """
        Takes the container circuit reference as an argument. 

        Called once the external terminals have been connected and the
        non-default parameters have been set. Make sanity checks
        here. Internal terminals/devices should also be defined here
        (use circuit reference for this).  Raise cir.CircuitError if a fatal
        error is found.
        """
        # Access to global variables is through the glVar 

        # Adjust according to temperature
        self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Absolute temperature (note temp is in deg. C)
        # T = const.T0 + temp
        deltaT = temp - self.tnom
        self._adjIdc = self.idc * (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)

    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

    def getDCSource(self):
        return self._adjIdc

    # def getTDSource(self, ctime):
    #     """
    #     Documentation (isTDSource = True)
    #     ctime is the current time
    #     """
    #     # return current at ctime
    
    # def getFDSource(self, fvec):
    #     """
    #     Documentation (isFDSource = True)
    #     """
    #     # should return a np.array with currents for each frequency

# Here you can add additional functions and classes that only are
# visible withing this module.

