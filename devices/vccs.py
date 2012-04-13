"""
:mod:`vccs` -- Voltage-controlled current source
------------------------------------------------

.. module:: vccs
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir

class Device(cir.Element):
    r"""
    Voltage-controlled current source
    ---------------------------------

    Schematic::

                  g Vcont
                   ,---,    
        0 o-------( --> )---------o 1
                   `---`     


        2 o      + Vcont -        o 3

    Temperature dependence:

    .. math::
        
      g(T) = g(T_{nom}) (1 + t_{c1} \Delta T + t_{c2} \Delta T^2)

      \Delta T = T - T_{nom}

    Netlist example::

        vccs:g1 gnd 4 3 gnd g=2mS

    """

    # devtype is the 'model' name
    devType = "vccs"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 4
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M, or linear controlled sources
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
        cir.Element.tempItem,
        g = ('Linear transconductance', 'S', float, 0.),
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


    def process_params(self):
        # Access to global variables is through the glVar 

        # Adjust according to temperature
        self.set_temp_vars(self.temp)
        self.linearVCCS = [((2,3), (0,1), self._gT)]
        

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Absolute temperature (note temp is in deg. C)
        # T = const.T0 + temp
        deltaT = temp - self.tnom
        self._gT = self.g * (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)




