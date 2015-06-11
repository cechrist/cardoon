"""
:mod:`capacitor` -- Linear capacitor
------------------------------------

.. module:: capacitor
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from cardoon.globalVars import glVar
import cardoon.circuit as cir

class Device(cir.Element):
    """
    Linear Capacitor
    ----------------

    Connection diagram::


           +----------/\/\/\/\/\-----------+
           |             rint              |
           |                               |
           |                               |
           |       || C         resr       |
      0 o--+-------||----+---/\/\/\/\/\----+--o 1
                   ||  Term:cr

    Netlist example::

        cap:c1 1 2 c=10uF
        cap:c1 1 2 c=10uF resr = 0.01Ohm
        cap:c1 1 2 c=10uF rint = 1MOhm
        cap:c1 1 2 c=10uF resr = 0.01Ohm rint = 1MOhm

    """
    # Device category
    category = "Basic components"

    # devtype is the 'model' name
    devType = "nidcap"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    paramDict = dict(
        c = ('Capacitance', 'F', float, 0.),
        resr = ('Equivalent series resistance', 'Ohm', float, 0.),
        rint = ('Internal or dielectric resistance', 'Ohm', float, 0.)
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


        # remove any existing internal connections
        self.clean_internal_terms()

        # Access to global variables is through the glVar 
        if not self.c:
            raise cir.CircuitError(self.instanceName 
                                   + ': Capacitance can not be zero')
        
        # Adjust according to temperature (not needed so far)
        # self.set_temp_vars(self.temp)

        # Set up model
        # Add internal terminal
        
        if self.resr != 0.:
            cr = self.add_internal_term('cr', 'A')
            self.gesr = 1/self.resr
            self.linearVCQS = [((0, cr), (0, cr), self.c)]
            self.linearVCCS = [((1, cr), (1, cr), self.gesr)]
            # Add Internal/Dielectric resistance only if one is enterred
            if self.rint != 0.:
                self.linearVCCS.append(((0, cr), (0, cr), 1/self.rint))
        else:
            self.linearVCQS = [((0, 1), (0, 1), self.c)]
            # Add Internal/Dielectric resistance only if one is enterred
            if self.rint != 0.:
                self.linearVCCS = [((0, 1), (0, 1), 1/self.rint)]



#    def set_temp_vars(self, temp):
#        """
#        Calculate temperature-dependent variables for temp given in C
#        """
#        # Absolute temperature (note temp is in deg. C)
#        # T = const.T0 + temp
#        deltaT = temp - self.tnom
#        self.g /= (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)


