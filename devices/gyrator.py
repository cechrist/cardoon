"""
:mod:`gyrator` -- Basic gyrator
-------------------------------

.. module:: gyrator
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    """
    Gyrator
    -------

    The gyrator converts Port 1 voltage into Port 2 current and
    *vice-versa*. Combined with the VCCS device it can be used to
    implement all the remaining controlled sources:

      * VCVS = VCCS + gyrator
      
      * CCCS = gyrator + VCCS

      * CCVS = gyrator + VCCS + gyrator

    Connection diagram::

            0  o---------+            +----------o 2
        +                |            |                +
                        /|\          /^\               
       Vin1     g Vin2 ( | )        ( | ) g Vin1      Vin2
                        \V/          \|/               
        -                |            |                -
            1  o---------+            +----------o 3

    Netlist example::

        gyr:gg 1 0 2 0 g=1m

    """
    # Device category
    category = "Controlled Sources"

    # devtype is the 'model' name
    devType = "gyr"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 4
    
    paramDict = dict(
        g = ('Gyrator gain', 'Ohms', float, 1e-3)
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

        # Setup gyrator
        # Access to global variables is through the glVar 
        self.linearVCCS = [((0, 1), (3, 2), self.g), 
                           ((2, 3), (0, 1), self.g)]



