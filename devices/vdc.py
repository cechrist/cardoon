"""
:mod:`vdc` -- DC voltage source
-------------------------------

.. module:: vdc
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import glVar
import circuit as cir

class Device(cir.Element):
    """
    DC voltage source
    -----------------

    Schematic::
                          
                   ,---,  vdc       Rint
       1 o--------( - + )---------/\/\/\/\--------o 0
                   '---'  
   
    Rint is independent of temperature. Teperature dependence of vdc
    is as follows:

    .. math::
        
      v_{DC}(T) = v_{DC}(T_{nom}) (1 + t_{c1} \Delta T + t_{c2} \Delta T^2)

      \Delta T = T - T_{nom}

    Netlist example::

        vdc:vdd vddnode gnd vdc=3V


    Internal Topology
    +++++++++++++++++

    A gyrator is used to convert a current source into a voltage
    source if Rint is zero::

                                       i/gyr      Term: i
        0  o---------+            +----------------+
                     | gyr V(i)   |                |
          +         /|\          /|\              /^\ 
        vin        ( | )        ( | ) gyr vin    ( | ) gyr vdc
          -         \V/          \V/              \|/  
                     |            |                |
        1  o---------+            +----------------+
                                          |
                                         --- tref
                                          V

    Otherwise a Norton equivalent circuit is used.

    """
    # Device category
    category = "Sources"

    # devtype is the 'model' name
    devType = "vdc"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    isDCSource = True
    # isTDSource = True
    # isFDSource = True


    paramDict = dict(
        cir.Element.tempItem,
        vdc = ('DC voltage', 'V', float, 0.),
        rint = ('Internal resistance', 'Ohms', float, 0.),
        tnom = ('Nominal temperature', 'C', float, 27.),
        tc1 = ('Voltage temperature coefficient 1', '1/C', float, 0.),
        tc2 = ('Voltage temperature coefficient 2', '1/C^2', float, 0.)
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
        if self.rint != 0.:
            # Independent source attribute: output port 
            self.sourceOutput = (1, 0)
            # Can use Norton equivalent, no need for internal terminals
            self._idc = self.vdc / self.rint
            self.linearVCCS = [((0,1), (0,1), 1. / self.rint)]
        else:
            # Connect internal terminal
            ti = self.add_internal_term('i', '{0} A'.format(glVar.gyr))
            tref = self.add_reference_term()
            # Independent source attribute: output port 
            self.sourceOutput = (tref, ti)
            # Setup gyrator
            # Access to global variables is through the glVar (e.g.,
            # glVar.u0)
            self.linearVCCS = [((0, 1), (ti, tref), glVar.gyr), 
                               ((ti, tref), (0, 1), glVar.gyr)]
            # sourceOutput should be already OK
            self._idc = glVar.gyr * self.vdc

        # Adjust according to temperature
        self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Absolute temperature (note temp is in deg. C)
        # T = const.T0 + temp
        deltaT = temp - self.tnom
        self._adjIdc = self._idc \
            * (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)

    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

    def get_DCsource(self):
        return self._adjIdc



