"""
:mod:`vpulse` -- Pulse voltage source
-------------------------------------

.. module:: vpulse
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from cardoon.globalVars import glVar
import cardoon.circuit as cir
import ipulse

class Device(cir.Element):
    """
    Pulse voltage source
    --------------------

    Connection diagram::
                          
                   ,---,  vout       Rint
       1 o--------( - + )---------/\/\/\/\--------o 0
                   '---'  
                 
           vout = vpulse(t)
   
    This source only works for time domain. It is equivalent to a
    short circuit (or rint) for DC or frequency-domain.

    Netlist example::

        vpulse:vin 4 0 v1=-1V v2=1V td=1ms pw=10ms per=20ms


    Internal Topology
    +++++++++++++++++

    Same as vdc.

    """
    # Device category
    category = "Sources"

    # devtype is the 'model' name
    devType = "vpulse"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    #isDCSource = True
    isTDSource = True
    #isFDSource = True

    paramDict = dict(
#        cir.Element.tempItem,
        v1 = ('Initial value', 'V', float, 0.),
        v2 = ('Pulsed value', 'V', float, 0.),
        td = ('Delay time', 's', float, 0.),
        tr = ('Rise time', 's', float, 0.),
        tf = ('Fall time', 's', float, 0.),
        pw = ('Pulse width', 's', float, np.inf),
        per = ('Period', 's', float, np.inf),
        rint = ('Internal resistance', 'Ohms', float, 0.)
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

        if self.rint != 0.:
            # Independent source attribute: output port 
            self.sourceOutput = (1, 0)
            # Can use Norton equivalent, no need for internal terminals
            self._i1 = self.v1 / self.rint
            self._i2 = self.v2 / self.rint

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
            self._i1 = glVar.gyr * self.v1
            self._i2 = glVar.gyr * self.v2

        self._deltai = self._i2 - self._i1
        if self.pw < np.inf:
            self._rem = self.per - self.tr - self.pw - self.tf 
        else:
            self._rem = 0.

#    def set_temp_vars(self, temp):
#        """
#        Calculate temperature-dependent variables for temp given in C
#        """
#        # Absolute temperature (note temp is in deg. C)
#        # T = const.T0 + temp
#        deltaT = temp - self.tnom
#        self._adjIdc = self._idc \
#            * (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)

    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

# Note: we could calculate some DC component but seems too much work
# for little use.
#
#    def get_DCsource(self):
#        return self._i1

    # Use current-source functions for easier maintenance
    get_TDsource = ipulse.get_TDsource
    get_next_event = ipulse.get_next_event
