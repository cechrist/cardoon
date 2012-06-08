"""
:mod:`ipulse` -- Pulse current source
-------------------------------------

.. module:: ipulse
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir

# The following functions are shared with vpulse

def get_TDsource(dev, time):
    """
    Returns source value at ctime
    """
    # There may be a more efficient way of implementing this
    # creating (in process_params()) a list of event times and
    # actions
    if time > dev.per:
        # Remove whole periods from time
        t = time % dev.per
    else:
        t = time

    if t <= dev.td:
        return dev._i1
    else:
        t -= dev.td
        
    if t < dev.tr:
        # Ramp up
        iout = dev._i1 + dev._deltai * t / dev.tr
        return iout
    else:
        t -= dev.tr

    if t < dev.pw:
        return dev._i2
    else:
        t -= dev.pw

    # import pdb; pdb.set_trace()
    if t < dev.tf:
        # Ramp down
        iout = dev._i2 - dev._deltai * t / dev.tf
        return iout

    return dev._i1


def get_next_event(dev, time):
    """
    Returns time of next event
    """
    if time > dev.per:
        # Remove whole periods from time
        t = time % dev.per
    else:
        t = time

    if t <= dev.td:
        return time -t + dev.td
    else:
        t -= dev.td
        
    if t < dev.tr:
        # Ramp up
        return time - t + dev.tr
    else:
        t -= dev.tr


class Device(cir.Element):
    """
    Pulse current source
    --------------------

    Connection diagram::
                           
                   ,---,  iout
        0 o-------( --> )---------o 1
                   '---'    

        iout = pulse(t)

    This source only works for time domain. It is equivalent to an
    open circuit for DC or frequency-domain.

    Netlist example::

        ipulse:i1 gnd 4 i1=-1V i2=1V td=1ms pw=10ms per=20ms

    """
    # Device category
    category = "Sources"

    # devtype is the 'model' name
    devType = "ipulse"

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

    # Independent source attribute: output port
    sourceOutput = (0, 1)

    paramDict = dict(
#        cir.Element.tempItem,
        i1 = ('Initial value', 'A', float, 0.),
        i2 = ('Pulsed value', 'A', float, 0.),
        td = ('Delay time', 's', float, 0.),
        tr = ('Rise time', 's', float, 0.),
        tf = ('Fall time', 's', float, 0.),
        pw = ('Pulse width', 's', float, np.inf),
        per = ('Period', 's', float, np.inf)
        )

    def __init__(self, instanceName):
        cir.Element.__init__(self, instanceName)

    def process_params(self):
        # Access to global variables is through the glVar 
        # Adjust according to temperature
        #self.set_temp_vars(self.temp)

        self._i1 = self.i1
        self._i2 = self.i2

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
#        self._adjIdc = self.idc * (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)

    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

#    def get_DCsource(self):
#        return self.idc

    get_TDsource = get_TDsource
    get_next_event = get_next_event


