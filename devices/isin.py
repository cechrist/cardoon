"""
:mod:`isin` -- Sinusoidal current source
----------------------------------------

.. module:: isin
.. moduleauthor:: Carlos Christoffersen
"""

import cmath as cm
import numpy as np
from globalVars import const, glVar
import circuit as cir

class Device(cir.Element):
    """
    (Co-)Sinusoidal current source::

                           
                   ,---,  iout
        0 o-------( --> )---------o 1
                   '---'    

        iout = idc + mag * cos(2 * pi * freq + phase)

    This source works for time and frequency domain. For AC analysis,
    the 'acmag' parameter is provided. By default acmag = mag.

    Netlist example::

        isin:vdd gnd 4 idc=2mA amp=2mA freq=1GHz phase=90 

    """

    # devtype is the 'model' name
    devType = "isin"

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
    isTDSource = True
    isFDSource = True

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
        mag = ('Amplitude', 'A', float, 0.),
        acmag = ('Amplitude for AC analysis only', 'A', float, None),
        phase = ('Phase', 'degrees', float, 0.),
        freq = ('Frequency', 'Hz', float, 1e3)
        )

    def __init__(self, instanceName):
        cir.Element.__init__(self, instanceName)

    def process_params(self):
        # Access to global variables is through the glVar 
        # Adjust according to temperature
        #self.set_temp_vars(self.temp)
        if self.acmag == None:
            # Use regular amplitude instead
            self._acmag = self.mag
        else:
            # Use AC amplitude
            self._acmag = self.acmag
            
        self._omega = 2. * np.pi * self.freq
        self._phase = self.phase * np.pi / 180.

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

    def get_DCsource(self):
        return self.idc

    def get_TDsource(self, ctime):
        """
        ctime is the current time
        """
        # used if isTDSource = True
        # return current at ctime
        return self.amp * np.cos(self._omega * ctime + self._phase)

#    def get_next_event(self, ctime):
#        """
#        Returns time of next discontinuity in function/derivative
#        """
#        # Sinewave is smooth. Not needed
#        pass
    
    def get_FDsource(self):
        """
        Returns a tuple with a frequency and a current phasor vectors
    
          (fvec, currentVec)
    
        """
        # used if isFDSource = True. fvec is defined by the source
        # parameters. 
        
        # Only one frequency component for sinusoidal wave
        fvec = np.array([self.freq])
        currentVec = np.array([cm.rect(self.mag, self._phase)])
        return (fvec, currentVec)

    def get_AC(self):
        """ 
        Returns AC magnitude and phase
        """
        return cm.rect(self._acmag, self._phase)
