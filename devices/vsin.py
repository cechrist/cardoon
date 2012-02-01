"""
:mod:`vsin` -- Sinusoidal voltage source
----------------------------------------

.. module:: vsin
.. moduleauthor:: Carlos Christoffersen
"""

import cmath as cm
import numpy as np
from globalVars import glVar
import circuit as cir

class Device(cir.Element):
    """
    (Co-)Sinusoidal voltage source
    ------------------------------

    Connection diagram::
                          
                   ,---,  vout       Rint
       0 o--------( - + )---------/\/\/\/\--------o 1
                   '---'  
                 
           vout = vdc + mag * cos(2 * pi * freq * t + phase)
   
    This source works for time and frequency domain. For AC analysis,
    the 'acmag' parameter is provided. By default acmag = mag.

    Netlist example::

        vsin:vin gnd 4 vdc=2V amp=1V freq=1GHz phase=90 


    Internal Topology
    +++++++++++++++++

    Implemented using a gyrator if Rint is zero::

                                       i/gyr       ti
        0  o---------+            +----------------+
                     | gyr V23    |                |
          +         /|\          /|\              /^\ 
        vin        | | |        | | | gyr vin    | | | gyr vout
          -         \V/          \V/              \|/  
                     |            |                |
        1  o---------+            +----------------+
                                          |
                                         --- tref
                                          V

    """

    # devtype is the 'model' name
    devType = "vsin"

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

    paramDict = dict(
        cir.Element.tempItem,
        vdc = ('DC voltage', 'V', float, 0.),
        rint = ('Internal resistance', 'Ohms', float, 0.),
        mag = ('Amplitude', 'V', float, 0.),
        acmag = ('Amplitude for AC analysis only', 'V', float, None),
        phase = ('Phase', 'degrees', float, 0.),
        freq = ('Frequency', 'Hz', float, 1e3)
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

        if self.acmag == None:
            # Use regular amplitude instead
            self._acmag = self.mag
        else:
            # Use AC amplitude
            self._acmag = self.acmag
            
        self._omega = 2. * np.pi * self.freq
        self._phase = self.phase * np.pi / 180.

        if self.rint:
            # Independent source attribute: output port 
            self.sourceOutput = (1, 0)
            # Can use Norton equivalent, no need for internal terminals
            self._idc = self.vdc / self.rint
            self._imag = self.mag / self.rint
            self._acimag = self._acmag / self.rint

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
            self._imag = glVar.gyr * self.mag
            self._acimag = glVar.gyr * self._acmag

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

    def get_DCsource(self):
        return self._idc

    def get_TDsource(self, time):
        """
        Returns source value at ctime
        """
        return self._imag * np.cos(self._omega * time + self._phase)


    def get_FDsource(self):
        """
        Returns a tuple with a frequency and a current phasor vectors
    
          (fvec, currentVec)
    
        """
        # used if isFDSource = True. fvec is defined by the source
        # parameters. 
        
        # Only one frequency component for sinusoidal wave
        fvec = np.array([self.freq])
        currentVec = np.array([cm.rect(self.imag, self._phase)])
        return (fvec, currentVec)

    def get_AC(self):
        """ 
        Returns AC magnitude and phase
        """
        return cm.rect(self._acimag, self._phase)



