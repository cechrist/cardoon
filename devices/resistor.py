"""
:mod:`resistor` -- Resistor model
---------------------------------

.. module:: resistor
.. moduleauthor:: Carlos Christoffersen

Linear resistor plus definitions for nonlinear electrothermal model
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    """
    Resistor
    --------

    Connection diagram::

                    R
      0 o--------/\/\/\/---------o 1

    Normally a linear device. If the electro-thermal version is used
    (res_t), the device is nonlinear.

    Netlist examples::

        # Linear resistor (2 terminals)
        res:r1 1 2 r=1e3 tc1=10e-3

        # Electro-thermal resistor (nonlinear, 4 terminals)
        res_t:r1 1 2 3 4 r=1e3 tc1=10e-3

    """
    # Device category
    category = "Basic components"

    # devtype is the 'model' name
    devType = "res"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    
    # Create electrothermal device
    makeAutoThermal = True

    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    # isNonlinear = True
    # needsDelays = True
    # isFreqDefined = True
    # isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Nonlinear device attributes (used for electrothermal model)
    csOutPorts = [(0, 1)]
    noisePorts = [(0, 1)]
    controlPorts = [(0, 1)]
    qsOutPorts = [ ]
    vPortGuess = np.array([0.])

    paramDict = dict(
        cir.Element.tempItem,
        r = ('Resistance', 'Ohms', float, 0.),
        rsh = ('Sheet resistance', 'Ohms', float, 0.),
        l = ('Lenght', 'm', float, 0.),
        w = ('Width', 'm', float, 0.),
        narrow = ('Narrowing due to side etching', 'm', float, 0.),
        tnom = ('Nominal temperature', 'C', float, 27.),
        tc1 = ('Temperature coefficient 1', '1/C', float, 0.),
        tc2 = ('Temperature coefficient 2', '1/C^2', float, 0.)
        )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)


    def process_params(self, thermal = False):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # Delete AD tape (if any)
        ad.delete_tape(self)
        # Access to global variables is through the glVar 
        if not self.r and not (self.rsh and self.l and self.w):
            raise cir.CircuitError(self.nodeName 
                                   + ': Resistance can not be zero')

        if self.r:
            # if R is given it overrides rsh et ale
            self.g = 1. / self.r
        else:
            self.g = (self.w-self.narrow) / (self.l-self.narrow) / self.rsh

        if not thermal:
            # Adjust according to temperature
            self.set_temp_vars(self.temp)

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Delete AD tape (if any)
        ad.delete_tape(self)
        # Absolute temperature (note temp is in deg. C)
        # T = const.T0 + temp
        deltaT = temp - self.tnom
        self.g /= (1. + (self.tc1 + self.tc2 * deltaT) * deltaT)
        self._St =  4. * const.k * (temp + const.T0) * self.g
        # Set this here to keep VCCS in sync with temperature 
        self.linearVCCS = [((0, 1), (0, 1), self.g)]


    # Implemented nonlinear functions for electrothermal model only:
    def eval_cqs(self, vPort):
        """
        Return current given port voltage. Charge vector is empty
        """
        return (self.g * vPort, np.array([]))


    def power(self, vPort, ioutV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages and currents from eval_cqs()
        """
        return vPort[0] * ioutV[0]

    def get_OP(self, vPort):
        """
        Calculates operating point information

        (Should add calculation of simple capacitive model)
        Input:  vPort = [vr] (may also have temperature)
        Output: dictionary with OP variables
        """
        # Get g at the requested temperature (in case of thermal resistor)
        self.eval_cqs(vPort)
        self.OP = {'R': 1./self.g, 'Sthermal': self._St}
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density (constant)
        
        Requires a previous call to get_OP() 
        """
        return np.array([self.OP['Sthermal']])
