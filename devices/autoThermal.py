"""
Generates an electrothermal device from a nonlinear device class

One assumption is that the base class defines numTerms directly in the
class definition and it is not changed in process_params().

-------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
#import circuit as cir
#from globalVars import const, glVar
# For automatic differentiation:
import cppaddev as ad


def thermal_device(nle):
    class ThermalDevice(nle):
        """
        Generic electrothermal nonlinear element

        Inherits from a regular nonlinear device class (nle). nle can
        be a linear device (like resistor.py) but it must implement
        all nonlinear functions to use this template. The
        electrothermal device is always nonlinear.

        Adds one thermal port (pair of terminals) connected after the
        regular terminals. Temperature in this port is in degree C.  A
        current source proportional to the instantaneous power
        dissipated in device is connected to the thermal port.
        """
    
        # devtype is the 'model' name
        devType = nle.devType + '_t'
    
        # Number of terminals. If numTerms is set here, the parser knows
        # in advance how many external terminals to expect. By default the
        # parser makes no assumptions and allows any number of connections
        #
        # Add two thermal terminals
        numTerms = nle.numTerms + 2
        
        # Force nonlinear behaviour (even if base class is linear, see
        # resistor.py)
        isNonlinear = True
        # Initial guess for input ports:
        vPortGuess = np.concatenate((nle.vPortGuess,[27.]), axis=0)
        
        def __init__(self, instanceName):
            """
            Here the base class constructor is called
            """
            nle.__init__(self, instanceName)
            self.__thermalFlag = False
    
        def process_params(self, circuit):
            """
            Takes the container circuit reference as an argument. 
    
            Basically do whatever the base class needs and add extra
            thermal terminals
            """
            # Let the base class do its stuff first
            nle.process_params(self, circuit)

            # Add thermal terminals to control and output tuples only
            # once
            if not self.__thermalFlag:
                thermalPort = (self.numTerms-1, self.numTerms-2)
                self.csOutPorts = self.csOutPorts + (thermalPort, )
                self.controlPorts = self.controlPorts + (thermalPort, )
                # Thermal output number
                self.__ton = len(self.csOutPorts) - 1
                # Thermal control port number
                self.__tpn = len(self.controlPorts) - 1
                # Mark this instance to avoid repeating this
                self.__thermalFlag = True

        def eval_cqs(self, vPort, saveOP = False):
            """
            vPort is a vector with control voltages (last port is
            thermal)
        
            Returns a numpy vector: currents first, thermal and then
            charges.
            """
            # set temperature in base class first
            nle.set_temp_vars(self, vPort[self.__tpn])

            # now calculate currents and charges
            if saveOP:
                (outV1, opV) = nle.eval_cqs(self, vPort, saveOP)
            else:
                outV1 = nle.eval_cqs(self, vPort)
            # Calculate instantaneous power 
            pout = nle.power(self, vPort, outV1)
            # Re-arrange output vector
            outV = np.concatenate((outV1[:self.__ton], [pout], 
                                   outV1[self.__ton:]), axis = 0)
            if saveOP:
                return (outV, opV)
            else:
                return outV
    
        # Create these using the AD facility
        eval_and_deriv = ad.eval_and_deriv
        eval = ad.eval
    
    # Return template class
    return ThermalDevice
