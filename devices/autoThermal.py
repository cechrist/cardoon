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
        
        def __init__(self, instanceName):
            nle.__init__(self, instanceName)
            self.__addThermalPorts = True
    
        def process_params(self):
            """
            Process parameters as in base class

            Add extra thermal terminals if __addThermalPorts is True
            """
            # Let the base class do its stuff first
            nle.process_params(self)

            # Add thermal terminals to control and output tuples only
            # if needed. Base class must reset __addThermalPorts to
            # True if needed.
            if self.__addThermalPorts:
                # Add units to thermal port
                self.neighbour[self.numTerms-1].unit = 'C'
                self.neighbour[self.numTerms-2].unit = 'C'
                self.csOutPorts.append((self.numTerms-1, self.numTerms-2))
                self.controlPorts.append((self.numTerms-2, self.numTerms-1))
                # Thermal output number
                self.__ton = len(self.csOutPorts) - 1
                # Thermal control port number
                self.__tpn = len(self.controlPorts) - 1
                self.__addThermalPorts = False

            # Initial guess for input ports: 
            try:
                if len(self.vPortGuess) < len(self.controlPorts):
                    self.vPortGuess = np.concatenate(
                        (self.vPortGuess,[27.]), axis=0)
            except AttributeError:
                # Ignore if vPortGuess not provided
                pass

        def eval_cqs(self, vPort, saveOP = False):
            """
            vPort is a vector with control voltages (last port is thermal)
            """
            # set temperature in base class first
            nle.set_temp_vars(self, vPort[self.__tpn])

            # now calculate currents and charges
            if saveOP:
                (iVec1, qVec, opV) = nle.eval_cqs(self, vPort, saveOP)
            else:
                (iVec1, qVec) = nle.eval_cqs(self, vPort)
            # Calculate instantaneous power 
            pout = nle.power(self, vPort, iVec1)
            # Re-arrange output vector
            iVec = np.concatenate((iVec1, [pout]), axis = 0)
            if saveOP:
                return (iVec, qVec, opV)
            else:
                return (iVec, qVec)
    
        # Create these using the AD facility
        eval_and_deriv = ad.eval_and_deriv
        eval = ad.eval
    
    # Return template class
    return ThermalDevice
