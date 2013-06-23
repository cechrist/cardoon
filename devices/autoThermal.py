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
from globalVars import glVar
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
        regular terminals. Temperature in this port is the difference
        with ambient temperature in degrees C.  A current source
        proportional to the instantaneous power dissipated in device
        is connected to the thermal port.
        """
    
        # devtype is the 'model' name
        devType = nle.devType + '_t'
            
        # Force nonlinear behaviour (even if base class is linear, see
        # resistor.py)
        isNonlinear = True
        
        def __init__(self, instanceName):
            nle.__init__(self, instanceName)
            self.__addThermalPorts = True
            if nle.numTerms:
                # Add two thermal terminals
                self.numTerms = nle.numTerms + 2
                self.__varTerms = False
            else:
                self.__varTerms = True
    
        def process_params(self):
            """
            Process parameters as in base class

            Add extra thermal terminals if __addThermalPorts is True
            """
            # Process parameters in base class
            nle.process_params(self, thermal = True)

            # Add thermal terminals to control and output tuples only
            # if needed. Base class must reset __addThermalPorts to
            # True if needed.
            if self.__addThermalPorts:
                # Add units to thermal port
                self.connection[self.numTerms-1].unit = \
                    '+{0} C'.format(glVar.temp)
                self.connection[self.numTerms-2].unit = \
                    '+{0} C'.format(glVar.temp)
                self.csOutPorts = self.csOutPorts + [(self.numTerms-1, 
                                                      self.numTerms-2)]
                self.controlPorts = self.controlPorts + [(self.numTerms-2, 
                                                          self.numTerms-1)]
                # Thermal output number
                self.__ton = len(self.csOutPorts) - 1
                # Thermal control port number
                self.__tpn = len(self.controlPorts) - 1
                self.__addThermalPorts = False

            # Initial guess for input ports: 
            try:
                # Consider time-delayed ports come after regular control ports
                if len(self.vPortGuess) < len(self.controlPorts) + self.nDelays:
                    self.vPortGuess = np.insert(self.vPortGuess, self.__tpn, 0.)
                        
            except AttributeError:
                # Ignore if vPortGuess not provided
                pass

        def eval_cqs(self, vPort, saveOP = False):
            """
            vPort is a vector with control voltages (last port is thermal)
            """
            # set temperature in base class first
            nle.set_temp_vars(self, vPort[self.__tpn] + glVar.temp)

            # Remove thermal port from vPort (needed in case
            # time-delayed ports follow regular control ports)
            vPort1 = np.delete(vPort, self.__tpn)
            # now calculate currents and charges
            if saveOP:
                (iVec1, qVec, opV) = nle.eval_cqs(self, vPort1, saveOP)
            else:
                (iVec1, qVec) = nle.eval_cqs(self, vPort1)
            # Calculate instantaneous power 
            pout = nle.power(self, vPort1, iVec1)
            # Re-arrange output vector
            iVec = np.append(iVec1, pout)
            if saveOP:
                return (iVec, qVec, opV)
            else:
                return (iVec, qVec)
    
        # Create these using the AD facility
        eval_and_deriv = ad.eval_and_deriv
        eval = ad.eval
    
    # Return template class
    return ThermalDevice
