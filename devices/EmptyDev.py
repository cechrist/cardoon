"""
Empty device code to use as a model to create new ones

-------------------------------------------------------------------

General Device Class Requirements
=================================

__init__(self, instanceName)

instanceName is a string 


Mandatory attributes
--------------------

Argument: devType = 'string' 

Specifies the name of the device model, for example 'res' for a
resistor model

List of all internal linear VCCS using the following format:

Argument: linearVCCS = [[(t0, t1), (t2, t3), g], ... ]


        0  o--------+      +------o 2
                           |      
          +               /|\       
        Vin              | | | g Vin     
          -               \V/       
                           |      
        1  o--------+      +------o 3

List of all internal linear VCQS (charge sources) using the following
format:

Argument: linearVCQS = [[(t0, t1), (t2, t3), c], ... ]

linearVCCS and linearVCQS may be empty lists and may be modified
by process_params() according to paramenter values.


-------------------------------------------------------------------

Nonlinear Devices
=================

Mandatory attributes
--------------------

isNonlinear = True

vPortGuess = numpy vector with a valid guess for the argument of
eval_cqs()

needsDelays = True or False

- Current source output ports (csOutPorts): for each current source
  in the device, list ports as (n1, n2). Current flows from n1 to n2.
  
  Example for a 3-terminal BJT with be and ce current sources,
  assuming teminals are connected c (0) - b (1) - e (2):
  
  csOutPorts = ((1, 2), (0, 2))

- Controlling ports (controlPorts): list here all ports whose voltages
  are needed to calculate the nonlinear currents / charges in same format.

  Example for BJT without intrinsic RC, RB and RE (vbc, vbe):
  controlPorts = ((1, 0), (1, 2))

  - Time-delayed port voltages (csDelayedContPorts): if needsDelays
    is True, list port voltages in triplet format (n1, n2,
    delay). Otherwise, this attribute is not required.

Similar vectors are required for output ports of charge sources:

qsOutPorts 

cs* and qs* (and noisePorts) can be modified by process_params()
according to parameter values.


-------------------------------------------------------------------

Independent Sources
===================

Must provide: 

1. At least one of the source flags set to True

2. A tuple with output port.
   Example: sourceOutput = (0, 1) for a current source.

3. Implement at least one of the source-related functions

-------------------------------------------------------------------

Noise current spectral density sources
======================================

Same format as csOutPorts (for nonlinear devices). Default is an empty
tuple.

Example:

noisePorts = ((1, 2), (0, 2))


-------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
# For automatic differentiation:
#import cppaddev as ad

class Device(cir.Element):
    """
    Class name must be "Device".
    Document here device type and terminal connection order. Example:

    This element does nothing. 

    Terminal order: 0 Drain, 1 Gate, 2 Source, 3 Bulk

             Drain 0
                     o
                     |
                     |
                 |---+
                 |
    Gate 1 o-----|<-----o 3 Bulk
                 |
                 |---+
                     |
                     |
                     o
            Source 2
    """

    # devtype is the 'model' name
    devType = "emptydev"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    # numTerms = 2
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    # isNonlinear = True
    # Initial guess for input ports:
    # vPortGuess = np.array([0., 0.])
    # needsDelays = True
    # isFreqDefined = True
    # isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Nonlinear device attributes
    # csOutPorts = ((0, 2), )
    # controlPorts = ((0, 3), (1, 3), (2, 3))
    # qsOutPorts = ( )
    # noisePorts = ( )
    # csDelayedContPorts = ( )

    # Independent source attribute: output port
    # sourceOutput = (0, 1)

    # Define parameters in a dictionary as follows: parameter name is
    # the key. Parameters are converted to class attributes after
    # circuit initialization.  If model is dependent on temperature,
    # include 'cir.Element.tempItem' ('temp' parameter description).
    paramDict = dict(
        cir.Element.tempItem,
        r = ('Resistance', 'Ohms', float, 0.),
        rsh = ('Sheet resistance', 'Ohms', float, 0.),
        l = ('Lenght', 'm', float, 0.),
        w = ('Width', 'm', float, 0.),
        )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed


    def process_params(self, circuit):
        """
        Takes the container circuit reference as an argument. 

        Called once the external terminals have been connected and the
        non-default parameters have been set. Make sanity checks
        here. Internal terminals/devices should also be defined here
        (use circuit reference for this).  Raise cir.CircuitError if a fatal
        error is found.
        """
        # Use the following to make sure connections to internal
        # terminals are not repeated if this process_params is called
        # many times. 
        self.clean_internal_terms(circuit)

        # Calculate temperature-dependent variables
        # self.set_temp_vars(self.temp)

        # Access to global variables is through const and glVar (e.g.,
        # const.u0)

        # if device is based on cppaddev, make sure tape is re-generated
        # ad.delete_tape(self)
        pass

#     For now do not implement this:
#
#     def newParams(self, circuit):
#         """
#         Recalculate any variable affected by new parameter values (or
#         even change the internal terminal structure if needed).
#         """
#         pass


#    def set_temp_vars(self, temp):
#        """
#        Calculate temperature-dependent variables for temp given in C
#        """
#        # Absolute temperature (note temp is in deg. C)
#        T = const.T0 + temp
#        # Thermal voltage
#        self.Vt = const.k * T / const.e


    #---------------------------------------------------------------
    # Nonlinear device methods (isNonlinear = True)
    #---------------------------------------------------------------


    # Returns a numpy vector: currents first and then charges. 
    # vPort is a vector (normally a numpy
    # array) with the values of controlling port voltages, including
    # delayed ports, if any.
    #
    # def eval_cqs(self, vPort, saveOP=False):
    #     """
    #     vPort is a vector with control voltages
    #
    #     Returns a numpy vector: currents first and then charges.
    #     If saveOP = True, return tuple with normal vector and OP 
    #     variables (only needed if ever saveOP is True, see resistor)
    #     """
    #     # should return np.array with current in each port
    #     #
    #     # To avoid automatic differentiation problems, use the 
    #     # condassign() function provided in cppaddev.py to replace 
    #     # if statements dependent on vPort
    #     pass
    #     if saveOP:
    #         return (outVec, opVars)
    #     else:
    #         return outVec

    # Required for automatic electrothermal device generation
    # def power(self, vPort, currV):
    #     """ 
    #     Calculate total instantaneous power 
    # 
    #     Input: control voltages and currents from eval_cqs()
    #     """
    #     vds = vPort[0] - vPort[2]
    #     # pout = vds*ids + vdb*idb + vsb*isb
    #     pout = vds*currV[0] + vPort[0] * currV[1] + vPort[2] * currV[2] 
    #     return pout

    # Operating point
    # def get_OP(self, vPort):
    #     """
    #     Calculates operating point information
    # 
    #     Input:  vPort = [vdb , vgb , vsb]
    #     Output: dictionary with OP variables
    #     """
    #     # First we need the Jacobian
    #     (outV, jac) = self.eval_and_deriv(vPort)
    #     # if this is not needed then saveOP flag does not have 
    #     # to be implemented
    #     opV = self.get_op_vars(vPort) 
    # 
    #     # Check things that change if the transistor is reversed
    #     if opV[11] > 0.:
    #         reversed = False
    #         gds = jac[0,0]
    #     else:
    #         reversed = True
    #         gds = jac[0,2]
    #         
    #     self.OP = {'VD': vPort[0],
    #                'VG': vPort[1],
    #                'VS': vPort[2],
    #                'IDS': outV[0]}

    #---------------------------------------------------------------------
    # Noise: in general requires a previous call to get_OP 
    #---------------------------------------------------------------------
    # 
    # def get_noise(self, f):
    #     """
    #     Return noise spectral density at frequency f
    #     
    #     Requires a previous call to get_OP() 
    #     """
    #     s1 = self.OP['Sthermal'] + self.OP['kSflicker'] / pow(f, self.af)
    #     s2 = something
    #     return np.array([s1, s2])

    # The following two functions should be present, normally
    # implemented by evaluating the AD tape (i.e. they run faster than
    # eval_cqs()). But we could also implement them manually by other
    # means:
    #
    # def eval(self, vPort): same as eval_cqs()
    # def eval_and_deriv(self, vPort): returns a tuple, (outVec, Jacobian)
    #
    # Use functions directly from cppaddev imported as ad:
    #eval_and_deriv = ad.eval_and_deriv
    #eval = ad.eval

    #---------------------------------------------------------------
    # Linear frequency-defined device methods (isFreqDefined = True)
    #---------------------------------------------------------------

    # def get_ymatrix(self, fvec)
    #     """
    #     Documentation 
    #     """
    #     # should return 3-D np.array with Y matrix for each frequency
    #     pass

    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

    # def get_DCsource(self):
    #     """
    #     Documentation (isDCSource = True)
    #     """
    #     # return current value

    # def get_TDsource(self, ctime):
    #     """
    #     Documentation (isTDSource = True)
    #     ctime is the current time
    #     """
    #     # return current at ctime
    
    # def get_FDsource(self, fvec):
    #     """
    #     Documentation (isFDSource = True)
    #     """
    #     # should return a np.array with currents for each frequency



# Here you can add additional functions and classes that only are
# visible withing this module.

