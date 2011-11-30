"""
:mod:`empty` -- Description goes here
-------------------------------------

.. module:: empty
.. moduleauthor:: <your name>

Empty device code to use as a model to create new ones
"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
# For automatic differentiation:
import cppaddev as ad

class Device(cir.Element):
    """
    Class name must be "Device".
    Document here device type and terminal connection order. Example:

    This element does nothing. 

    Terminal order: 0 Drain, 1 Gate, 2 Source, 3 Bulk::

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
    numTerms = 2
    
    # Flags: at least one should be enabled for other than
    # linear R, L, C, M
    #
    # isNonlinear = True
    # Initial guess for input ports:
    # vPortGuess = np.array([0., 0.])
    # needsDelays = True
    # isFreqDefined = True
    # fPortsDefinition = [(0, 1), (2, 3)]
    # isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Local reference for internal voltages (by default last external terminal)
    #localReference = 3

    # Nonlinear device attributes
    # csOutPorts = [(0, 2)]
    # controlPorts = [(0, 3), (1, 3), (2, 3)]
    # qsOutPorts = [ ] 
    # noisePorts = [ ]
    # csDelayedContPorts = [ ]

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
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed


    def process_params(self):
        """
        This should be called each time parameter values are changed
        """
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # Use the following to make sure connections to internal
        # terminals are not repeated if this process_params is called
        # many times. 
        self.clean_internal_terms()

        # Calculate temperature-dependent variables
        # self.set_temp_vars(self.temp)

        # Access to global variables is through const and glVar (e.g.,
        # const.u0)

        # if device is based on cppaddev, make sure tape is re-generated
        # ad.delete_tape(self)
        pass


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Absolute temperature (note temp is in deg. C)
        T = const.T0 + temp
        # Thermal voltage
        self.Vt = const.k * T / const.e


    #---------------------------------------------------------------
    # Nonlinear device methods (isNonlinear = True)
    #---------------------------------------------------------------

    def eval_cqs(self, vPort, saveOP=False):
        """
        vPort is a vector with control voltages
    
        Returns tuple with two numpy vectors: one for currents and
        another for charges.

        If saveOP = True, return tuple with normal vectors and OP 
        variables 
        """
        # Note that saveOP may be ommitted, see documentation

        # calculation here
        iVec = np.array([i1, i2])
        qVec = np.array([q1])
        if saveOP:
            # calculate opVars
            return (iVec, qVec, opVars)
        else:
            return (iVec, qVec)


    # Required for automatic electrothermal device generation
    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 
    
        Input: control voltages and currents from eval_cqs()
        """
        vds = vPort[0] - vPort[2]
        # pout = vds*ids + vdb*idb + vsb*isb
        pout = vds*currV[0] + vPort[0] * currV[1] + vPort[2] * currV[2] 
        return pout

    def get_OP(self, vPort):
        """
        Calculates operating point information
    
        Input:  vPort = [vdb , vgb , vsb]
        Output: dictionary with OP variables
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)
        # if this is not needed then saveOP flag does not have 
        # to be implemented
        opV = self.get_op_vars(vPort) 
    
        # Check things that change if the transistor is reversed
        if opV[11] > 0.:
            reversed = False
            gds = jac[0,0]
        else:
            reversed = True
            gds = jac[0,2]
            
        self.OP = {'VD': vPort[0],
                   'VG': vPort[1],
                   'VS': vPort[2],
                   'IDS': outV[0]}

    #---------------------------------------------------------------------
    # Noise: in general requires a previous call to get_OP 
    #---------------------------------------------------------------------
    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        May require a previous call to get_OP() 
        """
        s1 = self.OP['Sthermal'] + self.OP['kSflicker'] / pow(f, self.af)
        s2 = something
        return np.array([s1, s2])

    # The following two functions should be present, normally
    # implemented by evaluating the AD tape (i.e. they run faster than
    # eval_cqs()). But we could also implement them manually by other
    # means:
    #
    # def eval(self, vPort): same as eval_cqs()
    # def eval_and_deriv(self, vPort): returns a tuple, (outVec, Jacobian)
    #
    # Use functions directly from cppaddev imported as ad instead:
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval

    #---------------------------------------------------------------
    # Linear frequency-defined device methods (isFreqDefined = True)
    #---------------------------------------------------------------
    def get_Y_matrix(self, fvec):
        """
        Returns Y matrix for all frequencies
        fvec: vector/scalar
        Frequencies in fvec can not be zero
        """
        # should return 3-D np.array, one vector per element of Y matrix
        pass

    def get_G_matrix(self):
        """
        Returns a matrix with the DC G parameters
        """
        # Return DC Y matrix
        return np.array([[y11, y12],
                         [y12, y11]])


    #---------------------------------------------------------------
    # Sources: DC always is always added to TD or FD
    #---------------------------------------------------------------

    def get_DCsource(self):
        # used if isDCSource = True
        # return current value
        pass
    
    def get_TDsource(self, ctime):
        """
        ctime is the current time
        """
        # used if isTDSource = True
        # return current at ctime
        pass
    
    def get_FDsource(self, fvec):
        """
        fvec: frequency vector
        May not work for f=0
        """
        # used if isFDSource = True
        # should return a np.array with currents for each frequency
        pass



# Here you can add additional functions and classes that only are
# visible within this module.

