"""
:mod:`svdiode` -- State-variable-based diode model
--------------------------------------------------

.. module:: svdiode
.. moduleauthor:: Carlos Christoffersen

Diode model (from freeda/carrot source, in turn derived from spice
model)
"""

# use: import pdb; pdb.set_trace() for debuging

import numpy as np
from globalVars import const, glVar
import circuit as cir
from svjunction import Junction
import cppaddev as ad

class Device(cir.Element):
    """
    State-Variable-Based Diode device (based on Spice model)::

            o  1                           
            |                            
          --+--
           / \     
          '-+-'
            |                          
            o  0    	                  

    Internally represented as::

        0  o
           |
           \ 
           / Rs
           \ 
           / 
           |                                     2
        4  o---------+                  +----------------+
                     | i(x)+dq/dt       |                |
          +         /|\                /|\ gyr vin      /^\ 
        vin        | | |              | | |            | | | gyr v(x)
          -         \V/                \V/              \|/  
                     |                  |                |
        1  o---------+                  +------+---------+
                                               |
                                              --- (terminal 3 is gnd)
                                               V

    Terminal 4 not present if Rs = 0
    """

    # devtype is the 'model' name
    devType = "svdiode"

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
    isNonlinear = True
    # Initial guess for state variable
    vPortGuess = np.array([0.])
    # needsDelays = True
    # isFreqDefined = True
    # isDCSource = True
    # isTDSource = True
    # isFDSource = True

    # Nonlinear device attributes
    # csOutPorts = ((0, 2), )
    # csContPorts = ((0, 3), (1, 3), (2, 3))
    qsOutPorts = ( )
    # csDelayedContPorts = ( )

    # Independent source attribute: output port
    # sourceOutput = (0, 1)

    # Define parameters (note most parameters defined in svjunction.py)
    paramDict = dict(
        cir.Element.tempItem + Junction.paramList,
        ibv = ('Current at reverse breakdown voltage', 'A', float, 1e-10),
        bv = ('Breakdown voltage', 'V', float, 0.),
        area = ('Area multiplier', ' ', float, 1.),
        rs = ('Series resistance', 'Ohms', float, 0.),
        kf = ('Flicker noise coefficient', '', float, 0.),
        af = ('Flicker noise exponent', '', float, 1.),
       )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed
        self.jtn = Junction()

    def process_params(self, circuit):
        """
        Takes the container circuit reference as an argument. 
        """
        # Remove internal terminals
        self.clean_internal_terms(circuit)
        # Define topology first
        # Needs at least one internal terminal: 
        termList = [self.nodeName + ':n2', 'gnd']
        circuit.connect_internal(self, termList)
        self.linearVCCS = [[(0,1), (2,3), glVar.gyr]]
        # Nonlinear device attributes
        self.csOutPorts = ((0, 1), (3, 2))
        self.noisePorts = ((0, 1), )
        self.controlPorts = ((2, 3), )

        if self.rs:
            # Needs one more terminal
            circuit.connect_internal(self, [self.nodeName + ':n3,'])
            g = 1. / self.rs / self.area
            self.linearVCCS.append([(0,2), (0,2), g])
            # Nonlinear device outputs change
            self.csOutPorts = ((4, 1), (3, 2))
            self.noisePorts = ((0, 4), (4, 1))

        if self.tt or self.cj0:
            # Add charge source (otherwise the charge calculation is ignored)
            self.qsOutPorts = (self.csOutPorts[0], )

        # Calculate variables in junction
        self.jtn.process_params(self)
        # Calculate temperature-dependent variables
        self.set_temp_vars(self.temp)
        # Make sure tape is re-generated
        ad.delete_tape(self)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Everything is handled by the PN junction
        self.jtn.set_temp_vars(temp, self, ad)
        self._Sthermal = 4. * const.k * (temp + const.T0) * self.rs


    #---------------------------------------------------------------
    # Nonlinear device methods (isNonlinear = True)
    #---------------------------------------------------------------
    def eval_cqs(self, vPort):
        """
        Models intrinsic diode and charge.  

        vPort[0]: x as in Rizolli's equations

        Returns a vector with 3 elements: current, voltage and charge
        """
        # Calculate state-variable PN junction current, voltage and charge
        outV = self.jtn.eval_cqs(vPort[0], self, ad)

        # add breakdown current: here we need safe_exp() because state
        # variable only transforms the forward exponential
        if (self.bv > 0.):
            outV[0] -= self.ibv * \
                ad.safe_exp(-(outV[1] + self.bv) / self.n / self._vt)

        # area effect (voltage is not scaled by area)
        outV[0] *= self.area
        outV[2] *= self.area

        # Scale voltage to make it more similar to currents
        outV[1] *= glVar.gyr

        return outV

    def power(self, vPort, ioutV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages and currents from eval_cqs()
        """
        # ioutV[0]: current
        # ioutV[1]: gyr*voltage
        return ioutV[0] * ioutV[1] / glVar.gyr


    def get_OP(self, vPort):
        """
        Calculates operating point information

        Input:  vPort = [x]  (state variable)
        Output: dictionary with OP variables
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)
        # opV = self.get_op_vars(vPort) not needed for now 

        # Use negative indexing for charge in case power is inserted
        # between current and charge.
        # Note scaling factor in voltage must be considered
        self.OP = dict(
            x = vPort[0],
            VD = outV[1] / glVar.gyr,
            ID = outV[0],
            gd = glVar.gyr * jac[0,0] / jac[1,0],
            Cd = jac[-1,0],
            Sthermal = self._Sthermal,
            Sshot = 2. * const.q * outV[0],
            kFliker = self.kf * outV[0]
            )

        return self.OP
                   

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 
        """
        sj = self.OP['Sshot'] + self.OP['kSflicker'] / pow(f, self.af)
        if self.rs:
            sV = np.array([self.OP['Sthermal'], sj])
        else:
            sV = np.array([sj])
        return sV



    # Use automatic differentiation for eval and deriv functions
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval


# Here you can add additional functions and classes that only are
# visible withing this module.

