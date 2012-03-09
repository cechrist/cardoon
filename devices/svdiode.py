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
import cppaddev as ad

#-----------------------------------------------------------------------
class SVJunction:
    """
    P-N Juction model.
    
    Based on carrot source, in turn based on spice/freeda diode
    model. This is intended to model any regular P-N junction such as
    Drain-Bulk, diodes, Collector-Bulk, etc.
    
    Transit time, breakdown and area effects not included
    """
    def process_params(self, isat, n, fc, cj0, vj, m, xti, eg0, Tnomabs):
        """
        Calculate variables dependent on parameter values only
        """
        # Saturation current
        self.isat = isat 
        # Emission coefficient
        self.n = n
        # Coefficient for forward-bias depletion capacitance
        self.fc = fc
        # Zero-bias depletion capacitance
        self.cj0 = cj0
        # Built-in junction potential
        self.vj = vj
        # PN junction grading coefficient
        self.m = m
        # Set some handy variables
        self._k1 = xti / self.n
        self._k2 = const.q * eg0 / self.n / const.k / Tnomabs
        self._k3 = const.q * eg0 / self.n / const.k
        self._k4 = 1. - self.m


    def set_temp_vars(self, obj):
        """
        Calculate temperature-dependent variables for temp given in C

        obj is an object instance containing the following attributes:
        tnratio, Tabs, Tnomabs, vt, egapn, egap_t
        """
        self._t_is = self.isat * pow(obj.tnratio, self._k1) \
            * np.exp(self._k2 - self._k3 / obj.Tabs) 
        # Maximum argument in exponential function (no need to use
        # safe_exp() with this model)
        self._alpha = 1. / self.n / obj.vt
        max_exp_arg = 5e10
        self._svth = np.log(max_exp_arg / self._alpha) / self._alpha
        self._kexp = self.n * obj.vt * max_exp_arg
        # Capacitance
        if self.cj0:
            self._t_vj = self.vj * obj.tnratio \
                - 3. * obj.vt * np.log(obj.tnratio) \
                - obj.tnratio * obj.egapn + obj.egap_t
            self._t_cj0 = self.cj0 * (1. + self.m 
                                     * (.0004 * (obj.Tabs - obj.Tnomabs) 
                                         + 1. - self._t_vj / self.vj))
            self._k5 = self._t_vj * self._t_cj0 / self._k4
            self._k6 = self._t_cj0 * pow(1. - self.fc, - self.m - 1.)
            self._k7 = ((1. - self.fc * (1. + self.m)) * self._t_vj * self.fc 
                        + .5 * self.m * self._t_vj * self.fc * self.fc)

    def get_idvd(self, x):
        """
        Returns junction a tuple (current, voltage)

        x: state variable
        """
        # Static current
        b = self._svth - x
        c = self._t_is * (np.exp(self._alpha * x) - 1.)
        d = self._t_is * self._kexp * \
            (1. + self._alpha * (x - self._svth)) - self._t_is
        iD = ad.condassign(b, c, d)
        # Diode voltage
        d = self._svth + \
            np.log(1. + self._alpha * (x - self._svth)) / self._alpha
        vD = ad.condassign(b, x, d)
        return (iD, vD)


    def get_qd(self, vd):
        """
        Returns junction depletion charge

        vd: diode voltage
        """
        b = self.fc * self._t_vj - vd
        c = self._k5 * (1. - pow(1. - vd / self._t_vj, self._k4))
        d = self._k6 * ((1. - self.fc * (1. + self.m)) 
                        * vd + .5 * self.m * vd * vd / self._t_vj 
                        - self._k7) + self._k5 \
                        * (1. - pow(1. - self.fc, self._k4))
        return ad.condassign(b, c, d)

#-----------------------------------------------------------------------
class Device(cir.Element):
    """
    State-Variable-Based Diode
    --------------------------

    Based on spice model. Connection diagram::

            o  1                           
            |                            
          --+--
           / \     
          '-+-'
            |                          
            o  0    	                  

    This model has better convergence properties. Externally it
    behaves exactly like the regular diode device. 

    Implementation includes depletion and diffusion charges. 

    Netlist examples::

        svdiode:d1 1 0 isat=10fA cj0=20fF

        # Electrothermal device
        svdiode_t:d2 2 3 1000 gnd cj0=10pF tt=1e-12 rs=100 bv = 4.

        # Model statement
        .model dmodel1 svdiode (cj0 = 10pF tt=1ps)

    Internal Topology
    +++++++++++++++++

    The internal representation is the following::

        0  o
           |
           \ 
           / Rs
           \ 
           / 
           |  t2                                 tx
           o---------+                  +----------------+
                     | i(x)+dq/dt       |                |
          +         /|\                /|\ gyr vin      /^\ 
        vin        | | |              | | |            | | | gyr v(x)
          -         \V/                \V/              \|/  
                     |                  |                |
        1  o---------+                  +--------+-------+
                                                 |
                                                --- tref
                                                 V

    Terminal t2 not present if Rs = 0

    Important Note
    ++++++++++++++

    This implementation does not account for the power dissipation
    in Rs. Use an external thermal resistor if that is needed.
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

    # Independent source attribute: output port
    # sourceOutput = (0, 1)

    # Define parameters (note most parameters defined in svjunction.py)
    paramDict = dict(
        cir.Element.tempItem,
        isat = ('Saturation current', 'A', float, 1e-14),
        n = ('Emission coefficient', ' ', float, 1.),
        fc = ('Coefficient for forward-bias depletion capacitance', ' ', 
              float, .5),
        cj0 = ('Zero-bias depletion capacitance', 'F', float, 0.),
        vj = ('Built-in junction potential', 'V', float, 1.),
        m = ('PN junction grading coefficient', ' ', float, .5),
        tt = ('Transit time', 's', float, 0.),
        xti = ('Is temperature exponent', ' ', float, 3.),
        eg0 = ('Energy bandgap', 'eV', float, 1.11),
        tnom = ('Nominal temperature', 'C', float, 27.),
        ibv = ('Current at reverse breakdown voltage', 'A', float, 1e-10),
        bv = ('Breakdown voltage', 'V', float, np.inf),
        area = ('Area multiplier', ' ', float, 1.),
        rs = ('Series resistance', 'Ohms', float, 0.),
        kf = ('Flicker noise coefficient', '', float, 0.),
        af = ('Flicker noise exponent', '', float, 1.),
       )

    def __init__(self, instanceName):
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed
        self.jtn = SVJunction()

    def process_params(self, thermal = False):
        # Make sure tape is re-generated
        ad.delete_tape(self)
        # Remove internal terminals
        self.clean_internal_terms()

        # Set flag to add thermal ports if needed
        self.__addThermalPorts = True

        # Define topology first
        # Needs at least one internal terminal: 
        tx = self.add_internal_term('x', 's.v.')
        tref = self.add_reference_term() 
        self.linearVCCS = [((0, 1), (tx, tref), glVar.gyr)]
        # Nonlinear device attributes
        self.csOutPorts = [(0, 1), (tref, tx)]
        self.noisePorts = [(0, 1)]
        self.controlPorts = [(tx, tref)]

        if self.rs:
            # Needs one more terminal
            t2 = self.add_internal_term('Vd_int', 'V')
            g = self.area / self.rs
            self.linearVCCS = [((t2, 1), (tx, tref), glVar.gyr),
                               ((0, t2), (0, t2), g)]
            # Nonlinear device outputs change
            self.csOutPorts = [(t2, 1), (tref, tx)]
            self.noisePorts = [(t2, 1), (0, t2)]

        # Nonlinear device attributes (defined in process_params())
        self.qsOutPorts = [ ]
        self._qd = False
        if self.tt or self.cj0:
            # Add charge source (otherwise the charge calculation is ignored)
            self.qsOutPorts.append(self.csOutPorts[0])
            self._qd = True

        # Absolute nominal temperature
        self.Tnomabs = self.tnom + const.T0
        self.egapn = self.eg0 - .000702 * (self.Tnomabs**2) \
            / (self.Tnomabs + 1108.)

        # Calculate variables in junction
        self.jtn.process_params(self.isat, self.n, self.fc, self.cj0, self.vj, 
                                self.m, self.xti, self.eg0, self.Tnomabs)
        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C
        """
        # Make sure tape is re-generated
        ad.delete_tape(self)
        # Absolute temperature
        self.Tabs = temp + const.T0
        # Thermal voltage
        self.vt = const.k * self.Tabs / const.q
        # Temperature-adjusted egap
        self.egap_t = self.eg0 - .000702 * (self.Tabs**2) / (self.Tabs + 1108.)
        self._Sthermal = 4. * const.k * self.Tabs * self.rs
        # Normalized temp
        self.tnratio = self.Tabs / self.Tnomabs
        # Everything else is handled by the PN junction
        self.jtn.set_temp_vars(self)

    #---------------------------------------------------------------
    # Nonlinear device methods (isNonlinear = True)
    #---------------------------------------------------------------
    def eval_cqs(self, vPort):
        """
        Models intrinsic diode and charge.  

        vPort[0]: x as in Rizolli's equations

        Returns a tuple with 2 or 3 elements: current, voltage and
        charge. Charge is ommited if both cj0 and tt are zero.
        """
        # Calculate state-variable PN junction current, voltage and charge
        (iD, vD) = self.jtn.get_idvd(vPort[0])
        if self.cj0:
            qD = self.jtn.get_qd(vD)
        if self.tt:
            qD += self.tt * iD

        # add breakdown current
        if (self.bv < np.inf):
            iD -= self.ibv * \
                ad.safe_exp(-(vD + self.bv) / self.n / self.vt)

        # area effect
        iD *= self.area
        # Scale voltage (gyrator gain)
        vD *= glVar.gyr

        iVec = np.array([iD, vD])
        if self._qd:
            # area effect for charge
            qD *= self.area
            qVec = np.array([qD])
            return (iVec, qVec)
        else:
            return (iVec, np.array([]))


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
            Sthermal = self._Sthermal,
            Sshot = 2. * const.q * outV[0],
            kFliker = self.kf * outV[0]
            )
        # Add capacitor
        if self._qd:
            self.OP['Cd'] = jac[-1,0] / jac[1,0]

        return self.OP
                   

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 
        """
        sj = self.OP['Sshot'] + self.OP['kSflicker'] / pow(f, self.af)
        if self.rs:
            sV = np.array([sj, self.OP['Sthermal']])
        else:
            sV = np.array([sj])
        return sV



    # Use automatic differentiation for eval and deriv functions
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval



