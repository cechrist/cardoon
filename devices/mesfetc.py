"""
:mod:`mesfetc` -- Intrinsic MESFET using Curtice-Ettemberg cubic model
----------------------------------------------------------------------

.. module:: mesfetc
.. moduleauthor:: Carlos Christoffersen


"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad


class Device(cir.Element):
    """
    Intrinsic MESFET using Curtice-Ettemberg cubic model
    ----------------------------------------------------

    Model derived from fREEDA 1.4 MesfetCT device by Carlos
    E. Christoffersen and Hector Gutierrez.

    Terminal order: 0 Drain, 1 Gate, 2 Source::

               Drain 0
                       o
                       |
                       |
                   |---+
                   |
      Gate 1 o---->|
                   |
                   |---+
                       |
                       |
                       o
              Source 2

    Netlist example::

        mesfetc:m1 2 3 4 a0=0.09910 a1=0.08541 a2=-0.02030 a3=-0.01543

    """
    # Device category
    category = "Semiconductor devices"

    devType = "mesfetc"
    paramDict = dict(
        cir.Element.tempItem,
        a0 = ('Drain saturation current for Vgs=0', 'A', float, 0.1),
        a1 = ('Coefficient for V1', 'A/V', float, 0.05),
        a2 = ('Coefficient for V1^2', 'A/V^2', float, 0.),
        a3 = ('Coefficient for V1^3', 'A/V^3', float, 0.),
        beta = ('V1 dependance on Vds', '1/V', float, 0.),
        vds0 = ('Vds at which BETA was measured', 'V', float, 4.),
        gama = ('Slope of drain characteristic in the linear region', '1/V', 
                float, 1.5),
        vt0 = ('Voltage at which the channel current is forced to be zero\
 for Vgs<=Vto', 'V', float, -1e+10),
        cgs0 = ('Gate-source Schottky barrier capacitance for Vgs=0', 'F', 
                float, 0.),
        cgd0 = ('Gate-drain Schottky barrier capacitance for Vgd=0', 'F', 
                float, 0.),
        isat = ('Diode saturation current', 'A', float, 0.),
        n = ('Diode ideality factor', '', float, 1.),
        ib0 = ('Breakdown current parameter', 'A', float, 0.),
        nr = ('Breakdown ideality factor', '', float, 10.),
        t = ('Channel transit time', 's', float, 0.),
        vbi = ('Built-in potential of the Schottky junctions', 'V', float, 0.8),
        fcc = ('Forward-bias depletion capacitance coefficient', 'V', 
               float, 0.5),
        vbd = ('Breakdown voltage', 'V', float, 1e+10),
        tnom = ('Reference Temperature', 'K', float, 293),
        avt0 = ('Pinch-off voltage (VP0 or VT0) linear temp. coefficient',
                '1/K', float, 0.),
        bvt0 = ('Pinch-off voltage (VP0 or VT0) quadratic temp. coefficient',
                '1/K^2', float, 0),
        tbet = ('BETA power law temperature coefficient', '1/K', float, 0),
        tm = ('Ids linear temp. coeff.', '1/K', float, 0),
        tme = ('Ids power law temp. coeff.', '1/K^2', float, 0),
        eg = ('Barrier height at 0 K', 'eV', float, 0.8),
        m = ('Grading coefficient', '', float, 0.5),
        xti = ('Diode saturation current temperature exponent', '', float, 2.),
        b = ('Thermal conductivity temperature exponent', '', float, 1.22),
        area = ('Area multiplier', '', float, 1.),
        )
    
    numTerms = 3

    # Create electrothermal device
    makeAutoThermal = True

    isNonlinear = True

    # For now only one current source between D and S
    csOutPorts = [(0, 2)]
    # Controling voltages are DB, GB and SB
    controlPorts = [(0, 3), (1, 3), (2, 3)]
    vPortGuess = np.array([0., 0., 0.])
    # No charge sources defined by now
    qsOutPorts = [ ]
    # No time-delayed port voltages required

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
        ad.delete_tape(self)

        if self.type == 'n':
            self._tf = 1.
        elif self.type == 'p':
            self._tf = -1.
        else:
            raise cir.CircuitError(
                '{0}: unrecognized type: {1}. Valid types are "n" or "p"'.format(self.nodeName, self.type))

        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        ad.delete_tape(self)

        # Absolute temperature (note self.temp is in deg. C)
        T = const.T0 + temp
        # Thermal voltage
        self.Vt = const.k * T / const.q
        

    def eval_cqs(self, vPort, saveOP = False):
        """
        Calculates drain current. Input is a vector as follows:
        vPort = [vdb , vgb , vsb]

        If saveOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # Invert all voltages in case of a P-channel device
        vPort1 = self._tf * vPort
        # The following formula (11.2.10) seems to work better
        # but it is just an approximation
        vp = pow(np.sqrt(vPort1[1]-self.vt0 \
                             + pow(np.sqrt(2.*self.phi)+\
                                       .5*self.gamma, 2))\
                     - .5*self.gamma, 2) - 2.*self.phi
        
        # Normalized currents.  
        i_f = inv_f((vp - vPort1[2]) / self.Vt)
        i_r = inv_f((vp - vPort1[0]) / self.Vt)
        
        # Calculate IS at this point
        # All this not needed if we just want IS
        eox = 34.5e-12
        cox = eox / self.tox
        mu0 = self.kp / cox
        mus = mu0 / (1. + self.theta * np.sqrt(vp + 2.*self.phi))
        chi = self.Vt * mus / self.l / self.vsat
        
        # From Page 464, also (A.2.3.2d)
        n = 1. + self.gamma / 2. / np.sqrt(2.*self.phi + vp)
        
        # Here kp = u0 Cox
        # 1 / (1 + theta*sqrt(vp+2*phi)) accounts for effective mobility
        IS = self.kp / (1. + self.theta*np.sqrt(vp+2.*self.phi))\
            * n * self.Vt * self.Vt / 2. * self.w/self.l
        
        # Get drain current (including vsat effect)
        # 1/(1 + chi * abs(sqrt(1+i_f)-sqrt(1+i_r))) is the velocity
        # saturation factor. It needs absolute value to work with
        # reverse biasing
        idrain = IS * (i_f - i_r) \
            / (1. + chi * np.abs(np.sqrt(1.+i_f)-np.sqrt(1.+i_r))) * self._tf

        iVec, qVec = np.array([idrain]), np.array([])

        if saveOP:
            # Create operating point variables vector
            opV = np.array([vp, i_f, i_r, IS, n])
            return (iVec, qVec, opV)
        else:
            # Return numpy array with one element per current source.
            return (iVec, qVec)


    # Use AD for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    get_op_vars = ad.get_op_vars
    
    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages and currents from eval_cqs()
        """
        vds = vPort[0] - vPort[2]
        # pout = vds*ids + vdb*idb + vsb*isb
        pout = vds*currV[0] 
        # Not available yet:  vPort[0] * currV[1] + vPort[2] * currV[2] 
        return pout

    def get_OP(self, vPort):
        """
        Calculates operating point information

        For now it is quite incomplete
        Input:  vPort = [vdb , vgb , vsb]
        Output: dictionary with OP variables
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)
        opV = self.get_op_vars(vPort)

        self.OP = dict(
            VD = vPort[0],
            VG = vPort[1],
            VS = vPort[2],
            IDS = outV[0],
            vp = opV[0],
            i_f = opV[1],
            i_r = opV[2],
            IS = opV[3],
            n = opV[4]
            )
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None




