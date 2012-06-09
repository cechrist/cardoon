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
    E. Christoffersen and Hector Gutierrez. ---> modify this to carrot
    model: 2 diode juntions plus a current source for ids.

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

    Internal Topology::

                   ,----------------,-------------,--------,--o 0 (D)
                   |                |             |        |
                  /^\               |             |        |
                 ( | ) igd(Vgd)   ----- Cgd       |        |
                  \|/             -----           |        |
                   |                |             |       /|\
        (G) 1 o----+----------------,      Cds  -----    ( | ) ids(Vgs, Vgd)
                   |                |           -----     \V/               
                  /|\               |             |        |
                 ( | ) igs(Vgs)   ----- Cgs       |        |
                  \V/             -----           |        |
                   |                |             |        |
                   `----------------'-------------'--------'--o 2 (S)

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

    # igs, igd, ids
    csOutPorts = [(1, 2), (1, 0), (0, 2)]
    # Controling voltages are Vgs, Vgd
    controlPorts = [(1, 2), (1, 0)]
    # Time-delayed control port added later. Guess includes time-delayed port
    vPortGuess = np.array([0., 0., 0.])
    # charge sources
    qsOutPorts = [(1, 2), (1, 0), (0, 2)]

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
        # Time-delayed control port
        self.delayedContPorts = [(1, 2, self.t)]

        self._k2 = self.cgs0 / sqrt(1. - self.fcc);
        self._k3 = self.cgd0 / sqrt(1. - self.fcc);
        self._bm1 = self.b - 1.;

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
        self._Vt = const.k * T / const.q
        delta_T = temp - self.tnom
        tn = T / (self.tnom + const.T0)
        self._k5 = self.n  * self._Vt
        self._k6 = self.nr * self._Vt
        self._Vt0 = self.vt0 * (1. + (delta_T * self.avt0) 
                                + (delta_T**2 * self.bvt0))
        if (self.tbet):
            self._Beta = self.beta * pow(1.01, (delta_T * self.tbet))
        else:
            self._Beta = self.beta
        Ebarr = self.eg -.000702 * T * T / (T + 1108.)
        EbarrN = eg -.000702 * tnom * tnom / (tnom + 1108.)
        Nn = 1. / 38.696 / self._Vt
        self._Is = self.isat * np.exp((tn - 1.) * Ebarr / Nn / Vt)
        if (self.xti):
            self._Is *= pow(tn, xti / Nn)
        self._Vbi = self.vbi * tn -3. * Vt * np.log(tn) + tn * EbarrN - Ebarr
        self._k1 = self.fcc * self._Vbi
        self._k4 = 2. * Vbi * (one - fcc)
        

    def eval_cqs(self, vPort, saveOP = False):
        """
        Calculates gate and drain current. Input is a vector as follows:
        vPort = [vgs(t), vgd(t), vgs(t-td)]

        If saveOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # static igs current
        igs = self._Is * (np.exp(vPort[0] / self._k5) - 1.) \
            - self.ib0 * np.exp(-(vPort[0] + self.vbd) / self._k6)
      
        # # Calculate cgs, including temperature effect.
        # condassign(cgs, k1 - vPort[0],
      	# cgs0 / sqrt(one - vPort[0] / Vbi), k2 * (one + (vPort[0] - k1) / k4))
        # cgs *= (one + m * (0.0004 * delta_T + one - Vbi / vbi))
      
        # static igd current
        igd = self._Is * (np.exp(vPort[1] / self._k5) - 1.) \
            - self.ib0 * np.exp(-(vPort[1] + self.vbd) / self._k6)
      
        # Power dissipated on this junction
        ip[2] -= igd * vPort[1]
      
        # Calculate cgd, including temperature effect.
        condassign(cgd, k1 - vPort[1],
      	cgd0 / sqrt(one - vPort[1] / Vbi), k3 * (one + (vPort[1] - k1) / k4))
        cgd *= (one + m * (0.0004 * delta_T + one - Vbi / vbi))
      
      	# Calculate the total current igd = static + dq_dt
        igd += cgd * x[4]
      
        # Calculate ids. Include temperature effects.
        vx = x[5] * (one + Beta * (vds0 - vp[1]))
      
        itmp = (a0 + vx*(a1 + vx*(a2 + vx  * a3)))* tanh(gama * vp[1])
        condassign(ids, (itmp * vp[1]) * (vPort[0] - Vt0), itmp, zero)
      
        if (tme && tm)
          ids *= pow((1 + delta_T * tm), tme)




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




