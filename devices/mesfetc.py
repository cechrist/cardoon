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
from diode import Junction

class Device(cir.Element):
    """
    Cubic Curtice-Ettemberg Intrinsic MESFET Model
    ----------------------------------------------

    Model derived from fREEDA 1.4 MesfetCT model adapted to re-use
    junction code from ``diode.py``. Some parameter names have been
    changed: ``isat``, ``tau``. Uses symmetric diodes and
    capacitances. Works in reversed mode.

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

                   ,----------------,------------,--o 0 (D)
                   |                |            |
                  /^\               |            |
                 ( | ) igd(Vgd)   ----- Cgd      |
                  \|/             -----          |
                   |                |           /|\ 
        (G) 1 o----+----------------,          ( | ) ids(Vgs, Vgd)
                   |                |           \V/               
                  /|\               |            |
                 ( | ) igs(Vgs)   ----- Cgs      |
                  \V/             -----          |
                   |                |            |
                   `----------------'------------'--o 2 (S)

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
 for Vgs<=Vto', 'V', float, -np.inf),
        cgs0 = ('Gate-source Schottky barrier capacitance for Vgs=0', 'F', 
                float, 0.),
        cgd0 = ('Gate-drain Schottky barrier capacitance for Vgd=0', 'F', 
                float, 0.),
        isat = ('Diode saturation current', 'A', float, 0.),
        n = ('Diode ideality factor', '', float, 1.),
        ib0 = ('Breakdown current parameter', 'A', float, 0.),
        nr = ('Breakdown ideality factor', '', float, 10.),
        vbd = ('Breakdown voltage', 'V', float, np.inf),
        tau = ('Channel transit time', 's', float, 0.),
        vbi = ('Built-in potential of the Schottky junctions', 'V', float, 0.8),
        fcc = ('Forward-bias depletion capacitance coefficient', 'V', 
               float, 0.5),
        tnom = ('Nominal temperature', 'C', float, 27.),
        avt0 = ('Pinch-off voltage (VP0 or VT0) linear temp. coefficient',
                '1/K', float, 0.),
        bvt0 = ('Pinch-off voltage (VP0 or VT0) quadratic temp. coefficient',
                '1/K^2', float, 0.),
        tbet = ('BETA power law temperature coefficient', '1/K', float, 0),
        tm = ('Ids linear temp. coeff.', '1/K', float, 0.),
        tme = ('Ids power law temp. coeff.', '1/K^2', float, 0.),
        eg0 = ('Barrier height at 0 K', 'eV', float, 0.8),
        mgs = ('Gate-source grading coefficient', '', float, 0.5),
        mgd = ('Gate-drain grading coefficient', '', float, 0.5),
        xti = ('Diode saturation current temperature exponent', '', float, 2.),
        area = ('Area multiplier', '', float, 1.),
        )
    
    numTerms = 3

    # Create electrothermal device
    makeAutoThermal = True

    isNonlinear = True
    nDelays = 2

    # igs, igd, ids
    csOutPorts = [(1, 2), (1, 0), (0, 2)]
    # Controling voltages are Vgs, Vgd
    controlPorts = [(1, 2), (1, 0)]
    # Time-delayed control port added later. Guess includes time-delayed port
    vPortGuess = np.array([0., 0., 0., 0.])
    # charge sources
    qsOutPorts = [(1, 2), (1, 0)]

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        self.diogs = Junction()
        self.diogd = Junction()

    def process_params(self, thermal = False):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.
        ad.delete_tape(self)
        # Time-delayed control port
        self.delayedContPorts = [(1, 2, self.tau), (1, 0, self.tau)]

        # Absolute nominal temperature
        self.Tnomabs = self.tnom + const.T0
        self.egapn = self.eg0 - .000702 * (self.Tnomabs**2) \
            / (self.Tnomabs + 1108.)

        # Calculate variables in junctions
        self.diogs.process_params(self.isat, self.n, self.fcc, 
                                  self.cgs0, self.vbi, 
                                  self.mgs, self.xti, self.eg0, self.Tnomabs)
        self.diogd.process_params(self.isat, self.n, self.fcc, 
                                  self.cgd0, self.vbi, 
                                  self.mgd, self.xti, self.eg0, self.Tnomabs)
        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        ad.delete_tape(self)
        # Absolute temperature (note self.temp is in deg. C)
        self.Tabs = const.T0 + temp
        # Thermal voltage
        self.vt = const.k * self.Tabs / const.q
        # Temperature-adjusted egap
        self.egap_t = self.eg0 - .000702 * (self.Tabs**2) / (self.Tabs + 1108.)
        # Everything else is handled by junctions
        self.diogs.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                                 self.egapn, self.egap_t)
        self.diogd.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                                 self.egapn, self.egap_t)

        delta_T = temp - self.tnom
        self._k6 = self.nr * self.vt
        self._Vt0 = self.vt0 * (1. + (delta_T * self.avt0) 
                                + (delta_T**2 * self.bvt0))
        if (self.tbet):
            self._Beta = self.beta * pow(1.01, (delta_T * self.tbet))
        else:
            self._Beta = self.beta

        if self.tme * self.tm != 0.:
            self._idsFac = pow((1. + delta_T * self.tm), self.tme)
        else:
            self._idsFac = 1.
        

    def eval_cqs(self, vPort, getOP = False):
        """
        Calculates gate and drain current. Input is a vector as follows:
        vPort = [vgs(t), vgd(t), vgs(t-tau), vgd(t-tau)]

        getOP has no effect for now
        """
        # Calculate junction currents
        igs = self.diogs.get_id(vPort[0])
        qgs = self.diogs.get_qd(vPort[0])
        igd = self.diogd.get_id(vPort[1])
        qgd = self.diogd.get_qd(vPort[1])

        # Add breakdown current
        igs -= self.ib0 * np.exp(-(vPort[0] + self.vbd) / self._k6)
        igd -= self.ib0 * np.exp(-(vPort[1] + self.vbd) / self._k6)

        DtoSswap = ad.condassign(vPort[0] - vPort[1], 1., -1.)
        vds = DtoSswap * (vPort[0] - vPort[1])
        vgsi =  ad.condassign(DtoSswap, vPort[0], vPort[1])
        vgsiT =  ad.condassign(DtoSswap, vPort[2], vPort[3])
        # Calculate ids. 
        vx = vgsiT * (1. + self._Beta * (self.vds0 - vds))
        ids = (self.a0 + vx * (self.a1 + vx * (self.a2 + vx  * self.a3))) \
            * np.tanh(self.gama * vds) * self._idsFac
        # vgsiT makes more sense than vgsi below? (vgsi in original doc)
        ids = ad.condassign((vgsi - self._Vt0), ids, 0.)
        # Must ensure ids > 0 for power conservation
        ids = ad.condassign(ids, ids, 0.)
      
        # Return numpy array with one element per current source.
        iVec = np.array([igs, igd, ids * DtoSswap]) * self.area
        qVec = np.array([qgs, qgd]) * self.area

        return (iVec, qVec)


    # Use AD for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    
    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages and currents from eval_cqs()
        """
        vds = vPort[0] - vPort[1]
        # pout = vds*ids + vgs*igs + vgd*igd
        pout = vds*currV[2] + vPort[0] * currV[0]  + vPort[1] * currV[1]
        return pout

    def get_OP(self, vPort):
        """
        Calculates operating point information

        Input:  vPort = [vgs , vgd , vgs]
        Output: dictionary with OP variables
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)

        opDict = dict(
            VGS = vPort[0],
            VDS = vPort[0] - vPort[1],
            IDS = outV[2],
            IGS = outV[0],
            IGD = outV[1]
            )
        return opDict

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None




