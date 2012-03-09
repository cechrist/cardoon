"""
:mod:`mosEKV` -- MOS EKV 2.6 intrinsic model
--------------------------------------------

.. module:: mosEKV
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    """
    Intrinsic EPFL EKV 2.6 MOSFET
    -----------------------------

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

    Mostly based on [1], but some updates from a later revision (dated
    1999) are also included.
    
    [1] The EPFL-EKV MOSFET Model Equations for Simulation, Technical
    Report, Model Version 2.6, June, 1997, Revision I, September,
    1997, Revision II, July, 1998, Bucher, Christophe Lallement,
    Christian Enz, Fabien Theodoloz, Francois Krummenacher,
    Electronics Laboratories, Swiss Federal Institute of Technology
    (EPFL), Lausanne, Switzerland
    
    This implementation includes accurate current interpolation
    function (optional), works for negative VDS and includes
    electrothermal model, DC operating point paramenters and noise
    equations.
    
    Code originally based on fREEDA 1.4 implementation
    <http://www.freeda.org>::
    
        // Element information
        ItemInfo Mosnekv::einfo =
        {
          "mosnekv",
          "EPFL EKV MOSFET model",
          "Wonhoon Jang",
          DEFAULT_ADDRESS"transistor>mosfet",
          "2003_05_15"
        };
    
    Parameter limit checking, simple capacitance calculations for
    operating point are not yet implemented.

    Netlist examples::

        mosekv:m1 2 3 4 gnd w=30e-6 l=1e-6 type = n ekvint=0

        # Electro-thermal version
        mosekv_t:m1 2 3 4 gnd 1000 gnd w=30e-6 l=1e-6 type = n

        # Model statement
        .model ekvn mosekv (type = n kp = 200u theta = 0.6)

    Internal Topology
    +++++++++++++++++

    The internal topology is the following::

             ,----------------------------+-------------+--o 0 (D)
             |                            |             |
            /|\                           |             |
           ( | ) idb (Vds > 0)          -----           |
            \V/                         ----- qd        |       
             |             1 (G)          |            /|\       
             |               o            |           ( | ) ids    
             |               |            |            \V/      
             |               |            |             |       
             |             -----          |             |
             |             ----- qg       |      qs     |
             |               |            |      ||     |
     (B) 3 o-+---------------+------------+------||-----+--o 2 (S)
                                                 ||

    The impact ionization current (idb) is normally added to the drain
    current, but if the device is in reverse (Vds < 0 for N-channel)
    mode, it is added to the source current.
    """

    devType = "mosekv"
    
    # Create electrothermal device
    makeAutoThermal = True

    numTerms = 4

    isNonlinear = True

    paramDict = dict(
        cir.Element.tempItem,
        type = ('N- or P-channel MOS (n or p)', '', str, 'n'),
        l = ('Gate length', 'm', float, 1e-06),
        w = ('Gate width', 'm', float, 1e-06),
        np = ('Parallel multiple device number', '', float, 1.),
        ns = ('Serial multiple device number', '', float, 1.),
        cox = ('Gate oxide capacitance per area', 'F/m^2', float, 0.7e-3),
        xj = ('Junction depth', 'm', float, 1e-07),
        dw = ('Channel width correction', 'm', float, 0.),
        dl = ('Channel length correction', 'm', float, 0.),
        vt0 = ('Long_channel threshold voltage', 'V', float, 0.5),
        gamma = ('Body effect parameter', 'V^1/2', float, 1.),
        phi = ('Bulk Fermi potential', 'V', float, 0.7),
        kp = ('Transconductance parameter', 'A/V^2', float, 50e-6),
        e0 = ('Mobility reduction coefficient', 'V/m', float, 1e12),
        ucrit = ('Longitudinal critical field', 'V/m', float, 2e6),
        tox = ('Oxide thickness', 'm', float, None),
        nsub = ('Channel doping', '1/cm^3', float, None),
        vfb = ('Flat-band voltage', 'V', float, None),
        u0 = ('Low-field mobility', 'cm^2/(V.s)', float, None),
        vmax = ('Saturation velocity', 'm/s', float, None),
        theta = ('Mobility recuction coefficient', '1/V', float, 0.),
        Lambda = ('Channel-length modulation', '', float, 0.5),
        weta = ('Narrow-channel effect coefficient', '', float, 0.25),
        leta = ('Short-channel effect coefficient', '', float, 0.1),
        q0 = ('Reverse short channel effect peak charge density', 'A.s/m^2', 
              float, 0.),
        lk = ('Reverse short channel effect characteristic length', 'm', 
              float, 2.9e-07),
        iba = ('First impact ionization coefficient', '1/m', float, 0.),
        ibb = ('Second impact ionization coefficient', 'V/m', float, 3e+08),
        ibn = ('Saturation voltage factor for impact ionization', '', 
               float, 1.),
        tcv = ('Threshold voltage temperature coefficient', 'V/K', 
               float, 0.001),
        bex = ('Mobility temperature exponent', '', float, -1.5),
        ucex = ('Longitudinal critical field temperature exponent', '', 
                float, 0.8),
        ibbt = ('Temperature coefficient for IBB', '1/K', float, 0.0009),
        avto = ('Area related threshold voltage mismatch parameter', 'Vm', 
                float, 0.),
        akp = ('Area related gain mismatch parameter', 'm', float, 0.),
        agamma = ('Area related body effect mismatch parameter', 'V^(1/2)m', 
                  float, 0.),
        kf = ('Flicker noise coefficient', '', float, 0.),
        af = ('Flicker noise exponent', '', float, 1.),
        satlim = ('Ratio defining the saturation limit if/ir', '', 
                  float, 54.5982),
# The following commented parameters not implemented
#        nqs = ('Non-Quasi-Static operation switch', '', float, 0.),
#        xqc = ('Charge/capacitance model selector', '', float, 0.4),
        tnom = ('Nominal temperature of model parameters', 'C', float, 27.),
        ekvint = ('Interpolation function (0: accurate, 1: simple)', '', 
                  int, 0)
        )

    # Ids, Idb, Isb
    csOutPorts = [(0, 2), (0, 3), (2, 3)]
    # Controling voltages are DB, GB and SB
    controlPorts = [(0, 3), (1, 3), (2, 3)]
    vPortGuess = np.array([0., 0., 0.])
    # One charge source connected to each D, G, S
    qsOutPorts = [(0, 3), (1, 3), (2, 3)]
    # Noise sources: one between drain and source
    noisePorts = [(0, 2)]
    # No time-delayed port voltages required

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        # If qox code is disabled, qox is zero by default
        self.qox = 0.

    def process_params(self, thermal = False):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.
        # import pdb; pdb.set_trace()

        # Delete AD tape (if any)
        ad.delete_tape(self)
        # Any value other than zero uses the simple function
        if self.ekvint:
            self.interp = f_simple
        else:
            self.interp = f_accurate
        # Set Constants
        self.Tref = 300.15
        self._Tn = self.tnom + const.T0
        vtTref = (const.k * self.Tref) / const.q
        self.egTref = 1.16 - 0.000702 * self.Tref * self.Tref \
            / (self.Tref + 1108.)
        vtTnom = (const.k * self._Tn) / const.q
        self.egTnom = 1.16 - .000702 * self._Tn * self._Tn / (self._Tn + 1108.)
        niTnom = 1.45e10 * (self._Tn / self.Tref) \
            * np.exp(self.egTref / (2. * vtTref) - self.egTnom / (2. * vtTnom))
        # For N/P channel
        self._tcv = np.abs(self.tcv)
        if self.type == 'n':
            self.eta = 0.5
            self._tf = 1.
        elif self.type == 'p':
            self.eta = 1./3.
            self._tf = -1.
        else:
            raise cir.CircuitError(
                '{0}: unrecognized type: {1}. Valid types are "n" or "p"'.format(self.nodeName, self.type))
        #---------------------------------------------------------------
        # Calculate any missing parameters from user-defined settings
        # COX
        if (not self.is_set('tox')) and self.tox:
            self.cox = const.epOx / self.tox
        # GAMMA
        if (not self.is_set('gamma')) and self.nsub:
            self.gamma = np.sqrt(2. * const.q * const.epSi 
                                 * self.nsub * 1.e6) / self.cox
        # PHI
        if (not self.is_set('phi')) and self.nsub:
            self.phi = 2. * vtTnom * np.log(self.nsub / niTnom)
        # VT0: if specified, must be with the correct sign, otherwise
        # we have to adjust for P-channel
        if (not self.is_set('vt0')):
            if self.vfb:
                self.vt0 = self._tf * (self.vfb + self.phi 
                                       + self.gamma * np.sqrt(self.phi))
            else:
                self.vt0 *= self._tf 
        # Make sure vt0 is positive for calculations
        self._vt0 = abs(self.vt0)
        # KP
        if (not self.is_set('kp')) and self.u0:
            self.kp = self.u0 * 1.e-4 * self.cox    # /*(m^2/cm^2)*/
        # UCRIT
        if (not self.is_set('ucrit')) and (self.vmax > 0.) and (self.u0 > 0.):
            self.ucrit = self.vmax / (self.u0 * 1.e-4)
        # E0: no need for anything since theta != 0 triggers simple
        # mobility model

        # Initialize more variables
        self._weff = self.w + self.dw
        self._leff = self.l + self.dl
        # sqrt gate area
        self._sga = np.sqrt(self.np * self._weff * self.ns * self._leff)
        self._gammaa = self.gamma + self.agamma / self._sga
        # Clip to zero if negative
        if self._gammaa < 0.:
            self._gammaa = 0.
        cEps = 0.001936  #  was: 4. * pow(22.e-3, 2)
        cA = 0.028
        xi  = cA * (10. * self._leff / self.lk - 1.)
        self._deltavRSCE = 2. * self.q0 / \
            (self.cox * pow(1. + 0.5 * (xi + np.sqrt(xi*xi + cEps)), 2))

        # constants used in eval_c1qs()
        self._lc = np.sqrt(const.epSi * self.xj / self.cox)
        self._lmin = self.ns * self._leff / 10.
        self._Cox = self.cox * self.np * self._weff * self.ns * self._leff

        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        # Delete AD tape (if any)
        ad.delete_tape(self)
        # Absolute temperature (note self.temp is in deg. C)
        T = const.T0 + temp
        # Thermal voltage
        self._Vt = const.k * T / const.q
        # threshold voltage
        vt0T = self._vt0 - self._tcv * (T - self._Tn)
        self._vt0a = vt0T + self.avto / self._sga
        kpT = self.kp * pow(T / self._Tn, self.bex)
        self.kpa = kpT * (1 + self.akp / self._sga)
        # Clip to zero if negative
        self.kpa = ad.condassign(self.kpa, self.kpa, 0.)
        self._t_ucrit = self.ucrit * pow(T / self._Tn, self.ucex)
        # Energy gap
        egT = 1.16 - 0.000702 * T * T / (T + 1108.)		
        self.phiT = self.phi * T / self.Tref \
            - 3. * self._Vt * np.log(T / self.Tref) \
            - self.egTref * T / self.Tref + egT
        self._t_ibb = self.ibb * (1. + self.ibbt * (T - self.Tref))

        self._vc = self._t_ucrit * self._leff * self.ns
        self._qb0 = self._gammaa * np.sqrt(self.phiT)
        # Noise variables
        self._kSt = 4. * const.k * T

        # ----------------- self.qox calculation ------------------------
        #
        # Disabled since this is not strictly necessary if the
        # analysis code is written carefully. In particular, it should
        # not assume that all charges are initally set to zero.
        # 
        # Calculate gate charge with VG=VD=VS=0  (-qox)
        # # Effective gate voltage including reverse short channel effect
        # vgprime = - self._vt0a - self._deltavRSCE \
        #     + self.phiT + self._gammaa * np.sqrt(self.phiT)
        # # Effective substrate factor including charge-sharing for
        # # short and narrow channels.  Pinch-off voltage for
        # # narrow-channel effect
        # vp0 = vgprime - self.phiT - \
        #     self._gammaa * (np.sqrt(vgprime + self._gammaa*self._gammaa / 4.) 
        #                     - self._gammaa / 2.)
        # # Effective substrate factor accounting for charge-sharing
        # tmp = 16. * self._Vt * self._Vt
        # vxprime = 0.5 * (self.phiT + self.phiT + tmp)
        # tmp = self.leta / self._leff * 2. * np.sqrt(vxprime)
        # tmp -= 3. * self.weta * np.sqrt(vp0 + self.phiT) / self._weff
        # gammao = self._gammaa - const.epSi / self.cox * tmp
        # gammaprime = 0.5 * (gammao + np.sqrt(gammao * gammao + 0.1 * self._Vt))
        # # Pinch-off voltage including short- and narrow-channel effect
        # vp = vgprime - self.phiT 
        # vp -= gammaprime * (np.sqrt(vgprime + pow(.5 * gammaprime, 2)) 
        #                     - gammaprime / 2.)
        # sqvpphi = np.sqrt(vp + self.phiT + 1.e-6)
        # # 'oxide' charge
        # self.qox = self._gammaa * sqvpphi / self._Vt 


    def eval_cqs(self, vPort, saveOP = False):
        """
        Calculates Ids, Idb, Isb currents and D, G, S, charges. 

        Input:  vPort = [vdb , vgb , vsb]

        Output: vector with Ids, Idb, Isb currents and vector with D,
        G, S charges.
        
        If saveOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # Invert all voltages in case of a P-channel device
        vPort1 = self._tf * vPort
        # If vds is negative, swap Vd and Vs and (also swap charges
        # later)
        DtoSswap = ad.condassign(vPort1[0] - vPort1[2], 1., -1.)
        # perform the swap (need tmp variable)
        tmp = ad.condassign(DtoSswap, vPort1[0], vPort1[2])
        vPort1[2] = ad.condassign(DtoSswap, vPort1[2], vPort1[0])
        vPort1[0] = tmp

        # Effective gate voltage including reverse short channel effect
        vgprime =  vPort1[1] - self._vt0a - self._deltavRSCE \
            + self.phiT + self._gammaa * np.sqrt(self.phiT)
        # Effective substrate factor including charge-sharing for
        # short and narrow channels.  Pinch-off voltage for
        # narrow-channel effect
        vp0 = vgprime - self.phiT - \
            self._gammaa * (np.sqrt(vgprime + self._gammaa*self._gammaa / 4.) 
                            - self._gammaa / 2.)
        vp0 = ad.condassign(vgprime, vp0, -self.phiT)
        
        # Effective substrate factor accounting for charge-sharing
        tmp = 16. * self._Vt * self._Vt
        vsprime = 0.5 * (vPort1[2] + self.phiT + 
                         np.sqrt(pow(vPort1[2] + self.phiT, 2) + tmp))
        vdprime = 0.5 * (vPort1[0] + self.phiT + 
                         np.sqrt(pow(vPort1[0] + self.phiT, 2) + tmp))
        
        tmp = self.leta / self._leff * (np.sqrt(vsprime) + np.sqrt(vdprime))
        tmp -= 3. * self.weta * np.sqrt(vp0 + self.phiT) / self._weff
        gammao = self._gammaa - const.epSi / self.cox * tmp

        gammaprime = 0.5 * (gammao + np.sqrt(gammao * gammao + 0.1 * self._Vt))

        # Pinch-off voltage including short- and narrow-channel effect
        vp = vgprime - self.phiT 
        vp -= gammaprime * (np.sqrt(vgprime + pow(.5 * gammaprime, 2)) 
                            - gammaprime / 2.)
        vp = ad.condassign(vgprime, vp, -self.phiT)
        
        # Slope factor
        n = 1 + self._gammaa / (2. * np.sqrt(vp + self.phiT + 4. * self._Vt))
        
        # Forward normalized current
        i_f = (vp - vPort1[2]) / self._Vt
        i_f = self.interp(i_f)

        # Velocity saturation voltage
        vdss = np.sqrt(0.25 + self._Vt * np.sqrt(i_f) / self._vc) - 0.5
        vdss *= self._vc
        
        # Drain-to-source saturation voltage for reverse normalized current
        tmp = np.sqrt(i_f) - 0.75 * np.log(i_f)
        vdssprime = np.sqrt(0.25 + self._Vt * tmp / self._vc) - 0.5 
        vdssprime *= self._vc
        vdssprime += self._Vt * (np.log(self._vc / (2. * self._Vt)) - 0.6)

        # Channel-length modulation
        tmp = np.sqrt(i_f) - vdss / self._Vt
        deltav = np.sqrt(self.Lambda * tmp + 1./64)
        deltav *= 4. * self._Vt
        vds = .5 * (vPort1[0] - vPort1[2]) 
        vip = np.sqrt(vdss*vdss + deltav*deltav) 
        vip -= np.sqrt(pow(vds - vdss, 2) + deltav * deltav)
        deltal = np.log(1. + (vds - vip) / self._lc / self._t_ucrit)
        deltal *= self.Lambda * self._lc
        
        # Equivalent channel length including channel-length
        # modulation and velocity saturation
        lprime = self.ns * self._leff - deltal + (vds + vip) / self._t_ucrit
        leq = 0.5 * (lprime + np.sqrt(lprime*lprime + self._lmin*self._lmin))
        
        # Reverse normalized current
        tmp = vp - vds - vPort1[2] 
        tmp -= np.sqrt(vdssprime*vdssprime + deltav*deltav) 
        tmp += np.sqrt(pow(vds-vdssprime, 2) + deltav*deltav)
        irprime = self.interp(tmp / self._Vt)

        # Reverse normalized currect for mobility model, intrinsic
        # charges/capacitances, thermal noise model and NQS time-constant ???
        ir = self.interp((vp -  vPort1[0]) / self._Vt)
        
        # Quasi-static model equations
        # Dynamic model for the intrinsic node charges
        sqvpphi = np.sqrt(vp + self.phiT + 1.e-6)
        nq = 1. + self._gammaa / 2. / sqvpphi
        
        # Normalized intrinsic node charges
        xf = np.sqrt(0.25 + i_f)
        xr = np.sqrt(0.25 + ir)
        tmp = pow(xf + xr, 2)

        qd = 3. * xr**3 + 6. * xr*xr*xf + 4. * xr*xf*xf + 2. * xf**3
        qd *= 4. / 15. / tmp
        qd = - nq * (qd - .5)

        qs = 3. * xf**3 + 6. * xf*xf*xr + 4. * xf*xr*xr + 2. * xr**3
        qs *= 4. / 15. / tmp
        qs = - nq * (qs - .5)

        qi = qs + qd

        qb1 = - self._gammaa * sqvpphi / self._Vt 
        qb1 -= (nq - 1.) * qi / nq
        qb = ad.condassign(vgprime, qb1, -vgprime / self._Vt)
        # qg = -qi - qox - qb, but qox == 0, so:
        qg = -qi - qb - self.qox
        
        # Transconductance factor and mobility reduction due to vertical field
        betao = self.kpa * self.np * self._weff / leq 
        if self.theta:
            # Simple mobility reduction model
            vpprime = 0.5 * (vp + np.sqrt(vp * vp + 2. * self._Vt * self._Vt))
            beta = betao / (1. + theta * vpprime)
        else:
            # Rigorous mobility reduction model
            betaoprime = 1. + self.cox * self._qb0 / self.e0 / const.epSi
            betaoprime *= betao
            tmp = np.abs(qb + self.eta * qi)
            tmp = 1. + self.cox * self._Vt * tmp / self.e0 / const.epSi
            beta = betaoprime / tmp

        # Specific current
        IS = 2. * n * beta * self._Vt * self._Vt
        
        # Drain-to-source current
        ids = IS * (i_f - irprime)

        # import pdb; pdb.set_trace()        
        
        # Impact ionization current
        vib =  vPort1[0] - vPort1[2] - 2. * self.ibn * vdss
        idb1 = ids * self.iba * vib / self._t_ibb
        idb1 *= ad.safe_exp(-self._t_ibb * self._lc / vib) 
        idb = ad.condassign(vib, idb1, 0.)
        
        # -------------------------------------------------------------
        # Create output vectors
        qVec = np.array([0., qg, 0.], dtype=type(idb))
        # have to switch charges if Drain and Source voltages switched
        qVec[0] = ad.condassign(DtoSswap, qd, qs)
        qVec[2] = ad.condassign(DtoSswap, qs, qd)
        # De-normalize charge and invert if needed
        qVec *= self._Cox * self._Vt * self._tf

        iVec = np.array([0., 0., 0.], dtype=type(idb))
        iVec[0] = DtoSswap * ids
        iVec[1] = ad.condassign(DtoSswap, idb, 0.)
        iVec[2] = ad.condassign(DtoSswap, 0., idb)
        # Revert currents if needed
        iVec *= self._tf
        
        #--------------------------------------------------------------
        # Operating point information
        if saveOP:
            # Vth
            Vth = self._vt0a + self._deltavRSCE \
                + gammaprime * np.sqrt(vsprime) \
                - self._gammaa * np.sqrt(self.phiT)
            # Non quasi-static equations
            tau0 = .5 * self._Cox / self._Vt / beta
            tmp = (xf**2 + 3. * xf*xr + xr**2) / pow(xf + xr, 3)
            tau = tau0 * 4. / 15. * tmp 
            # Create operating point variables vector
            opV = np.array([vp, n, beta, IS, i_f, ir, irprime, 
                            Vth, tau0, tau, qi, DtoSswap])
            return (iVec, qVec, opV)
        else:
            return (iVec, qVec)


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

        (Should add calculation of simple capacitive model)
        Input:  vPort = [vdb , vgb , vsb]
        Output: dictionary with OP variables
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)
        opV = self.get_op_vars(vPort)

        # Check things that change if the transistor is reversed
        if opV[11] > 0.:
            reversed = False
            gds = jac[0,0]
        else:
            reversed = True
            gds = jac[0,2]
            
        # Use negative index for charges as power may be inserted in
        # between currents and charges by electrothermal model
        self.OP = dict(
            VD = vPort[0],
            VG = vPort[1],
            VS = vPort[2],
            IDS = outV[0],
            IDB = outV[1],
            ISB = outV[2],
            QD = outV[-3],
            QG = outV[-2],
            QS = outV[-1],
            Vp = opV[0],
            n = opV[1],
            Beta = opV[2],
            IS = opV[3],
            IF = opV[4],
            IR = opV[5],
            IRprime = opV[6],
            tef = 1. / (np.sqrt(.25 + opV[4]) + .5),
            Vth = opV[7],
            Vov = self._tf * opV[1] * (opV[0] - vPort[2]),
            Vdsat = self._tf * self._Vt * (2. * np.sqrt(opV[4]) + 4.),
            gm = jac[0,1], 
            gmbs = - jac[0,2] - jac[0,1] - jac[0,0],
            gds = gds,
            tau0 = opV[8],
            tau = opV[9],
            Sthermal = self._kSt * opV[2] * np.abs(opV[10]),
            kSfliker = self.kf * pow(jac[0,1], 2) / self._Cox,
            Reversed = reversed
            )
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 
        """
        s = self.OP['Sthermal'] + self.OP['kSflicker'] / pow(f, self.af)
        return np.array([s])


    # Use AD for eval and deriv function
    eval = ad.eval
    eval_and_deriv = ad.eval_and_deriv
    get_op_vars = ad.get_op_vars
    
#-------------------------------------------------------------------------
# Helper functions

def f_simple(v):
    """
    Simple interpolation function for F(v)

    Not accurate for moderate inversion
    """
    # Have to treat the function for large negative v specially
    # otherwise exp(.5*v) << 1 and we get log(1) = 0
    b = v + 20.
    d = np.exp(.5 * v)
    c = np.log(1. + d)
    # f = (b > -20) ? c : d
    f = ad.condassign(b, c, d)
    return f*f


def f_accurate(v):
    """
    Calculate a more accurate value of F(v) 

    Eq. (41) and (42) in [1]
    
    Refine f_simple with a few Newton iterations. This function
    performs the Newton iterations even if not needed on purpose to
    store the operations in the AD tape.
    """
    # get initial estimate
    i = f_simple(v)
    # Apply 3 Newton iterations to refine approximation
    for j in range(3):
        k1 = np.sqrt(0.25 + i)
        f = v - (2.*k1 - 1. + np.log(k1 - .5))
        vprime = -1./(2.*k1*(0.5 - k1)) + 1./k1
        i += f / vprime
    b = v + 20.
    # If b is very negative we need this to avoid the singularity that
    # happens when i is very small: k1 = 0.5 and the log goes to infinity
    i = ad.condassign(b, i, np.exp(v))
    #print v, i, np.exp(v)
    return i
