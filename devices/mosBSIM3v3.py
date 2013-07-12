"""
:mod:`mosBSIM3v3` -- Intrinsic BSIM3v3 MOSFET model
---------------------------------------------------

.. module:: mosBSIM3v3
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad
import mosExt

# KboQ = 8.617087e-5  #/* Kb / q  where q = 1.60219e-19 */
# EPSOX = 3.453133e-11
# EPSSI = 1.03594e-10 #this is in meters --> 1.03594e-8 F/cm
MAX_EXP = 5.834617425e14
MIN_EXP = 1.713908431e-15
MIN_1 = MIN_EXP * (1. + 2. * MIN_EXP)
EXP_THRESHOLD = 34.0

def log1pexp(x):
    """
    Safe calculation of log(1 + exp(x))
    """
    y1 = ad.condassign(x - EXP_THRESHOLD,
                       x,
                       np.log(1.+np.exp(x)))
    return ad.condassign(-x -  EXP_THRESHOLD,
                          np.exp(x),
                          y1)

class BSIM3(cir.Element):
    """
    Intrinsic BSIM3 MOSFET Model (version 3.2.4)
    --------------------------------------------
    
    This model mainly converted from fREEDA 2.0 mosnbsim3 model
    written by Ramya Mohan (http://www.freeda.org/) with some
    improvements. Also includes some code taken from ngspice
    (http://ngspice.sourceforge.net/) and pyEDA EDA Framework
    (https://github.com/cogenda/pyEDA).  *Results are reasonable but
    requires more testing*

    Default parameters listed for NMOS type. Default values for some
    parameters such as u0 and vth0 are different for PMOS type.

    Notes:

       * Most parameters are not checked for valid values

       * According to ngspice documentation, temperature specification
         is not functional (probably the same applies here)

       * Parameter descriptions need reviewing

       * The code to internally calculate k1 and k2 is disabled by
         default because using default values seems to give more
         reasonable results (use ``k1enable`` to enable).

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

    Netlist examples::

        bsim3_i:m1 2 3 4 gnd w=10e-6 l=1e-6 type = n 
        bsim3_i:m2 4 5 6 6 w=30e-6 l=1e-6 type = p 

    Internal topology
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


    """
    # Device category
    category = "Semiconductor devices"

    devType = "bsim3_i"
    paramDict = dict(
        cir.Element.tempItem,
        type = ('N- or P-channel MOS (n or p)', '', str, 'n'),
        k1enable = ('Enable k1, k2 internal calculation', '', bool, False),
        vth0 = (
            'Threshold voltage of long channel device at Vbs=0 and small Vds',
            'V', float, 0.7),
        l = ('Length', 'm', float, 1e-06),
        w = ('Width', 'm', float, 1e-06),
        tox = ('Gate oxide thickness', 'm', float, 1.5e-08),
        toxm = ('Gate oxide thickness used in extraction', '', float, 1.5e-08),
        cdsc = ('Drain/Source and channel coupling capacitance', 'F/m^2', 
                float, 0.00024),
        cdscb = ('Body-bias dependence of cdsc', 'F/V/m^2', float, 0),
        cdscd = ('Drain-bias dependence of cdsc', 'F/V/m^2', float, 0),
        cit = ('Interface state capacitance', '', float, 0),
        nfactor = ('Subthreshold swing coefficient', '', float, 1),
        xj = ('Junction depth', 'm', float, 1.5e-07),
        vsat = ('Saturationvelocity at tnom', 'm/s', float, 80000),
        at = ('Temperature coefficient of vsat', 'm/s', float, 33000),
        a0 = ('Non-uniform depletion width effect coefficient', '', float, 1),
        ags = ('Gate bias coefficient of Abulk', '', float, 0),
        a1 = ('Non-saturation effect coefficient', '', float, 0),
        a2 = ('Non-saturation effect coefficient', '', float, 1),
        keta = ('Body-bias coefficient of non-uniform depletion width effect', 
                '', float, -0.047),
        nsub = ('Substrate doping concentration', 'cm^{-3}', float, 6e+16),
        nch = ('Channel doping concentration', 'cm^{-3}', float, 1.7e+17),
        ngate = ('Poly-gate doping concentration', 'cm^{-3}', float, 0),
        vbm = ('Maximum body voltage', 'V', float, -3),
        xt = ('Doping depth', 'm', float, 1.55e-07),
        kt1 = ('Temperature coefficient of Vth', 'V', float, -0.11),
        kt1l = ('Temperature coefficient of Vth', 'V m', float, 0),
        kt2 = ('Body-coefficient of kt1', '', float, 0.022),
        k3 = ('Narrow width effect coefficient', '', float, 80),
        k3b = ('Body effect coefficient of k3', '', float, 0),
        w0 = ('Narrow width effect parameter', 'm', float, 2.5e-06),
        nlx = ('Lateral non-uniform doping effect', 'm', float, 1.74e-07),
        dvt0 = ('Short channel effect coefficient 0', '', float, 2.2),
        dvt1 = ('Short channel effect coefficient 1', '', float, 0.53),
        dvt2 = ('Short channel effect coefficient 2', 'V^{-1}', float, -0.032),
        dvt0w = ('Narrow width effect coefficient 0', 'm^{-1}', float, 0),
        dvt1w = ('Narrow width effect coefficient 1', 'm^{-1}', float, 5.3e+06),
        dvt2w = ('Narrow width effect coefficient 2', 'V^{-1}', float, -0.032),
        drout = ('DIBL coefficient of output resistance', '', float, 0.56),
        dsub = ('DIBL coefficient in the subthreshold region', '', float, 0.56),
        ua = ('Linear gate dependence of mobility', 'm/V', float, 2.25e-09),
        ub = ('Quadratic gate dependence of mobility', '(m/V)^2', 
              float, 5.87e-19),
        uc = ('Body-bias dependence of mobility', 'm/V^2', float, -4.65e-11),
        u0 = ('Low-field mobility at Tnom', 'cm^2/V/s', float, 670),
        voff = ('Threshold voltage offset', 'V', float, -0.08),
        tnom = ('Nominal temperature', 'C', float, 27.),
        elm = ('Non-quasi-static Elmore Constant Parameter', '', float, 5),
        delta = ('Effective Vds parameter', 'V', float, 0.01),
        rdsw = ('Sorce-drain resistance per width', '', float, 0),
        prwg = ('Gate-bias effect on parasitic resistance', '', float, 0),
        prwb = ('Body-effect on parasitic resistance', '', float, 0),
        prt = ('Temperature coefficient of parasitic resistance', '', float, 0),
        eta0 = ('Subthreshold region DIBL coefficeint', '', float, 0.08),
        etab = ('Subthreshold region DIBL coefficeint', '', float, -0.07),
        pclm = ('Channel length modulation coefficient', '', float, 1.3),
        pdibl1 = ('Drain-induced barrier lowering oefficient', '', float, 0.39),
        pdibl2 = ('Drain-induced barrier lowering oefficient', '', 
                  float, 0.0086),
        pdiblb = ('Body-effect on drain induced barrier lowering', '', 
                  float, 0),
        pscbe1 = ('Substrate current body-effect coeffiecient 1', 'V/m', 
                  float, 4.24e+08),
        pscbe2 = ('Substrate current body-effect coeffiecient 2', 'V/m', 
                  float, 1e-05),
        pvag = ('Gate dependence of output resistance parameter', '', float, 0),
        vfb = ('Flat band voltage', 'V', float, -1),
        acde = ('Exponential coefficient for finite charge thickness', '', 
                float, 1),
        moin = ('Coefficient for gate-bias dependent surface potential', '', 
                float, 15),
        noff = ('C-V turn-on/off parameter', '', float, 1),
        voffcv = ('C-V lateral shift parameter', '', float, 0),
        lint = ('Length reduction parameter', 'm', float, 0),
        ll = ('Length reduction parameter', '', float, 0),
        llc = ('Length reduction parameter for CV', '', float, 0),
        lln = ('Length reduction parameter', '', float, 1),
        lw = ('Length reduction parameter', '', float, 0),
        lwc = ('Length reduction parameter for CV', '', float, 0),
        lwn = ('Length reduction parameter', '', float, 1),
        lwl = ('Length reduction parameter', '', float, 0),
        lwlc = ('Length reduction parameter for CV', '', float, 0),
        wr = ('Width dependence of rds', '', float, 1),
        wint = ('Width reduction parameter', 'm', float, 0),
        dwg = ('Width reduction parameter', 'm/V', float, 0),
        dwb = ('Width reduction parameter', 'm/V', float, 0),
        wl = ('Width reduction parameter', '', float, 0),
        wlc = ('Width reduction parameter for CV', '', float, 0),
        wln = ('Width reduction parameter', '', float, 1),
        ww = ('Width reduction parameter', '', float, 0),
        wwc = ('Width reduction parameter for CV', '', float, 0),
        wwn = ('Width reduction parameter', '', float, 1),
        wwl = ('Width reduction parameter', '', float, 0),
        wwlc = ('Width reduction parameter for CV', '', float, 0),
        b0 = ('Abulk narrow width parameter', '', float, 0),
        b1 = ('Abulk narrow width parameter', '', float, 0),
        clc = ('Vdsat paramater for C-V model', '', float, 1e-07),
        cle = ('Vdsat paramater for C-V model', '', float, 0.6),
        alpha0 = ('Substrate current model parameter', 'm/V', float, 0),
        alpha1 = ('Substrate current model parameter', 'V^{-1}', float, 0),
        beta0 = ('Diode limiting current', 'V', float, 30),
        ute = ('Temperature coefficient of mobility', '', float, -1.5),
        k1 = ('First order body effect coefficient', 'V^{0.5}', float, 0.53),
        k2 = ('Second order body effect coefficient', '', float, -0.0186),
        ua1 = ('Temperature coefficient for ua', 'm/V', float, 4.31e-09),
        ub1 = ('Temperature coefficient for ub', '(m/V)^2', float, -7.61e-18),
        uc1 = ('Temperature coefficient for uc', 'm/V^2', float, -5.6e-11),
        )
    
    numTerms = 4

    # Create electrothermal device (not for now as temperature
    # dependence not reliable)
    # makeAutoThermal = True

    isNonlinear = True

    # Ids, Idb, Isb
    csOutPorts = [(0, 2), (0, 3), (2, 3)]
    # Controling voltages are DB, GB and SB
    controlPorts = [(0, 3), (1, 3), (2, 3)]
    vPortGuess = np.array([0., 0., 0.])
    # One charge source connected to each D, G, S
    qsOutPorts = [(0, 3), (1, 3), (2, 3)]

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
            # Change parameter default values
            if not self.is_set('u0'):
                self.u0 = 250
        else:
            raise cir.CircuitError(
                '{0}: unrecognized type: {1}. Valid types are "n" or "p"'.format(self.instanceName, self.type))

        # Nominal abs temperature
        self._Tn = const.T0 + self.tnom
        # Nominal Thermal voltage
        self._Vtn = const.k * self._Tn / const.q

        self.factor1 = np.sqrt(const.epSi / const.epOx * self.tox)
        Eg0 = 1.16 - 7.02e-4 * (self._Tn**2) / (self._Tn + 1108.0)
        ni = 1.45e10 * (self._Tn / 300.15) * np.sqrt(self._Tn / 300.15) \
            * np.exp(21.5565981 - Eg0 / (2. * self._Vtn))
        #self.esi = 11.7 * const.epsilon0 (replaced by const.epSi)
        self.ldrn = self.l
        self.wdrn = self.w
        
        t0 = pow(self.ldrn, self.lln)
        t1 = pow(self.wdrn, self.lwn)
        tmp1 = self.ll / t0 + self.lw / t1 + self.lwl / (t0 * t1)
        self.dl = self.lint + tmp1
        #tmp2 = llc / t0 + lwc / t1 + lwlc / (t0 * t1)
        #self.dlc = dlc + tmp2  # ???
        
        t2 = pow(self.ldrn, self.wln)
        t3 = pow(self.wdrn, self.wwn)
        tmp3 = self.wl / t2 + self.ww / t3 + self.wwl / (t2 * t3)
        self.dw = self.wint + tmp3
        #tmp4 = wlc / t2 + wwc / t3 + wwlc / (t2 * t3)
        #self.dwc = dwc + tmp4 # ???
        
        self.leff = self.l - 2.0 * self.dl
        self.weff = self.w - 2.0 * self.dw

        self.AbulkCVfactor = (1. + pow(self.clc/self.leff, self.cle))

        #self.leffCV = l - 2.0 * dlc
        #self.t11 = leffCV * leffCV

        # was epOx = 3.453133e-11
        self.cox = const.epOx / self.tox

        self.phi = 2. * self._Vtn * np.log(self.nch / ni)
        self.sqrtPhi = np.sqrt(self.phi)
        self.phis3 = self.sqrtPhi * self.phi

        self.Xdep0 = np.sqrt(2. * const.epSi / 
                          (const.q * self.nch * 1e6)) * self.sqrtPhi
        self.litl = np.sqrt(3. * self.xj * self.tox)
        self.vbi = self._Vtn * np.log(1.0e20 * self.nch / (ni**2))
        self.cdep0 = np.sqrt(const.q * const.epSi * self.nch * 1e6 
                          / 2. / self.phi)

        if not self.is_set('toxm'):
            self.toxm = self.tox
        if not self.is_set('dsub'):
            self.dsub = self.drout

        self.ldeb = np.sqrt(const.epSi * self._Vtn
                           / (const.q * self.nch * 1e6)) / 3.
        
        #import pdb; pdb.set_trace()
        if (self.k1enable != 0.) and \
                not (self.is_set('k1') or self.is_set('k2')):
            vbx = self.phi - 7.7348e-4  * self.nch * self.xt**2
            # From ngspice
            vbx = -abs(vbx)
            Vbm = -abs(self.vbm)
            gamma1 = 5.753e-12 * np.sqrt(self.nch) / self.cox
            gamma2 = 5.753e-12 * np.sqrt(self.nsub) / self.cox
            T0 = gamma1 - gamma2
            T1 = np.sqrt(self.phi - vbx) - self.sqrtPhi
            T2 = np.sqrt(self.phi * (self.phi - Vbm)) - self.phi
            self._k2 = T0 * T1 / (2. * T2 + Vbm)
            self._k1 = gamma2 - 2. * self._k2 * np.sqrt(self.phi - Vbm)
            # print self._k1, self._k2
        else:
            self._k1 = self.k1
            self._k2 = self.k2
        
        if not self.is_set('vth0'):
            self._vth0 = self.vfb + self.phi + self._k1 * self.sqrtPhi
        else:
            self._vth0 = abs(self.vth0)

        self.k1ox = self._k1 * self.tox / self.toxm
        self.k2ox = self._k2 * self.tox / self.toxm

        t1 = np.sqrt(const.epSi / const.epOx * self.tox * self.Xdep0)
        t0 = ad.safe_exp(-0.5 * self.dsub * self.leff / t1)
        self.theta0vb0 = (t0 + 2.0 * t0**2)
        # From freeda, ngspice:
        self._Tox = 1e8 * self.tox
        
        #Calculation of vbsc(Vbc) and Vbseff
        if self._k2 < 0.:
            self.vbsc = .9 * (self.phi - (.5 * self._k1 / self._k2)**2)
            if self.vbsc > -3.:
                self.vbsc = -3.
            elif self.vbsc < -30.:
                self.vbsc = -30.
        else:
            self.vbsc = -30.

        if self.u0 > 1.:
            self._u0 = self.u0 * 1e-4
        else:
            self._u0 = self.u0 

        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        ad.delete_tape(self)
        # Absolute temperature (note self.temp is in deg. C)
        absTemp = const.T0 + temp
        # Thermal voltage
        self._Vt = const.k * absTemp / const.q
        # -------------------------------------------------------

        self._ToTnm1 = absTemp / self._Tn - 1.

        self.vsattemp = self.vsat - self.at * self._ToTnm1
        self.rds0 = (self.rdsw + self.prt * self._ToTnm1) \
            / pow(self.weff * 1e6, self.wr)

        self.V0 = self.vbi - self.phi

        #Mobility calculation
        self._ua = self.ua + self.ua1 * self._ToTnm1
        self._ub = self.ub + self.ub1 * self._ToTnm1
        self._uc = self.uc + self.uc1 * self._ToTnm1

        self.u0temp = self._u0 * pow(absTemp / self._Tn, self.ute) 
        

    def eval_cqs(self, vPort, getOP = False):
        """
        Calculates currents and charges

        Input: vPort = [vdb , vgb , vsb]
        Output: iVec = [ids, idb, isb], qVec = [qd, qg, qs]

        If getOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # import pdb; pdb.set_trace()
        # Invert all voltages in case of a P-channel device
        vPort1 = self._tf * vPort
        # If vds is negative, swap Vd and Vs and (also swap charges
        # later)
        DtoSswap = ad.condassign(vPort1[0] - vPort1[2], 1., -1.)
        # perform the swap (need tmp variable)
        tmp = ad.condassign(DtoSswap, vPort1[0], vPort1[2])
        vPort1[2] = ad.condassign(DtoSswap, vPort1[2], vPort1[0])
        vPort1[0] = tmp

        # Calculate VDS, VGS and VBS for bsim model
        VDS = vPort1[0] - vPort1[2]
        VGS = vPort1[1] - vPort1[2]
        VBS = -vPort1[2]

        # ----------------------------------------------------------------
        T0 = VBS - self.vbsc - 0.001 
        T1 = np.sqrt(T0 * T0 - 0.004 * self.vbsc)
        Vbseff = self.vbsc + .5 * (T0 + T1)
        Vbseff = ad.condassign(-Vbseff + VBS,
                                VBS,
                                Vbseff)
        
        #Calculation of Phis, sqrtPhis and Xdep
        Phis = ad.condassign(Vbseff,
                             self.phi**2 / (self.phi + Vbseff),
                             self.phi - Vbseff)
        
        sqrtPhis = ad.condassign(Vbseff,
                                 self.phis3 / (self.phi + 0.5 * Vbseff),
                                 np.sqrt(Phis))
        
        Xdep = self.Xdep0 * sqrtPhis / self.sqrtPhi
        
        #Calculation of Threshold voltage-vth
        T3 = np.sqrt(Xdep)
        T1 = ad.condassign(self.dvt2 * Vbseff + .5,
                           1. + self.dvt2 * Vbseff,
                           (1. + 3. * self.dvt2 * Vbseff) 
                           / (3. + 8. * self.dvt2 * Vbseff))

        ltl = self.factor1 * T3 * T1

        T1 = ad.condassign(self.dvt2w * Vbseff + .5,
                           1. + self.dvt2w * Vbseff,
                           (1. + 3. * self.dvt2w * Vbseff) 
                           / (3. + 8. * self.dvt2w * Vbseff))

        ltw = self.factor1 * T3 * T1
        
        # Alternative to prevent overflow (apparently not needed)
        #T2 = ad.safe_exp(-.5 * self.dvt1 * self.leff / ltl)
        k_temp = -.5 * self.dvt1 * self.leff / ltl
        T2 = ad.condassign(k_temp + EXP_THRESHOLD,
                           ad.safe_exp(k_temp),
                           MIN_EXP)
        Theta0 =  T2 * (1. + 2. * T2)

        thetavth = self.dvt0 * Theta0
        Delt_vth = thetavth * self.V0

        # Alternative to prevent overflow (apparently not needed)
        #T2 = ad.safe_exp(-.5 * self.dvt1w * self.weff * self.leff / ltw)
        # Modified from freeda's source to prevent using uninitialized
        # variable
        k_temp = -.5 * self.dvt1w * self.weff * self.leff / ltw
        T2 = ad.condassign(k_temp + EXP_THRESHOLD,
                           ad.safe_exp(k_temp),
                           MIN_EXP)
        T2 *= (1. + 2. * T2)
                
        T0 = self.dvt0w * T2
        T2 = T0 * self.V0

        T0 = np.sqrt(1. + self.nlx / self.leff)
        T1 = self.k1ox * (T0 - 1.) * self.sqrtPhi \
            + (self.kt1 + self.kt1l / self.leff + self.kt2 * Vbseff) \
            * self._ToTnm1
        TMP2 = self.tox * self.phi / (self.weff + self.w0)
        
        T3 = self.eta0 + self.etab * Vbseff
        T4 = ad.condassign(-T3 + 1.0e-4,
                            1. / (3. - 2e4 * self.eta0 + self.etab * Vbseff),
                            1.)

        dDIBL_Sft_dVd = T3 * self.theta0vb0
        DIBL_Sft = dDIBL_Sft_dVd * VDS
        
        Vth = self._vth0 - self._k1 * self.sqrtPhi + self.k1ox * sqrtPhis \
            - self.k2ox * Vbseff - Delt_vth - T2 \
            + (self.k3 + self.k3b * Vbseff) * TMP2 + T1 - DIBL_Sft
        
        #Calculate n
        temp_tmp2 = self.nfactor * const.epSi / Xdep
        temp_tmp3 = self.cdsc + self.cdscb * Vbseff + self.cdscd * VDS
        temp_tmp4 = (temp_tmp2 + temp_tmp3 * Theta0 + self.cit) / self.cox
        
        n = ad.condassign(temp_tmp4 + .5,
                          1. + temp_tmp4,
                          (1. + 3. * temp_tmp4) \
                              * (1. / (3. + 8. * temp_tmp4)))
                          
        #Poly Gate Si Depletion Effect
        Vgs_eff = VGS
        
        Vgst = Vgs_eff - Vth # not in Nikhil's code
        
        #Effective Vgst (Vgsteff) Calculation
        T10 = 2. * n * self._Vt
        VgstNVt = Vgst / T10
        ExpArg = -(2. * 0.08 + Vgst) / T10
        
        T1 = T10 * log1pexp(VgstNVt)
        T2 = 1. + T10 * self.cox * ad.safe_exp(ExpArg) / self._Vt / self.cdep0
        
        Vgsteff = ad.condassign(VgstNVt - EXP_THRESHOLD,
                                Vgst,
                                T1 / T2)
        Vgsteff = ad.condassign(
            ExpArg - EXP_THRESHOLD,
            self._Vt * self.cdep0 \
                / self.cox / ad.safe_exp((Vgst + 0.08) / n / self._Vt),
            T1 / T2)
        
        T3 = T2 * T2
        
        # Calculate Effective Channel Geometry
        T9 = sqrtPhis - self.sqrtPhi
        k_temp = self.weff - 2. * (self.dwg * Vgsteff + self.dwb * T9)
        
        Weff = ad.condassign(-k_temp + 2.0e-8,
                              2e-8 * (4.0e-8 - k_temp) * T0,
                              k_temp)

        T0 = self.prwg * Vgsteff + self.prwb * (sqrtPhis - self.sqrtPhi)
        Rds = ad.condassign(T0 + 0.9,
                            self.rds0 * (1. + T0),
                            self.rds0 * (.8 +T0) / (17. + 20. * T0))

        #Calculate Abulk
        T1 = 0.5 * self.k1ox / sqrtPhis
        T9 = np.sqrt(self.xj * Xdep)
        T5 = self.leff / (self.leff + 2. * T9)
        T2 = (self.a0 * T5) + self.b0 / (self.weff + self.b1)
        T6 = T5 * T5
        T7 = T5 * T6
        
        Abulk0 = 1. + T1 * T2
        
        T8 = self.ags * self.a0 * T7
        Abulk = Abulk0 - T1 * T8 * Vgsteff

        Abulk0 = ad.condassign(-Abulk0 + .1,
                                (.2 - Abulk0) / (3. - 20. * Abulk0),
                                Abulk0)

        Abulk = ad.condassign(-Abulk + .1,
                               (.2 - Abulk) / (3. - 20. * Abulk),
                               Abulk)
        
        T2 = self.keta * Vbseff
        
        T0 = ad.condassign(T2 + 0.9,
                           1. / (1. + T2),
                           (17. + 20. * T2) / (0.8 + T2))

        Abulk *= T0
        Abulk0 *= T0
        
        T0 = Vgsteff + 2. * Vth
        T2 = self._ua + self._uc * Vbseff
        T3 = T0 / self.tox
        T5 = T3 * (T2 + self._ub * T3)
        
        Denomi =  ad.condassign(T5 + .8,
                                1. + T5,
                                (.6 + T5) / (7. + 10. * T5))

        ueff = self.u0temp / Denomi
        
        Esat = 2. * self.vsattemp / ueff
        
        # Saturation Drain Voltage Vdsat
        WVCox = Weff * self.vsattemp * self.cox
        WVCoxRds = WVCox * Rds
        
        EsatL = Esat * self.leff
        Lambda = self.a2
        
        Vgst2Vtm = Vgsteff + 2. * self._Vt
        
        T0 = 1. / (Abulk * EsatL + Vgst2Vtm)
        T3 = EsatL * Vgst2Vtm
        Vdsat = T3 * T0
        
        # Effective Vds(Vdseff) Calculation
        T1 = Vdsat - VDS - self.delta
        T2 = np.sqrt(T1**2 + 4. * self.delta * Vdsat)
        #T0 = T1 / T2
        #T3 = 2. * self.delta / T2
        k_temp = Vdsat - .5 * (T1 + T2)
        Vdseff = ad.condassign(k_temp - VDS,
                               VDS,
                               k_temp)

        # The following seems unnnecessary:
        # Added to eliminate non-zero Vdseff at Vds=0.0
        #Vdseff = ad.condassign(abs(VDS),
        #                       Vdseff,
        #                       0.)
        
        # Calculate Vasat
        T6 = 1. - .5 * Abulk * Vdsat / Vgst2Vtm #T6=tmp4
        T9 = WVCoxRds * Vgsteff #expanded
        T0 = EsatL + Vdsat + 2. * T9 * T6
        
        T9 = WVCoxRds * Abulk
        T1 = 2. / Lambda - 1. + T9
        
        Vasat = T0 / T1
        
        diffVds = VDS - Vdseff
        
        #Calculate VACLM
        VACLM = ad.condassign(
            (diffVds - 1.0e-10) * self.pclm,
            self.leff * (Abulk + Vgsteff / EsatL) * diffVds \
                / (self.pclm * Abulk * self.litl),
            MAX_EXP)

        #Calculate VADIBL
        T1 = np.sqrt(const.epSi / const.epOx * self.tox * self.Xdep0)
        
        T0 = ad.safe_exp(-.5 * self.drout * self.leff / T1)
        T2 = T0 + 2. * T0**2
        thetaRout = self.pdibl1 * T2 + self.pdibl2 #drout, pdibl1,
                                                   #pdibl2 are given
        VADIBL = ad.condassign(
            thetaRout,
            (Vgst2Vtm - Vgst2Vtm * Abulk * Vdsat / (Vgst2Vtm + Abulk * Vdsat)) \
                / thetaRout,
            MAX_EXP)
            
        VADIBL = ad.condassign(self.pdiblb * Vbseff + 0.9,
                               VADIBL / (1. + self.pdiblb * Vbseff),
                               VADIBL * (17. + 20. * self.pdiblb * Vbseff) \
                                   / (.8 + self.pdiblb * Vbseff))
                               
        #Calculate Va
        T8 = self.pvag / EsatL
        T9 = T8 * Vgsteff
        
        T0 = ad.condassign(T9 + 0.9,
                           1. + T9,
                           (.8 + T9) / (17. + 20. * T9))
        
        T3 = VACLM + VADIBL #tmp3 = T3
        T1 = VACLM * VADIBL / T3
        Va = Vasat + T0 * T1
        
        #Calculate VASCBE
        if self.pscbe1 != 0.:
            rcpVASCBE = ad.condassign(
                abs(diffVds),
                self.pscbe2 * ad.safe_exp(-self.pscbe1 * self.litl/diffVds) \
                    / self.leff,
                0.)
        else:
            rcpVASCBE = self.pscbe2 / self.leff

        # Original:
        #VASCBE = ad.condassign(
        #    diffVds - self.pscbe1 * self.litl / EXP_THRESHOLD,
        #    self.leff * np.exp(self.pscbe1 * self.litl/diffVds) / self.pscbe2,
        #    MAX_EXP * self.leff / self.pscbe2)

        #Calculate Ids
        CoxWovL = self.cox * Weff / self.leff
        beta = ueff * CoxWovL
        
        T0 = 1. - .5 * Abulk * Vdseff / Vgst2Vtm
        
        fgche1 = Vgsteff * T0
        
        T9 = Vdseff / EsatL
        fgche2 = 1. + T9
        
        gche = beta * fgche1 / fgche2
        
        T0 = 1. + gche * Rds
        T9 = Vdseff / T0
        Idl = gche * T9
        
        T9 = diffVds / Va
        T0 = 1. + T9
        Idsa = Idl * T0
        
        T9 = diffVds * rcpVASCBE
        T0 = 1. + T9
        Ids = Idsa * T0
        
        #Substrate current begins
        T1 = self.alpha0 + self.alpha1 * self.leff
        k_temp = T1 / self.leff * diffVds
        #T2 = k_temp * ad.safe_exp(-self.beta0 / diffVds)
        # Original T2:
        ad.condassign(
            diffVds - self.beta0 / EXP_THRESHOLD,
            k_temp * ad.safe_exp(-self.beta0 / diffVds),
            k_temp * MIN_EXP)

        k_temp = ad.condassign(T1, T2 * Idsa, 0.)
        Isub = ad.condassign(self.beta0, k_temp, 0.)
        
        # ********* Output current vector **************
        iVec = np.array([0., 0., 0.], dtype=type(Ids))
        iVec[0] = DtoSswap * Ids
        iVec[1] = ad.condassign(DtoSswap, Isub, 0.)
        iVec[2] = ad.condassign(DtoSswap, 0., Isub)
        # Revert currents if needed
        iVec *= self._tf
        # **********************************************

        #--------------------------------------------------------------
        # Charge calculation follows (does not work):

        #calculation of vfbzb
        T0 = -.5 * self.dvt1w * self.weff * self.leff \
            / self.factor1 / np.sqrt(self.Xdep0)
        
        #T2 = ad.safe_exp(T0) * (1. + 2. * ad.safe_exp(T0))
        T2 = ad.condassign(T0 + EXP_THRESHOLD,
                           ad.safe_exp(T0) * (1. + 2. * ad.safe_exp(T0)),
                           MIN_1)
        T0 = self.dvt0w * T2
        T2 = T0 * (self.vbi - self.phi)
        
        T0 = -.5 * self.dvt1 * self.leff / (self.factor1 * np.sqrt(self.Xdep0))

        #T3 = ad.safe_exp(T0) * (1. + 2. * ad.safe_exp(T0))
        T3 = ad.condassign(T0 + EXP_THRESHOLD,
                           ad.safe_exp(T0) * (1. + 2. * ad.safe_exp(T0)),
                           MIN_1)
        
        T3 = self.dvt0 * T3 * (self.vbi - self.phi)
        
        T4 = self.tox * self.phi / (self.weff + self.w0)
        T5 = self.k1ox * (T0 - 1.) * self.sqrtPhi \
            + (self.kt1 + self.kt1l / self.leff) * self._ToTnm1

        T0 = np.sqrt(1. + self.nlx / self.leff)
        T6 = self._vth0 - T2 -T3 + self.k3 * T4 + T5
        vfbzb = T6 - self.phi - self._k1 * self.sqrtPhi
        
        #Calculation for VbseffCV
        VbseffCV = ad.condassign(Vbseff,
                                 self.phi - Phis,
                                 Vbseff)
        
        #Calculation for VgsteffCV
        T0 = n * self.noff * self._Vt
        T1 = (Vgs_eff - Vth) / T0

        Vgsteff = ad.condassign(T1 - EXP_THRESHOLD,
                                Vgs_eff - Vth - self.voffcv,
                                T0 * log1pexp(T1))
                
        # This (after the previous) can not be right:
        #Vgsteff = ad.condassign(-T1 - EXP_THRESHOLD,
        #                         T0 * log(1. + MIN_EXP),
        #                         T0 * log(1. + np.exp(T1)))

        #Calculation for Vfbeff
        V3 = vfbzb - Vgs_eff + VbseffCV - .02
        T0 = np.sqrt(V3**2 + .08 * abs(vfbzb))
        Vfbeff = vfbzb - 0.5 * (V3 + T0)
        
        T0 = (Vgs_eff - VbseffCV - vfbzb) / self._Tox
        
        #Calculation for Tcen
        T1 = T0 * self.acde
        
        Tcen = ad.condassign(EXP_THRESHOLD + T1,
                             self.ldeb * np.exp(T1),
                             self.ldeb * MIN_EXP)
        Tcen = ad.condassign(-EXP_THRESHOLD + T1,
                             self.ldeb * MAX_EXP,
                             Tcen)
        
        V3 = self.ldeb - Tcen - 1e-3 * self.tox
        V4 = np.sqrt(V3**2 + 4e-3 * self.tox * self.ldeb)
        Tcen = self.ldeb - .5 * (V3 + V4)
        Ccen = const.epSi / Tcen
        Coxeff = Ccen * self.cox / (Ccen + self.cox)
        
        #Calculation for QoverlapCox
        CoxWLcen = Weff * self.leff * Coxeff 
        Qac0 = CoxWLcen * (Vfbeff - vfbzb)
        # QovCox = Qac0 / Coxeff
        
        T0 = .5 * self.k1ox
        T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff
        T1 = np.sqrt(T0 * T0 + T3)
        T2 = CoxWLcen * T0 / T1

        T1 = ad.condassign(-T3,
                            T0 + T3 / self.k1ox,
                            np.sqrt(T0**2 + T3))

        T2 = ad.condassign(-T3,
                            CoxWLcen,
                            CoxWLcen * T0 / T1)
        
        Qsub0 = CoxWLcen * self.k1ox * (T1 - T0)
        # QovCox = Qsub0 / Coxeff
        
        #Calculation for Delta_phis
        T2 = ad.condassign(self.k1ox,
                           self.moin * self._Vt * self.k1ox**2,
                           .25 * self.moin * self._Vt)

        T0 = ad.condassign(self.k1ox,
                           self.k1ox * np.sqrt(self.phi),
                           .5 * np.sqrt(self.phi))

        T1 = 2. * T0 + Vgsteff
        DeltaPhi = self._Vt * np.log(1. + T1 * Vgsteff / T2)
        
        #The calculation for Tcen must be done once more
        # fREEDA:
        #T0 = (Vgsteff + 4.*(self._vth0 - self.vfb - self.phi))/ (2. * self._Tox)
        # ngspice:
        T3 =  4. * (Vth - vfbzb - self.phi)
        T0 = ad.condassign(T3,
                           .5 * (Vgsteff + T3) / self._Tox,
                           .5 * (Vgsteff + 1.0e-20) / self._Tox)

        k_temp = 2.01375270747048 * T0  # was np.exp(.7 * log(T0))
        T1 = 1. + k_temp
        T2 = 0.35 * k_temp / (T0 * self._Tox)
        Tcen = 1.9e-9 / T1
        
        Ccen = const.epSi / Tcen
        Coxeff = Ccen * self.cox / (Ccen + self.cox)
        CoxWLcen = Weff * self.leff * Coxeff 

        AbulkCV = Abulk0 * self.AbulkCVfactor
        VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV

        T0 = VdsatCV - VDS - .02
        T1 = np.sqrt(T0**2 + .08 * VdsatCV)

        # From ngspice: internal version BSIM3v32V32
        VdseffCV = VdsatCV - .5 * (T0 + T1)
        # From freeda:
#        VdseffCV = ad.condassign(T0,
#                                 VdsatCV - .5 * (T0 + T1),
#                                 VdsatCV * (1. - .04/(T1-T0)))
        # Need this to prevent NaN in charge Jacobian
        VdseffCV += 1e-200
                                 
        # Seems not needed
        #VdseffCV = ad.condassign(abs(VDS),
        #                         Vdseff,
        #                         0.)
        
        T0 = AbulkCV * VdseffCV
        T1 = Vgsteff - DeltaPhi
        T2 = 12. * (T1 - .5 * T0 + 1e-20)
        T3 = T0 / T2
        T4 = 1. - 12. * T3**2
        T5 = AbulkCV * (6. * T0 * (4. * T1 - T0) / (T2**2) - .5)
        T6 = T5 * VdseffCV / AbulkCV

        qgate = CoxWLcen * (T1 - T0 * (.5 - T3))
        qbulk = CoxWLcen * (1. - AbulkCV) * (.5*VdseffCV - T0*VdseffCV/T2)
        
        #QovCox = qbulk / Coxeff
        
        T2 = T2 / 12.
        T3 = .5 * CoxWLcen / (T2**2)
        T4 = T1 * (2. * T0**2 / 3. + T1*(T1 - 4. * T0 / 3.)) - 2. * T0**3 / 15.
        
        qsrc = -T3 * T4
        qgate += Qac0 + Qsub0 - qbulk
        qbulk -= (Qac0 + Qsub0)
        qdrn = -(qbulk + qgate + qsrc)
        
        # ************ Output charge vector ****************
        qVec = np.array([0., qgate, 0.], dtype=type(qsrc))
        # have to switch charges if Drain and Source voltages switched
        qVec[0] = ad.condassign(DtoSswap, qdrn, qsrc)
        qVec[2] = ad.condassign(DtoSswap, qsrc, qdrn)
        # invert sign if needed
        qVec *= self._tf

        # Return numpy array with one element per current source.
        return (iVec, qVec)


    # Use AD for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    
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

        For now it is quite incomplete
        Input:  vPort = [vdb , vgb , vsb]
        Output: dictionary with OP variables

        (This needs more work)
        """
        outV = self.eval(vPort)

        opDict = dict(
            temp = self.temp,
            VD = vPort[0],
            VG = vPort[1],
            VS = vPort[2],
            IDS = outV[0]
            )
        return opDict


# Define extrinsic model
ExtBSIM3 = mosExt.extrinsic_mos(BSIM3)

devList = [BSIM3, ExtBSIM3]
