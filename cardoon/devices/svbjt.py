"""
:mod:`svbjt` -- State-Variable-Based Bipolar Transistor
-------------------------------------------------------

.. module:: svbjt
.. moduleauthor:: Carlos Christoffersen

"""

from bjt import extrinsic_bjt

import numpy as np
import cardoon.circuit as cir
from cardoon.globalVars import const, glVar
import cppaddev as ad
from diode import Junction
from svdiode import SVJunction

class SVBJTi(cir.Element):
    """
    State-variable-based Gummel-Poon intrinsic BJT model based

    This implementation based mainly on previous implementation in
    carrot and some equations from Pspice manual, with the addition of
    the state-variable definitions.
    
    Terminal order: 0 Collector, 1 Base, 2 Emitter, (3 Bulk, not included)::

                      
      C (0) o----,         4----o  E (2)
                  \       /
                   \     /
                  ---------
                      |
                      o 
   
                      B (1)

    Can be used for NPN or PNP transistors.

    Intrinsic Internal Topology
    +++++++++++++++++++++++++++

    The state variable formulation is achieved by replacing the BE and
    BC diodes (Ibf, Ibr) with state-variable based diodes. This
    requires two additional variables (nodes) but eliminates large
    positive exponentials from the model::

                                  Term : x2 
                      +--------------------------+
                      |                          |
                     /|\                        /^\ 
                    ( | ) gyr v2               ( | ) gyr vbc(x)
                     \V/                        \|/  
             tref     |                          |
                 ,----+--------------------------+ 
                 |    |                          |               
                 |   /^\                        /|\              
                 |  ( | ) gyr v1               ( | ) gyr vbe(x)  
                ---  \|/                        \V/  
                 V    |                          |
                      +--------------------------+
                                  Term : x1                
                                                  
    All currents/charges in the model are functions of voltages v3
    (x2) and v4 (x1). Note that vbc and vbe are now also functions of
    x1, x2.

    In addition we may need 2 additional nodes (plus reference) if rb
    is not zero: Bi for the internal base node and tib to measure the
    internal base current and calculate Rb(ib).

    1. If RB == 0::

                           +----------------+--o 0 (C)
                    -      |                |
                          /^\               |
                   v2    ( | ) ibc(x2)      |
                          \|/               |       
                    +      |               /|\       
           (B) 1 o---------+              ( | ) ice(x1,x2)
                    +      |               \V/      
                          /|\               |       
                   v1    ( | ) ibe(x1)      |
                          \V/               |
                    -      |                |
                           +----------------+--o 2 (E)

    2. If RB != 0::

                                     +----------------+--o 0 (C)
                                -    |                |
                                    /^\               |
                  gyr tib      v2  ( | ) ibc(x2)      |
                                    \|/               |       
                     ,---,      +    |               /|\       
         (B) 1 o----( --> )----------+ Term : Bi    ( | ) ice(x1,x2)
                     `---`      +    |               \V/      
                                    /|\               |       
                               v1  ( | ) ibe(x1)      |
                                    \V/               |
                                -    |                |
                   gyr v(1,Bi)       +----------------+--o 2 (E)
                                  
                      ,---,       
                 +---( <-- ) -----+
                 |    `---`       |
          tref   |                | ib/gyr
              ,--+                |
              |  |    ,---,       | Term : ib
              |  +---( --> )------+
              |       `---`       
             --- 
              V     gyr ib Rb(ib)
                                           
    Charge sources are connected between internal nodes defined
    above. If xcjc is not 1 but RB is zero, xcjc is ignored.
    """

    devType = "svbjt"
    
    paramDict = dict(
        cir.Element.tempItem,
        type = ('Type (npn or pnp)', '', str, 'npn'),
        isat = ('Transport saturation current', 'A', float, 1e-16),
        bf = ('Ideal maximum forward beta', '', float, 100.),
        nf = ('Forward current emission coefficient', '', float, 1.),
        vaf = ('Forward early voltage', 'V', float, 0.),
        ikf = ('Forward-beta high current roll-off knee current', 'A', 
               float, 0.),
        ise = ('Base-emitter leakage saturation current', 'A', float, 0.),
        ne = ('Base-emitter leakage emission coefficient', '', float, 1.5),
        br = ('Ideal maximum reverse beta', '', float, 1.),
        nr = ('Reverse current emission coefficient', '', float, 1.),
        var = ('Reverse early voltage', 'V', float, 0.),
        ikr = ('Corner for reverse-beta high current roll off', 'A', float, 0.),
        isc = ('Base collector leakage saturation current', 'A', float, 0.),
        nc = ('Base-collector leakage emission coefficient', '', float, 2.),
        rb = ('Zero bias base resistance', 'Ohm', float, 0.),
        rbm = ('Minimum base resistance', 'Ohm', float, 0.),
        irb = ('Current at which rb falls to half of rbm', 'A', float, 0.),
        eg = ('Badgap voltage', 'eV', float, 1.11),
        cje = ('Base emitter zero bias p-n capacitance', 'F', float, 0.),
        vje = ('Base emitter built in potential', 'V', float, 0.75),
        mje = ('Base emitter p-n grading factor', '', float, 0.33),
        cjc = ('Base collector zero bias p-n capacitance', 'F', float, 0.),
        vjc = ('Base collector built in potential', 'V', float, 0.75),
        mjc = ('Base collector p-n grading factor', '', float, 0.33),
        xcjc = ('Fraction of cbc connected internal to rb', '', float, 1.),
        fc = ('Forward bias depletion capacitor coefficient', '', float, 0.5),
        tf = ('Ideal forward transit time', 's', float, 0.),
        xtf = ('Transit time bias dependence coefficient', '', float, 0.),
        vtf = ('Transit time dependency on vbc', 'V', float, 0.),
        itf = ('Transit time dependency on ic', 'A', float, 0.),
        tr = ('Ideal reverse transit time', 's', float, 0.),
        xtb = ('Forward and reverse beta temperature coefficient', '', 
               float, 0.),
        xti = ('IS temperature effect exponent', '', float, 3.),
        tnom = ('Nominal temperature', 'C', float, 27.),
        area = ('Current multiplier', '', float, 1.)
        )

    def __init__(self, instanceName):
        cir.Element.__init__(self, instanceName)
        # Use junctions to model diodes and capacitors
        self.jif = SVJunction()
        self.jir = SVJunction()
        self.jile = Junction()
        self.jilc = Junction()
        # collector/emitter terminal numbers: may be re-mapped by
        # extrinsic device
        self._ct = 0
        self._et = 2

    def process_params(self):
        """
        Adjusts internal topology and makes preliminary calculations
        according to parameters.
        """
        # Define topology first. Add state variable nodes
        x2 = self.add_internal_term('x2','s.v.')
        x1 = self.add_internal_term('x1','s.v.')
        tref = self.add_reference_term()
        # Default configuration assumes rb == 0
        # ibe, vbe, ibc, vbc, ice 
        self.csOutPorts = [(1, self._et), (tref, x1), (1, self._ct), 
                           (tref, x2), (self._ct, self._et)]
        # Controling voltages are x1, x2
        self.controlPorts = [(x1, tref), (x2, tref)]
        # qbe, qbc
        self.qsOutPorts = [(1, self._et), (1, self._ct)]

        # Flag to signal if the extra charge Qbx is needed or not
        self._qbx = False
        # Default state-variable VCCSs
        self.linearVCCS = [((1, self._et), (x1, tref), glVar.gyr),
                           ((1, self._ct), (x2, tref), glVar.gyr)]
        if self.rb != 0.:
            # rb is not zero: add internal terminals
            tBi = self.add_internal_term('Bi', 'V') 
            tib = self.add_internal_term('ib', '{0} A'.format(glVar.gyr)) 
            # Add Linear VCCS for gyrator(s)
            self.linearVCCS = [((tBi, self._et), (x1, tref), glVar.gyr),
                               ((tBi, 0), (x2, tref), glVar.gyr),
                               ((1, tBi), (tib, tref), glVar.gyr),
                               ((tib, tref), (1, tBi), glVar.gyr)]
            # ibe, vbe, ibc, vbc, ice, Rb(ib) * ib
            self.csOutPorts = [(tBi, self._et), (tref, x1), (tBi, self._ct), 
                               (tref, x2), (self._ct, self._et), (tref, tib)]
            # Controling voltages are x1, x2 and gyrator port
            self.controlPorts = [(x1, tref), (x2, tref), (tib, tref)]
            # qbie, qbic
            self.qsOutPorts = [(tBi, self._et), (tBi, self._ct)]
            # Now check if Cjbc must be splitted (since rb != 0)
            if (self.cjc != 0.) and (self.xcjc < 1.):
                # add extra charge source
                self.qsOutPorts.append((1, self._ct))
                self._qbx = True
        
        # Make sure the guess is consistent
        self.vPortGuess = np.zeros(len(self.controlPorts))
#        # Try guess in active region
#        self.vPortGuess[0] = 100. # x1
#        self.vPortGuess[1] = -1.  # x2
#        if self.rb != 0.:
#            self.vPortGuess[2] = 1e-6 / glVar.gyr  # ib 
        # In principle we may not need any charge
        keepPorts = [ ]
        if self.cje + self.tf != 0.:
            # keep qbe
            keepPorts.append(self.qsOutPorts[0])
        if self.cjc + self.tr != 0.:
            # keep qbc, qbx (if any)
            if self._qbx:
                keepPorts += self.qsOutPorts[-2:]
            else:
                keepPorts.append(self.qsOutPorts[-1])
        self.qsOutPorts = keepPorts

        # keep track of how many output variables are needed
        self.ncurrents = len(self.csOutPorts)
        self.ncharges = len(self.qsOutPorts)

        # NPN or PNP 
        if self.type == 'pnp':
            self._typef = -1.
        else:
            self._typef = 1.

        # Calculate common variables
        # Absolute nominal temperature
        self.Tnomabs = self.tnom + const.T0
        self.egapn = self.eg - .000702 * (self.Tnomabs**2) \
            / (self.Tnomabs + 1108.)
        # jif produces if, cje
        self.jif.process_params(self.isat, self.nf, self.fc, self.cje, 
                                self.vje, self.mje, self.xti, self.eg, 
                                self.Tnomabs)
        # jir produces ir, cjc
        self.jir.process_params(self.isat, self.nr, self.fc, self.cjc, 
                                self.vjc, self.mjc, self.xti, self.eg, 
                                self.Tnomabs)
        if self.ise != 0.:
            # jile produces ile
            self.jile.process_params(self.ise, self.ne, 0, 0, 0, 0, 
                                     self.xti, self.eg, self.Tnomabs)
        if self.isc != 0.:
            # jilc produces ilc
            self.jilc.process_params(self.isc, self.nc, 0, 0, 0, 0, 
                                     self.xti, self.eg, self.Tnomabs)
        # Constants needed for rb(ib) calculation
        if self.irb != 0.:
            self._ck1 = 144. / self.irb / self.area /np.pi/np.pi
            self._ck2 = np.pi*np.pi * np.sqrt(self.irb * self.area) / 24.

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        # Absolute temperature (note self.temp is in deg. C)
        self.Tabs = const.T0 + temp
        # Normalized temp
        self.tnratio = self.Tabs / self.Tnomabs
        tnXTB = pow(self.tnratio, self.xtb)
        # Thermal voltage
        self.vt = const.k * self.Tabs / const.q
        # Temperature-adjusted egap
        self.egap_t = self.eg - .000702 * (self.Tabs**2) / (self.Tabs + 1108.)
        # set temperature in juctions
        self.jif.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                               self.egapn, self.egap_t)
        self.jir.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                               self.egapn, self.egap_t)
        # Adjust ise and isc (which have different temperature variation)
        if self.ise != 0.:
            self.jile.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                                    self.egapn, self.egap_t)
            self.jile._t_is /= tnXTB
        if self.isc != 0.:
            self.jilc.set_temp_vars(self.Tabs, self.Tnomabs, self.vt, 
                                    self.egapn, self.egap_t)
            self.jilc._t_is /= tnXTB
        # Now some BJT-only variables
        self._bf_t = self.bf * tnXTB
        self._br_t = self.br * tnXTB
        

    def eval_cqs(self, vPort):
        """
        Calculates currents/charges

        Input is a vector may be one of the following, depending on
        parameter values::

          vPort = [xbe, xbc]
          vPort = [xbe, xbc, v4_i] (gyrator voltage, irb != 0)

        Output also depends on parameter values. Charges only present
        if parameters make them different than 0 (i.e., cje, tf, cjc,
        etc. are set to nonzero values)::
        
          iVec = [ibe, vbe, ibc, vbc, ice]
          iVec = [ibe, vbe, ibc, vbc, ice, gyr*ib*Rb] (rb != 0)

          qVec = [qbe, qbc]
          qVec = [qbe, qbc, qbx] (rb != 0 and cjc != 1)

        """
        # Invert state variables if needed
        vPort1 = self._typef * vPort

        # Calculate junctions currents and voltages
        (ibf, vbe) = self.jif.get_idvd(vPort1[0])
        (ibr, vbc) = self.jir.get_idvd(vPort1[1])
        if self.ise != 0.:
            ile = self.jile.get_id(vbe)
        else:
            ile = 0.
        if self.isc != 0.:
            ilc = self.jilc.get_id(vbc)
        else:
            ilc = 0.
        # Kqb
        q1m1 = 1.
        if self.var != 0.:
            q1m1 -= vbe / self.var
        if self.vaf != 0.:
            q1m1 -= vbc / self.vaf
        kqb = 1. / q1m1
        # We need extra checking to consider the following
        # possibilities to create the AD tape:
        # 
        # 1. both ikf and ikr are zero -> no tape generated
        # 2. One of them is nonzero but both ibf and ibr are zero -> want tape
        #    but only for the nonzero parameter
        if self.ikf + self.ikr != 0.:
            q2 = 0.
            if self.ikf != 0.:
                q2 += ibf / self.ikf 
            if self.ikr != 0.:
                q2 += ibr / self.ikr
            kqb *= .5 * (1. + np.sqrt(1. + 4. * q2))

        # Create output vector [ibe, ibc, ice, ...]
        iVec = np.zeros(self.ncurrents, dtype = type(ibf))
        qVec = np.zeros(self.ncharges, dtype = type(ibf))
        # ibe, vbe
        iVec[0] = ibf / self._bf_t + ile
        iVec[1] = glVar.gyr * vbe
        # ibc, vbc
        iVec[2] = ibr / self._br_t + ilc
        iVec[3] = glVar.gyr * vbc
        # ice
        iVec[4] = (ibf - ibr) / kqb

        # RB
        if self.rb != 0.:
            # Using gyrator
            # vPort1[2] not defined if rb == 0
            # ib has area effect included (removed by _ck1 and _ck2)
            ib = vPort1[2] * glVar.gyr
            if self.irb != 0.:
                ib1 = np.abs(ib)
                x = np.sqrt(1. + self._ck1 * ib1) - 1.
                x *= self._ck2 / np.sqrt(ib1)
                tx = np.tan(x)
                c = self.rbm + 3. * (self.rb - self.rbm) \
                    * (tx - x) / (x * tx * tx)
                rb = ad.condassign(ib1, c, self.rb)
            else:
                rb = self.rbm + (self.rb - self.rbm) / kqb
            # Output is gyr * ib * rb.  It is divided by area^2 to
            # compensate that the whole vector is multiplied by area
            # at the end
            iVec[5] = glVar.gyr * ib * rb / pow(self.area, 2)
            vbcx = ib * rb / self.area + vbc

        # Charges ----------------------------------------------- 

        # Note that if tf == 0 and cje == 0, nothing is calculated and
        # nothing is assigned to the output vector.

        # qbe is the first charge (0)
        if self.tf != 0.:
            # Effective tf
            tfeff = self.tf
            if self.vtf != 0.:
                x = ibf / (ibf + self.itf)
                # safe_exp() not needed since positive vbc grows
                # logarithmically
                tfeff *= (1. + self.xtf * x*x * 
                          np.exp(vbc /1.44 /self.vtf))
            qVec[0] = tfeff * ibf
        if self.cje != 0.:
            qVec[0] += self.jif.get_qd(vbe) 

        # qbc 
        if self._qbx:
            if self.tr != 0.:
                qVec[-2] = self.tr * ibr
            if self.cjc != 0.:
                qVec[-2] += self.jir.get_qd(vbc) * self.xcjc 
                # qbx
                qVec[-1] = self.jir.get_qd(vbcx) * (1. - self.xcjc)
        else:
            if self.tr != 0.:
                qVec[-1] = self.tr * ibr
            if self.cjc != 0.:
                qVec[-1] += self.jir.get_qd(vbc) 

        # Consider area effect and invert currents if needed
        iVec *= self.area * self._typef
        qVec *= self.area * self._typef

        return (iVec, qVec)


    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages as in eval_cqs() and currents from
        returned by eval_cqs()
        """
        # vce = vbe - vbc
        gyrvce = currV[1] - currV[3]
        if self.rb != 0.:
            # currV[5] = ib * Rb * gyr
            # vPort[2] = ib / gyr
            pRb = currV[5] * vPort[2]
        else:
            pRb = 0.

        # pout = ibe * vbie + ibc * vbic + vce * ice + pRb
        pout = (currV[0] * currV[1] + currV[2] * currV[3] 
                + currV[4] * gyrvce) / glVar.gyr + pRb

        return pout

    def get_OP(self, vPort):
        """
        Calculates operating point information

        Input: same as eval_cqs
        Output: dictionary with OP variables
        For now it is quite incomplete
        """
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)
        power = self.power(vPort, outV)

        # calculate gm, etc. in terms od jac for state-variable
        # formulation
        opDict = dict(
            VBE = outV[1] / glVar.gyr,
            VCE = (outV[1] - outV[3]) / glVar.gyr,
            IB = outV[0] + outV[2],
            IC = outV[4] - outV[2],
            IE = - outV[4] - outV[0],
            Temp = self.temp,
            Power = power,
            )
        return opDict

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None

#--------------------------------------------------------------------------

Device = extrinsic_bjt(SVBJTi)
