"""
:mod:`bjt` -- Bipolar Transistor
--------------------------------

.. module:: bjt
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
import cppaddev as ad
from diode import Junction

def extrinsic_bjt(IBJT):
    class BJT(IBJT):
        """
    Bipolar Junction Transistor
    ---------------------------

    This device accepts 3 or 4 terminal connections.

    Netlist examples::
    
        bjt:q1 2 3 4 1 model = mypnp isat=4e-17 bf=147 iss=10fA
        bjt:q2 2 3 4  model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m
        svbjt:q3 2 3 4 1 model = mypnp vaf=80 ikf=4m iss=15fA
    
        # Electro-thermal versions
        bjt_t:q2 2 3 5 1 pout gnd model = mypnp
        svbjt_t:q3 2 3 5 1 pout gnd model = mypnp
    
        # Model statement
        .model mypnp bjt_t (type=pnp isat=5e-17 cje=60fF vje=0.83 mje=0.35)

    Extrinsic Internal Topology
    +++++++++++++++++++++++++++

    RC, RE and a Collector-Bulk connection are added to intrinsic
    BJT models::

                  RC      ct             et    RE
      C (0) o---/\/\/\/--+-----,         4----/\/\/\/----o  E (2)
                         |      \       /
                         |       \     /     
                       -----    ---------
                        / \         |
                       /   \        o 
                       -----
                         |          B (1)
                         o Bulk (3)

    If RE or RC are zero the internal nodes (ct, et) are not
    created. If only 3 connections are specified then the
    Bulk-Collector junction is not connected.

    Important Note
    ++++++++++++++

    This implementation does not account for the power dissipation
    in RE, RC. Use external thermal resistors if that is needed.

    Intrinsic Model Information
    +++++++++++++++++++++++++++

        """ 
        # Additional documentation
        extraDoc = IBJT.__doc__

        # Create electrothermal device
        makeAutoThermal = True
        
        isNonlinear = True
    
        # devtype is the 'model' name: remove the 'i' from intrinsic name
        devType = IBJT.devType 
    
        # Do not set numTerms to allow 3 or 4 terminals to be connected
        
        paramDict = dict(
            IBJT.paramDict.items(),
            re = ('Emitter ohmic resistance', 'W', float, 0.),
            rc = ('Collector ohmic resistance', 'W', float, 0.),
            cjs = ('Collector substrate capacitance', 'F', float, 0.),
            mjs = ('substrate junction exponential factor', '', float, 0.),
            vjs = ('substrate junction built in potential', 'V', float, 0.75),
            ns = ('substrate p-n coefficient', '', float, 1.),
            iss = ('Substrate saturation current', 'A', float, 0.)
            )

        def __init__(self, instanceName):
            IBJT.__init__(self, instanceName)
            # Collector-bulk junction
            self.cbjtn = Junction()
            self.__doc__ += IBJT.__doc__

        def process_params(self, thermal = False):
            if thermal:
                extraTerms = 2
            else:
                extraTerms = 0

            # First check external connections
            if self.numTerms == 3 + extraTerms:
                self.__addCBjtn = False
            elif self.numTerms == 4 + extraTerms:
                self.__addCBjtn = True
            else:
                raise cir.CircuitError( 
                    '{0}: Wrong number of connections. \
Can only be {1} or {2}, {3} found.'.format(self.nodeName, 
                                           3 + extraTerms,
                                           4 + extraTerms,
                                           self.numTerms))

            # Remove internal terminals
            self.clean_internal_terms()
            # Remove tape if present
            ad.delete_tape(self)
            # Tell autothermal to re-generate thermal ports
            self.__addThermalPorts = True

            extraVCCS = list()
            if self.re:
                # Add ei node and translate port descriptions
                self._et = self.add_internal_term('et', 'V')
                extraVCCS += [((2, self._et), (2, self._et), 1./self.re)]
            if self.rc:
                # Add ci node and translate port descriptions
                self._ct = self.add_internal_term('ci', 'V')
                extraVCCS += [((0, self._ct), (0, self._ct), 1./self.rc)]

            # Process parameters from intrinsic device: emitter and
            # collector terminals are already substituted.
            IBJT.process_params(self)
            # Calculate variables in junction
            self.cbjtn.process_params(self.iss, self.ns, self.fc, self.cjs, 
                                      self.vjs, self.mjs, self.xti, self.eg, 
                                      self.Tnomabs)
            # Add RE, RC resistors (if any)
            self.linearVCCS += extraVCCS
            if self.__addCBjtn:
                # Add bulk-collector junction
                self.__bccn = len(self.controlPorts)
                self.__bcon = len(self.csOutPorts) 
                self.controlPorts.append((3, self._ct))
                self.csOutPorts.append((3, self._ct))
                if self.cjs:
                    self.qsOutPorts.append((3, self._ct))

                # Initial guess for input ports: 
                try:
                    if len(self.vPortGuess) < len(self.controlPorts):
                        self.vPortGuess = np.concatenate(
                            (self.vPortGuess, [1]), axis=0)
                except AttributeError:
                    # Ignore if vPortGuess not provided
                    pass

            # Adjust temperature
            self.set_temp_vars(self.temp)
            
        def set_temp_vars(self, temp):
            """
            Calculate temperature-dependent variables, given temp in deg. C
            """
            # First calculate variables from base class
            IBJT.set_temp_vars(self, temp)
            # Adjust collector-bulk junction temperature
            self.cbjtn.set_temp_vars(self)


        def eval_cqs(self, vPort, saveOP = False):
            """
            vPort is a vector with control voltages 
            """
            if self.__addCBjtn:
                # calculate currents and charges in base class
                (iVec1, qVec1) = IBJT.eval_cqs(self, vPort)
                # Add contribution of substrate. Ignore extra control
                # ports only used for power calculation
                v1 = vPort[self.__bccn] * self._typef
                isub = self.cbjtn.get_id(v1) * self.area * self._typef
                iVec = np.concatenate((iVec1, [isub]), axis = 0)
                if self.cjs:
                    qsub = self.cbjtn.get_qd(v1) * self.area * self._typef
                    qVec = np.concatenate((qVec1, [qsub]), axis = 0)
                else:
                    qVec = qVec1
                return (iVec, qVec)
            else:
                return IBJT.eval_cqs(self, vPort)
    
        # Create these using the AD facility
        eval_and_deriv = ad.eval_and_deriv
        eval = ad.eval
    
        def power(self, vPort, currV):
            """ 
            Calculate total instantaneous power 

            Power in RE, RC not considered.

            Input: control voltages as in eval_cqs() and currents
            returned by eval_cqs()
            """
            pout = IBJT.power(self, vPort, currV)
            if self.__addCBjtn:
                # add power from substrate junction and RE, RC
                pout += vPort[self.__bccn] * currV[self.__bcon]

            return pout


    # Return template class
    return BJT

#----------------------------------------------------------------------------

class BJTi(cir.Element):
    """
    Gummel-Poon intrinsic BJT model

    This implementation based mainly on previous implementation in
    carrot and some equations from Pspice manual.
    
    Terminal order: 0 Collector, 1 Base, 2 Emitter::

                      
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

    Internally may add 2 additional nodes (plus reference) if rb is
    not zero: Bi for the internal base node and tib to measure the
    internal base current and calculate Rb(ib). The possible
    configurations are described here.

    1. If RB == 0::

                         +----------------+--o 0 (C)
                         |                |
                        /^\               |
                       ( | ) ibc(vbc)     |
                        \|/               |       
                         |               /|\       
         (B) 1 o---------+              ( | ) ice    
                         |               \V/      
                        /|\               |       
                       ( | ) ibe(vbe)     |
                        \V/               |
                         |                |
                         +----------------+--o 2 (E)

    2. If RB != 0::

                                     +----------------+--o 0 (C)
                                     |                |
                                    /^\               |
                                   ( | ) ibc(vbc)     |
                    gyr * tib       \|/               |       
                     ,---,           |               /|\       
         (B) 1 o----( --> )----------+ Bi           ( | ) ice    
                     `---`           |               \V/      
                                    /|\               |       
                                   ( | ) ibe(vbe)     |
                                    \V/               |
                                     |                |
                                     +----------------+--o 2 (E)
                     gyr v(1,Bi)  
                      ,---,       
                 +---( <-- )------+
                 |    `---`       |
          tref   |                | voltage: ib/gyr
             ,---+                |
             |   |    ,---,       |         
             |   +---( --> )------+ tib
             |        `---`       
            ---     gyr ib Rb(ib)
             V      
                                           
    Charge sources are connected between internal nodes defined
    above. If xcjc is not 1 but RB is zero, xcjc is ignored.
   
    """

    devType = "bjt"
    
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
        rb = ('Zero bias base resistance', 'W', float, 0.),
        rbm = ('Minimum base resistance', 'W', float, 0.),
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
        tf = ('Ideal forward transit time', 'S', float, 0.),
        xtf = ('Transit time bias dependence coefficient', '', float, 0.),
        vtf = ('Transit time dependency on vbc', 'V', float, 0.),
        itf = ('Transit time dependency on ic', 'A', float, 0.),
        tr = ('Ideal reverse transit time', 'S', float, 0.),
        xtb = ('Forward and reverse beta temperature coefficient', '', 
               float, 0.),
        xti = ('IS temperature effect exponent', '', float, 3.),
        tnom = ('Nominal temperature', 'C', float, 27.),
        area = ('Current multiplier', '', float, 1.)
        )

    def __init__(self, instanceName):
        cir.Element.__init__(self, instanceName)
        # Use junctions to model diodes and capacitors
        self.jif = Junction()
        self.jir = Junction()
        self.jile = Junction()
        self.jilc = Junction()
        # collector/emitter terminal numbers: may be re-mapped by
        # extrinsic device
        self._ct = 0
        self._et = 2

    def process_params(self):
        """
        Adjusts internal topology and makes preliminary calculations
        according to parameters
        """
        # Default configuration assumes rb == 0
        # ibe, ibc, ice 
        self.csOutPorts = [(1, self._et), (1, self._ct), 
                           (self._ct, self._et)]
        # Controling voltages are vbe, vbc 
        self.controlPorts = [(1, self._et), (1, self._ct)]
        self.vPortGuess = np.array([0., 0.])
        # qbe, qbc
        self.qsOutPorts = [(1, self._et), (1, self._ct)]

        # Define topology first
        # Flag to signal if the extra charge Qbx is needed or not
        self._qbx = False
        if self.rb:
            # rb is not zero: add internal terminals
            tBi = self.add_internal_term('Bi', 'V')
            tib = self.add_internal_term('ib', '{0} A'.format(glVar.gyr)) 
            tref = self.add_reference_term()
            # Linear VCCS for gyrator(s)
            self.linearVCCS = [((1, tBi), (tib, tref), glVar.gyr),
                               ((tib, tref), (1, tBi), glVar.gyr)]
            # ibe, ibc, ice, Rb(ib) * ib
            self.csOutPorts = [(tBi, self._et), (tBi, self._ct), 
                               (self._ct, self._et), (tref, tib)]
            # Controling voltages are vbie, vbic and gyrator port
            self.controlPorts = [(tBi, self._et), (tBi, self._ct), 
                                 (tib, tref)]
            self.vPortGuess = np.array([0., 0., 0.])
            # qbie, qbic
            self.qsOutPorts = [(tBi, self._et), (tBi, self._ct)]
            # Now check if Cjbc must be splitted (since rb != 0)
            if self.cjc and (self.xcjc < 1.):
                self.qsOutPorts.append((1, self._ct))
                self._qbx = True
            
        # In principle we may not need any charge
        keepPorts = [ ]
        if self.cje + self.tf:
            # keep qbe
            keepPorts.append(self.qsOutPorts[0])
        if self.cjc + self.tr:
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
        if self.ise:
            # jile produces ile
            self.jile.process_params(self.ise, self.ne, 0, 0, 0, 0, 
                                     self.xti, self.eg, self.Tnomabs)
        if self.isc:
            # jilc produces ilc
            self.jilc.process_params(self.isc, self.nc, 0, 0, 0, 0, 
                                     self.xti, self.eg, self.Tnomabs)
        # Constants needed for rb(ib) calculation
        if self.irb:
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
        self.jif.set_temp_vars(self)
        self.jir.set_temp_vars(self)
        # Adjust ise and isc (which have different temperature variation)
        if self.ise:
            self.jile.set_temp_vars(self)
            self.jile._t_is /= tnXTB
        if self.isc:
            self.jilc.set_temp_vars(self)
            self.jilc._t_is /= tnXTB
        # Now some BJT-only variables
        self._bf_t = self.bf * tnXTB
        self._br_t = self.br * tnXTB
        

    def eval_cqs(self, vPort):
        """
        Calculates currents/charges

        Input is a vector may be one of the following, depending on
        parameter values::

          vPort = [vbe, vbc]
          vPort = [vbie, vbic, v1_i] (gyrator voltage, rb != 0)

        Output also depends on parameter values. Charges only present
        if parameters make them different than 0 (i.e., cje, tf, cjc,
        etc. are set to nonzero values)::
        
          iVec = [ibe, ibc, ice]
          iVec = [ibe, ibc, ice, gyr*ib*Rb] (rb != 0)

          qVec = [qbe, qbc]
          qVec = [qbe, qbc, qbx] (rb != 0 and cjc != 1)
        """
        # Invert control voltages if needed
        vPort1 = self._typef * vPort

        # Calculate regular PN junctions currents and charges
        ibf = self.jif.get_id(vPort1[0])
        ibr = self.jif.get_id(vPort1[1])
        if self.ise:
            ile = self.jile.get_id(vPort1[0])
        else:
            ile = 0.
        if self.isc:
            ilc = self.jilc.get_id(vPort1[1])
        else:
            ilc = 0.
        # Kqb
        q1m1 = 1.
        if self.var:
            q1m1 -= vPort1[0] / self.var
        if self.vaf:
            q1m1 -= vPort1[1] / self.vaf
        kqb = 1. / q1m1
        q2 = 0.
        if self.ikf:
            q2 += ibf / self.ikf 
        if self.ikr:
            q2 += ibr / self.ikr
        if q2:
            kqb *= .5 * (1. + np.sqrt(1. + 4. * q2))

        # Create output vector [ibe, ibc, ice, ...]
        iVec = np.zeros(self.ncurrents, dtype = type(ibf))
        qVec = np.zeros(self.ncharges, dtype = type(ibf))
        # ibe
        iVec[0] = ibf / self._bf_t + ile
        # ibc
        iVec[1] = ibr / self._br_t + ilc
        # ice
        iVec[2] = (ibf - ibr) / kqb

        # RB
        if self.rb:
            # Using gyrator
            # vPort1[2] not defined if rb == 0
            # ib has area effect included (removed by _ck1 and _ck2)
            ib = vPort1[2] * glVar.gyr
            if self.irb:
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
            iVec[3] = glVar.gyr * ib * rb / pow(self.area, 2)
            vbcx = ib * rb / self.area + vPort1[1]

        # Charges ----------------------------------------------- 

        # Note that if tf == 0 and cje == 0, nothing is calculated and
        # nothing is assigned to the output vector.

        # qbe is the first charge (0)
        if self.tf:
            # Effective tf
            tfeff = self.tf
            if self.vtf:
                x = ibf / (ibf + self.itf)
                tfeff *= (1. + self.xtf * x*x * 
                          ad.safe_exp(vPort1[1] /1.44 /self.vtf))
            qVec[0] = tfeff * ibf
        if self.cje:
            qVec[0] += self.jif.get_qd(vPort1[0]) 

        # qbc 
        if self._qbx:
            if self.tr:
                qVec[-2] = self.tr * ibr
            if self.cjc:
                qVec[-2] += self.jir.get_qd(vPort1[1]) * self.xcjc 
                # qbx
                qVec[-1] = self.jir.get_qd(vbcx) * (1. - self.xcjc)
        else:
            if self.tr:
                qVec[-1] = self.tr * ibr
            if self.cjc:
                qVec[-1] += self.jir.get_qd(vPort1[1]) 

        # Consider area effect and invert currents if needed
        iVec *= self.area * self._typef
        qVec *= self.area * self._typef

        return (iVec, qVec)


    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages as in eval_cqs() and currents returned
        by eval_cqs()
        """
        vce = vPort[0] - vPort[1]
        if self.rb:
            # currV[3] = ib * Rb * gyr
            # vPort[2] = ib / gyr
            pRb = currV[3] * vPort[2]
        else:
            pRb = 0.

        # pout = ibe * vbie + ibc * vbic + vce * ice + pRb
        pout = currV[0] * vPort[0] + currV[1] * vPort[1] \
            + currV[2] * vce + pRb

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

        self.OP = dict(
            VBE = vPort[0],
            VCE = vPort[0] - vPort[1],
            IB = outV[0] + outV[1],
            IC = outV[2] - outV[1],
            IE = - outV[2] - outV[0],
            gm = jac[2,0] - jac[1,0],
            rpi = 1./(jac[0,0] + jac[1,0]),
            )
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None


#------------------------------------------------------------------------

Device = extrinsic_bjt(BJTi)
