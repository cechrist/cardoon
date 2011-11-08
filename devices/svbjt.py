"""
:mod:`svbjt` -- State-Variable Intrinsic Bipolar Transistor
-----------------------------------------------------------

.. module:: svbjt
.. moduleauthor:: Carlos Christoffersen

################################# This is not ready yet!
"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
import cppaddev as ad
from diode import Junction
from svdiode import SVJunction

class Device(cir.Element):
    """
    State-variable-based Gummel-Poon intrinsic BJT model based

    This implementation based mainly on previous implementation in
    carrot and some equations from Pspice manual.
    
    Terminal order: 0 Collector, 1 Base, 2 Emitter, (3 Bulk, not included)::

                      
      C (0) o----,         4----o  E (2)
                  \       /
                   \     /
                  ---------
                      |
                      o 
   
                      B (1)

    Can be used for NPN or PNP transistors.

    The state variable formulation is achieved by replacing the BE and
    BC diodes (Ibf, Ibr) with state-variable based diodes. This
    requires two additional variables (nodes) but eliminates large
    positive exponentials from the model.  In addition we may need up
    to 2 additional nodes (plus gnd) if rb is not zero: Bi(3) for the
    internal base node and, if rbm is specified, ib(4) to measure the
    internal base current and calculate Rb(ib).

    Bulk connection, RC, RE are not included for now.

    Netlist examples::

        bjt:q1 2 3 4 model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m

        # Electro-thermal version
        bjt_t:q2 2 3 5 pout gnd model = mypnp

        # Model statement
        .model mypnp bjt_t (type=pnp isat=5e-17 cje=60fF vje=0.83 mje=0.35)

    """

    devType = "bjt"
    
    numTerms = 3  # for now

    # Create electrothermal device
    makeAutoThermal = True

    isNonlinear = True

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
#        re = ('Emitter ohmic resistance', 'W', float, 0.),
        rb = ('Zero bias base resistance', 'W', float, 0.),
        rbm = ('Minimum base resistance', 'W', float, 0.),
        irb = ('Current at which rb falls to half of rbm', 'A', float, 0.),
#        rc = ('Collector ohmic resistance', 'W', float, 0.),
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

    # Default configuration assumes rb == 0
    # ibe, ibc, ice 
    csOutPorts = ((1, 2), (1, 0), (0, 2))
    # Controling voltages are vbe, vbc 
    controlPorts = ((1, 2), (1, 0))
    vPortGuess = np.array([0., 0.])
    # qbe, qbc
    qsOutPorts = ((3, 2), (3, 0))

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        # Use junctions to model diodes and capacitors
        self.jif = SVJunction()
        self.jir = SVJunction()
        self.jile = Junction()
        self.jilc = Junction()

    def process_params(self, circuit):
        """
        Adjusts internal topology and makes preliminary calculations
        according to parameters.
        """
        # Remove internal terminals
        self.clean_internal_terms(circuit)
        # Remove tape if present
        ad.delete_tape(self)

        # Define topology first
        # Flag to signal if the extra charge Qbx is needed or not
        self._qbx = False
        if self.rb:
            # rb is not zero: add internal terminal
            termList = [self.nodeName + ':Bi']
            if self.irb:
                # in addition we need gyrator
                termList += [self.nodeName + ':ib', 'gnd']
                # Linear VCCS for gyrator(s)
                linearVCCS = [[(1, 3), (4, 5), glVar.gyr],
                              [(4, 5), (1, 3), glVar.gyr]]
                # ibe, ibc, ice, Rb(ib) * ib
                self.csOutPorts = ((3, 2), (3, 0), (0, 2), (5, 4))
                # Controling voltages are vbie, vbic and gyrator port
                self.controlPorts = ((3, 2), (3, 0), (4, 5))
                # qbie, qbic
                self.qsOutPorts = ((3, 2), (3, 0))
            else:
                # We can use current source
                # ibe, ibc, ice
                self.csOutPorts = ((3, 2), (3, 0), (0, 2), (1, 3))
                # Controling voltages are vbie, vbic
                self.controlPorts = ((3, 2), (3, 0), (1, 3))
                # qbie, qbic
                self.qsOutPorts = ((3, 2), (3, 0))

            # Connect required nodes
            circuit.connect_internal(self, termList)
            # Now check if Cjbc must be splitted (since rb != 0)
            if self.cjc and (self.xcjc < 1.):
                # add extra charge source and control voltage
                self.controlPorts += ((1, 0), )
                self.qsOutPorts += ((1, 2), )
                self._qbx = True

        # Initially we may not need any charge
        keepPorts = ( )
        if self.cje + self.tf:
            # keep qbe
            keepPorts = self.qsOutPorts[1:]
        if self.cjc + self.tr:
            # keep qbc, qbx (if any)
            if self._qbx:
                keepPorts += self.qsOutPorts[-2:]
            else:
                keepPorts += self.qsOutPorts[-1:]
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
        # Calculate temperature-dependent variables
        self.set_temp_vars(self.temp)


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

          vPort = [xbe, xbc]
          vPort = [xbe, xbc, vbbi]  (rb != 0, irb == 0)
          vPort = [xbe, xbc, v4gnd] (gyrator voltage, irb != 0)
          vPort = [xbe, xbc, vbbi, vbc] (xcjc < 1)
          vPort = [xbe, xbc, v4gnd, vbc] (xcjc < 1)

        Output also depends on parameter values. Charges only present
        if parameters make them different than 0 (i.e., cje, tf, cjc,
        etc. are set to nonzero values)::
        
          outV = [ibe, vbe, ibc, vbc, ice, qbe, qbc]
          outV = [ibe, vbe, ibc, vbc, ice, vbbi/Rb, qbe, qbc] 
                 (rb != 0, irb == 0)
          outV = [ibe, vbe, ibc, vbc, ice, gyr*ib*Rb, qbe, qbc] (irb != 0)
          outV = [ibe, vbe, ibc, vbc, ice, vbbi/Rb, qbe, qbc, qbx] 
                 (rb != 0, irb == 0)
          outV = [ibe, vbe, ibc, vbc, ice, gyr*ib*Rb, qbe, qbc, qbx] (irb != 0)

        """
        # Invert state variables if needed
        vPort1 = self._typef * vPort

        # Calculate junctions currents and voltages
        (ibf, vbe) = self.jif.get_idvd(vbe)
        (ibr, vbc) = self.jif.get_idvd(vbc)
        if self.ise:
            ile = self.jile.get_id(vbe)
        else:
            ile = 0.
        if self.isc:
            ilc = self.jilc.get_id(vbc)
        else:
            ilc = 0.
        # Kqb
        q1m1 = 1.
        if self.var:
            q1m1 -= vbe / self.var
        if self.vaf:
            q1m1 -= vbc / self.vaf
        kqb = 1. / q1m1
        q2 = 0.
        if self.ikf:
            q2 += ibf / self.ikf 
        if self.ikr:
            q2 += ibr / self.ikr
        if q2:
            kqb *= .5 * (1. + np.sqrt(1. + 4. * q2))

        # Create output vector [ibe, ibc, ice, ...]
        outV = np.zeros((self.ncurrents + self.ncharges), dtype = type(ibf))
        # ibe, vbe
        outV[0] = ibf / self._bf_t + ile
        outV[1] = glVar.gyr * vbe
        # ibc, vbc
        outV[2] = ibr / self._br_t + ilc
        outV[3] = glVar.gyr * vbc
        # ice
        outV[4] = (ibf - ibr) / kqb

        # RB
        if self.irb:
            # Using gyrator
            # vPort1[2] not defined if irb == 0
            # ib has area effect included (removed by _ck1 and _ck2)
            ib = vPort1[2] * glVar.gyr
            ib1 = np.abs(ib)
            x = np.sqrt(1. + self._ck1 * ib1) - 1.
            x *= self._ck2 / np.sqrt(ib1)
            tx = np.tan(x)
            c = self.rbm + 3. * (self.rb - self.rbm) \
                * (tx - x) / (x * tx * tx)
            rb = ad.condassign(ib1, c, self.rb)
            # Output is gyr * ib * rb.  It is divided by area^2 to
            # compensate that the whole vector is multiplied by area
            # at the end
            outV[5] = glVar.gyr * ib * rb / pow(self.area, 2)
        elif self.rb:
            # Using current source = vPort1[2] / rb
            rb = self.rbm + (self.rb - self.rbm) / kqb
            # Output is vbib / Rb
            outV[5] = vPort1[2] / rb 

        # Charges ----------------------------------------------- 

        # Note that if tf == 0 and cje == 0, nothing is calculated and
        # nothing is assigned to the output vector.

        # qbe is the first charge (0)
        if self.tf:
            # Effective tf
            tfeff = self.tf
            if self.vtf:
                x = ibf / (ibf + self.itf)
                # safe_exp() not needed since positive vbc grows
                # logarithmically
                tfeff *= (1. + self.xtf * x*x * 
                          np.exp(vbc /1.44 /self.vtf))
            outV[self.ncurrents] = tfeff * ibf
        if self.cje:
            outV[self.ncurrents] += self.jif.get_qd(vbe) 

        # qbc 
        if self._qbx:
            if self.tr:
                outV[-2] = self.tr * ibr
            if self.cjc:
                outV[-2] += self.jir.get_qd(vbc) * self.xcjc 
                # qbx
                outV[-1] = self.jir.get_qd(vPort1[-1]) * (1. - self.xcjc)
        else:
            if self.tr:
                outV[-1] = self.tr * ibr
            if self.cjc:
                outV[-1] += self.jir.get_qd(vbc) 

        # Consider area effect and invert currents if needed
        outV *= self.area * self._typef

        return outV


    # Use AD for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    
    def power(self, vPort, currV):
        """ 
        Calculate total instantaneous power 

        Input: control voltages as in eval_cqs() and currents from
        returned by eval_cqs()
        """
        # vce = vbe - vbc
        gyrvce = currV[1] - currV[3]
        if self.rb or self.irb:
            # currV[5] = ib * Rb * gyr
            # vPort[2] = ib / gyr
            # or
            # currV[5] = VRb / Rb
            # vPort[2] = VRb
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

        # calculate gm, etc. in terms od jac for state-variable
        # formulation
        self.OP = dict(
            VBE = outV[1] / glVar.gyr,
            VCE = (currV[1] - currV[3]) / glVar.gyr,
            IB = outV[0] + outV[2],
            IC = outV[4] - outV[2],
            IE = - outV[4] - outV[0],
            )
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None



