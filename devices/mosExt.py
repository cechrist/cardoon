"""
:mod:`mosExt` -- Extrinsic MOSFET model
---------------------------------------

.. module:: mosExt
.. moduleauthor:: Carlos Christoffersen

This module add the extrinsic part to an intrinsic mosfet model.
Requirements for intrinsic model classes:

 1. Control voltages: VDB, VGB, VSB

 2. Output currents: IDS, IDB, ISB

 3. Output charges: QD, QG, QS

 4. csOutPorts, controlPorts and qsOutPorts may be overwritten here

 5. linearVCCS and linearVCQS handled by this module

 6. Intrinsic device name must not include ``mos`` and must end with
    ``_i``. Extrinsic device name is: ``mos<intrinsic minus _i>``

 7. Define absolute nominal temperature in param_init(): self._Tn

 8. Noise ports/model not implemented here yet.

Usage::

    import mosExt

    EKV = mosExt.extrinsic_mos(IntEKV)

"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
import cppaddev as ad
from diode import Junction

def extrinsic_mos(IMOS):
    class ExtMOS(IMOS):
        """
    Extrinsic Silicon MOSFET 
    ------------------------

    Netlist examples::
    
        mos<type>:m1 2 3 4 gnd w=10e-6 l=1e-6 type = n as=4e-12 ps=8e=12
        mos<type>:m2 4 5 6 6 w=30e-6 l=1e-6 type = p pd=8u ps=16u 
    
    Extrinsic Internal Topology
    +++++++++++++++++++++++++++

    The model adds the following to the intrinsic model (for NMOS)::

                                     o D (0)
                                     |
                                     \ 
                      Cgdo           / Rd       Drain/source area plus
                                     \          sidewall model
                       ||            |-----------,-----,
                ,------||------------|           |     |   
                |      ||            |         ----- ----- 
                |                ||---         -----  / \  
                |                ||              |   -----
      G (1) o---+----------------||<-------------+-----+------o B (3)
                |                ||              |   -----
                |                ||---         -----  \ / 
                |      ||            |         ----- -----
                `------||------------|           |     |
                       ||            |-----------'-----'
                                     \ 
                      Cgso           / Rs 
                                     \ 
                                     |
                                     o S (2)

        """ 
        # Additional documentation
        extraDoc = """
    Intrinsic model
    +++++++++++++++

    See {0} intrinsic model documentation.
        """.format(IMOS.devType)

        # devtype is the 'model' name: remove the '_i' from intrinsic name
        devType = 'mos' + IMOS.devType.split('_i')[0]

        paramDict = dict(
            IMOS.paramDict.items(),
            cgdo = ('Gate-drain overlap capacitance per meter channel width',
                    'F/m', float, 0.),
            cgso = ('Gate-source overlap capacitance per meter channel width',
                    'F/m', float, 0.),
            cgbo = ('Gate-bulk overlap capacitance per meter channel length',
                    'F/m', float, 0.),
            rsh = ('Drain and source diffusion sheet resistance', 
                   'Ohm/square', float, 0.),
            js = ('Source drain junction current density', 'A/m^2', float, 0.),
            pb = ('Built in potential of source drain junction', 
                  'V', float, .8),
            mj = ('Grading coefficient of source drain junction', 
                  '', float, .5),
            pbsw = ('Built in potential of source, drain junction sidewall',
                    'V', float, .8),
            mjsw = ('Grading coefficient of source drain junction sidewall',
                    '', float, .33),
            cj = ('Source drain junction capacitance per unit area',
                  'F/m^2', float, 0.),
            cjsw = (
                'Source drain junction sidewall capacitance per unit length', 
                'F/m', float, 0.),
            ad = ('Drain area', 'm^2', float, 0.),
            asrc = ('Source area', 'm^2', float, 0.),
            pd = ('Drain perimeter', 'm', float, 0.),
            ps = ('Source perimeter', 'm', float, 0.),
            nrd = ('Number of squares in drain', 'squares', float, 1.),
            nrs = ('Number of squares in source', 'squares', float, 1.),
            fc = ('Coefficient for forward-bias depletion capacitances', ' ', 
                  float, .5),
            xti = ('Junction saturation current temperature exponent', '', 
                   float, 3.),
            eg0 = ('Energy bandgap', 'eV', float, 1.11)
            )

        def __init__(self, instanceName):
            IMOS.__init__(self, instanceName)
            # Create one junction for each (nonlinear) capacitance
            self.dj = Junction()
            self.djsw = Junction()
            self.sj = Junction()
            self.sjsw = Junction()
            self.__doc__ += IMOS.__doc__

        def process_params(self, thermal = False):
            # Remove tape if present
            ad.delete_tape(self)

            # Remove internal terminals (there should be none created
            # by intrinsic model)
            self.clean_internal_terms()
            # Tell autothermal to re-generate thermal ports
            self.__addThermalPorts = True

            # By default drain and source are terminals 0 and 2
            self.__di = 0
            self.__si = 2

            # Resistances
            extraVCCS = list()
            if self.rsh:
                if self.nrd:
                    # Drain resistor
                    self.__di = self.add_internal_term('di', 'V')
                    extraVCCS += [((0, self.__di), (0, self.__di), 
                                   1. / self.rsh / self.nrd)]
                if self.nrs:
                    # Source resistor
                    self.__si = self.add_internal_term('si', 'V')
                    extraVCCS += [((2, self.__si), (2, self.__si), 
                                   1. / self.rsh / self.nrs)]
            # Linear capacitances
            extraVCQS = list()
            if self.cgdo:
                # Gate-drain ovelrlap cap
                extraVCQS += [((1, self.__di), (1, self.__di),
                               self.cgdo * self.w)]
            if self.cgso:
                # Gate-source ovelrlap cap
                extraVCQS += [((1, self.__si), (1, self.__si),
                               self.cgso * self.w)]
            if self.cgbo:
                # Gate-bulk ovelrlap cap
                extraVCQS += [((1, 3), (1, 3),
                               self.cgdo * self.l)]

            # Add extra linear resistors/caps (if any)
            self.linearVCCS = extraVCCS
            self.linearVCQS = extraVCQS

            # Override nonlinear port specs if needed
            if extraVCCS:
                # Ids, Idb, Isb
                self.csOutPorts = [(self.__di, self.__si), (self.__di, 3), 
                                   (self.__si, 3)]
                # Controling voltages are DB, GB and SB
                self.controlPorts = [(self.__di, 3), (1, 3), (self.__si, 3)]
                # One charge source connected to each D, G, S
                self.qsOutPorts = [(self.__di, 3), (1, 3), (self.__si, 3)]

            # Process parameters from intrinsic device: drain and
            # source terminals are already substituted.
            IMOS.process_params(self)

            self.egapn = self.eg0 - .000702 * (self._Tn**2) \
                / (self._Tn + 1108.)
            # Initialize variables in junctions
            if self.ad:
                self.dj.process_params(isat = self.js * self.ad, 
                                       cj0 = self.cj * self.ad, 
                                       vj = self.pb, m = self.mj, 
                                       n = 1., fc = self.fc, 
                                       xti = self.xti, eg0 = self.eg0, 
                                       TnomAbs = self._Tn)
            if self.asrc:
                self.sj.process_params(isat = self.js * self.asrc, 
                                       cj0 = self.cj * self.asrc, 
                                       vj = self.pb, m = self.mj, 
                                       n = 1., fc = self.fc, 
                                       xti = self.xti, eg0 = self.eg0, 
                                       TnomAbs = self._Tn)
            if self.pd:
                self.djsw.process_params(isat = 0., 
                                         cj0 = self.cjsw * self.pd, 
                                         vj = self.pbsw, m = self.mjsw, 
                                         n = 1., fc = self.fc, 
                                         xti = self.xti, eg0 = self.eg0, 
                                         TnomAbs = self._Tn)
            if self.ps:
                self.sjsw.process_params(isat = 0., 
                                         cj0 = self.cjsw * self.ps, 
                                         vj = self.pbsw, m = self.mjsw, 
                                         n = 1., fc = self.fc, 
                                         xti = self.xti, eg0 = self.eg0, 
                                         TnomAbs = self._Tn)
            # Adjust temperature
            self.set_temp_vars(self.temp)
            
        def set_temp_vars(self, temp):
            """
            Calculate temperature-dependent variables, given temp in deg. C
            """
            # Remove tape if present
            ad.delete_tape(self)
            # First calculate variables from base class
            IMOS.set_temp_vars(self, temp)
            # Absolute temperature
            Tabs = temp + const.T0
            # Temperature-adjusted egap
            egap_t = self.eg0 - .000702 * (Tabs**2) / (Tabs + 1108.)
            # Adjust junction temperatures
            if self.ad:
                self.dj.set_temp_vars(Tabs, self._Tn, self._Vt, 
                                      self.egapn, egap_t)
            if self.pd:
                self.djsw.set_temp_vars(Tabs, self._Tn, self._Vt, 
                                        self.egapn, egap_t)
            if self.asrc:
                self.sj.set_temp_vars(Tabs, self._Tn, self._Vt, 
                                      self.egapn, egap_t)
            if self.ps:
                self.sjsw.set_temp_vars(Tabs, self._Tn, self._Vt, 
                                        self.egapn, egap_t)
            

        def eval_cqs(self, vPort, saveOP = False):
            """
            vPort is a vector with control voltages 
            """
            # calculate currents and charges in base class
            if saveOP:
                (iVec, qVec, opVec) = IMOS.eval_cqs(self, vPort, saveOP)
            else:
                (iVec, qVec) = IMOS.eval_cqs(self, vPort)
            # Add contribution drain diode
            v1 = -vPort[0] * self._tf
            if self.ad:
                # substract to idb
                iVec[1] -= self.dj.get_id(v1) * self._tf
                if self.cjs:
                    # substract to qd
                    qVec[0] -= self.dj.get_qd(v1) * self._tf
            if self.pd and self.cjsw:
                qVec[0] -= self.djsw.get_qd(v1) * self._tf
            # Add contribution source diode
            v1 = -vPort[2] * self._tf
            if self.asrc:
                # substract to isb
                iVec[2] -= self.sj.get_id(v1) * self._tf
                if self.cjs:
                    # substract to qs
                    qVec[2] -= self.sj.get_qd(v1) * self._tf
            if self.ps and self.cjsw:
                qVec[2] -= self.sjsw.get_qd(v1) * self._tf

            if saveOP:
                return (iVec, qVec, opVec)
            else:
                return (iVec, qVec)

    
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
            pout = IMOS.power(self, vPort, currV)
            if self.__addCBjtn:
                # add power from substrate junction and RE, RC
                pout += vPort[self.__bccn] * currV[self.__bcon]

            return pout


    # Return template class
    return ExtMOS

#------------------------------------------------------------------------


