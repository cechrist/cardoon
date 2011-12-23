"""
:mod:`extbjt` -- Extrinsing Bipolar Transistor Template Class
-------------------------------------------------------------

.. module:: extbjt
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from diode import Junction
# For automatic differentiation:
import cppaddev as ad


def extrinsic_bjt(IBJT):
    class BJT(IBJT):
        """
        Generic extrinsic bipolar transistor (BJT)

        Netlist examples::
        
            bjt:q1 2 3 4 1 model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m
            svbjt:q1 2 3 4 1 model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m
        
            # Electro-thermal versions
            bjt_t:q2 2 3 5 1 pout gnd model = mypnp
            svbjt_t:q2 2 3 5 1 pout gnd model = mypnp
        
            # Model statement
            .model mypnp bjt_t (type=pnp isat=5e-17 cje=60fF vje=0.83 mje=0.35)

        Internal Topology
        +++++++++++++++++

        Adds RC, RE and a C-Bulk connection to intrinsic bjt models
        (bjti or svbjti)::

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

        If RE or RC are zero the internal nodes (ct, et) are not created.

        Important Note
        ++++++++++++++

        This implementation does not account for the power dissipation
        in RE, RC. Use external thermal resistors if that is needed.
        """
    
        # devtype is the 'model' name: remove the 'i' from intrinsic name
        devType = IBJT.devType[:-1] 
    
        # Number of terminals. If numTerms is set here, the parser knows
        # in advance how many external terminals to expect. By default the
        # parser makes no assumptions and allows any number of connections
        #
        # Add bulk terminal
        numTerms = IBJT.numTerms + 1
        
        # Create electrothermal device
        makeAutoThermal = True

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
    
        def process_params(self):
            # Remove internal terminals
            self.clean_internal_terms()

            extraVCCS = list()
            # Extra control ports needed for power calculation
            extraPorts = list()
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
            # add power from substrate junction and RE, RC
            pout += vPort[self.__bccn] * currV[self.__bcon]
    
            return pout


    # Return template class
    return BJT
