"""
:mod:`mosACM` -- Simplified intrinsic ACM MOSFET model
------------------------------------------------------

.. module:: mosACM
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad
from mosACM import inv_f, inv_f1

class Device(cir.Element):
    """
    Simplified ACM Intrinsic MOSFET
    -------------------------------

    This model uses the simple equations for hand analysis. Only DC
    equations (with temperature dependence) included for now. 

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

        acms_i:m1 2 3 4 gnd w=10e-6 l=1e-6 type = n 
        acms_i:m2 4 5 6 6 w=30e-6 l=1e-6 type = p 

    Internal topology
    +++++++++++++++++

    Only ids is implemented. In the future charges will be added::

                           ,--o 0 (D)
                           |
                           |
                           |
                           |       
                          /|\       
          (G) 1 o-       ( | ) ids(VD, VG, VS, VB)
                          \V/      
                           |       
                           |
                           |
                           |
          (B) 3 o-         `--o 2 (S)
                  


    """
    # Device category
    category = "Semiconductor devices"

    devType = "acms_i"
    paramDict = dict(
        cir.Element.tempItem,
        type = ('N- or P-channel MOS (n or p)', '', str, 'n'),
        w = ('Channel width', 'm', float, 10e-6),
        l = ('Channel length', 'm', float, 10e-6),
        vth = ('Threshold Voltage', 'V', float, 0.5),
        isq = ('Sheet normalization current', 'A/V^2', float, 100e-9),
        n = ('Subthreshold slope factor', 'F/m^2', float, 1.3),
        cox = ('Gate oxide capacitance per area', 'F/m^2', float, 0.7e-3),
        tcv = ('Threshold voltage temperature coefficient', 'V/K', 
               float, 0.001),
        bex = ('Mobility temperature exponent', '', float, -1.5),
        tnom = ('Nominal temperature of model parameters', 'C', float, 27.)
        )
    
    numTerms = 4

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
                '{0}: unrecognized type: {1}. Valid types are "n" or "p"'.format(self.instanceName, self.type))

        # Make sure vth is positive for calculations
        self._vth = abs(self.vth)
        self._tcv = np.abs(self.tcv)

        # Nominal abs temperature
        self._Tn = const.T0 + self.tnom
        # Nominal Thermal voltage
        self._Vtn = const.k * self._Tn / const.q
        self._mu = self.isq * 2. / (self.n * self.cox * self._Vtn**2)
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
        # threshold voltage
        self._vthT = self._vth - self._tcv * (T - self._Tn)
        # Mobility
        muT = self._mu * (T / self._Tn)**self.bex
        # IS (specific current)
        s = self.w / self.l
        self._IS = .5 * s * self.n * muT * self.cox * self._Vt**2


    def eval_cqs(self, vPort, getOP = False):
        """
        Calculates drain current. Input is a vector as follows:
        vPort = [vdb , vgb , vsb]

        If getOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # Invert all voltages in case of a P-channel device
        vPort1 = self._tf * vPort

        vp = (vPort1[1] - self._vthT) / self.n
        
        # Normalized currents
        i_f = inv_f((vp - vPort1[2]) / self._Vt)
        i_r = inv_f((vp - vPort1[0]) / self._Vt)
        
        idrain = self._IS * (i_f - i_r) 

        iVec, qVec = np.array([idrain]), np.array([])

        if getOP:
            # Create operating point variables dictionary
            return {'VD': vPort[0],
                    'VG': vPort[1],
                    'VS': vPort[2],
                    'IDS': idrain,
                    'Vth': self._vthT,
                    'IS': self._IS, 
                    'vp': vp,
                    'i_f': i_f,
                    'i_r': i_r}
        else:
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
        # In the future add transconductances, etc.
        return self.eval_cqs(vPort, True)

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None



#     # Try fixed-point approach for now
#     while deltaif > glVar.abstol:
#         if (f > 0.):
#             i_fnew = pow(f + 2. - np.log(sqi1-1.), 2.) - 1.
#         else:
#             i_fnew = pow(np.exp(f-sqi1+2.) + 1., 2.) - 1.
#         deltaif = np.abs(i_fnew - i_f)
#         i_f = .5 * i_f + .5 * i_fnew
#         # counter += 1
#     # print counter
    

# Saved here in case some of the expressions are needed later
#
#   def inv_f(f):
#       """
#       Solve f^(-1)(i_f(r)) to get i_f and i_r using Newton's method
#   
#       Uses inv_f1() to get a good starting guess. Unfortunately not used
#       as it seems to cause convergence problems.
#       """
#       i_f = inv_f1(f)
#       # Use fixed number of iterations for easy taping. Solution has
#       # good accuracy for all values.
#       for counter in xrange(4):
#           sqi1 = np.sqrt(1. + i_f)
#           sqi11 = abs(sqi1 - 1.)
#           # f > 0
#           fodfp = (sqi1 - 2. + np.log(sqi11) - f) * 2. * sqi11
#   #        fodfp = ad.condassign(
#   #            f + 4., 
#   #            (sqi1 - 2. + np.log(sqi11) - f) * 2. * sqi11,
#   #            (-1. + .5 * i_f + np.log(.5 * i_f) - f) * i_f)
#           # f < 0
#           fodfn = (sqi1 - 1. - np.exp(f - sqi1 + 2.)) * 2. * sqi1 \
#               / (np.exp(f - sqi1 + 2.) + 1.)
#   #        fodfn = ad.condassign(
#   #            80 - f,
#   #            (sqi1 - 1. - np.exp(f - sqi1 + 2.)) * 2. * sqi1 /
#   #            (np.exp(f - sqi1 + 2.) + 1.),
#   #            2. * sqi1)
#           # Uses different formulation for f>0 and f<0 to ensure
#           # convergence in few iterations
#           ifnew = ad.condassign(f, i_f - fodfp, i_f - fodfn)
#   #        ifnew = i_f - fodfn
#   #        deltai = abs(ifnew - i_f)
#   #        print(deltai)
#           i_f = ad.condassign(ifnew, ifnew, 1e-300)
#   
#       return i_f
#   

