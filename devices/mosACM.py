"""
:mod:`mosACM` -- Incomplete intrinsic ACM MOSFET model
------------------------------------------------------

.. module:: mosACM
.. moduleauthor:: Carlos Christoffersen

Needs more work: temperature dependence, charge calculation, noise
equations, etc.

This module also includes ACM utility functions
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad
from mosEKV import f_simple

# Utility functions ----------------------------------------------------

def inv_f1(f):
    """
    Approximately solve f^(-1)(i_f(r)) to get i_f and i_r using relaxation 
    """
    # Obtain a good guess first using EKV's simple interpolation function
    i_f = 4. * f_simple(f)
    # Use fixed number of iterations for easy taping
    for i in xrange(30):
        sqi1 = np.sqrt(1. + i_f + 1e-15)
        # Uses different formulation for f>0 and f<0 for convergence
        i_fnew = ad.condassign(f, 
                               (f + 2. - np.log(sqi1 - 1.))**2 - 1.,
                               (np.exp(f - sqi1 + 2.) + 1.)**2 - 1.)
        i_f = .5 * i_f + .5 * i_fnew
    return i_f

def inv_f(f):
    """
    Solve f^(-1)(i_f(r)) to get i_f and i_r using Newton's method

    Uses a fixed number of iterations for easy taping and a different
    formulation for f>0 and f<0 to ensure convergence in few
    iterations. Solution has good accuracy for all values.
    """
    # Obtain a good guess first using EKV's simple interpolation function
    i1 = 4. * f_simple(f)
    i2 = i1
    for counter in xrange(4):
        # (f > 0) => (i1 >= 3.) The 1e-15 term needed to prevent AD
        # library from choking as without it we have log(0) later
        sqi1p = np.sqrt(1. + i1 + 1e-15) 
        sqi11 = sqi1p - 1. 
        fodfp = (sqi1p - 2. + np.log(sqi11) - f) * 2. * sqi11
        i1 = abs(i1 - fodfp)
        # f < 0
        sqi1 = np.sqrt(1. + i2) 
        fodfn = (sqi1 - 1. - np.exp(f - sqi1 + 2.)) * 2. * sqi1 \
            / (np.exp(f - sqi1 + 2.) + 1.)
        i2 = i2 - fodfn

    return ad.condassign(f, i1, i2)


def fifr(i):
    """
    f(i_f(r)) from ACM model
    """
    if i < 1e-8:
        # sqi1 = 1 + i/2
        return -1. + .5 * i + np.log(.5 * i)
    else:
        sqi1 = np.sqrt(1. + i)    
        return sqi1 - 2. + np.log(sqi1 - 1.)

# ------------------------------------------------------------------------

class Device(cir.Element):
    """
    Incomplete Intrinsic ACM MOSFET
    -------------------------------

    Only (some) DC equations are implemented for now. Temperature
    dependence is not complete.  Terminal order: 0 Drain, 1 Gate, 2
    Source, 3 Bulk::

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

        acm_i:m1 2 3 4 gnd w=10e-6 l=1e-6 type = n 
        acm_i:m2 4 5 6 6 w=30e-6 l=1e-6 type = p 

    Internal topology
    +++++++++++++++++

    For now only ids is implemented::

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

    devType = "acm_i"
    paramDict = dict(
        cir.Element.tempItem,
        type = ('N- or P-channel MOS (n or p)', '', str, 'n'),
        w = ('Channel width', 'm', float, 10e-6),
        l = ('Channel length', 'm', float, 10e-6),
        vt0 = ('Threshold Voltage', 'V', float, 0.532),
        phi = ('Surface Potential', 'V', float, 0.55),
        gamma = ('Bulk Threshold Parameter', 'V^(1/2)', float, 0.631),
        kp = ('Transconductance Parameter', 'A/V^2', float, 510.6e-6),
        theta = ('Mobility Saturation Parameter', '1/V', float, 0.814),
        vsat = ('Saturation Velocity', 'm/s', float, 8e4),
        tox = ('Oxide Thickness', 'm', float, 7.5e-9),
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
                '{0}: unrecognized type: {1}. Valid types are "n" or "p"'.format(self.nodeName, self.type))

        # Nominal abs temperature
        self._Tn = const.T0 + self.tnom

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
        self.Vt = const.k * T / const.q
        

    def eval_cqs(self, vPort, saveOP = False):
        """
        Calculates drain current. Input is a vector as follows:
        vPort = [vdb , vgb , vsb]

        If saveOP = True, return normal output vector plus operating
        point variables in tuple: (iVec, qVec, opV)
        """
        # Invert all voltages in case of a P-channel device
        vPort1 = self._tf * vPort
        # The following formula (11.2.10) seems to work better
        # but it is just an approximation
        vp = pow(np.sqrt(vPort1[1]-self.vt0 \
                             + pow(np.sqrt(2.*self.phi)+\
                                       .5*self.gamma, 2))\
                     - .5*self.gamma, 2) - 2.*self.phi
        
        # Normalized currents.  
        i_f = inv_f((vp - vPort1[2]) / self.Vt)
        i_r = inv_f((vp - vPort1[0]) / self.Vt)
        
        # Calculate IS at this point
        # All this not needed if we just want IS
        eox = 34.5e-12
        cox = eox / self.tox
        mu0 = self.kp / cox
        mus = mu0 / (1. + self.theta * np.sqrt(vp + 2.*self.phi))
        chi = self.Vt * mus / self.l / self.vsat
        
        # From Page 464, also (A.2.3.2d)
        n = 1. + self.gamma / 2. / np.sqrt(2.*self.phi + vp)
        
        # Here kp = u0 Cox
        # 1 / (1 + theta*sqrt(vp+2*phi)) accounts for effective mobility
        IS = self.kp / (1. + self.theta*np.sqrt(vp+2.*self.phi))\
            * n * self.Vt * self.Vt / 2. * self.w/self.l
        
        # Get drain current (including vsat effect)
        # 1/(1 + chi * abs(sqrt(1+i_f)-sqrt(1+i_r))) is the velocity
        # saturation factor. It needs absolute value to work with
        # reverse biasing
        idrain = IS * (i_f - i_r) \
            / (1. + chi * np.abs(np.sqrt(1.+i_f)-np.sqrt(1.+i_r))) * self._tf

        iVec, qVec = np.array([idrain]), np.array([])

        if saveOP:
            # Create operating point variables vector
            opV = np.array([vp, i_f, i_r, IS, n])
            return (iVec, qVec, opV)
        else:
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




