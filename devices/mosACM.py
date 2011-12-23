"""
Simplified intrinsic ACM MOSFET model

Needs more work: charge calculation, noise equations, etc.

---------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    """
    Implements a simplified ACM MOSFET model. 

    Only (some) DC equations are considered for now.
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
    """

    devType = "mosacm"
    paramDict = dict(
        cir.Element.tempItem,
        w = ('Channel width', 'm', float, 10e-6),
        l = ('Channel length', 'm', float, 10e-6),
        vt0 = ('Threshold Voltage', 'V', float, 0.532),
        phi = ('Surface Potential', 'V', float, 0.55),
        gamma = ('Bulk Threshold Parameter', 'V^(1/2)', float, 0.631),
        kp = ('Transconductance Parameter', 'A/V^2', float, 510.6e-6),
        theta = ('Mobility Saturation Parameter', '1/V', float, 0.814),
        vsat = ('Saturation Velocity', 'm/s', float, 8e4),
        tox = ('Oxide Thickness', 'm', float, 7.5e-9)
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

        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)

        ad.delete_tape(self)


    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables, given temp in deg. C
        """
        # Absolute temperature (note self.temp is in deg. C)
        T = const.T0 + temp
        # Thermal voltage
        self.Vt = const.k * T / const.q
        

    def eval_cqs(self, vPort, saveOP = False):
        """
        Calculates drain current. Input is a vector as follows:
        vPort = [vdb , vgb , vsb]

        saveOP is ignored for now but it will be needed
        """
        # The following formula (11.2.10) seems to work better
        # but it is just an approximation
        vp = pow(np.sqrt(vPort[1]-self.vt0 \
                             + pow(np.sqrt(2.*self.phi)+\
                                       .5*self.gamma, 2))\
                     - .5*self.gamma, 2) - 2.*self.phi
        
        i_f = solve_ifr(vPort[2],vp, self.Vt)
        i_r = solve_ifr(vPort[0],vp, self.Vt)
        
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
        idrain = IS*(i_f - i_r) \
            / (1. + chi * np.abs(np.sqrt(1.+i_f)-np.sqrt(1.+i_r)))
        # Return numpy array with one element per current source.
        return (np.array([idrain]), np.array([]))


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
        # First we need the Jacobian
        (outV, jac) = self.eval_and_deriv(vPort)

        self.OP = dict(
            VD = vPort[0],
            VG = vPort[1],
            VS = vPort[2],
            IDS = outV[0]
            )
        return self.OP

    def get_noise(self, f):
        """
        Return noise spectral density at frequency f
        
        Requires a previous call to get_OP() 

        Not implemented yet
        """
        return None



def solve_ifr(vs, vp, Vt):
    """
    Solve equations to get i_f and i_r using relaxation for now
    """
    i_f = 1.
    deltaif = 10.
    c1 = (vp - vs) / Vt
    # counter = 0
    # Try fixed-point with fixed number of iterations for easy taping
    for i in range(30):
        i_fnew = ad.condassign(c1,
                    pow(c1 + 2. - np.log(np.sqrt(1.+i_f)-1.), 2.) - 1.,
                    pow(np.exp(c1-np.sqrt(1.+i_f)+2.) + 1., 2.) - 1.)
        i_f = .5 * i_f + .5 * i_fnew

#     # Try fixed-point approach for now
#     while deltaif > glVar.abstol:
#         if (c1 > 0.):
#             i_fnew = pow(c1 + 2. - np.log(np.sqrt(1.+i_f)-1.), 2.) - 1.
#         else:
#             i_fnew = pow(np.exp(c1-np.sqrt(1.+i_f)+2.) + 1., 2.) - 1.
#         deltaif = np.abs(i_fnew - i_f)
#         i_f = .5 * i_f + .5 * i_fnew
#         # counter += 1
#     # print counter
    return i_f
    

