"""
State-Variable-Based P-N Juction model.

Based on carrot source, in turn based on spice/freeda diode
model. This is intended to model any regular P-N junction such as
Drain-Bulk, diodes, Collector-Bulk, etc. (if a state-variable approach
is desired)

-------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
from globalVars import const, glVar

class Junction:
    """
    State-Variable PN Junction model based on spice equations

    Breakdown and area effects not included
    """
    paramList = (
        ('isat', ('Saturation current', 'A', float, 1e-14)),
        ('n', ('Emission coefficient', ' ', float, 1.)),
        ('fc', ('Coefficient for forward-bias depletion capacitance', ' ', 
                float, .5)),
        ('cj0', ('Zero-bias depletion capacitance', 'F', float, 0.)),
        ('vj', ('Built-in junction potential', 'V', float, 1.)),
        ('m', ('PN junction grading coefficient', ' ', float, .5)),
        ('tt', ('Transit time', 's', float, 0.)),
        ('xti', ('Is temperature exponent', ' ', float, 3.)),
        ('eg0', ('Energy bandgap', 'eV', float, 1.11)),
        ('tnom', ('Nominal temperature', 'C', float, 27.))
    )

    def process_params(self, obj):
        """
        Calculate variables dependent on parameter values only

        obj: instance with parameter values as attributes
        """
        # Set some handy variables
        self._k1 = obj.xti / obj.n
        tnomabs = obj.tnom + const.T0
        self._k2 = const.q * obj.eg0 / obj.n / const.k / tnomabs
        self._k3 = const.q * obj.eg0 / obj.n / const.k
        self._k4 = 1. - obj.m
        self._egapn = obj.eg0 \
            - .000702 * tnomabs * tnomabs / (tnomabs + 1108.)

    def set_temp_vars(self, temp, obj, ad):
        """
        Calculate temperature-dependent variables for temp given in C

        obj: instance with parameter values as attributes
        ad: AD module (or could be helperf module)
        """
        # Absolute temperature (note temp is in deg. C)
        Tabs = const.T0 + temp
        # Thermal voltage
        self._vt = const.k * Tabs / const.q
        # Normalized temp
        tn = Tabs / (obj.tnom + const.T0)
        self._t_is = obj.isat * pow(tn, self._k1) \
            * np.exp(self._k2 - self._k3 / Tabs) 
        self._t_egap = obj.eg0 - .000702 * Tabs * Tabs / (Tabs + 1108.)
        self._t_vj = obj.vj * tn - 3. * self._vt * np.log(tn) \
            - tn * self._egapn + self._t_egap
        self._t_cj0 = obj.cj0 * (1. + obj.m \
                      * (.0004 * (temp - obj.tnom) \
                              + 1. - self._t_vj / obj.vj))
        self._k5 = self._t_vj * self._t_cj0 / self._k4
        self._k6 = self._t_cj0 * pow(1. - obj.fc, - obj.m - 1.)
        self._k7 = ((1. - obj.fc * (1. + obj.m)) * self._t_vj * obj.fc \
                        + .5 * obj.m * self._t_vj * obj.fc * obj.fc)
        # Maximum argument in exponential function (no need to use
        # safe_exp() with this model)
        self._alpha = 1. / obj.n / self._vt
        max_exp_arg = 5e8
        self._svth = np.log(max_exp_arg / self._alpha) / self._alpha
        self._kexp = obj.n * self._vt * max_exp_arg


    def eval_cqs(self, x, obj, ad):
        """
        Models intrinsic junction currrent and charge.  

        x: state variable
        obj: instance with parameter values as attributes
        ad: AD module (or could be helperf module)

        Returns a vector with two elements, one for current and the
        other for charge.
        """
        # Static current
        b = self._svth - x
        c = self._t_is * (np.exp(self._alpha * x) - 1.)
        d = self._t_is * self._kexp * \
            (1. + self._alpha * (x - self._svth)) - self._t_is
        # Define output vector: [id, vd, qd]
        outV = np.zeros((3), dtype=type(b))
        outV[0] = ad.condassign(b, c, d)
        # Diode voltage
        d = self._svth + \
            np.log(1. + self._alpha * (x - self._svth)) / self._alpha
        vj = ad.condassign(b, x, d)
        outV[1] = vj
        # Charges
        if obj.cj0:
            b = obj.fc * self._t_vj - vj
            c = self._k5 * (1. - pow(1. - vj / self._t_vj, self._k4))
            d = self._k6 * ((1. - obj.fc * (1. + obj.m)) 
                            * vj + .5 * obj.m * vj 
                            * vj / self._t_vj - self._k7) + self._k5 \
                            * (1. - pow(1. - obj.fc, self._k4))
            outV[2] = ad.condassign(b, c, d)

        if obj.tt:
            outV[2] += obj.tt * outV[0]

        return outV
