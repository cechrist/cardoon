"""
Do not use: needs additional work to provide same functionality as
cppaddev.py

Adol-C device interface. Provides functions to implement evaluation of
nonlinear equations and derivatives using the pyadolc library:

https://github.com/b45ch1/pyadolc

cppaddev.py is preferred as the cppad library is generally faster.

Usage: same as cppaddev.py

-------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
import adolc as ad

def condassign(b, c, d):
    """
    Returns the result of (b>0)? c : d 

    Use instead of if statements to force taping both c and d
    expressions.
    """
    if type(b) == ad._adolc.adub:
        ad_a = ad.adouble(0.)
        ad.condassign(ad_a, ad.adouble(b), ad.adouble(c), ad.adouble(d))
        return ad_a
    else:
        if b > 0.:
            return c
        else:
            return d


def delete_tape(dev):
    """
    Forces deletion of Adol-C tape

    Call this function (in process_params) to make sure tape is
    re-generated
    """
    try:
        del(dev._tag)
    except AttributeError:
        # do nothing if tape does not exist
        pass

def create_tape(dev, vPort):
    """
    Generate Adol-C tape

    Normally there is no need to manually call this.
    """
    assert dev.isNonlinear
    try:
        tag = dev._tag
    except AttributeError:
        tag = dev.adolcID
        dev._tag = tag
    ad.trace_on(tag)
    # Create derivative vector
    a_vPort = ad.adouble(vPort)
    ad.independent(a_vPort)
    # perform actual calculation (for now re-tape every time)
    a_out = dev.eval_cqs(a_vPort)
    ad.dependent(a_out)
    ad.trace_off()

def eval_and_deriv(dev, vPort):
    """
    Evaluates current and charge sources of a nonlinear device. 

    vPort is a numpy vector with input voltages
    
    Returns a tuple with one vector for currents and charges and
    another for the jacobian.
    """
    try:
        tag = dev._tag
    except AttributeError:
        create_tape(dev, vPort)
        tag = dev._tag
    
    fout = ad.function(tag, vPort)
    jac = ad.jacobian(tag, vPort)
    
    return (fout, jac)

def eval(dev, vPort):
    """
    Evaluates current and charge sources of a nonlinear device. 

    vPort is a numpy vector with input voltages
    """
    try:
        tag = dev._tag
    except AttributeError:
        create_tape(dev, vPort)
        tag = dev._tag

    return ad.function(tag, vPort)

