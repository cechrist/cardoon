"""
CppAD device interface. Provides functions to implement evaluation of
nonlinear equations and derivatives using the pycppad library:

http://www.seanet.com/~bradbell/pycppad/index.xml

Usage:
=====

import cppaddev as ad
...

    def process_params(self):
        ...
        # Add the following at the end to make sure tape is re-generated
        ad.delete_tape(self)

    def set_temp_vars(self, temp):
        ...
        # Add the following at the end to make sure tape is re-generated
        ad.delete_tape(self)

    # Use automatic differentiation for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    
Important: use the condassign() function below to replace any if
statements dependent on AD variables in eval_cqs() and
set_temp_vars(). Otherwise, you have to re-tape at each iteration
(call delete_tape() at the end of eval() and eval_and_deriv() in that
device.

If statements dependent on device parameters are OK (if we are not
doing sensitivities).

-------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
import pycppad as ad

def safe_exp(x):
    """
    Same as exp() except when x is greater than threshold. It has
    continuous derivatives.
    """    
    threshold = 50.
    c = np.exp(x)
    d = np.exp(threshold) * (x - threshold + 1.)
    if type(x) == ad.cppad_.a_float:
        return ad.condexp_lt(x, ad.ad(threshold), c, d)
    else:
        if x < threshold:
            return c
        else:
            return d

def condassign(b, c, d):
    """
    Returns the result of (b>0)? c : d 

    Use instead of if statements to force taping both c and d
    expressions. Uses type(b) to decide if ad.condexp_gt() should be
    called, so make sure this an ad type if you are generating a tape.
    """
    if type(b) == ad.cppad_.a_float:
        try:
            return ad.condexp_gt(b, ad.ad(0.), c, d)
        except:
            # Had to catch all because could not find exception name
            if type(c) != ad.cppad_.a_float:
                c = ad.ad(c)
            if type(d) != ad.cppad_.a_float:
                d = ad.ad(d)
            return ad.condexp_gt(b, ad.ad(0.), c, d)
    else:
        if b > 0.:
            return c
        else:
            return d


def delete_tape(dev):
    """
    Delete CppAD tape, if any.

    Call this function (in process_params) to make sure tape is
    re-generated
    """
    try:
        del(dev._func)
        del(dev._opfunc)
    except AttributeError:
        # do nothing if tape does not exist
        pass

def create_tape(dev, vPort):
    """
    Generate main CppAD tape

    Normally there is no need to call this function manually as tapes
    are generated as needed.
    """
    #import pdb; pdb.set_trace()
    assert dev.isNonlinear
    # Create derivative vector
    a_vPort = ad.independent(vPort)
    # perform actual calculation 
    (i_out, q_out) = dev.eval_cqs(a_vPort)
    # Concatenate vectors as we want only one tape to be generated
    a_out = np.concatenate((i_out, q_out), axis=0)
    # Save main function tape
    dev._func = ad.adfun(a_vPort, a_out)
    # optimize main function tape
    dev._func.optimize()


def create_OP_tape(dev, vPort):
    """
    Generate operating point CppAD tape

    Normally there is no need to call this function manually as tapes
    are generated as needed.
    """
    assert dev.isNonlinear
    a_vPort = ad.independent(vPort)
    (i_out, q_out, a_opvars) = dev.eval_cqs(a_vPort, saveOP=True)
    # Save operating point variable tape
    dev._opfunc = ad.adfun(a_vPort, a_opvars)
    # optimize tape (if needed uncomment)
    # dev._opfunc.optimize()


def eval_and_deriv(dev, vPort):
    """
    Evaluates current and charge sources of a nonlinear device. 

    vPort is a numpy vector with input voltages
    
    Returns a tuple with one vector for currents and charges and
    another for the jacobian.
    """
    #import pdb; pdb.set_trace()
    try:
        fout = dev._func.forward(0, vPort)
    except AttributeError:
        create_tape(dev, vPort)
        fout = dev._func.forward(0, vPort)
    jac = dev._func.jacobian(vPort)

    return (fout, jac)


def get_op_vars(dev, vPort):
    """
    Evaluates OP variables of a nonlinear device.

    vPort is a numpy vector with input voltages
    
    Returns OP variables
    """
    #import pdb; pdb.set_trace()
    try:
        opvars = dev._opfunc.forward(0, vPort)
    except AttributeError:
        create_OP_tape(dev, vPort)
        opvars = dev._opfunc.forward(0, vPort)

    return opvars


def eval(dev, vPort):
    """
    Evaluates current and charge sources of a nonlinear device. 

    vPort is a numpy vector with input voltages
    """
    try:
        return dev._func.forward(0, vPort)
    except AttributeError:
        create_tape(dev, vPort)
        return dev._func.forward(0, vPort)

