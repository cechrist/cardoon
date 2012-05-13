"""
:mod:`fsolve` -- Nonlinear equation solve functions
---------------------------------------------------

.. module:: fsolve
.. moduleauthor:: Carlos Christoffersen and others

This module provides the main function that tries different equation
solving strategies provided by an object passed as an argument.

A simple Newton method routine is also provided here. In the future we
could implement here additional nonlinear equation solving routines.

"""

from __future__ import print_function
import numpy as np
from globalVars import glVar

class NoConvergenceError(Exception):
    pass

def solve(x0, sV, convergence_helpers):
    """
    Attempt solving circuit equations using several strategies

    x0: initial guess

    sV: source vector

    convergence_helpers: list of functions that can be used to solve equations

    Example of helper functions::

        solve_simple(x0, sV)
        solve_homotopy_gmin(x0, sV)
        solve_homotopy_source(x0, sV)

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    for algorithm in convergence_helpers:
        if algorithm == None:
            raise NoConvergenceError(
                'Giving up. No convergence with any method')
        else:
            if algorithm.__doc__:
                print('\nTrying ' + algorithm.__doc__)
            try:
                (x, res, iterations) = algorithm(x0, sV)
            except NoConvergenceError as ce:
                print(str(ce))
            else:
                break

    return (x, res, iterations)


def fsolve_Newton(x0, get_deltax, f_eval):
    r"""
    Solves a multidimensional non-linear equation with
    Newton-Raphson's method.  In each iteration the linear system:
    
    .. math::

            J(x_n)(x_{n+1}-x_n) + F(x_n) = 0 \; ,

    is solved and a new value for :math:`x_{n+1}` is obtained. Arguments:

    x0: initial guess

    get_deltax(x): function that returns :math:`J(x_n)^{-1} F(x_n)`

    f_eval(x): function that returns error function

    Returns (x, res, niter)

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    
    ier = 2
    # This overwrites input vector (could use copy(x0))
    x = x0
    for i in xrange(glVar.maxiter):

        deltax = get_deltax(x)
        
        # Do not allow updates greater than glVar.maxdelta
        maxDelta = max(abs(deltax))
        if maxDelta > glVar.maxdelta:
            deltax *= glVar.maxdelta/maxDelta
        xnew = x - deltax

        if glVar.errfunc:
            # Check if error function is small
            errFunc = abs(f_eval(xnew))
            n1 = np.all(errFunc < (glVar.reltol * max(errFunc) + glVar.abstol))
            res1 = np.linalg.norm(errFunc)
        else:
            # Do not check error function to save time
            n1 = True
            res1 = 0.
        # Check if deltax is small
        n2 = np.all(abs(deltax) < (abs(glVar.reltol * np.maximum(x, xnew))
                                   + glVar.abstol))
        x = xnew
        if n1 and n2:
            ier = 1
            break

    res = max(res1, np.linalg.norm(deltax))

    if ier == 2:
        success = False
    else:
        success = True

    return (x, res, i+1, success)


