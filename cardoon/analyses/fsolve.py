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
from cardoon.globalVars import glVar

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

      ``x0``: initial guess vector

      ``get_deltax(x)``: function that returns the :math:`-J(x_n)^{-1}
      F(x_n)` vector

      ``f_eval(x)``: function that returns error function vector

    Convergence parameters are taken from the global options
    (:doc:`global_vars`). Relative and absolute tolerances are checked
    independently for each variable. 

    Returns the following tuple: ``(x, res, niter, success)``

      ``x``: estimated solution vector at the last iteration
      
      ``res``: maximum residual from all variables

      ``niter``: number of iterations 

      ``success``: True if method converged to specified tolerance,
      False otherwise

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    
    success = False
    # This overwrites input vector (could use copy(x0))
    x = x0
    for i in xrange(glVar.maxiter):

        # Actually this is -deltax as it must be substracted from x
        deltax = get_deltax(x)
        
        # Do not allow updates greater than glVar.maxdelta
        maxDelta = max(abs(deltax))
        if maxDelta > glVar.maxdelta:
            deltax *= glVar.maxdelta/maxDelta
        xnew = x + deltax

        # Check if deltax is small
        n1 = np.all(abs(deltax) < (glVar.reltol * np.maximum(abs(x), abs(xnew))
                                   + glVar.abstol))
        res = max(abs(deltax))
        if n1:
            if glVar.errfunc:
                # Optional: check if error function is small. Only
                # check for absolute error since nominal value is zero
                errFunc = max(abs(f_eval(xnew)))
                n2 = errFunc < glVar.abstol
                res = max(errFunc, res)
            else:
                # Do not check error function to save time
                n2 = True
        x = xnew
        if n1 and n2:
            success = True
            break


    return (x, res, i+1, success)


