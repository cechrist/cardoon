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

def solve(x0, sV, obj):
    """
    Attempt solving circuit equations using several strategies

    x0: initial guess
    sV: source vector
    obj: object that provides the following attribute::

        obj.convergence_helpers          # list of functions that can be used
                                         # to solve equations

    Example of helper functions:

        obj.solve_simple(x0, sV)
        obj.solve_homotopy_gmin(x0, sV)
        obj.solve_homotopy_source(x0, sV)

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    for algorithm in obj.convergence_helpers:
        if algorithm == None:
            raise NoConvergenceError(
                'Giving up. No convergence with any method')
        else:
            if algorithm.__doc__:
                print('Trying ' + algorithm.__doc__)
            try:
                (x, res, iterations) = algorithm(x0, sV)
            except NoConvergenceError as ce:
                print(str(ce))
            else:
                break

    return (x, res, iterations)



def fsolve_Newton(x0, f_Jac_eval, f_eval):
    """
    Solve a multidimensional non-linear equation with Newton-Raphson's method

    In each iteration the linear system:

            J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained (x_{n+1})

    x0: initial guess
    sV: source vector
    f_Jac_eval(x): function that returns (f, Jac) (error function and Jacobian)
    f_eval(x): function that returns error function

    Returns: (x, res, niter)

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    
    ier = 2
    # This overwrites input vector (could use copy(x0))
    x = x0
    for i in xrange(glVar.maxiter):
        (errFunc, Jac) = f_Jac_eval(x)
        try:
            deltax = np.linalg.solve(Jac, errFunc)
        except:
            print('Singular Jacobian')
            # Use pseudo-inverse
            deltax = np.dot(np.linalg.pinv(Jac), errFunc)
        # Do not allow updates greater than 10
        maxDelta = max(abs(deltax))
        # print(maxDelta)
        if maxDelta > glVar.maxdelta:
            deltax *= glVar.maxdelta/maxDelta
        xnew = x - deltax

        # Verify both error function and deltax
        errFunc = abs(f_eval(xnew))
        n1 = np.all(errFunc < (glVar.reltol * max(errFunc) + glVar.abstol)) 
        n2 = np.all(abs(deltax) < (abs(glVar.reltol * np.maximum(x, xnew))
                                   + glVar.abstol))
        if n1 and n2:
            ier = 1
            break

        x = xnew

    res = max(np.linalg.norm(deltax), np.linalg.norm(errFunc))

    if ier == 2:
        raise NoConvergenceError('No convergence. iter = {0} res = {1}'.format(
                i, res))

    return (x, res, i)


