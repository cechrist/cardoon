"""
:mod:`fsolve` -- Nonlinear equation solve functions
---------------------------------------------------

.. module:: fsolve
.. moduleauthor:: Carlos Christoffersen and others
"""

from __future__ import print_function
import logging
import numpy as np
from globalVars import glVar
from analysis import AnalysisError

class NoConvergenceError(Exception):
    pass


def fsolve(x0, f_Jac_eval, f_eval):
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
            logging.warning('Singular Jacobian')
            # Use least-squares approach
            deltax = np.linalg.lstsq(Jac, errFunc)
        # Do not allow updates greater than 10
        maxDelta = max(abs(deltax))
        maxAllowed = 200.
        # print(maxDelta)
        if maxDelta > maxAllowed:
            deltax *= maxAllowed/maxDelta
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


def solve(x0, sV, obj):
    """
    Attempt solving circuit equations using several strategies

    x0: initial guess
    sV: source vector
    obj: object that provides the following methods::

        (iVec, Jac) = obj.get_i_Jac(x)   # Returns current and Jacobian
        iVec = obj.get_i(x)              # Returns current

    This function originally adapted from pycircuit
    (https://github.com/henjo/pycircuit)
    """
    convergence_helpers = [solve_simple, solve_homotopy_gmin, 
                           solve_homotopy_source, 
                           None]

    for algorithm in convergence_helpers:
        if algorithm == None:
            raise AnalysisError('No convergence!')
        else:
            if algorithm.__doc__:
                print('Trying ' + algorithm.__doc__)
            try:
                (x, res, iter) = algorithm(x0, sV, obj)
            except NoConvergenceError as ce:
                logging.warning('Problems encountered: ' + str(ce))
            else:
                break

    return (x, res, iter)


def solve_simple(x0, sV, obj):
    """Simple Newton's method"""
    def f_Jac_eval(x):
        (iVec, Jac) = obj.get_i_Jac(x) 
        return (iVec - sV, Jac)

    def f_eval(x):
        iVec = obj.get_i(x) 
        return iVec - sV

    return fsolve(x0, f_Jac_eval, f_eval)


def solve_homotopy_gmin(x0, sV, obj):
    """Newton's method with gmin stepping"""
    x = x0
    for gmin in (1, 1e-1, 1e-2, 0):
        n_nodes = len(x0)
        Ggmin = gmin * np.eye(n_nodes)

        def f_Jac_eval(x):
            (iVec, Jac) = obj.get_i_Jac(x) 
            return (iVec - sV, Jac + Ggmin)
        
        def f_eval(x):
            iVec = obj.get_i(x) 
            return iVec - sV

    return fsolve(x0, f_Jac_eval, f_eval)


def solve_homotopy_source(x0, sV, obj):
    """Newton's method with source stepping"""
    x = x0
    for lambda_ in (0, 1e-2, 1e-1, 1):
        def f_Jac_eval(x):
            (iVec, Jac) = obj.get_i_Jac(x) 
            return (iVec - lambda_ * sV, Jac)
        
        def f_eval(x):
            iVec = obj.get_i(x) 
            return iVec - lambda_ * sV

    return fsolve(x0, f_Jac_eval, f_eval)


