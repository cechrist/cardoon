"""
:mod:`helperfunc` -- General helper function collection
-------------------------------------------------------

.. module:: helperfunc
.. moduleauthor:: Carlos Christoffersen and others
"""

from __future__ import print_function
import logging
import numpy as np
from globalVars import glVar

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

    if ier == 2:
        raise NoConvergenceError('No convergence: norm(deltax) = {0}'.format(
                np.linalg.norm(deltax)))

    return x


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
                x = algorithm(x0, sV, obj)
            except NoConvergenceError as ce:
                logging.warning('Problems encountered: ' + str(ce))
            else:
                break

    return x


def solve_simple(x0, sV, obj):
    """Simple Newton's method"""
    def f_Jac_eval(x):
        (iVec, Jac) = obj.get_i_Jac(x) 
        return (iVec - sV, Jac)

    def f_eval(x):
        iVec = obj.get_i(x) 
        return iVec - sV

    return fsolve(x0, f_Jac_eval, f_eval)


def solve_homotopy_gmin(obj, x0):
    """Newton's method with gmin stepping"""
    x = x0
    for gmin in (1, 1e-1, 1e-2, 0):
        n_nodes = len(x0)
        Ggmin = gmin * np.eye((n_nodes, n_nodes))

        def f_Jac_eval(x):
            (iVec, Jac) = obj.get_i_Jac(x) 
            return (iVec - sV, Jac + Ggmin)
        
        def f_eval(x):
            iVec = obj.get_i(x) 
            return iVec - sV

    return fsolve(x0, f_Jac_eval, f_eval)


def solve_homotopy_source(obj, x0):
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


