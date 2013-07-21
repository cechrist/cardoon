"""
Miscellany helper functions, mostly needed if an AD package is not
being used (normally, use cppaddev.py instead)
"""

import numpy as np

def safe_exp(x):
    """
    Same as exp() except when x is greater than threshold. It has
    continuous derivatives.
    """    
    threshold = 50.
    if x < threshold:
        return np.exp(x)
    else:
        return np.exp(threshold) * (x - threshold + 1.)


def condassign(a, b, c, d):
    """
    a = (b>0)? c : d

    Used in case that we want to evaluate code written for AD packages
    such as Adol-C or CPPAD.
    """
    if b > 0.:
        return c
    else:
        return d



