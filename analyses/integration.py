"""
:mod:`integration` -- Integration methods for ODEs
--------------------------------------------------

.. moduleauthor:: Carlos Christoffersen

This module implement different integration methods to be used in
transient analysis.
"""

import numpy as np

class BEuler:
    r"""
    Implements Backwards Euler method:

    .. math::

        \dot{q_{n+1}} = (q_{n+1} - q_n) / h

    """
    def set_h(self, h):
        """
        Set time step size to h
        """
        self.a0 = 1. / h

    def f_n1(self, q):
        r"""
        Returns :math:`f_{n-1}(q)`
        """
        return self.a0 * q


class Trapezoidal:
    r"""
    Implements Trapezoidal method:

    .. math::

        \dot{q_{n+1}} = \frac{2}{h} (q_{n+1} - q_n) - \dot{q_n}

    """
    def __init__(self):
        """
        For now assume we always start from equilibrium
        """
        self.qnm1 = None

    def set_h(self, h):
        """
        Set time step size to h
        """
        self.a0 = 2. / h    

    def f_n1(self, q):
        r"""
        Returns :math:`f_{n-1}(q)`
        
        Also stores q vector and calculates :math:`\dot{q}_{n-1}`
        """
        if self.qnm1 == None:
            self.dqnm1 = np.zeros(len(q))
        else:
            # dqnm1 = \dot{q_{n-1}}
            self.dqnm1 = self.a0 * (q - self.qnm1) - self.dqnm1
        # qnm1 = q_{n-1}
        self.qnm1 = np.copy(q)
        return self.a0 * q + self.dqnm1


