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
    def init(self, h, q):
        """
        Initialize for integration

        Set time step size to h. q is ignored as it is not needed
        """
        self.h = h
        self.a0 = 1. / h
        self.qn1 = np.zeros_like(q)

    def set_h(self, h):
        """
        Change time step size to h
        """
        self.h = h
        self.a0 = 1. / h

    def accept(self, q):
        """
        Accept q as last valid value
        """
        self.qn1[:] = q

    def f_n1(self):
        r"""
        Returns :math:`f_{n-1}(q)`
        """
        return self.a0 * self.qn1


class Trapezoidal:
    r"""
    Implements Trapezoidal method:

    .. math::

        \dot{q_{n+1}} = \frac{2}{h} (q_{n+1} - q_n) - \dot{q_n}

    """
    def init(self, h, q, dq = None):
        """
        Initialize for integration

        Set time step size to h, previous charge to q previous
        derivative to dq
        """
        self.h = h
        self.a0 = 2. / h    
        self.qnm1 = np.copy(q)
        if dq == None:
            self.dqnm1 = np.zeros_like(q)
        else:
            self.dqnm1 = np.copy(dq)

    def set_h(self, h):
        """
        Change time step size to h in next call to f_n1
        """
        self.h = h
        self.a0 = 1. / h

    def accept(self, q):
        """
        Accept q as last valid value
        """
        # dqnm1 = \dot{q_{n-1}}
        self.dqnm1 *= -1.
        self.dqnm1 += self.a0 * (q - self.qnm1)
        # qnm1 = q_{n-1}
        self.qnm1[:] = q

    def f_n1(self):
        r"""
        Returns :math:`f_{n-1}(q)`
        """
        return self.a0 * self.qnm1 + self.dqnm1


