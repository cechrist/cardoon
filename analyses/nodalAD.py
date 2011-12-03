"""
:mod:`nodalAD` -- Nodal Approach: full AD 
-----------------------------------------

.. module:: nodal
.. moduleauthor:: Carlos Christoffersen

************** This is experimental/incomplete ****************

This module generates a CPPAD tape for the whole circuit. Use as an
experimental replacement for nodal.py for illustration purposes
only. In principle the regular approach seems faster.

"""

import numpy as np
import pycppad as ad
from nodal import DCNodal, set_xin, set_i

#------------------------------------------------------------------------

class DCNodalAD(DCNodal):
    """
    Replacement for DCNodal using CPPAD to evaluate currents and Jacobian
    """

    def __init__(self, ckt):

        DCNodal.__init__(self, ckt)
            
    def _create_tape(self, xVec):
        """
        Generate main CppAD tape
    
        Normally there is no need to call this function manually as tapes
        are generated as needed.
        """
        # Create derivative vector
        a_xVec = ad.independent(xVec)
        # perform actual calculation 
        # Linear contribution 
        a_out = np.dot(self.G, a_xVec)
        xin = np.zeros(len(xVec), dtype = type(a_out[0]))
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from a_xVec
            xin[:len(elem.controlPorts)] = 0.
            #import pdb; pdb.set_trace()
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, a_xVec)
            (outV, qVec) = elem.eval_cqs(xin)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(a_out, elem.nD_cpos, elem.nD_cneg, outV)
        # Save main function tape
        self._func = ad.adfun(a_xVec, a_out)
        # optimize main function tape
        self._func.optimize()


    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G xVec + i(xVec)

        xVec: input vector of nodal voltages. 
        iVec: output vector of currents
        """
        # Erase arrays
        self.iVec[:] = 0.
        # Linear contribution
#        self.iVec += np.dot(self.G, xVec)
        try:
            self.iVec += self._func.forward(0, xVec)
        except AttributeError:
            self._create_tape(xVec)
            self.iVec += self._func.forward(0, xVec)
        return self.iVec


    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns: (iVec, Jac)

            iVec = G xVec + i(xVec)
            Jac = G + Ji(xVec)

        xVec: input vector of nodal voltages. 
        iVec: output vector of currents
        Jac: system Jacobian
        """
        # Erase arrays
        self.iVec[:] = 0.
        self.Jac[:] = 0.
        # Linear contribution
#        self.iVec += np.dot(self.G, xVec)
#        self.Jac += self.G
        try:
            self.iVec = self._func.forward(0, xVec)
        except AttributeError:
            self._create_tape(xVec)
            self.iVec = self._func.forward(0, xVec)

        self.Jac += self._func.jacobian(xVec)

        return (self.iVec, self.Jac)

