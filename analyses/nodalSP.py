"""
:mod:`nodalSP` -- Nodal Analysis using pysparse
-----------------------------------------------

.. module:: nodalSP
.. moduleauthor:: Carlos Christoffersen

This module contains basic classes/functions for nodal analysis. These
are part of the standard netlist analyses but they can also be used
independently.

This implementation uses UMFPACK
(http://www.cise.ufl.edu/research/sparse/umfpack/) to solve sparse
matrix linear systems with the pysparse interface
(http://pysparse.sourceforge.net/). For medium-size and large circuits
it is much more efficient than the dense implementation but still a
lot could be gained with finer control over matrix factorization. At
this time the main Jacobian is (almost) created and factored from
scratch at every iteration.

Many of the functions are directly imported from the ``nodal`` module
to avoid redundancy. Some of these should be optimized for better
performance.

"""

from __future__ import print_function
import numpy as np
import pysparse
from fsolve import fsolve_Newton, NoConvergenceError
from integration import BEuler
from nodal import set_quad, set_xin, set_i, set_Jac, make_nodal_circuit,\
    restore_RCnumbers, process_nodal_element,_NLFunction

#-------------------------------------------------------------------------
# ****************************** Classes *********************************
#-------------------------------------------------------------------------

class _NLFunctionSP(_NLFunction):
    """
    Nonlinear function interface class using sparse matrices

    This is an abstract class to be use as a base class by other
    nodal-based classes such as DCNodal. Only methods that differ from
    nodal._NLFunction are defined here.
    """

    def __init__(self):
        # List here the functions that can be used to solve equations
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_gmin, 
                                    self.solve_homotopy_source, 
                                    None]

    def _get_deltax(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc
        """
        umf = pysparse.umfpack.factorize(
            Jac, strategy="UMFPACK_STRATEGY_UNSYMMETRIC")
        umf.solve(errFunc, self.deltaxVec)
        return self.deltaxVec


    def solve_homotopy_gmin(self, x0, sV):
        """Newton's method with gmin stepping"""
        idx = np.arange(self.ckt.nD_nterms)
        def f(_lambda):
            gbase = 1e-5
            self.gmin = gbase / _lambda**3 - gbase
            self.val = self.gmin * np.ones(self.ckt.nD_nterms)
        # Add gmin from ground to all external nodes. Assumes all
        # external nodes are sorted first in the vector. This will
        # not work if the terminal order is changed.
        def get_deltax(xvec):
            (iVec, Jac) = self.get_i_Jac(xvec) 
            iVec[idx] += self.gmin * xvec[idx]
            Jac.update_add_at(self.val, idx, idx)
            return self._get_deltax(iVec - sV, Jac)
        def f_eval(xvec):
            iVec = self.get_i(xvec)
            iVec[idx] += self.gmin * xvec[idx]
            return iVec - sV
        (x, res, iterations, success) = \
            self._homotopy(0.5, f, x0, get_deltax, f_eval)
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError('gmin stepping did not converge')


#---------------------------------------------------------------------------

class DCNodal(_NLFunctionSP):
    """
    Calculates the DC part of currents and Jacobian

    Matrices and vectors (G, Jac, s) are allocated here. This is to
    centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt):
        # Init base class
        super(DCNodal, self).__init__()

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        # Allocate matrices/vectors
        # G here is G1 = G + G0 in documentation (allocated in refresh())

        # Jac is (G1 + di/dv) in doc
        self.Jac = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.iVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)
        if hasattr(self.ckt, 'nD_namRClist'):
            # Allocate external currents vector
            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parameter sweeps
        """
        self.Gll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.Gll, *vccs)
        # Frequency-defined elements
        for elem in self.ckt.nD_freqDefinedElem:
            set_Jac(self.Gll, elem.nD_fpos, elem.nD_fneg, 
                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        # Free unused memory
        self.Gll.compress()        
        # Create G in csr form for efficient matrix-vector multiplication
        self.G = self.Gll.to_csr()

            
    def set_ext_currents(self, extIvec):
        """
        Set external currents applied to subcircuit

        extIvec: vector of external currents. Length of this vector
        should be equal to the number of external connections. The sum
        of all currents must be equal to zero (KCL)

        This will fail if not a subcircuit
        """
        # This idea still needs some testing
        assert sum(extIvec[:ncurrents]) == 0
        ncurrents = len(self.ckt.nD_namRClist)
        # Must do the loop in case there are repeated connections
        for val,rcnum in zip(extIvec, self.ckt.nD_namRClist):
            self.extSVec[rcnum] = val

    def get_guess(self):
        """
        Retrieve guesses from vPortGuess in each nonlinear device
        
        Returns a guess vector
        """
        x0 = np.zeros(self.ckt.nD_dimension)
        for elem in self.ckt.nD_nlinElem:
            try:
                # Only add to positive side. This is not the only way
                # and may not work well in some cases but this is a
                # guess anyway
                for i,j in elem.nD_vpos:
                    x0[j] += elem.vPortGuess[i]
            except AttributeError:
                # if vPortGuess not given just leave things unchanged
                pass
        return x0

    def get_source(self):
        """
        Get the source vector considering only the DC source components
        """
        # Erase vector first. 
        try:
            # If subcircuit add external currents
            self.sVec[:] = self.extSVec
        except AttributeError:
            # Not a subcircuit
            self.sVec.fill(0.)
        for elem in self.ckt.nD_sourceDCElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_DCsource()
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        return self.sVec
        

    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G xVec + i(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents
        """
        # Linear contribution
        self.G.matvec(xVec, self.iVec)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            outV = elem.eval(xin)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
        return self.iVec

    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)::

            iVec = G xVec + i(xVec)
            Jac = G + (di/dx)(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        # Erase sparse matrix
        self.Jac.scale(0.)
        # Linear contribution
        self.G.matvec(xVec, self.iVec)
        self.Jac.shift(1., self.Gll)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_Jac(self.Jac, elem.nD_cpos, elem.nD_cneg, 
                    elem.nD_vpos, elem.nD_vneg, outJac)

        return (self.iVec, self.Jac)

    def save_OP(self, xVec):
        """
        Save nodal voltages in terminals and set OP in elements

        The following information is saved:

          * The nodal voltage vector for the circuit (self.xop)

          * The nodal voltage in each terminal (term.nD_vOP)

          * The port voltages in each nonlinear device (elem.nD_xOP)

          * The operating point (OP) information in nonlinear devices

        """
        # Set nodal voltage of reference to zero
        self.ckt.nD_ref.nD_vOP = 0.
        # Save nodal vector
        self.xop = xVec
        for v,term in zip(xVec, self.ckt.nD_termList):
            term.nD_vOP = v
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            # Set OP in element (discard return value)
            elem.nD_xOP = xin
            elem.get_OP(xin)

    
#----------------------------------------------------------------------

class TransientNodal(_NLFunctionSP):
    """
    Keeps track of transient analysis equations. 

    This class only sets up transient equations. Equations are solved
    elsewhere. Circuit history required for numerical integration is
    kept by an integration method class. 

    Matrices and vectors (G, C, JacI, JacQ, s, etc.) are allocated
    here. This is to centralize all allocations and avoid repetitions.

    Requires a nodal-ready Circuit instance (ckt) instance (see
    make_nodal_circuit())
    """

    def __init__(self, ckt, im):
        """
        Arguments:

        ckt: circuit instance

        im: Integration method instance. 

        """
        # Init base class
        super(TransientNodal, self).__init__()

        # Save ckt instance
        self.ckt = ckt
        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref

        self.im = im

        # Allocate matrices/vectors
        # G, C and G'
        self.Gll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        self.Cll = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # iVec = G' x + i'(x)   total current
        self.iVec = np.empty(self.ckt.nD_dimension)
        # System Jacobian: G' + di'/dx
        self.Jac = pysparse.spmatrix.ll_mat(self.ckt.nD_dimension, 
                                            self.ckt.nD_dimension)
        # Total charge: C x + q(x)
        self.qVec = np.empty(self.ckt.nD_dimension)
        # Source vector at current time s(t) 
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.deltaxVec = np.empty(self.ckt.nD_dimension)

#        if hasattr(self.ckt, 'nD_namRClist'):
#            # Allocate external currents vector
#            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        self.Gll.scale(0.)
        self.Cll.scale(0.)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.Gll, *vccs)
            for vccs in elem.nD_linVCQS:
                set_quad(self.Cll, *vccs)
#        # Frequency-defined elements not included for now
#        for elem in self.ckt.nD_freqDefinedElem:
#            set_Jac(self.G, elem.nD_fpos, elem.nD_fneg, 
#                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
        self.Gll.compress()
        self.Cll.compress()
        self.C = self.Cll.to_csr()

    def set_IC(self, h):
        """
        Set initial conditions (ICs)

        h: (initial) time step

        Retrieves ICs from DC operating point info in elements
        """
        # Nodal variables
        xVec = np.zeros(self.ckt.nD_dimension)
        # Get nodal voltages for tstart
        for i,term in enumerate(self.ckt.nD_termList):
            xVec[i] = term.nD_vOP
        # Calculate total charges
        self.update_q(xVec)
        # initialize integration element
        self.im.init(h, self.qVec)
        # Generate Gp 
        self.update_Gp()


    def update_Gp(self):
        """
        Recalculate Gp from im information
        """
        self.Gpll = self.Gll.copy()
        self.Gpll.shift(self.im.a0, self.Cll)
        self.Gp = self.Gpll.to_csr()
        return self.Gp


    def update_q(self, xVec):
        """ 
        Recalculate qVec for a given value of xVec
        """
        # Calculate total q vector
        # Calculate linear charges first
        self.C.matvec(xVec, self.qVec)
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)  
            outV = elem.eval(xin)
            set_i(self.qVec, elem.nD_qpos, elem.nD_qneg, 
                  outV[len(elem.csOutPorts):])
        return self.qVec
            

    def get_source(self, ctime):
        """
        Get the source vector considering DC and TD source components

        ctime: current time.
        """
        # Erase vector first. 
        try:
            # If subcircuit add external currents
            self.sVec[:] = self.extSVec
        except AttributeError:
            # Not a subcircuit
            self.sVec.fill(0.)
        for elem in self.ckt.nD_sourceDCElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_DCsource()
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        for elem in self.ckt.nD_sourceTDElem:
            # first get the destination row/columns 
            outTerm = elem.nD_sourceOut
            current = elem.get_TDsource(ctime)
            # This may not need optimization because we usually do not
            # have too many independent sources
            # import pdb; pdb.set_trace()
            if outTerm[0] >= 0:
                self.sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                self.sVec[outTerm[1]] += current
        return self.sVec
        

    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G' xVec + i'(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents
        """
        # Linear contribution
        self.Gp.matvec(xVec, self.iVec)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            outV = elem.eval(xin)
            # Update iVec. outV may have extra charge elements but
            # they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_i(self.iVec, elem.nD_qpos, elem.nD_qneg, 
                  self.im.a0 * outV[len(elem.csOutPorts):])
        return self.iVec

    def get_i_Jac(self, xVec):
        """
        Calculate total current and Jacobian

        Returns (iVec, Jac)::

            iVec = G' xVec + i'(xVec)
            Jac = G' + (di'/dx)(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents

        Jac: system Jacobian
        """
        # Erase sparse matrix
        self.Jac.scale(0.)
        # Linear contribution
        self.Gp.matvec(xVec, self.iVec)
        self.Jac.shift(1., self.Gpll)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(len(elem.controlPorts))
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            (outV, outJac) = elem.eval_and_deriv(xin)
            # Update iVec and Jacobian now. outV may have extra charge
            # elements but they are not used in the following
            set_i(self.iVec, elem.nD_cpos, elem.nD_cneg, outV)
            set_i(self.iVec, elem.nD_qpos, elem.nD_qneg, 
                  self.im.a0 * outV[len(elem.csOutPorts):])
            set_Jac(self.Jac, elem.nD_cpos, elem.nD_cneg, 
                    elem.nD_vpos, elem.nD_vneg, outJac)
            qJac = self.im.a0 * outJac[len(elem.csOutPorts):,:]
            set_Jac(self.Jac, elem.nD_qpos, elem.nD_qneg, 
                    elem.nD_vpos, elem.nD_vneg, qJac)

        return (self.iVec, self.Jac)


