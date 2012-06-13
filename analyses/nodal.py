"""
:mod:`nodal` -- Nodal Analysis
------------------------------

.. module:: nodal
.. moduleauthor:: Carlos Christoffersen

This module contains basic classes/functions for nodal analysis. These
are part of the standard netlist analyses but they can also be used
independently.

"""

from __future__ import print_function
import numpy as np
from fsolve import fsolve_Newton, NoConvergenceError
from integration import BEuler

# ****************** Stand-alone functions to be optimized ****************

def set_quad(G, row1, col1, row2, col2, g):
    """
    Set transconductance/transcapacitance quad

    G: target matrix
    """
    # good cython candidate
    #import pdb; pdb.set_trace()
    if col1 >= 0:
        if row1 >= 0:
            G[row1, col1] += g
        if row2 >= 0:
            G[row2, col1] -= g
    if col2 >= 0:
        if row1 >= 0:
            G[row1, col2] -= g
        if row2 >= 0:
            G[row2, col2] += g


def set_xin(xin, posCols, negCols, xVec):
    """
    Calculate input port voltage xin
    """
    # good candidate for cython
    for i,j in posCols:
        xin[i] += xVec[j]
    for i,j in negCols:
        xin[i] -= xVec[j]


def set_i(iVec, posRows, negRows, current):
    """
    Set current contributions of current into iVec
    """
    # good candidate for cython
    for i,j in posRows:
        iVec[j] += current[i]
    for i,j in negRows:
        iVec[j] -= current[i]


def set_Jac(M, posRows, negRows, posCols, negCols, Jac):
    """
    Set current contributions of Jac into M
    """
    # good candidate for cython
    for i1, i in posRows:
        for j1, j in posCols:
            M[i,j] += Jac[i1,j1]
        for j1, j in negCols:
            M[i,j] -= Jac[i1,j1]
    for i1, i in negRows:
        for j1, j in negCols:
            M[i,j] += Jac[i1,j1]
        for j1, j in posCols:
            M[i,j] -= Jac[i1,j1]

def delay_interp(td, vp, h, tsV, vpV):
    """
    Use linear interpolation to find v(t-td)

    td: time delay
    vp: current port voltage
    h: current time step
    tsV: vector with past time steps
    vpV: vector with past voltage samples

    returns (v(t-td), deriv_coeff)
    """
    #import pdb; pdb.set_trace()
    if td <= h:
        deriv_coeff = (h - td) / h
        vtd = vpV[-1] + (vp - vpV[-1]) * deriv_coeff
    else:
        deriv_coeff = 0.
        timeacc = h
        for i in xrange(-1, -len(tsV), -1):
            timeacc += tsV[i]
            if td < timeacc:
                # Found position: do interpolation
                vtd = vpV[i-1] + (vpV[i] - vpV[i-1]) \
                    * (timeacc - td) / tsV[i]
                break
    return (vtd, deriv_coeff) 


# ********************** Regular functions *****************************

def make_nodal_circuit(ckt, reference='gnd'):
    """
    Add attributes to Circuit/Elements/Terminals for nodal analysis

    This functionality should be useful for any kind of nodal-based
    analysis (DC, AC, TRAN, HB, etc.)

    Takes a Circuit instance (ckt) as an argument. If the circuit
    contains the 'gnd' node, it is used as the reference. Otherwise a
    reference node must be indicated.

    New attributes are added in Circuit/Element/Terminal
    instances. All new attributes start with ``nD_``

    Works with subcircuits too (see ``nD_namRClist`` attribute)
    """
    # get ground node
    ckt.nD_ref = ckt.get_term(reference)

    # make a list of all non-reference terminals in circuit 
    ckt.nD_termList = ckt.termDict.values() + ckt.get_internal_terms()
    # remove ground node from terminal list
    ckt.nD_termList.remove(ckt.nD_ref)
    # Assign a number (0-inf) to all nodes. For reference nodes
    # assign -1 
    ckt.nD_ref.nD_namRC = -1
    # Make a list of all elements
    ckt.nD_elemList = ckt.elemDict.values()
    # Set RC number of reference terminals to -1
    for elem in ckt.nD_elemList:
        if elem.localReference:
            elem.neighbour[elem.localReference].nD_namRC = -1
    # For the future: use graph techniques to find the optimum
    # terminal order
    for i, term in enumerate(ckt.nD_termList):
        term.nD_namRC = i

    # Store internal RC numbers for later use
    for elem in ckt.nD_elemList:
        elem.nD_intRC = [term.nD_namRC for term in 
                         elem.neighbour[elem.numTerms:]]

    # Dimension is the number of unknowns to solve for
    ckt.nD_dimension = len(ckt.nD_termList)
    # Number of external terminals excluding reference
    ckt.nD_nterms = len(ckt.termDict.values()) - 1
    # Maximum delay in circuit
    ckt.nD_maxDelay = 0.

    # Create specialized element lists
    ckt.nD_nlinElem = filter(lambda x: x.isNonlinear, ckt.nD_elemList)
    ckt.nD_delayElem = filter(lambda x: x.nDelays, ckt.nD_nlinElem)
    ckt.nD_freqDefinedElem = filter(lambda x: x.isFreqDefined, ckt.nD_elemList)
    ckt.nD_sourceDCElem = filter(lambda x: x.isDCSource, ckt.nD_elemList)
    ckt.nD_sourceTDElem = filter(lambda x: x.isTDSource, ckt.nD_elemList)
    ckt.nD_sourceFDElem = filter(lambda x: x.isFDSource, ckt.nD_elemList)

    # Map row/column numbers directly into VC*S descriptions
    for elem in ckt.nD_elemList:
        process_nodal_element(elem, ckt)

    # Subcircuit-connection processing
    try:
        connectTerms = ckt.get_connections()
        # List of RC numbers of external connections
        ckt.nD_namRClist = [term.nD_namRC for term in connectTerms]
    except AttributeError:
        # Not a subcircuit
        pass

def restore_RCnumbers(elem):
    """
    Restore RC numbers in internal terminals

    Assumption is number of internal terminals is the same
    """
    for term, i in zip(elem.neighbour[elem.numTerms:], elem.nD_intRC):
        term.nD_namRC = i

def process_nodal_element(elem, ckt):
    """
    Process element for nodal analysis
    """
    # Create list with RC numbers (choose one)
    # rcList = map(lambda x: x.nD_namRC, elem.neighbour)
    rcList = [x.nD_namRC for x in elem.neighbour]

    # Translate linear VCCS/VCQS format
    def convert_vcs(x):
        """
        Converts format of VC*S 

        input: [(contn1, contn2), (outn1, outn2), g]

        output: [row1, col1, row2, col2, g]
        """
        col1 = rcList[x[0][0]]
        col2 = rcList[x[0][1]]
        row1 = rcList[x[1][0]]
        row2 = rcList[x[1][1]]
        return [row1, col1, row2, col2, x[2]]
    elem.nD_linVCCS = map(convert_vcs, elem.linearVCCS)
    elem.nD_linVCQS = map(convert_vcs, elem.linearVCQS)

    # Translate positive and negative terminal numbers
    def create_list(portlist):
        """
        Converts an internal port list into 2 lists with (+-) nodes

        The format of each list is::

            (internal term number, namRC number)
        """
        tmp0 = [rcList[x1[0]] for x1 in portlist]
        tmp1 = [rcList[x1[1]] for x1 in portlist]
        return ([(i, j) for i,j in enumerate(tmp0) if j > -1],
                [(i, j) for i,j in enumerate(tmp1) if j > -1])

    # Convert nonlinear device descriptions to a format more
    # readily usable for the NA approach
    if elem.isNonlinear:
        if elem.nDelays:
            # Delayed control voltages (not needed in this form)
            # (elem.nD_dpos, elem.nD_dneg) = create_list(elem.delayedContPorts)
            # Store delays
            elem.nD_delay = np.array([x1[2] for x1 in elem.delayedContPorts])
            # Keep track of longest delay
            maxD = np.amax(elem.nD_delay)
            if maxD > ckt.nD_maxDelay:
                ckt.nD_maxDelay = maxD
            # Control voltages plus delayed ports
            allPorts = elem.controlPorts + elem.delayedContPorts
            (elem.nD_vpos, elem.nD_vneg) = create_list(allPorts)
            elem.nD_nxin = len(allPorts)
        else:
            # Control voltages
            (elem.nD_vpos, elem.nD_vneg) = create_list(elem.controlPorts)
            elem.nD_nxin = len(elem.controlPorts)
        # Current source terminals
        (elem.nD_cpos, elem.nD_cneg) = create_list(elem.csOutPorts)
        # Charge source terminals
        (elem.nD_qpos, elem.nD_qneg) = create_list(elem.qsOutPorts)
        # Create list with external ports for gmin calculation
        nports = elem.numTerms - 1
        elem.nD_extPorts = [(i, nports) for i in range(nports)]
        (elem.nD_epos, elem.nD_eneg) = create_list(elem.nD_extPorts)

    # Convert frequency-defined elements
    if elem.isFreqDefined:
        (elem.nD_fpos, elem.nD_fneg) =  create_list(elem.fPortsDefinition)

    # Translate source output terms
    if elem.isDCSource or elem.isTDSource or elem.isFDSource:
        # first get the destination row/columns 
        n1 = rcList[elem.sourceOutput[0]]
        n2 = rcList[elem.sourceOutput[1]]
        elem.nD_sourceOut = (n1, n2)

#----------------------------------------------------------------------    

def run_AC(ckt, fvec):
    """
    Set up and solve AC equations

    ckt: nodal-ready Circuit instance (ckt) instance with OP already set

    fvec: frequency vector. All frequencies must be different from zero

    DC solution not calculated here as it needs special treatment and
    that solution is already obtained by the OP analysis.

    Set results in terminals in vectors named ``aC_V``

    Returns a matrix with results. Dimension: (nvars x nfreq)

    """
    # Number of frequencies
    nfreq = np.shape(fvec)[0]
    # Allocate matrices/vectors. 
    # G is used for G + dI/dv
    G = np.zeros((ckt.nD_dimension, ckt.nD_dimension))
    # C is used for C + dQ/dv
    C = np.zeros((ckt.nD_dimension, ckt.nD_dimension))

    # Linear contributions G and C in AC formulation
    for elem in ckt.nD_elemList:
        # All elements have nD_linVCCS (perhaps empty)
        for vccs in elem.nD_linVCCS:
            set_quad(G, *vccs)
        # All elements have nD_linVCQS (perhaps empty)
        for vcqs in elem.nD_linVCQS:
            set_quad(C, *vcqs)
    # Calculate dI/dx and dQ/dx and add to G and C
    for elem in ckt.nD_nlinElem:
        (outV, outJac) = elem.eval_and_deriv(elem.nD_xOP)
        set_Jac(G, elem.nD_cpos, elem.nD_cneg, 
                elem.nD_vpos, elem.nD_vneg, outJac)
        if len(elem.qsOutPorts):
            # import pdb; pdb.set_trace()
            # Get charge derivatives from outJac
            qJac = outJac[len(elem.csOutPorts):,:]
            set_Jac(C, elem.nD_qpos, elem.nD_qneg, 
                    elem.nD_vpos, elem.nD_vneg, qJac)

    # Frequency-dependent matrices: a matrix of complex vectors, each
    # vector contains the values for all frequencies. This format is
    # compatible with the format returned by elem.get_Y_matrix(fvec)
    Y = np.zeros((ckt.nD_dimension, ckt.nD_dimension, nfreq), dtype=complex)
    # Frequency-defined elements
    for elem in ckt.nD_freqDefinedElem:
        # get first Y matrix for all frequencies. 
        set_Jac(Y, elem.nD_fpos, elem.nD_fneg, 
                elem.nD_fpos, elem.nD_fneg, elem.get_Y_matrix(fvec))

    # Sources
    sVec = np.zeros(ckt.nD_dimension, dtype = complex)
    for source in ckt.nD_sourceFDElem:
        try:
            current = source.get_AC()
            outTerm = source.nD_sourceOut
            # This may not need optimization because we usually do not
            # have too many independent sources
            if outTerm[0] >= 0:
                sVec[outTerm[0]] -= current
            if outTerm[1] >= 0:
                sVec[outTerm[1]] += current
        except AttributeError:
            # AC routine not implemented, ignore
            pass

#    if hasattr(ckt, 'nD_namRClist'):
#        # Allocate external currents vector
#        extSVec = np.zeros(ckt.nD_dimension)

    # Loop for each frequency: create and solve linear system
    j = complex(0., 1.)
    omegaVec = 2. * np.pi * fvec
    xVec = np.zeros((ckt.nD_dimension, nfreq), dtype=complex)
    for k, omega in enumerate(omegaVec):
        N = G + j * omega * C + Y[:,:,k] 
        xVec[:,k] = np.linalg.solve(N, sVec)
        
    # Save results in circuit
    # Set nodal voltage of reference to zero
    ckt.nD_ref.aC_V = 0.
    for k,term in enumerate(ckt.nD_termList):
        term.aC_V = xVec[k, :]

    return xVec


#-------------------------------------------------------------------------
# ****************************** Classes *********************************
#-------------------------------------------------------------------------

class _NLFunction(object):
    """
    Nonlinear function interface class

    This is an abstract class to be use as a base class by other
    nodal-based classes such as DCNodal
    """

    def __init__(self):
        # List here the functions that can be used to solve equations
        self.convergence_helpers = [self.solve_simple, 
                                    self.solve_homotopy_gmin2, 
                                    self.solve_homotopy_source, 
                                    self.solve_homotopy_gmin, 
                                    None]

    def _get_deltax(self, errFunc, Jac):
        """
        Solves linear system: Jac deltax = errFunc
        """
        try:
            deltax = np.linalg.solve(Jac, errFunc)
        except:
            print('Singular Jacobian')
            # Use pseudo-inverse
            deltax = np.dot(np.linalg.pinv(Jac), errFunc)
        return deltax

    def _set_gmin(self, _lambda):
        """
        Sets gmin value given lambda (used for homotopy)

        Range of lambda: [1e-4, 1]
        Range of gmin: [10, 0]
        """
        gbase = 1e-3
        self.gmin = gbase / _lambda - gbase

    def _homotopy(self, _lambda, f, x0, get_deltax, f_eval):
        """
        Controls _lambda and the homotopy flow
        
        The lambda parameter is varied from 0 to 1. 0 corresponds to a
        problem easy to solve and 1 correstponds to the original problem.
        Uses bisection to find lambda step step size to reach 1.

        Inputs:

            _lambda: Initial value for lambda
            f: f(lambda) (for example gmin(lambda))        
            x0: initial guess (modified on output to best approximation)
            get_deltax and f_eval: functions passed to Newton's solver

        Output:

          (x, res, iterations, success)

        """
        stack = [0., 1.]
        step = _lambda
        small = 1e-4
        x = np.copy(x0)
        success = True
        totIter = 0
        sepline = '===================================================='
        print('    lambda      |   Iterations    |   Residual')
        print(sepline)
        while stack:
            # Update parameter in nonlinear functions
            f(_lambda)
            (x, res, iterations, success1) = \
                fsolve_Newton(x, get_deltax, f_eval)
            print('{0:15} | {1:15} | {2:15}'.format(
                    _lambda, iterations, res), end='')
            totIter += iterations
            if success1:
                print('')
                # Save result
                x0[:] = x
                # Recover value of lambda_ from stack
                step = stack[-1] - _lambda
                _lambda = stack.pop()
            else:
                print('  <--- Backtracking')
                # Check if residual was big
                #if res > 1e-3: (not reliable)
                # Restore previous better guess
                x[:] = x0
                # push _lambda into stack
                stack.append(_lambda)
                step *= .5
                _lambda -= step
                if (_lambda < small) or (step < small):
                    success = False
                    break
        print(sepline)
        print('Total iterations: ', totIter)
        return (x, res, totIter, success)

    # The following functions used to solve equations, originally from
    # pycircuit but since they have evolved quite a bit
    def solve_simple(self, x0, sV):
        #"""Simple Newton's method"""
        # Docstring removed to avoid printing this all the time
        def get_deltax(x):
            (iVec, Jac) = self.get_i_Jac(x) 
            return self._get_deltax(iVec - sV, Jac)
        def f_eval(x):
            iVec = self.get_i(x) 
            return iVec - sV
        (x, res, iterations, success) = \
            fsolve_Newton(x0, get_deltax, f_eval)
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError(
                'No convergence. iter = {0} res = {1}'.format(iterations, res))
        
    def solve_homotopy_gmin(self, x0, sV):
        """Newton's method with gmin stepping"""
        idx = np.arange(self.ckt.nD_nterms)
        # Add gmin from ground to all external nodes. Assumes all
        # external nodes are sorted first in the vector. This will
        # not work if the terminal order is changed.
        def get_deltax(xVec):
            (iVec, Jac) = self.get_i_Jac(xVec) 
            iVec[idx] += self.gmin * xVec[idx]
            Jac[idx, idx] += self.gmin
            return self._get_deltax(iVec - sV, Jac)
        def f_eval(xVec):
            iVec = self.get_i(xVec)
            iVec[idx] += self.gmin * xVec[idx]
            return iVec - sV
        (x, res, iterations, success) = \
            self._homotopy(0.5, self._set_gmin, x0, get_deltax, f_eval)
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError('gmin stepping did not converge')

    def solve_homotopy_gmin2(self, x0, sV):
        """Newton's method with gmin stepping (nonlinear ports)"""
        # Adds gmin in parallel with nonlinear element (external) ports
        # Create Gmin matrix (structure is fixed)
        Gones = np.zeros((self.ckt.nD_dimension, self.ckt.nD_dimension))
        for elem in self.ckt.nD_nlinElem:
            Gadd = np.eye(len(elem.nD_extPorts))
            # import pdb; pdb.set_trace()
            set_Jac(Gones, elem.nD_epos, elem.nD_eneg, 
                    elem.nD_epos, elem.nD_eneg, Gadd)
        def get_deltax(xVec):
            (iVec, Jac) = self.get_i_Jac(xVec) 
            iVec += self.gmin * np.dot(Gones, xVec)
            Jac += self.gmin * Gones
            return self._get_deltax(iVec - sV, Jac)
        def f_eval(xVec):
            iVec = self.get_i(xVec)
            iVec += self.gmin * np.dot(Gones, xVec)
            return iVec - sV
        (x, res, iterations, success) = \
            self._homotopy(0.5, self._set_gmin, x0, get_deltax, f_eval)
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError('gmin stepping did not converge')

    def solve_homotopy_source(self, x0, sV):
        """Newton's method with source stepping"""
        def f(_lambda):
            self._lambda = _lambda
        def get_deltax(x):
            (iVec, Jac) = self.get_i_Jac(x) 
            return self._get_deltax(iVec - self._lambda * sV, Jac)
        def f_eval(x):
            return self.get_i(x) - self._lambda * sV
        (x, res, iterations, success) = \
            self._homotopy(0.5, f, x0, get_deltax, f_eval)
        if success:
            return (x, res, iterations)
        else:
            raise NoConvergenceError('Source stepping did not converge')

#---------------------------------------------------------------------------

class DCNodal(_NLFunction):
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
        # G here is G1 = G + G0 in documentation
        self.G = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        # Jac is (G1 + di/dv) in doc
        self.Jac = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        self.sVec = np.empty(self.ckt.nD_dimension)
        self.iVec = np.empty(self.ckt.nD_dimension)
        if hasattr(self.ckt, 'nD_namRClist'):
            # Allocate external currents vector
            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parameter sweeps
        """
        self.G.fill(0.)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.G, *vccs)
        # Frequency-defined elements
        for elem in self.ckt.nD_freqDefinedElem:
            set_Jac(self.G, elem.nD_fpos, elem.nD_fneg, 
                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())
            
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
                # Only use positive side. This is not the only way and
                # may not work well in some cases but this is a guess
                # anyway
                for i,j in elem.nD_vpos:
                    x0[j] = elem.vPortGuess[i]
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
        # Erase arrays
        self.iVec.fill(0.)
        # Linear contribution
        self.iVec += np.dot(self.G, xVec)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
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
        # Erase arrays
        self.iVec.fill(0.)
        self.Jac.fill(0.)
        # Linear contribution
        self.iVec += np.dot(self.G, xVec)
        self.Jac += self.G
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
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
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            # Set OP in element (discard return value)
            elem.nD_xOP = xin
            elem.get_OP(xin)

    
#----------------------------------------------------------------------

class TransientNodal(_NLFunction):
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

        # Make sure circuit is ready (analysis should take care)
        assert ckt.nD_ref
        # Save circuit and integration method instance
        self.ckt = ckt
        self.im = im

        # Allocate matrices/vectors
        # G, C and G'
        self.G = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        self.C = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        self.Gp = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        # iVec = G' x + i'(x)   total current
        self.iVec = np.empty(self.ckt.nD_dimension)
        # System Jacobian: G' + di'/dx
        self.Jac = np.empty((self.ckt.nD_dimension, self.ckt.nD_dimension))
        # Total charge: C x + q(x)
        self.qVec = np.empty(self.ckt.nD_dimension)
        # Source vector at current time s(t) 
        self.sVec = np.empty(self.ckt.nD_dimension)

#        if hasattr(self.ckt, 'nD_namRClist'):
#            # Allocate external currents vector
#            self.extSVec = np.empty(self.ckt.nD_dimension)
        self.refresh()

    def refresh(self):
        """
        Re-generate linear matrices

        Used for parametric sweeps
        """
        self.G.fill(0.)
        self.C.fill(0.)
        # Generate G matrix (never changes)
        for elem in self.ckt.nD_elemList:
            # All elements have nD_linVCCS (perhaps empty)
            for vccs in elem.nD_linVCCS:
                set_quad(self.G, *vccs)
            for vccs in elem.nD_linVCQS:
                set_quad(self.C, *vccs)
#        # Frequency-defined elements not included for now
#        for elem in self.ckt.nD_freqDefinedElem:
#            set_Jac(self.G, elem.nD_fpos, elem.nD_fneg, 
#                    elem.nD_fpos, elem.nD_fneg, elem.get_G_matrix())


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

        # Initialize time delay structures
        nsteps = np.ceil(self.ckt.nD_maxDelay / h)
        self._dpsteps = nsteps * [h]
        for elem in self.ckt.nD_nlinElem:
            if elem.nDelays:
                elem.tran_dpvolt = []
                # first have to retrieve port voltages from xVec
                xin = np.zeros(elem.nD_nxin)
                set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec) 
                # Extract delayed ports and append to element list
                for i in xrange(-elem.nDelays, 0):
                    elem.tran_dpvolt.append(nsteps * [xin[i]])


    def update_Gp(self):
        """
        Recalculate Gp from im information
        """
        self.Gp[:,:] = self.G + self.im.a0 * self.C
        return self.Gp


    def update_q(self, xVec):
        """ 
        Recalculate qVec for a given value of xVec
        """
        # Add linear charges
        self.qVec[:] = np.dot(self.C, xVec)
        # Nonlinear charges
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)  
            outV = elem.eval(xin)
            set_i(self.qVec, elem.nD_qpos, elem.nD_qneg, 
                  outV[len(elem.csOutPorts):])
        return self.qVec
            

    def accept(self, xVec):
        """
        Accept xVec for current time step and store state
        """
        self.im.accept(self.update_q(xVec))
        # Append current time step to list
        self._dpsteps.append(self.im.h)
        # Store here xVec components. Optimally this would be
        # performed in the same loop as update_q, but kept here for
        # clarity.
        for elem in self.ckt.nD_delayElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec) 
            # Extract delayed ports and append to element lists
            for i in xrange(-elem.nDelays, 0):
                elem.tran_dpvolt[i].append(xin[i])


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

    def get_rhs(self, ctime):
        """
        Returns system rhs vector: s' 

        s' includes source vector, charge history and convolution
        contributions

        ctime: current time
        """
        # TODO: add convolution
        return self.im.f_n1() + self.get_source(ctime)

    def get_i(self, xVec):
        """
        Calculate total current

        returns iVec = G' xVec + i'(xVec)

        xVec: input vector of nodal voltages. 

        iVec: output vector of currents
        """
        # Linear contribution
        self.iVec[:] = np.dot(self.Gp, xVec)
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            if elem.nDelays:
                # Apply delay to port voltages
                for i in xrange(-elem.nDelays, 0):
                    (xin[i], dc) = delay_interp(elem.nD_delay[i], 
                                                xin[i], self.im.h, 
                                                self._dpsteps, 
                                                elem.tran_dpvolt[i])
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
        # Linear contribution
        self.iVec[:] = np.dot(self.Gp, xVec)
        self.Jac[:,:] = self.Gp
        # Nonlinear contribution
        for elem in self.ckt.nD_nlinElem:
            # first have to retrieve port voltages from xVec
            xin = np.zeros(elem.nD_nxin)
            set_xin(xin, elem.nD_vpos, elem.nD_vneg, xVec)
            if elem.nDelays:
                # Derivative coefficients
                dc = np.empty_like(elem.nD_delay)
                # Apply delay to port voltages
                for i in xrange(-elem.nDelays, 0):
                    (xin[i], dc[i]) = delay_interp(elem.nD_delay[i], 
                                                   xin[i], self.im.h, 
                                                   self._dpsteps, 
                                                   elem.tran_dpvolt[i])
                (outV, outJac) = elem.eval_and_deriv(xin)
                # Multiply Jacobian columns by derivative factors
                for i in xrange(-elem.nDelays, 0):
                    outJac[:,i] *= dc[i]
            else:
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


