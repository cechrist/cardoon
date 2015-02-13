"""
:mod:`waveOP` -- Experimental Operating Point Analysis Based on Waves
-----------------------------------------------------

.. module:: waveOP
.. moduleauthor:: Tapan Savalia, Carlos Christoffersen

TODO: 

 * Make an analysis-independent block-creation and Newton solving
   module to be used for DC or transient analysis
 
 * Return results in a a proper format (easy to check results)

 * Implement parallel code

 * Add an automatic circuit partitioning algorithm

"""


from __future__ import print_function
import numpy as np
import scipy.linalg as linalg

from cardoon.paramset import ParamSet
from cardoon.globalVars import glVar
from analysis import AnalysisError, ipython_drop
import nodalSP
import nodal
from fsolve import solve, NoConvergenceError

class DCOP(ParamSet):            
    """
    Experimental DC Operating Point Based on Waves
    ----------------------------------------------

    This analysis is *experimental* and not yet ready for
    demonstration.  Circuit equations are formulated for parallel
    computation, but for now parallel processing not
    implemented. Assumptions:

      * Circuit must be given already divided in subcircuit blocks:
        can not be already flattened.  This analysis must be performed
        before any other analysis (since all other analyses flatten
        the circuit hierarchy).

      * Main circuit only contains subcircuit instances (no elements)

    Calculates the DC operating point of all subcircuits using a
    modification of nodal approach that uses wave variables at the
    subcircuit interfaces. After the analysis is complete, nodal
    voltages are saved in subcircuit block and terminals with the
    ``wND_`` prefix. 

    By default the voltage at all external voltages is printed after
    the analysis is complete. Optionally the operating points of
    nonlinear elements can be printed. 

    Example::

        .analysis waveop 

    """
    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "waveop"

    # Define parameters as follows
    paramDict = dict(
#        intvars = ('Print internal element nodal variables', '', bool, False),
#        elemop = ('Print element operating points', '', bool, False),
        maxiter = ('Maximum number of Newton iterations', '', int, 100),
        maxdelta = ('Maximum change in deltax for Newton iterations', '', 
                    float, 50.),
        z0 = ('Reference Impedance', 'Ohms', float, 50.),
#        rs = ('Series resistor between subcircuits', 'Ohms', float, 0.),
        nss = ('Number of source steps', '', int, 1)
        )

    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)
        

    def run(self, circuit):
        """
        Calculates the operating point by solving wave/nodal equations

        The state of the entire circuit is determined by the values of
        waves incident to subcircuit ports.
        """
        # for now just print some fixed stuff
        print('************************************************************')
        print('       Experimental Wave-Based Operating Point Analysis')
        print('************************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        # For now force dense matrices:
        if glVar.sparse:
            self.nd = nodalSP
        else:
            self.nd = nodal
            print('Using dense matrices\n')
        	
        # Initialize data structures for block calculations        
        self.init_blocks(circuit)
        # Set initial guess
        xnVec = np.zeros(self.systemSize, dtype=float)
        #xnVec[-self.circuit.nD_dimension:] = 0. <------------------------------
        rowAcc = 0
        for dcno in self.dcList:
            dim = dcno.ckt.nD_dimension
            xnVec[rowAcc:rowAcc + dim] = dcno.get_guess()
            rowAcc += dim

        # Source steppting loop goes from 1 to nss
        for i in range(1, self.nss+1):
            print('\nSource step: ', i)
            # Source stepping factor
            self.ssf = float(i) / self.nss

            success = False
            # Begin Newton iterations
            for nIter in xrange(self.maxiter):
                deltax = self.get_deltax(xnVec)
                #print(deltax)
                # Do not allow updates greater than maxdelta
                delta = max(abs(deltax))
                print('delta = ', delta)
                if delta > self.maxdelta:
                    limit = 1e5
                    if delta > limit:
                        mask = abs(deltax) > limit
                        deltax[mask] = limit
                        mask = abs(deltax) < -limit
                        deltax[mask] = -limit
                        delta = limit
                    deltax *= self.maxdelta/delta
                # Update x
                xnewVec = xnVec + deltax
    
                # Check if deltax is small enough
                n1 = np.all(abs(deltax) < (glVar.reltol *
                                           np.maximum(abs(xnVec), abs(xnewVec))
                                           + glVar.abstol))
                # Save new result for next iteration
                xnVec[:] = xnewVec

                res = max(abs(deltax))
                if n1:
                    success = True
                    break
            if not success:
                break
            
        # Print results
        if success:
            print('Success!')
        else:
            print('No convergence')
        print('iterations = ', nIter+1)

        # Must convert back to voltages to check results!
        # Currents are *not* converted
        aVec = xnVec[-self.nCurrents:] # incident waves
        rowAcc1 = 0
        for subckt, y0V, Nj in zip(self.subcktList,
                                   self.refCondList, self.NList):
            #import pdb; pdb.set_trace()
            dim = subckt.nD_dimension
            nports = len(y0V)
            xnVec[rowAcc1 : rowAcc1+nports] += np.dot(abs(Nj), aVec)[:nports]
            rowAcc1 += dim

        if xnVec.size < 20:
            print('xnewVec = ', xnVec)
        else:
            print('first part of xnewVec = ', xnVec[0:10])
              

    def init_blocks(self, circuit):
        """ 
        Initialize class attributes that are needed for subcircuit
        decomposition.

        List of attributes defined here:
        
        systemSize: full system dimension
        NList: list of Nj matrices
        dcList: list of dc objects for each subcircuit
        refCondList: list of reference conductance vectors (in S)
        nCurrents: number of extra interconnect current variables 
                   (from voltage sources)
        auxColList: list of columns to be extracted from Aj
        
        """
        # We want subcircuits arranged in alphabetical order
        nameList = [xsubckt.instanceName 
                    for xsubckt in circuit.subcktDict.itervalues()]
        nameList.sort()
        # To get subckt connections use circuit.subcktDict[name].connection
 
        # Before initializing the main circuit, flatten hierarchy of
        # subcircuits only. Also, create a list with all top-level
        # subcircuits (alhpabetical order).
        self.subcktList = []
        nCurrents = 0
        for j, name in enumerate(nameList):
	    # j: subcircuit number
            # get subcircuit definition
            subckt = circuit.cktDict[circuit.subcktDict[name].cktName] 
            subckt.flatten()
            # Append a *copy* of subckt definition with the same name as
            # the subcircuit instance
            defCopy = subckt.copy(name)
            defCopy.init()
            # print(defCopy.netlist_string())
            self.subcktList.append(defCopy)
            # Assign subcircuit owner to interconnect nodes (first
            # pass, Nj are created in second pass)
            for term in circuit.subcktDict[name].connection:
                try: 
                    if term.pop_owner != None:
                        # Terminal already in use by another
                        # subcircuit connection. Insert a 0V voltage
                        # source to separate terminals and create
                        # additional current variable.
                        term.pop_NjList.append((nCurrents, j))
                        nCurrents += 1
                except AttributeError:
                    # No one is using this terminal so far
                    term.pop_owner = j
                    term.pop_NjList = []
        # Initialize main circuit 
        circuit.init() 
        # systemSize: total dimension including subcircuit variables,
        # initially set to the number of interconnect currents in
        # voltage sources
        self.systemSize = nCurrents
        self.nCurrents = nCurrents
        print('ncurrents = ', nCurrents)
        # Create list of nodal object for subcircuits
        self.dcList = []
        # shape interconnect blocks Nj
        self.NList=[]
        # List of Aj columns to be extracted and copied in Nj (all)
        self.auxColList=[]
        # Reference conductance vectors list
        self.refCondList = []
        # Total number of interconnects in circuit
        self.nInterconnect = 0
        self.y0 = 1. / self.z0 # Conductance in S
        # Second pass: create DC objects and Nj
        for j, subckt in enumerate(self.subcktList):
	    # j: subcircuit number
            self.nd.make_nodal_circuit(subckt, subckt.get_connections())
            self.dcList.append(self.nd.DCNodal(subckt))
            y0V = self.y0 * np.ones_like(subckt.nD_extRClist)
            self.refCondList.append(y0V)
            # Create Nj
            Nj = np.zeros((subckt.nD_dimension, nCurrents), dtype=float)
            auxCol1 = []
            self.auxColList.append(auxCol1)
            for row, term in enumerate(
                circuit.subcktDict[subckt.name].connection):
                if term.pop_owner == j:
                    # Our terminal: assign ones
                    for col, subcktNum in term.pop_NjList:
                        Nj[row, col] = 1.
                        auxCol1.append(col)
                else:
                    # Not ours, assign -1 only in one element of Nj
                    for col, subcktNum in term.pop_NjList:
                        if subcktNum == j:
                            Nj[row, col] = -1.
                            auxCol1.append(col)
                    y0V[row] *= -1.
	    self.NList.append(Nj)
            self.systemSize += subckt.nD_dimension


    def get_deltax(self, xnVec):
        """
        Return deltax from:

        A deltax = sv - ivec

        System to be solved:

                 |A1          N1|   | delta_x1  |   |  b1  |
                 |    A2      N2|   | delta_x2  |   |  b2  |
                 |       .    : |   | :         | = |   :  |   (11)
                 |         .  : |   | delta_xj  |   |  bj  |
                 |M1  M2  ..  0 |   | delta_iic |   |  bC  |

 	Aj: jth subcircuit MNAM
	Nj: incidence matrix for voltage sources connecting subcircuit
	Mj = Nj^T
	delta_xj: update for nodal variables
	delta_iic: update for interconnect currents
                   (through voltage sources connecting subcircuits)
	bj = sj - iVecj (Newton's Method rhs vector)
	bC = -sum Mj xj

        """
	# j: subcircuit index
	# n: Newton iteration index        

        # iicVec -> delta_iic
        iicVec = xnVec[-self.nCurrents:]
	bC = np.zeros(self.nCurrents, dtype=float)

        # Node numbering for circuit blocks: External nodes first, then internal

        # Get matrix blocks and source vectors from DCNodal objects

        # Add voltage sources
        AccC = np.zeros((self.nCurrents, self.nCurrents), dtype=float)
        # Add resistors between subcircuits (in Ohms): doesn't work?
        # AccC = np.diag(self.rs * np.ones(self.nCurrents, dtype=float))
        
        AList = []  # AList = [A1, A2,.....Aj] : Subcircuit List
        AnList=[]
        bList = []  # BList = [b1, b2......bj] : Source vector of
        rowAcc = 0 
        NewNList = []
	for dcno, Nj, y0V, auxCol1 in zip(self.dcList, 
                                          self.NList, 
                                          self.refCondList, 
                                          self.auxColList):
            #import pdb; pdb.set_trace()
            # Equation__(14)
            dim = dcno.ckt.nD_dimension
            # xjVec : subcircuit source vectors as shown in Eq.(11) : x1 to xj
            xjVec = xnVec[rowAcc : rowAcc + dim] 
            (iVec, Aj) = dcno.get_i_Jac(xjVec)
            NewNj = np.zeros_like(Nj)
            AjBlock = self.nd.get_submatrix(Aj, len(y0V))
	    NewNj[:, auxCol1] = AjBlock
            NewNj -= self.y0 * Nj
            NewNList.append(NewNj)
            Aj = self.nd.add_to_diagonal(Aj, y0V)
            AList.append(Aj)
            bj = self.ssf * dcno.get_source() - iVec - np.dot(NewNj, iicVec)
            bList.append(bj)             
            # Build AccC and bC without the inverting Aj
            deltaXjStar = dcno.factor_and_solve(bj, Aj)
            # Added 1e-3 factor: currents in mA
            AccC -= np.dot(Nj.T, dcno.solve_linear_system(NewNj))
            bC -= np.dot(Nj.T, xjVec + deltaXjStar) 
            rowAcc += dim    

        deltaX = np.empty(self.systemSize, dtype=float)
        # Solve Eq. (14) to get nodal voltages at interconnect and
        # subcircuit external currents
        try:
            deltaIicVec = linalg.solve(AccC, bC)
        except linalg.LinAlgError:
            # Try solving least squares problem
            deltaIicVec = linalg.lstsq(AccC, bC)[0]
        except ValueError:
            print(AccC)
            print(bC)
            exit()
        # Write currents to output vector!
        deltaX[-self.nCurrents:] = deltaIicVec
        rowAcc = 0
        # Fill contributions from subcircuits
        for dcno, bj, NewNj in zip(self.dcList, bList, NewNList):
            bj = bj - np.dot(NewNj, deltaIicVec)   
            deltaXjVec = dcno.solve_linear_system(bj)
            dim = bj.size 				#|  Equation_____(12)
            deltaX[rowAcc:rowAcc + dim] = deltaXjVec
            rowAcc += dim
        return deltaX                                
        
aClass = DCOP

