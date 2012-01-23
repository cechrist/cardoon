"""
:mod:`op` -- Operating Point Analysis
-------------------------------------

.. module:: op
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import numpy as np

from paramset import ParamSet
from analysis import AnalysisError, ipython_drop
import nodal as nd
from nodalAD import DCNodalAD
from fsolve import solve, NoConvergenceError

class Analysis(ParamSet):
    r"""
    DC Operating Point
    ------------------

    Calculates the DC operating point of a circuit using the nodal
    approach. After the analysis is complete, nodal voltages are saved
    in circuit and terminals with the ``nD_`` prefix.  After this the
    analysis drops to an interactive shell if the ``shell`` global
    variable is set to ``True``.

    By default the voltage at all external voltages is printed after
    the analysis is complete. Optionally the operating points of
    nonlinear elements can be printed. 

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``.

    Example::

        .analysis op intvars=1 shell=1

    OP/DC formulation
    +++++++++++++++++

    In the following discussion we assume that :math:`v` is the vector
    of nodal variables. There are 4 types of devices to consider for
    DC operating point:
    
      1. Linear VCCS: considered in the :math:`G` matrix.
      
      2. Nonlinear VCCS: considered in the :math:`i(v)` vector. This and
         its Jacobian (:math:`di/dv`) are returned by
         ``eval_and_deriv()``.
      
      3. Frequency-defined devices: their DC contribution is added to
         the :math:`G` matrix.
      
      4. Sources: contribute the source vector, :math:`s`. 
    
    The analysis solves the following nonlinear equation iteratively
    using Newton's method:
    
    .. math::
    
        G v + i(v) - s = 0
    
    The iteration is defined by linearizing :math:`i(v)` as follows:

    .. math::

        G v_{k+1} + i_k + \frac{di_k}{dv} \, (v_{k+1} - v_k) - s = 0 \; ,

    where the :math:`k` suffix indicates the iteration
    number. :math:`v_k` is assumed to be known and :math:`v_{k+1}` is
    the unknown to solve for. The initial guess (:math:`v_0`) is set
    to the values suggested by the nonlinear devices, if any, or
    otherwise to zero. The previous equation can be re-arranged as as
    the following system of linear equations:

    .. math::

         (G + \frac{di_k}{dv}) \, v_{k+1} = 
                s - i_k + \frac{di_k}{dv} \, v_k \; ,

    This equation can be seen as the nodal equation of a circuit
    obtained by substituting the nonlinear devices by current sources
    and transcunductances that are dependent of the current
    approximation for the nodal voltages (:math:`v_k`).

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "op"

    # Define parameters as follows
    paramDict = dict(
        intvars = ('Print internal element nodal variables', '', bool, False),
        elemop = ('Print element operating points', '', bool, False),
        fullAD = ('Use CPPAD for entire nonlinear part', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)


    def run(self, circuit):
        """
        Calculates the operating point by solving nodal equations

        The state of all devices is determined by the values of the
        voltages at the controlling ports.
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('             Operating point analysis')
        print('******************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # Create nodal object
        nd.make_nodal_circuit(circuit)
        if self.fullAD:
            dc = DCNodalAD(circuit)
        else:
            dc = nd.DCNodal(circuit)
        x0 = dc.get_guess()
        sV = dc.get_source()
        # solve equations
        try: 
            (x, res, iterations) = solve(x0, sV, dc)
        except NoConvergenceError as ce:
            print(ce)
            return
        dc.save_OP(x)

        print('Number of iterations = ', iterations)
        print('Residual = ', res)

        print('\n Node      |  Value               | Unit ')
        print('----------------------------------------')
        for key in sorted(circuit.termDict.iterkeys()):
            term = circuit.termDict[key]
            print('{0:10} | {1:20} | {2}'.format(key, term.nD_vOP, term.unit))

        if self.intvars or self.elemop:
            for elem in circuit.nD_nlinElem:
                print('\nElement: ', elem.nodeName)
                if self.intvars:
                    print('\n    Internal nodal variables:\n')
                    for term in elem.get_internal_terms():
                        print('    {0:10} : {1:20} {2}'.format(
                                term.nodeName, term.nD_vOP, term.unit))
                if self.elemop:
                    print('\n    Operating point info:\n')
                    for line in elem.format_OP().splitlines():
                        print('    ' + line.replace('|',':'))
        print('\n')

        def getvar(termname):
            return circuit.termDict[termname].nD_vOP

        def getterm(termname):
            return circuit.termDict[termname]

        def getdev(elemname):
            return circuit.elemDict[elemname]

        if self.shell:
            ipython_drop("""
Available commands:
    getvar(<terminal name>) returns variable at <terminal name>
    getterm(<terminal name>) returns terminal reference
    getdev(<device name>) returns device reference
""", globals(), locals())






