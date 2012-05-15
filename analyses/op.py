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
from nodalAD import DCNodalAD
from globalVars import glVar
from fsolve import solve, NoConvergenceError
import nodalSP
import nodal

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
    the global variables in ``.options``. The type of matrix used in
    this analysis is controlled by the ``sparse`` option. Global
    options are documented in :doc:`global_vars`. 

    OP analysis formulation is documented in :doc:`analysis`, and
    internal classes and functions used in this analysis are
    documented in :doc:`analyses_classes`.

    Example::

        .analysis op intvars=1 shell=1

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

        if glVar.sparse:
            nd = nodalSP
        else:
            nd = nodal
            print('Using dense matrices\n')

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
            (x, res, iterations) = solve(x0, sV, dc.convergence_helpers)
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






