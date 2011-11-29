"""
:mod:`op` -- Operating Point Analysis
-------------------------------------

.. module:: op
.. moduleauthor:: Carlos Christoffersen and others

"""

from __future__ import print_function
import numpy as np

from paramset import ParamSet
from analysis import AnalysisError
import nodal as nd
from fsolve import solve

class Analysis(ParamSet):
    """
    DC Operating Point Calculation

    Calculates the DC operating point of a circuit. This is also
    useful for many other analyses such as DC, AC, Transient, HB,
    etc. Charge sources are ignored in this calculation.
    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "op"

    # Define parameters as follows
    paramDict = dict(
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
        # Create nodal object
        nd.make_nodal_circuit(circuit)
        dc = nd.DCNodal(circuit)
        x0 = dc.get_guess()
        sV = dc.get_source()
        # solve equations
        (x, res, iter) = solve(x0, sV, dc)
        dc.save_OP(x)

        # for now just print some fixed stuff
        print('******************************************************')
        print('             Operating point analysis')
        print('******************************************************')
        print('Number of iterations = ', iter)
        print('Residual = ', res)

        print('\n Node      |  Value               | Unit ')
        print('----------------------------------------')
        for key in sorted(circuit.termDict.iterkeys()):
            term = circuit.termDict[key]
            print('{0:10} | {1:20} | {2}'.format(key, term.nD_v, term.unit))

        for elem in circuit.nD_nlinElements:
            print('\nElement: ', elem.nodeName)
            print(' Variable  |  Value ')
            print('-------------------------')
            print(elem.format_OP())

        if self.shell: 
            from IPython.Shell import IPShellEmbed
            args = ['-pi1','In <\\#>: ','-pi2','   .\\D.: ',
                    '-po','Out<\\#>: ','-nosep']
            ipshell = IPShellEmbed(args, 
                                   banner = 'Dropping into IPython, type CTR-D to exit',
                                   exit_msg = 'Leaving Interpreter, back to program.')
            ipshell()






