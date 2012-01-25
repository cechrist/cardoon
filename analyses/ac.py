"""
:mod:`ac` -- AC sweep Analysis
------------------------------

.. module:: ac
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import numpy as np

from paramset import ParamSet
from analysis import AnalysisError, ipython_drop
import nodal as nd
from fsolve import solve, NoConvergenceError
import matplotlib.pyplot as plt


class Analysis(ParamSet):
    r"""
    AC Sweep
    --------

    Calculates a AC sweep of a circuit using the nodal approach. After
    the analysis is complete, nodal voltages are saved in circuit and
    terminals with the ``aC_`` prefix.  After this the analysis drops
    to an interactive shell if the ``shell`` global variable is set to
    ``True``.

    An OP analysis is performed first to obtain the Jacobian from
    nonlinear devices. Convergence parameters for the Newton method
    are controlled using the global variables in ``.options``.

    One plot window is generated for each ``.plot`` statement. Use the
    following request types for this analysis: ``ac_mag``,
    ``ac_phase`` or ``ac_dB``.

    AC formulation documented in :doc:`analysis`

    Example::

        .analysis ac start=100. stop=1MEG num=100 log=True
        .plot ac_dB 153 151 23
        .plot ac_phase 23

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "ac"

    # Define parameters as follows
    paramDict = dict(
        start = ('Frequency sweep start value', 'Hz', float, 1.),
        stop = ('Frequency sweep stop value', 'Hz', float, 10.),
        log = ('Use logarithmic scale', '', bool, False),
        num = ('Number of points in sweep', '', int, 50),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)


    def run(self, circuit):
        """
        Calculates a AC sweep by solving nodal equations around

        The parameter to be swept is specified in the analysis options
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('                 AC sweep analysis')
        print('******************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        nd.make_nodal_circuit(circuit)
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

        # Create frequency vector
        if self.log:
            fvec = np.logspace(start = np.log10(self.start), 
                               stop = np.log10(self.stop), 
                               num = self.num)
            pltfunc = plt.semilogx
        else:
            fvec = np.linspace(start = self.start, 
                               stop = self.stop, 
                               num = self.num)
            pltfunc = plt.plot

        # Perform analysis
        nd.run_AC(circuit, fvec)

        # Process output requests.  In the future this should be moved
        # to a common module that processes output requests such as
        # plot, print, save, etc.
        flag = False
        for outreq in circuit.outReqList:
            if outreq.type[:2] == 'ac':
                flag = True
                plt.figure()
                plt.grid(True)
                plt.xlabel('Frequency [Hz]')

                # First item may be output format
                pkey = outreq.type[3:]
                if pkey == 'mag':
                    f1 = np.abs
                    unit = ''
                elif pkey == 'phase':
                    f1 = lambda v: 180./np.pi * np.angle(v)
                    unit = 'degrees, '
                elif pkey == 'dB':
                    f1 = lambda v: 20. * np.log10(v)
                    unit = 'dB, '
                else:
                    # should never happen
                    assert False
                    
                for termname in outreq.varlist:
                    term = circuit.termDict[termname]
                    pltfunc(fvec, f1(term.aC_V), 
                            label = 'Term: {0} [{1}]'.format(
                            term.nodeName, unit + term.unit)) 
                if len(outreq.varlist) == 1:
                    plt.ylabel(
                        'Term: {0} [{1}]'.format(term.nodeName, term.unit))
                else:
                    plt.legend()
        if flag:
            plt.show()

        def getvec(termname):
            return circuit.termDict[termname].aC_V

        if self.shell:
            ipython_drop("""
Available commands:
    fvec: frequency vector
    getvec(<terminal>) to retrieve AC result vector
    plt.* to access pyplot commands (plt.plot(x,y), plt.show(), etc.)
""", globals(), locals())




