"""
:mod:`dc` -- DC sweep Analysis
------------------------------

.. module:: dc
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import numpy as np

from paramset import ParamSet
from analysis import AnalysisError, ipython_drop
from fsolve import solve, NoConvergenceError
from globalVars import glVar
import matplotlib.pyplot as plt

class Analysis(ParamSet):
    """
    DC Sweep
    --------

    Calculates a DC sweep of a circuit using the nodal approach. After
    the analysis is complete, nodal voltages are saved in circuit and
    terminals with the ``dC_`` prefix.  After this the analysis drops
    to an interactive shell if the ``shell`` global variable is set to
    ``True``.

    The following parameters can be swept: 

      * Any device parameter of ``float`` type (device name must be
        specified in this case)

      * Global temperature (no device specified): sweep temperature of
        all devices that do not explicitly have ``temp`` set.

    Convergence parameters for the Newton method are controlled using
    the global variables in ``.options``.

    One plot window is generated for each ``.plot`` statement. Use
    ``dc`` request type for this analysis.

    DC formulation documented in :doc:`analysis`

    Examples::

        # Device parameter sweep
        .analysis dc device=vsin:v1 param=vdc start=-2. stop=2. num=50 

        # Global temperature sweep
        .analysis dc param=temp start=-20C stop=80C 

        # Some options that affect convergence properties
        .options maxiter=300 gyr=1e-5 maxdelta=5.
        
        .plot dc 153 151 23

    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "dc"

    # Define parameters as follows
    paramDict = dict(
        device = ('Instance name of device to sweep variable', '', str, ''),
        param = ('Parameter to sweep', '', str, ''),
        start = ('Sweep start value', '(variable)', float, 0.),
        stop = ('Sweep stop value', '(variable)', float, 0.),
        num = ('Number of points in sweep', '', int, 50),
        verbose = ('Show iterations for each point', '', bool, False),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)


    def run(self, circuit):
        """
        Calculates a DC sweep by solving nodal equations

        The parameter to be swept is specified in the analysis options
        """
        # for now just print some fixed stuff
        print('******************************************************')
        print('                 DC sweep analysis')
        print('******************************************************')
        if hasattr(circuit, 'title'):
            print('\n', circuit.title, '\n')

        if glVar.sparse:
            import nodalSP as nd
        else:
            import nodal as nd
            print('Using dense matrices\n')

        # Only works with flattened circuits
        if not circuit._flattened:
            circuit.flatten()
            circuit.init()

        # get device (if any)
        paramunit = None
        # tempFlag indicates if we are doing a global temperature sweep
        tempFlag = False
        if self.device:
            # Device specified, try to find it
            try:
                dev = circuit.elemDict[self.device]
            except KeyError: 
                raise AnalysisError('Could not find: {0}'.format(self.device))
            if self.param:
                try:
                    pinfo = dev.paramDict[self.param]
                except KeyError:
                    raise AnalysisError('Unrecognized parameter: ' 
                                        + self.param)
                else:
                    if not pinfo[2] == float:
                        raise AnalysisError('Parameter must be float: ' 
                                            + self.param)
            else:
                raise AnalysisError("Don't know what parameter to sweep!")
            paramunit = pinfo[1]
        else:
            # No device, check if temperature sweep
            if self.param != 'temp':
                raise AnalysisError(
                    'Only temperature sweep supported if no device specified')
            paramunit = 'C'
            tempFlag = True

        # Create nodal object: for now assume devices do not change
        # topology during sweep
        nd.make_nodal_circuit(circuit)
        dc = nd.DCNodal(circuit)
        x = dc.get_guess()

        sweepvar = np.linspace(start = self.start, stop = self.stop, 
                               num = self.num)
        circuit.dC_sweep = sweepvar
        if tempFlag:
            circuit.dC_var = 'Global temperature sweep: temp'
        else:
            circuit.dC_var = 'Device: ' + dev.nodeName \
                + '  Parameter: ' + self.param
        circuit.dC_unit = paramunit
        print('System dimension: {0}'.format(circuit.nD_dimension))
        print('Sweep: ', circuit.dC_var)
        
        xVec = np.zeros((circuit.nD_dimension, self.num))
        tIter = 0
        tRes = 0.
        for i, value in enumerate(sweepvar):
            if tempFlag:
                for elem in circuit.nD_elemList:
                    # Only sweep temperature if not explicitly given
                    # for a device
                    if not elem.is_set('temp'):
                        try:
                            elem.set_temp_vars(value)
                        except AttributeError:
                            # It is OK if element independent of temperature
                            pass
            else:
                setattr(dev, self.param, value)
                # re-process parameters (topology must not change, for
                # now at least)
                dev.process_params()
                # Re-generate nodal attributes. 
                nd.restore_RCnumbers(dev)
                nd.process_nodal_element(dev)

            # re-process linear matrix
            dc.refresh()
            sV = dc.get_source()

            # solve equations
            try: 
                (x, res, iterations) = solve(x, sV, dc.convergence_helpers)
            except NoConvergenceError as ce:
                print(ce)
                return
            # Save result
            xVec[:,i] = x
            # Keep some info about iterations
            tIter += iterations
            tRes += res
            if self.verbose:
                print('{0} = {1}'.format(self.param , value))
                print('Number of iterations = ', iterations)
                print('Residual = ', res)

        # Calculate average residual and iterations
        avei = tIter / len(sweepvar)
        aver = tRes / len(sweepvar)
        print('Average iterations: {0}'.format(avei))
        print('Average residual: {0}\n'.format(aver))

        # Save results in nodes
        circuit.nD_ref.dC_v = np.zeros(self.num)
        for i,term in enumerate(circuit.nD_termList):
            term.dC_v = xVec[i,:]

        # Restore original attribute values ------------------------
        if tempFlag:
            for elem in circuit.nD_elemList:
                if not elem.is_set('temp'):
                    elem.set_attributes()
                    try:
                        elem.set_temp_vars(value)
                    except AttributeError:
                        # It is OK if element independent of temperature
                        pass
        else:
            dev.set_attributes()
            # re-process parameters (topology must not change, for
            # now at least)
            dev.process_params()
            # Re-generate nodal attributes. 
            nd.restore_RCnumbers(dev)
            nd.process_nodal_element(dev)
        # -----------------------------------------------------------

        # Process output requests.  In the future this should be moved
        # to a common module that processes output requests such as
        # plot, print, save, etc.
        flag = False
        for outreq in circuit.outReqList:
            if outreq.type == 'dc':
                flag = True
                plt.figure()
                plt.grid(True)
                plt.xlabel('{0} [{1}]'.format(circuit.dC_var, circuit.dC_unit))
                for termname in outreq.varlist:
                    term = circuit.termDict[termname]
                    plt.plot(sweepvar, term.dC_v, 
                             label = 'Term: {0} [{1}]'.format(
                            term.nodeName, term.unit)) 
                if len(outreq.varlist) == 1:
                    plt.ylabel(
                        'Term: {0} [{1}]'.format(term.nodeName, term.unit))
                else:
                    plt.legend()
        if flag:
            plt.show()

        def getvec(termname):
            return circuit.termDict[termname].dC_v

        if self.shell:
            ipython_drop("""
Available commands:
    sweepvar: vector with swept parameter
    getvec(<terminal>) to retrieve results
    plt.* to access pyplot commands (plt.plot(x,y), plt.show(), etc.)
""", globals(), locals())

