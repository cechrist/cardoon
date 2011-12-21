"""
:mod:`testdev` -- Test equations of one nonlinear device
--------------------------------------------------------

.. module:: testdev
.. moduleauthor:: Carlos Christoffersen

"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from paramset import ParamSet
from analysis import AnalysisError, ipython_drop

class Analysis(ParamSet):
    """
    Test equations of one nonlinear device

    One advantage of using this method over a DC sweep is that no
    Newton iterations are needed. The following internal functions are
    tested here:

    * process_params()
    * set_temp_vars()
    * eval_cqs()
    * eval()
    * get_OP()
    * power() (for electrothermal models)

    After completion the analysis drops to an interactive shell if the
    ``shell`` global variable is set to ``True``
    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "testdev"

    # Define parameters as follows
    paramDict = dict(
        ports_bias = ('Vector with default values of port voltages', 'V', 
                      list, []),
        sweep_port = ('Port number to be swept, starting from zero', '', int, 
                      0),
        start = ('Sweep start value', 'V', float, 0.),
        stop = ('Sweep stop value', 'V', float, 0.),
        sweep_num = ('Number of points in sweep', '', int, 0),
        device = ('Instance name of device to test', '', str, ''),
        plot = ('Auto-plot currents and charges', '', bool, True),
        useAD = ('Use automatic differentiation', '', bool, True),
        param = ('Parameter for outer sweep', '', str, ''),
        param_val = ('Vector with parameter values to sweep', '', list, []),
        shell = ('Drop to ipython shell after calculation', '', bool, False)
        )


    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)


    def run(self, circuit):
        """
        The selected device must be nonlinear. Retrieves the device
        instance from circuit and applies the voltages given by
        port_bias and the combination of sweep_port and
        sweep_val. Then plots currents/charges as a function of swept
        voltage.
        """
        # get device 
        try:
            dev = circuit.elemDict[self.device]
        except KeyError: 
            raise AnalysisError('Could not find: {0}'.format(self.device))
            return

        # Initialize device
        dev.init()

        nports = len(self.ports_bias)
        # Check that parameters make some sense
        if not dev.isNonlinear:
            raise AnalysisError(self.device + ' is linear')
        if nports != len(dev.controlPorts):
            raise AnalysisError('ports_bias for ' + self.device 
                                + ': wrong number of control ports given')
        if self.sweep_num < 1:
            raise AnalysisError('sweep_num = ' + str(self.num)
                                + ': wrong number')
        if self.sweep_port >= nports:
            raise AnalysisError('sweep_port = ' + str(self.port)
                                + ': wrong number (starts from zero)')
        param = False
        paramunit = None
        if self.param:
            try:
                pinfo = dev.paramDict[self.param]
            except KeyError:
                raise AnalysisError('Unrecognized parameter: ' 
                                    + self.param)
            else:
                if not ((pinfo[2] == float) or (pinfo[2] == int)):
                    raise AnalysisError('Parameter must be numeric: ' 
                                        + self.param)
                paramunit = ' ' + pinfo[1]
                param = True
            
        # Print some info about what is being tested
        vports  = np.array(self.ports_bias)
        # Generate operating point info
        if self.useAD:
            dev.get_OP(vports)
        print('******************************************************')
        print('Nonlinear device internal source test analysis')
        print('******************************************************')
        print(dev)
        print('\nPort bias voltages:',  self.ports_bias, 'V')
        print('Inner sweep port:', self.sweep_port, ', start:', 
              self.start, ', stop:', self.stop)

        if param:
            npsweep = len(self.param_val)
            print('Parameter sweep: {0} = {1}'.format(self.param, 
                                                      self.param_val))
        else:
            npsweep = 1
        # Print element variables
        dev.print_vars()

        vsweep = np.linspace(self.start, self.stop, self.sweep_num)
        nsamples = np.shape(vsweep)[0]
        ncurrents = len(dev.csOutPorts)
        iout = None
        if ncurrents: 
            iout = np.zeros((npsweep, nsamples, ncurrents))
        ncharges = len(dev.qsOutPorts)
        qout = None
        if ncharges:
            qout = np.zeros((npsweep, nsamples, ncharges))
        for j in range(npsweep):
            if param:
                setattr(dev, self.param, self.param_val[j])
                # re-process parameters
                dev.process_params()
            for i in range(np.shape(vsweep)[0]):
                vports[self.sweep_port] = vsweep[i]
                if self.useAD:
                    # Use AD tape to evaluate function
                    outvars = dev.eval(vports)
                else:
                    # The one below is slower as it does not use tapes:
                    outvars = np.concatenate(dev.eval_cqs(vports), axis=0)
                # The function below in addition calculates derivatives
                #(outvars, Jac) = dev.eval_and_deriv(vports)
                # Convert current, charge to vectors
                if ncurrents: 
                    iout[j,i,:] = outvars[0:ncurrents]
                if ncharges:
                    qout[j,i,:] = outvars[ncurrents:]

        # Reset device attributes
        dev.clean_attributes()
        dev.set_attributes()
        # Plot currents and voltages if requested
        if self.plot:
            self.plot_all(vsweep, iout, qout, param, npsweep, paramunit)

        if self.shell:
            ipython_drop('', globals(), locals())


    def plot_all(self, vsweep, iout, qout, param, npsweep, paramunit):
        """
        Plot all currents and charges
        """
        label = ''
        var = []
        varlegend = []
        if iout != None:
            ncurrents = iout.shape[0]
        else:
            ncurrents = 0
        if qout != None:
            ncharges = qout.shape[0]
        else:
            ncharges = 0
        if ncurrents:
            var.append(iout)
            varlegend.append(' Current (A)')
        if ncharges:
            var.append(qout)
            varlegend.append(' Charge (C)')
            
        for k in range(len(var)):
            plt.ioff()
            #plt.figure(k+1)
            # sweep for each output var first
            nvars = np.shape(var[k])[2]
            for i in range(nvars):
                #plt.subplot(nvars, 1, i+1)
                plt.figure()
                plt.grid(True)
                # Now go though samples
                for j in range(npsweep):
                    if param:
                        label = self.param + ':'\
                            + str(self.param_val[j]) + paramunit
                    plt.plot(vsweep, var[k][j, :, i], linewidth=1, 
                             label=label)
                plt.ylabel('Out ' + str(i) + varlegend[k])
                plt.xlabel('In ' + str(self.sweep_port) 
                           + ' Voltage (V)')
                if param and not i:
                    plt.legend(loc = 'upper left')
                #if i == 0:
                plt.title('Device ' + self.device )

        plt.show()

# Here you can add additional functions and classes that only are
# visible withing this module.

