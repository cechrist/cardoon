"""
:mod:`vccs` -- Voltage-controlled current source
------------------------------------------------

.. module:: vccs
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    r"""
    Voltage-controlled current source
    ---------------------------------

    Schematic::

                   g Vc   (or if nonlinear, i(vc))
                   ,---,    
        0 o-------( --> )---------o 1
                   `---`     


        2 o      +  Vc   -        o 3

    By default the source is linear. If a nonlinear function is
    provided, the linear gain (g) can not be specified and is not
    used.

    Netlist examples::

        vccs:g1 gnd 4 3 gnd g=2mS
        vccs:iout 0 cout 1 0 f='1e-3 * np.tanh(vc)' 

    """
    # Device category
    category = "Controlled sources"

    # devtype is the 'model' name
    devType = "vccs"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 4
    
    paramDict = dict(
#        cir.Element.tempItem,
        g = ('Linear transconductance', 'S', float, 1e-3),
        f = ('Nonlinear function i(vc)', 'A', str, '')
        )

    def __init__(self, instanceName):
        """
        Here the Element constructor must be called. Do not connect
        internal nodes here.
        """
        cir.Element.__init__(self, instanceName)
        

    def process_params(self):

        if self.is_set('f'):
            # Nonlinear source
            if self.is_set('g'):
                raise cir.CircuitError(
                    '{0}: can not specify both g and f'.format(self.nodeName))
            # test f expression to make sure it is valid
            try:
                vc = .5
                result = eval(self.f)
            except Exception as e:
                raise cir.CircuitError(
                    '{0}: Invalid expression: {1} ({2})'.format(
                        self.nodeName, self.f, e))
            try:
                abs(result)
            except TypeError:
                raise cir.CircuitError(
                    '{0}: Invalid expression: {1} result not a number)'.format(
                        self.nodeName, self.m))
            # Set nonlinear attributes
            self.isNonlinear = True
            self.controlPorts = [(2, 3)]
            self.csOutPorts = [(0, 1)]
            self.qsOutPorts = []

        else:
            # linear source
            self.linearVCCS = [((2,3), (0,1), self.g)]
        

    def eval_cqs(self, vPort):
        """
        Returns current given control voltage. Charge vector is empty

        vPort[0] = control voltage
        iout[0] = output current
        """
        vc = vPort[0]
        iout = np.array([eval(self.f)])
        return (iout, np.array([]))

    # Use automatic differentiation for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval

    def get_OP(self, vPort):
        """
        Calculates operating point information

        vPort[0] = control voltage
        """
        # Get g at the requested temperature (in case of thermal resistor)
        (iout, qout) = self.eval_cqs(vPort)
        self.OP = {'vc': vPort[0], 
                   'i': iout[0]}
        return self.OP




