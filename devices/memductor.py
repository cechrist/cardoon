"""
:mod:`memductor` -- Basic (nonlinear) memductor
-----------------------------------------------

.. module:: memductor
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    r"""
    Memductor
    ---------

    Connection diagram::


                  +  Vin   -      Iin
                _______________  ---->
        0      |_   _   _   _| |         1
          o----| |_| |_| |_| | |-------o    External view
               |_____________|_|
                                
    Device equation:  

    .. math::    

         q(t) = q(\varphi(t))

         \frac{dq}{dt} = \frac{dq}{d\varphi} \frac{d\varphi}{dt}

         I_{in} = W(\varphi) V_{in}

    :math:`W(\varphi)` is the memductance function.

    Netlist example::

        memd:m1 1 0 w = '1e3 * (np.cosh(1e6 * phi)-1.)' 

    Notes: 

      * the memductance function (``W(phi)``) is given as an
        expression in the ``w`` parameter. The independent variable is
        the memductor flux: ``phi``. Constants and mathematical
        functions can also be used in the definition of ``w``.

      * The initial flux can be adjusted with the ``phi0`` parameter

      * the memductor loses its memory as the capacitor discharges
        through Rleak (Rleak is necessary to ensure a unique DC
        solution). The values of C and Rleak can be adjusted to change
        the time constant

      * The capacitor value has no effect on the memductance, but has
        an effect in the internal model: a larger capacitor will
        produce lower voltages at vc.

    Internal Topology
    +++++++++++++++++

    The internal implementation uses a gyrator and adds one internal
    node: ``vc``. The voltage at ``vc`` is equal to ``(gyr/C) * phi``,
    where ``gyr`` is a global variable that can be changed with the
    ``.options`` keyword::


              --> Iin      
        0  o---------+     
                     | 
          +         /|\  i = w(phi) * Vin     
        Vin        ( | ) 
          -         \V/  phi = (C/gyr) * vc
                     |     
        1  o---------+     
                           
                                     Term: vc                  
        +       +----------------+--------+---------,
                |                |        |         |  
               /^\             -----      /        /^\       
        vc    ( | ) gyr Vin    ----- C    \ Rleak ( | ) phi0 * gyr / C / Rleak
               \|/               |        /        \|/     
                |                |        |         |       
        -       +----------------+--------+---------'     
                                 |                                 
                                --- tref                           
                                 -            

    """
    # Device category
    category = "Basic components"

    # devtype is the 'model' name
    devType = "memd"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    isNonlinear = True
    
    paramDict = dict(
        w = ('Memductance function W(phi)', 'Siemens', str, 'abs(1e-3*phi)'),
        phi0 = ('Initial flux', 'Vs', float, 0.),
        c = ('Auxiliary capacitance', 'F', float, 10e-6),
        rleak = ('Leackage resistance', 'Ohms', float, 1e9)
        )

    def __init__(self, instanceName):
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)


    def process_params(self):
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # remove any existing internal connections
        self.clean_internal_terms()

        # Test parameters
        if not self.rleak:
            raise cir.CircuitError(self.nodeName 
                                   + ': leackage resistance can not be zero')
        if not self.c:
            raise cir.CircuitError(self.nodeName 
                                   + ': capacitance can not be zero')
        # test m expression to make sure it is valid
        try:
            phi = .5
            result = eval(self.w)
        except Exception as e:
            raise cir.CircuitError(
                '{0}: Invalid expression: {1} ({2})'.format(self.nodeName, 
                                                            self.w, e))
        try:
            abs(result)
        except TypeError:
            raise cir.CircuitError(
                '{0}: Invalid expression: {1} (result not a number)'.format(
                    self.nodeName, self.w))

        # Connect internal terminal
        tvc = self.add_internal_term('vc', 'V') 
        tref = self.add_reference_term() 
        # Set up source if phi0 is given
        if self.phi0:
            self.isDCSource = True
            self.sourceOutput = (tref, tvc)
            self._i0 = self.phi0 * glVar.gyr / self.c / self.rleak

        # Setup gyrator
        # Access to global variables is through the glVar 
        self.linearVCCS = [((0,1), (tref, tvc), glVar.gyr), 
                           ((tvc, tref), (tvc, tref), 1./self.rleak)]
        self.linearVCQS = [((tvc, tref), (tvc, tref), self.c)]

        self.controlPorts = [(0,1), (tvc, tref)]
        self.csOutPorts = [(0,1)]
        self.qsOutPorts = []


    def eval_cqs(self, vPort, saveOP = False):
        """
        Returns memductor current given input voltage. Charge vector
        is empty

        vPort[0] = memductor voltage
        vPort[1] = internal cap voltage
        iout[0] = memductor current
        """
        phi = self.c * vPort[1] / glVar.gyr
        W = eval(self.w)
        iout = np.array([W * vPort[0]])
        if saveOP:
            opVars = np.array([W])
            return (iout, np.array([]), opVars)
        else:
            return (iout, np.array([]))

    # Use automatic differentiation for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval
    get_op_vars = ad.get_op_vars

    def get_OP(self, vPort):
        """
        Calculates operating point information

        vPort[0] = memductor voltage
        vPort[1] = internal cap voltage
        """
        (iout, qout) = self.eval_cqs(vPort)
        opV = self.get_op_vars(vPort)
        self.OP = {'v': vPort[0],
                   'i': iout[0],
                   'phi': self.c * vPort[1] / glVar.gyr, 
                   'W': opV[0]}
        return self.OP

    def get_DCsource(self):
        return self._i0
