"""
:mod:`memristor` -- Linear memristor
------------------------------------

.. module:: memristor
.. moduleauthor:: Carlos Christoffersen
"""

import numpy as np
from globalVars import const, glVar
import circuit as cir
import cppaddev as ad

class Device(cir.Element):
    """
    Memristor
    ---------

    Connection diagram::

                _______________
        0      |_   _   _   _  |         1
          o----+ |_| |_| |_| |_+-------o    External view
               |_______________|
                                
    Netlist example::

        mem:m1 1 0 m = '1e3 * (np.cosh(1e6 * q)-1.)' 

    Notes: 

      * The memristance function is given as an expression in the
        ``m`` parameter. Constants and mathematical functions can be
        used. The independent variable is ``q``.

      * the memristor loses its memory as the capacitor discharges
        through Rleak (Rleak is necessary to ensure a unique DC
        solution). The values of C and Rleak can be adjusted to change
        the time constant.

    Internal Topology
    +++++++++++++++++

    Internal implementation uses gyrators (adds 2 internal nodes)::

                                        im/gyr    Term: im
        0  o---------+            +----------------+
                     | gyr V(im)  |                |
          +         /|\          /^\              /|\ 
        Vin        ( | )        ( | ) gyr Vin    ( | ) gyr^2 * M(q) * V(im)
          -         \V/          \|/              \V/ 
                     |            |                |   q = C * vc 
        1  o---------+            +----------------+
                                          |
                                         --- tref 
                                          V 

                                        vc      Term: vc                  
                                  +----------------+--------,
                                  |                |        |                
                                 /^\             -----      /                
                                ( | ) gyr V(im)  ----- C    \ Rleak
                                 \|/               |        /                
                                  |                |        |                
                                  +----------------+--------'               
                                          |                                 
                                         --- tref                           
                                          V                                 

    """
    # Device category
    category = "Basic Components"

    # devtype is the 'model' name
    devType = "mem"

    # Number of terminals. If numTerms is set here, the parser knows
    # in advance how many external terminals to expect. By default the
    # parser makes no assumptions and allows any number of connections
    #
    numTerms = 2
    isNonlinear = True
    
    paramDict = dict(
        m = ('Memristance function M(q)', 'Ohms', str, 'abs(5e9*q)'),
        c = ('Auxiliary capacitance', 'F', float, 1e-6),
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
            q = 1.
            result = eval(self.m)
        except Exception as e:
            raise cir.CircuitError(
                '{0}: Invalid expression: {1} ({2})'.format(self.nodeName, 
                                                            self.m, e))
        try:
            abs(result)
        except TypeError:
            raise cir.CircuitError(
                '{0}: Invalid expression: {1} (result not a float)'.format(
                    self.nodeName, self.m))

        # Connect internal terminal
        tim = self.add_internal_term('im', '{0} A'.format(glVar.gyr)) 
        tvc = self.add_internal_term('vc', 'V') 
        tref = self.add_reference_term() 
        # Setup gyrator
        # Access to global variables is through the glVar 
        self.linearVCCS = [((0,1), (tref, tim), glVar.gyr), 
                           ((tim, tref), (0,1), glVar.gyr),
                           ((tim, tref), (tref, tvc), glVar.gyr),
                           ((tvc, tref), (tvc, tref), 1./self.rleak)]
        self.linearVCQS = [((tvc, tref), (tvc, tref), self.c)]

        self.controlPorts = [(tim, tref), (tvc, tref)]
        self.csOutPorts = [(tim, tref)]
        self.qsOutPorts = []


    def eval_cqs(self, vPort):
        """
        Returns memristor voltage given current. Charge vector is empty

        vPort[0] = memristor current / gyr
        vPort[1] = internal cap voltage
        iout[0] = gyr * memristor voltage
        """
        q = self.c * vPort[1]
        normM = glVar.gyr**2 * eval(self.m)
        iout = np.array([normM * vPort[0]])
        return (iout, np.array([]))

    # Use automatic differentiation for eval and deriv function
    eval_and_deriv = ad.eval_and_deriv
    eval = ad.eval

    def get_OP(self, vPort):
        """
        Calculates operating point information

        vPort[0] = memristor current / gyr
        vPort[1] = internal cap voltage
        """
        # Get g at the requested temperature (in case of thermal resistor)
        (iout, qout) = self.eval_cqs(vPort)
        self.OP = {'q': self.c * vPort[1], 
                   'i': iout[0]}
        return self.OP
