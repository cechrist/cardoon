"""
:mod:`1port` -- 1-Port network (S-parameter based)
--------------------------------------------------

.. module:: 1port
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
from scipy.interpolate import interp1d
import touchstone as ts
import cardoon.circuit as cir
from cardoon.globalVars import const, glVar
# For automatic differentiation:
#import cppaddev as ad

class Device(cir.Element):
    """1-Port network (S-parameter based)
    ----------------------------------

    This is a proof of concept model that should be generalized for
    n-ports. Parameters are read from a touchstone file (.s1p).

    By default, network is considered an open circuit outside the
    frequency range defined in the touchstone file.  Optionally, the
    network can be extended as a series or parallel RLC circuit.

    Netlist Example::

      1port:t1 vout gnd s1p=trans1.s1p cap=1pF

    Internal Topology
    +++++++++++++++++

    The internal schematic is the following::

               I                               v+ + v-    Term:   v-
              --->                               ---->      vp   ---->
          0 o--------,                          ,------------+----------,  2
       +             |                          |            |          |  
                     |                          |           ,-,        ,-, 
       V            /|\ vp (1-s11)/Z0          /^\          | |        | | 
                   ( | )                      ( | )       1 | |     s11| | 
       -            \V/                    V1  \|/          '-'        '-' 
                     |                          |            |          |  
          1 o--------+                          +---------+--+----------'   
                                                          |
                                                         --- lref    
                                                          V

    Internal terminal name: vp (keeps track of :math:`v^+`) 

    Note: Current in internal loop gets scaled by the global gyr
    variable.

    """
    # Device category
    category = "Distributed components" # Nport would go here

    # devtype is the 'model' name
    devType = "1port"
    numTerms = 2
    isFreqDefined = True
    fPortsDefinition = [(0, 1), (2, 3)]

    # Define parameters in a dictionary as follows: parameter name is
    # the key. Parameters are converted to class attributes after
    # circuit initialization.  If model is dependent on temperature,
    # include 'cir.Element.tempItem' ('temp' parameter description).
    paramDict = dict(
        s1p = ('Touchstone (s1p) file', '', str, ''),
        cap = ('Capacitor value (for extended frequency)', 'F', float, 0.),
        ind = ('Inductor value (for extended frequency)', 'H', float, 0.),
        res = ('Resistor value (for extended frequency)', 'Ohms', float, 0.),
        series = ('Extend as a series RLC', '', bool, True)
#        dcg = ('DC conductance', 'S', float, 0.)
        )

    def __init__(self, instanceName):
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed


    def process_params(self):
        """
        This should be called each time parameter values are changed.
        """
        # Called once the external terminals have been connected and
        # the non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined
        # here.  Raise cir.CircuitError if a fatal error is found.

        # Use the following to make sure connections to internal
        # terminals are not repeated if this process_params is called
        # many times. 
        self.clean_internal_terms()

        # Re-read touchstone file
        data = ts.touchstone(self.s1p)
        fvec, svec = data.get_sparameter_arrays()
        svec = svec[:,0,0]
        self.z0 = float(data.reference[0])

        # Create interpolation function
        self.sfunc = interp1d(fvec, svec)
        self.fmin = fvec[0]
        self.fmax = fvec[-1]
        
        # Add internal terminals
        self.add_internal_term('vp', 'V') # 3
        self.add_reference_term()         # 4

        def get_s_series(f):
            omega = 2. * np.pi * f
            j = complex(0., 1.)
            real1 = self.cap * (self.res - self.z0) * omega
            real2 = self.cap * (self.res + self.z0) * omega
            imag = omega**2 * self.ind * self.cap - 1.
            return (real1 + j * imag) / (real2 + j * imag)

        def get_s_parallel(f):
            omega = 2. * np.pi * f
            jZ0R = complex(0., 1.) * self.z0 * self.res
            real1 = self.ind * (self.res - self.z0) * omega
            real2 = self.ind * (self.res + self.z0) * omega
            imag = omega**2 * self.ind * self.cap - 1.
            return (real1 + jZ0R * imag) / (real2 + jZ0R * imag)
        
        if self.series:
            # Extend frequency using series RLC
            self.extend_s = get_s_series
        else:
            # Extend frequency using parallel RLC
            self.extend_s = get_s_parallel

        #debug:
        #fv=np.array([1e6,10e6,20e6,21e6])
        #print(self.get_Y_matrix(fv))

        
        # Calculate temperature-dependent variables
        # self.set_temp_vars(self.temp)


    def get_OP(self, vPort):
        """
        Calculates operating point information
    
        Input:  vPort = [v]

        Output: dictionary with OP variables

        The frequency-domain model is always used for this calculation.
        """
        ydc = self.get_dc_ymatrix()
        iout = np.dot(ydc, vPort)
        opDict = dict(
            V1 = vPort[0],
            I1 = iout[0],
            )
        return opDict

        

    #---------------------------------------------------------------------
    # Noise: in general requires a previous call to get_OP 
    #---------------------------------------------------------------------
    def get_noise(self, freq):
        """
        Return noise spectral density at frequency freq
        
        freq may be a scalar/vector
        Not implemented
        """
        pass

    #---------------------------------------------------------------
    # Linear frequency-defined device methods (isFreqDefined = True)
    #---------------------------------------------------------------
    def get_Y_matrix(self, freq):
        """
        Returns a 3-D matrix with Y parameters

        freq: frequency vector/scalar. **Frequency can not be zero**
        """
        # The following calculation is vectorized.
        # This is required to handle both scalar and vector fvec:
        try:
            nfreqs = np.shape(freq)[0]
            select = (freq >= self.fmin) * (freq <= self.fmax)
            # Only interpolate for selected frequencies
            s11 = np.zeros(freq.shape, dtype = 'complex128')
            s11[select] = self.sfunc(freq[select])
            if freq[-1] > self.fmax:
                nsel = freq > self.fmax
                s11[nsel] = self.extend_s(freq[nsel])
            y1 = (1. - s11) / self.z0
            y = np.zeros((2, 2, nfreqs), dtype = type(y1[0]))
        except IndexError:
            #import pdb; pdb.set_trace()
            if (freq < self.fmin) or (freq > self.fmax):
                s11 =  self.extend_s(freq)
            else:
                s11 = self.sfunc(freq)
            y1 = (1. - s11) / self.z0
            y = np.zeros((2, 2), dtype = type(y1))

        # external terminals block
        y[0,1] = y1
        # Internal subcircuit block
        y[1,0] = -glVar.gyr
        y[1,1] = glVar.gyr * (1. + s11)

        # return 3-D np.array
        return y


    def get_G_matrix(self):
        """
        Returns a matrix with the DC G parameters
        """
        y = np.zeros((2, 2), dtype = float)
        if self.series:
            s11 = 1.
        else:
            s11 = -1.

        #s11 = (1. + self.dcg * self.z0) / (1. - self.dcg * self.z0)
        # external terminals block
        y[0,1] = (1. - s11) / self.z0
        # Internal subcircuit block
        y[1,0] = -glVar.gyr
        y[1,1] = glVar.gyr * (1. + s11)

        # Return DC Y matrix
        return y




