"""
:mod:`tlinps4` -- Physical 4-terminal transmission line (S-matrix based)
------------------------------------------------------------------------

.. module:: diode
.. moduleauthor:: Carlos Christoffersen

"""

import numpy as np
import circuit as cir
from globalVars import const, glVar
# For automatic differentiation:
#import cppaddev as ad

class Device(cir.Element):
    """
    4-terminal physical transmission line model using scattering
    parameters::

             0 o===================================o 2
                               Z0
             1 o===================================o 3


    This model is similar to tlinpy4, but it is more robust and can
    handle lossless lines, even at DC, but internally requires 2
    additional ports to keep track of v1+ and v2+. This model is more
    suitable for convolution as the S parameters are better behaved
    than the Y parameters.

    Netlist Examples::

      tlinps4:tl1 in gnd out gnd z0mag=100. length=0.3m
      .model c_line tlins4 (z0mag=75.00 k=7 fscale=1.e10 alpha = 59.9)

    Internal Topology
    +++++++++++++++++

    The model is symmetric. The schematic for Port 1 is shown here::

               I1                              v1+ + v1-          v1-
              --->                               ---->     v1+   ---->
          0 o--------+                          +------------+----------+  4
       +             |                          |            |          |  
                     |                          |           +-+  s12 v2+|  
      V1            -|- (v1+ - s12 v2+)/Z0     ---          | |        -|- 
                   ( V )                      ( ^ )       1 | |       ( V )
       -            ---                    V1  -|-          +-+        --- 
                     |                          |            |          |  
          1 o--------+                          +------------+----------+  6
                                                             |
                                                            ---
                                                             -

    Note: for a matched transmission line, s11 = s22 = 0 and s12 =
    s21. The equivalent 'Y' matrix is::

               [              1/Z0    -s12/Z0 ]
               [                              ]
               [             -s21/Z0    1/Z0  ]           
           Y = [                              ]
               [ -1            1        s12   ]
               [                              ]
               [        -1    s21        1    ]

    """

    # devtype is the 'model' name
    devType = "tlinps4"
    numTerms = 4
    isFreqDefined = True
    fPortsDefinition = [(0, 1), (2, 3), (4, 6), (5, 6)]

    # Define parameters in a dictionary as follows: parameter name is
    # the key. Parameters are converted to class attributes after
    # circuit initialization.  If model is dependent on temperature,
    # include 'cir.Element.tempItem' ('temp' parameter description).
    paramDict = dict(
        k = ('Effective relative dielectric constant', '', float, 1.),
        alpha = ('Attenuation', 'dB/m', float, 0.1),
        z0mag = ('Magnitude of characteristic impedance', 'Ohms', float, 50.),
        fscale = ('Scaling frequency for attenuation', 'Hz', float, 0.),
        tand = ('Loss tangent', '', float, 0.),
        length = ('Line length', 'm', float, 0.1)
        )

    def __init__(self, instanceName):
        # Here the Element constructor must be called. Do not connect
        # internal nodes here.
        cir.Element.__init__(self, instanceName)
        # Ambient temperature (temp) by default set to 27 C 
        # Add statements as needed


    def process_params(self, circuit):
        """
        Process parameters

        circuit: circuit instance that contains the element

        This should be called each time parameter values are changed.
        """
        # Called once the external terminals have been connected and the
        # non-default parameters have been set. Make sanity checks
        # here. Internal terminals/devices should also be defined here
        # (use circuit reference for this).  Raise cir.CircuitError if a fatal
        # error is found.

        # Use the following to make sure connections to internal
        # terminals are not repeated if this process_params is called
        # many times. 
        self.clean_internal_terms(circuit)
        # Calculate the capacitance per unit length
        # k = epsilon_eff
        self.c = np.sqrt(self.k) / (self.z0mag * const.c0)

        # Calculate the inductance per unit length
        self.l = (self.z0mag * self.z0mag) * self.c

        # convert alpha from db/m to nepers/m
        self.alpha_nepers = self.alpha / const.Np2dB

        # Add internal terminals
        termlist = [self.nodeName + ':n4', self.nodeName + ':n5', 'gnd']
        circuit.connect_internal(self, termlist)

        # Calculate temperature-dependent variables
        # self.set_temp_vars(self.temp)
        # if device is based on cppaddev, make sure tape is re-generated
        # ad.delete_tape(self)


    def get_OP(self, vPort):
        """
        Calculates operating point information
    
        Input:  vPort = [v1, v2]

        Output: dictionary with OP variables

        The frequency-domain model is always used for this calculation.
        """
        ydc = self.get_dc_ymatrix()
        iout = np.dot(ydc, vPort)
        self.OP = dict(
            V1 = vPort[0],
            V2 = vPort[1],
            I1 = iout[0],
            I2 = iout[1]
            )
        return self.OP

        

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
    def get_ymatrix(self, freq):
        """
        Returns a 3-D matrix with Y parameters

        freq: frequency vector/scalar. **Frequency can not be zero**
        """
        # The following calculation is vectorized
        omega = 2. * np.pi * freq

        if self.fscale:
            alphaf = self.alpha_nepers * np.sqrt(freq / self.fscale)
        else:
            alphaf = self.alpha_nepers
        r = 2.0 * alphaf * self.z0mag
        g = self.tand * omega * self.c
        # define imaginary unit
        j = complex(0., 1.)
        z1 = r + j * self.l * omega
        y1 = g + j * self.c * omega

        # characterestic impedance
        # Zo = sqrt( (R + jwL)/(G + jwC)) (calculate inverse here)
        invZ0 = np.sqrt(y1 / z1)
        # attenuation
        # gamma = sqrt( (R + jwL)*(G + jwC))
        gamma = np.sqrt(z1 * y1)
         
        s12 = np.exp(-gamma * self.length)
        # This is required to handle both scalar and vector fvec:
        try:
            nfreqs = np.shape(freq)[0]
            y = np.zeros((4, 4, nfreqs), dtype = type(invZ0[0]))
        except IndexError:
            nfreqs = 1
            y = np.zeros((4, 4), dtype = type(invZ0))
        # proper Y-part block
        y[0,2] = invZ0
        y[0,3] = -invZ0*s12
        y[1,2] = -invZ0*s12
        y[1,3] = invZ0
        # Scattering circuit block
        y[2,0] = -1.
        y[2,2] = 1.
        y[2,3] = s12
        y[3,1] = -1.
        y[3,2] = s12
        y[3,3] = 1.

        # return 3-D np.array
        return y


    def get_dc_ymatrix(self):
        """
        Returns a matrix with the DC Y parameters

        """
        if self.fscale:
            s12 = 1.
        else:
            s12 = np.exp(-self.alpha_nepers * self.length)
        # Use impedance magnitude as reference
        invZ0 = 1. / self.z0mag

        y = np.zeros((4, 4), dtype = type(invZ0))
        # proper Y-part block
        y[0,2] = invZ0
        y[0,3] = -invZ0*s12
        y[1,2] = -invZ0*s12
        y[1,3] = invZ0
        # Scattering circuit block
        y[2,0] = -1.
        y[2,2] = 1.
        y[2,3] = s12
        y[3,1] = -1.
        y[3,2] = s12
        y[3,3] = 1.

        # Return DC Y matrix
        return y



# Here you can add additional functions and classes that only are
# visible withing this module.

