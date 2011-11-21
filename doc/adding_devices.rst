
Implementing a new device model
===============================

A current-source approach is adopted in this simulator. Therefore all
models must be implemented using independent and voltage-controlled
sources. This is not a limitation any circuit component can be
described in that way, including ideal voltage sources and
state-variable defined nonlinear models.

To implement a new model, simply create a Python file in the
``devices`` directory. That file constitutes a module to be imported
into the device library. The model itself must be implemented in a
class defined as follows::

    import numpy as np
    import circuit as cir
    # Physical constants, global variables
    from globalVars import const, glVar
    # Automatic differentiation
    import cppaddev as ad
    
    class Device(cir.Element):
        """
	Document new model here
        """

        # devtype is the 'model' name
        devType = "emptydev"
    
        # Number of terminals. 
        numTerms = 2
            
        # Model parameters
        paramDict = dict(
            cir.Element.tempItem,
            r = ('Resistance', 'Ohms', float, 0.),
            rsh = ('Sheet resistance', 'Ohms', float, 0.),
            l = ('Lenght', 'm', float, 0.),
            w = ('Width', 'm', float, 0.),
            )
    
        def __init__(self, instanceName):
            """
	    Initialization code
            """
	    # Do not include here parameter-dependent code
            cir.Element.__init__(self, instanceName)
    
    
        def process_params(self):
            """
	    Prepares the device for simulation
	    
            Raises cir.CircuitError if a fatal error is found
            """
            # if device is based on cppaddev, make sure tape is re-generated
            # ad.delete_tape(self)

            # Use the following to make sure connections to internal
            # terminals are not repeated if this process_params is called
            # many times. 
            self.clean_internal_terms()

	    # This adds 2 internal terminals (in addition to any
	    # existing ones): 
	    self.add_internal_terms(2)
    
            # Ambient temperature (temp) by default set to 27 C 
            # Calculate temperature-dependent variables (if any)
            # self.set_temp_vars(self.temp)

Module documentation
--------------------

Documentation for the device library catalog goes in the ``Device``
class docstring in reStructuredText (reST) format. If the device
contains internal nodes, document the internal topology here.  Example
for diode device::

        Diode device (based on spice model)::
        
                   o  1                           
                   |                            
                 --+--
                  / \     
                 '-+-' 
                   |                          
                   o  0 
    
        Includes depletion and diffusion charges.
    
        Netlist examples::
    
            diode:d1 1 0 isat=10fA cj0=20fF
    
            # Electrothermal device
            diode_t:d2 2 3 1000 gnd cj0=10pF tt=1e-12 rs=100 bv = 4.
    
            # Model statement
            .model dmodel1 diode (cj0 = 10pF tt=1ps)
    

Attributes and functions description
------------------------------------

Mandatory argument: ``devType = 'string'``. Specifies the netlist name
of the device model, for example 'res' for a resistor model

The following attributes are not mandatory and defaut to empty lists
if not specified. They can be used with any type of device model:
linear, nonlinear, frequency-defined.

* If ``numTerms`` is set, the parser knows in advance how many
  external terminals to expect. By default ``numTerms = 0`` and the
  program makes no assumptions and allows any number of connections.

* If internal linear VCCS are needed, they are specified using the
  following format::

    linearVCCS = [((t0, t1), (t2, t3), g), ... ]
  
  
    0  o--------+      +------o 2
                       |      
      +               /|\       
    Vin              | | | g Vin     
      -               \V/       
                       |      
    1  o--------+      +------o 3

  The format consists on a list of tuples, one per voltage-controlled
  current source (VCCS). Each tuple has 2 tuples for the control and
  output ports, respectively and the transconductance goes at the end.

* The same format is used for linear charge sources (VCQS)::

    linearVCQS = [((t0, t1), (t2, t3), c), ... ]

Both ``linearVCCS`` and ``linearVCQS`` may be empty lists and may be
modified by ``process_params()`` according to paramenter
values. Inductors are represented by a combination of VCCS and VCQS
(see inductor model as an example).

* Parameters are listed in a dictionary named ``paramDict`` as shown
  in the sample code. The parameter name is the key. The fields in the
  description tuple are: long description, unit, type, default
  value. The default value can be ``None``. Parameters are converted
  to class attributes after circuit initialization. For this reason
  parameter names can not be Python keywords (unfortunately ``is`` is
  a keyword). If model is dependent on temperature, the first item
  should be ``cir.Element.tempItem``, which contains the description
  for the device temperature parameter (temp).

* The ``process_params(self, circuit)`` function is called once the
  external terminals have been connected and the non-default
  parameters have been set. This function may be called multiple times
  for example for paramter sweeps or parameter sensitivity. Make
  sanity checks here. Internal terminals/devices must also be
  connected here.

Temperature Dependence
----------------------

As previously described, it should have a "temp" parameter.  Compared
with regular parameters, temperature is specially treated: by default
all devices take the global temperature defined in the ".options"
card. This can be overriden by the device ".model" line. In turn that
is overriden by the temperature specified in the element line
itself. For electrothermal devices, this parameter is ignored and the
temperature at the thermal port is used. All temperatures are
specified in degrees C.

Temperature-related code is included in the following (optional)
function::

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C

	temp: temperature in degree C
        """
        # Absolute temperature 
        T = const.T0 + temp
        # Thermal voltage
        self.Vt = const.k * T / const.e

Note that linear devices may be temperature-dependent. In that case
this function would modify the conductances and capacitances in
``linearVCCS`` and ``linearVCQS`` lists.
This function may be called multiple times and may be used to
auto-generate an electrothermal device (described in next section).

Operating Point
---------------

The following function generates a dictionary with operating point
variables should be implemented by all devices. For
frequency-dependent devices, f is assumed to be zero. Variable names
are arbitrary::

   def get_OP(self, vPort):
       """
       Calculates operating point information
   
       Input:  vPort = [vdb , vgb , vsb]
       Output: dictionary with OP variables
       """
       # First we need the Jacobian
       (outV, jac) = self.eval_and_deriv(vPort)
       # if this is not needed then saveOP flag does not have 
       # to be implemented
       opV = self.get_op_vars(vPort) 
   
       # Check things that change if the transistor is reversed
       if opV[11] > 0.:
           reversed = False
           gds = jac[0,0]
       else:
           reversed = True
           gds = jac[0,2]
           
       self.OP = {'VD': vPort[0],
                  'VG': vPort[1],
                  'VS': vPort[2],
                  'IDS': outV[0]}

If the model noise model is dependent on the operating point, this is
the place to calculate the corresponding variables. 


Noise current spectral density sources
--------------------------------------

Same format as ``csOutPorts`` (for nonlinear devices). Default is an
empty tuple.

Example::

  noisePorts = [(1, 2), (0, 2)]

The ``get_noise()`` function in general requires a previous call to
get_OP()::

     def get_noise(self, f):
         """
         Return noise spectral density at frequency f
         
	 f may be a scalar/vector
         Requires a previous call to get_OP() 
         """
         s1 = self.OP['Sthermal'] + self.OP['kSflicker'] / pow(f, self.af)
         s2 = something
         return np.array([s1, s2])

This interface is still experimental and may change.

Nonlinear models
----------------

The following attributes are required for nonlinear models::

  isNonlinear = True
  needsDelays = True or False

An optional attribute, ``vPortGuess`` is a numpy vector with a valid
set of controlling voltages to be used as an initial guess. If this is
not specified, the initial guess is set to zero.

* Current source output ports (``csOutPorts``): for each current
  source in the device, list ports as follows: ``(n1, n2)``. Current
  flows from ``n1`` to ``n2``.
  
  Example for a 3-terminal BJT with BE and CE current sources,
  assuming teminals are connected C (0) - B (1) - E (2)::
  
    csOutPorts = [(1, 2), (0, 2)]

* Controlling ports (``controlPorts``): list here all ports whose
  voltages are needed to calculate the nonlinear currents / charges in
  same format.

  Example for BJT without intrinsic RC, RB and RE (vbc, vbe)::

    controlPorts = [(1, 0), (1, 2)]

* Time-delayed port voltages (``csDelayedContPorts``): optional, if
  ``needsDelays`` is ``True``, list port voltages in triplet format::

    csDelayedContPorts = [(n1, n2, delay1), (n3, n4, delay2)]

Similar vectors are required for output ports of charge sources
(``qsOutPorts``). Some of these could be empty and can be modified by
``process_params()`` according to parameter values.

* The nonlinear model equations that are dependent on the control port
  voltages are implemented in the following function::

      def eval_cqs(self, vPort, saveOP=False):
          """
          vPort is a vector with control voltages
      
          Returns a numpy vector: currents first and then charges.
          If saveOP = True, return tuple with normal vector and OP 
          variables (only needed if ever saveOP is True, see resistor)
          """
          # calculation here
          outVec = np.array([var1, var2])
          if saveOP:
              # calculate opVars
              return (outVec, opVars)
          else:
              return outVec

  The ``saveOP`` argument is optional and may be ommitted if not
  needed. ``vPort`` contains control port voltages (or state
  variables) in the order defined by ``controlPorts``, followed by any
  voltages defined in ``csDelayedContPorts``.

  The variables in ``outVec`` are first currents following the order
  defined in ``csOutPorts``, followed by any charges defined in
  ``csOutPorts``.

  To avoid automatic differentiation problems, use the
  ``ad.condassign()`` function provided in cppaddev.py to replace
  ``if`` statements dependent on variables related to ``vPort``.

* The following two functions should be present, normally implemented
  by evaluating the AD tape (i.e. they run *much* faster than
  ``eval_cqs()``). But we could also implement them manually by other
  means::

     def eval(self, vPort): same as eval_cqs()
     def eval_and_deriv(self, vPort): returns a tuple, (outVec, Jacobian)

  To have those automatically implemented using cppad, add the
  following to the ``Device`` class::

     # Use functions directly from cppaddev (imported as ad)
     eval_and_deriv = ad.eval_and_deriv
     eval = ad.eval

* Automatic electrothermal model generation allows to implement one
  nonlinear model with two different netlist names: the normal one
  with electrical terminals only (e.g., "bjt") and an electrothermal
  model that has an additional pair of thermal terminals. The voltage
  in this thermal port is the temperature and the current is
  proportional to the power dissipated in the device. The netlist name
  for the electrothermal model is formed by adding "_t" to the
  original name (e.g., "bjt_t").

  To implement an automatic electrothermal model, set the following
  attribute::

      makeAutoThermal = True

  In addition, the following function must be implemented::

     def power(self, vPort, currV):
         """ 
         Returns total instantaneous power 
     
         Input: input (vPort) and output vectors in the format from 
	 eval_cqs()
         """
         vds = vPort[0] - vPort[2]
         # pout = vds*ids + vdb*idb + vsb*isb
         pout = vds*currV[0] + vPort[0] * currV[1] + vPort[2] * currV[2] 
         return pout
   
   This function takes the input vector and the results from
   ``eval_cqs()`` and returns the total power dissipated at the
   nonlinear current sources.


Independent Sources
-------------------

Must provide: 

1. At least one (perhaps more) of the source flags set to ``True``::

        # isDCSource = True
        # isTDSource = True
        # isFDSource = True

2. A tuple with output port. Voltage sources are implemented using a
   gyrator and a current source. Example::

     sourceOutput = (0, 1) # for a current source

3. Implement at least one of the source-related functions::

       def get_DCsource(self):
           """
           Documentation (isDCSource = True)
           """
           # return current value
    
       def get_TDsource(self, ctime):
           """
           Documentation (isTDSource = True)
           ctime is the current time
           """
           # return current at ctime
      
       def get_FDsource(self, fvec):
           """
           Documentation (isFDSource = True)
           """
           # should return a np.array with currents for each frequency


Linear frequency-defined 
------------------------

If the attribute ``isFreqDefined = True``, then the model must
also include the following attribute with the port definitions for the
frequency-domain part of the device::

  fPortsDefinition = [(0, 1), (2, 3)]

The format of this list is one tuple per port. In the example above,
there are two ports. The positive terminals are 0 and 2. The other
terminals, 1 and 3 are (local) references.

The Y parameters are calculated in the following functions::

    def get_ymatrix(self, fvec):
        """
        Documentation 
	fvec is a frequency vector/scalar, but frequency can not be zero
        """
        # should return 3-D np.array with Y matrix for each frequency
        return ymatrix
    

    def get_dc_ymatrix(self):
        """
        Returns a matrix with the DC Y parameters
	"""	
	return ymatrix

This interface is still experimental and may change. 

    




