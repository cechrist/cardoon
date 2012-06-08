
Implementing a new device model
===============================

A current-source approach is adopted in this simulator. Therefore all
models must be implemented using independent and voltage-controlled
current sources. This is not a limitation any circuit component can be
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
	# Device category
	category = "Basic Components"

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

	    # This adds one internal terminal (in addition to any
	    # existing ones). First argument is the internal variable
	    # name and second is the variable unit. Returns terminal index
	    ti1 = self.add_internal_term('i1', 'A')
    
            # Ambient temperature (temp) by default set to 27 C 
            # Calculate temperature-dependent variables (if any)
            # self.set_temp_vars(self.temp)

It is recommended to copy one of the existing device files to a new
file name and use that as a starting point to create a new device.


Module documentation
--------------------

Documentation for the device library catalog goes in the ``Device``
class docstring in reStructuredText (reST) format. Use an underlined
main title, to be included as an entry in the :doc:`device_library`.
If the device contains internal nodes, document the internal topology
here.  Example for diode device::

    """
    Junction Diode
    --------------

    Model based on spice model. Connection diagram::
    
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
    """


Attributes and functions description
------------------------------------

Mandatory attribute: ``devType = 'string'``. Specifies the netlist name
of the device model, for example 'res' for a resistor model. 

Another mandatory attribute is ``category = 'string'``. This is the
broad category to classify the current model in the
:doc:`device_library` (for example *Basic Components*). You can use
one of the existing categories or create a new one.

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
  for the device temperature parameter (``temp``).

* The ``process_params(self)`` function is called once the
  external terminals have been connected and the non-default
  parameters have been set. This function may be called multiple times
  for example for paramter sweeps or parameter sensitivity. Make
  sanity checks here. Internal terminals/devices must also be
  connected here (see next section).


Internal Terminals, Local References and Terminal Units
-------------------------------------------------------

Some models in addition to the external port voltages require
additional independent variables that can be be obtained by defining
internal terminals. For example, an inductor can be implemented using
current sources as shown below::

        0  o---------+            +----------------+ til
                     | til-tref   |                |
          +         /|\          /^\               |
        Vin        ( | )        ( | ) Vin        ----- L
          -         \V/          \|/             -----
                     |            |                |
        1  o---------+            +----------------+
                                          |
                                         --- tref
                                          V

The additional variable is the inductor current, which in this circuit
can be obtained as ``til - tref``. Here Node ``tref`` is used as a
local reference. Internal references are merged with the global
reference in nodal analysis and so do not add additional unknowns.
Both nodes ``til`` and ``tref`` are implemented using internal
terminals. Note that terminals in a device are internally numbered
consecutively after external terminals. If a model has 2 external
terminals (i.e., 0 to 1), the first internal terminal would be 3.
Internal terminals are normally created in ``process_params()`` as
follows::

	# This adds one internal terminal. Assume only 2 external
	# terminals are connected so far
	til = self.add_internal_term('i1', 'A') # til = 2
	# Add local reference terminal
	tref = self.add_reference_term() # tref = 3

The first argument in ``add_internal_terms()`` is the internal
variable name and second is the variable unit. Internal terminals can
be directly accessed from the terminal list of the device
(``self.neighbour``). The return value is the internal terminal index.
For models that are used as a base class for other devices such as
electrothermal models or extrinsic models, the number of external
terminals may change. For that reason it is *strongly recommended* to
use the return value from ``add_internal_terms()`` and
``add_reference_term()`` instead of fixed numbers. Example from BJT
model::

       # rb is not zero: add internal terminals
       tBi = self.add_internal_term('Bi', 'V')
       tib = self.add_internal_term('ib', '{0} A'.format(glVar.gyr)) 
       tref = self.add_reference_term()
       # Linear VCCS for gyrator(s)
       self.linearVCCS = [((1, tBi), (tib, tref), glVar.gyr),
                          ((tib, tref), (1, tBi), glVar.gyr)]


Terminals have an attribute called ``unit``.  The unit of any existing
terminal variable can be manually changed as follows::

        # Set unit for terminal 6
        self.neighbour[6].unit = 'C'


Temperature Dependence
----------------------

As previously described, most devices should have a ``temp`` parameter.
Compared with regular parameters, temperature is specially treated: by
default all devices take the global temperature defined in the
".options" card. This can be overriden by the device ".model" line. In
turn that is overriden by the temperature specified in the element
line itself. For electrothermal devices, this parameter is ignored and
the temperature at the thermal port is used. All temperatures are
specified in degrees C.

Temperature-related code is included in the following (optional)
function::

    def set_temp_vars(self, temp):
        """
        Calculate temperature-dependent variables for temp given in C

	temp: temperature in degree C
        """
        # if device is based on cppaddev, make sure tape is re-generated
        # ad.delete_tape(self)

        # Absolute temperature 
        T = const.T0 + temp
        # Thermal voltage
        self.Vt = const.k * T / const.e

Note that linear devices may be temperature-dependent. In that case
this function would modify the conductances and capacitances in
``linearVCCS`` and ``linearVCQS`` lists.  This function may be called
multiple times and is used to auto-generate electrothermal models.

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

This function should work when given for both scalar and vector
frequencies. It should take advantage of the vectorization
facilities in numpy.  This interface is still experimental and may
change.

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

* Time-delayed port voltages (``delayedContPorts``): optional, if
  ``needsDelays`` is ``True``, list port voltages in triplet format::

    delayedContPorts = [(n1, n2, delay1), (n3, n4, delay2)]

Similar vectors are required for output ports of charge sources
(``qsOutPorts``). Some of these could be empty and can be modified by
``process_params()`` according to parameter values.

* The nonlinear model equations that are dependent on the control port
  voltages are implemented in the following function::

      def eval_cqs(self, vPort, saveOP=False):
          """
          vPort is a vector with control voltages
      
          Returns tuple with two numpy vectors: one for currents and
          another for charges.

          If saveOP = True, return tuple with normal vectors and OP 
          variables 
          """
          # calculation here
          iVec = np.array([i1, i2])
	  qVec = np.array([q1])
          if saveOP:
              # calculate opVars
              return (iVec, qVec, opVars)
          else:
              return (iVec, qVec)

  The ``saveOP`` argument is optional and may be ommitted if it is
  never needed. ``vPort`` contains control port voltages (or state
  variables) in the order defined by ``controlPorts``, followed by any
  voltages defined in ``csDelayedContPorts``.

  The variables in ``iVec`` are first currents following the order
  defined in ``csOutPorts``, in ``qVec`` are the charges defined in
  ``csOutPorts``. If there are no currents/charges, return an empty
  vector.

  To avoid automatic differentiation problems, use the
  ``ad.condassign()`` function provided in cppaddev.py to replace
  conditional (``if``) statements dependent on variables related to
  ``vPort``.

* The following two functions should be present, normally implemented
  by evaluating the AD tape (they run *much* faster than
  ``eval_cqs()``). But we could also implement them manually by other
  means::

     def eval(self, vPort): same as eval_cqs()
     def eval_and_deriv(self, vPort): returns a tuple, (outVec, Jacobian)

  To have those automatically implemented using cppad, add the
  following to the ``Device`` class::

     # Use functions directly from cppaddev (imported as ad)
     eval_and_deriv = ad.eval_and_deriv
     eval = ad.eval

Automatic Eletro-Thermal Models
+++++++++++++++++++++++++++++++

Automatic electrothermal model generation allows to implement one
nonlinear model with two different netlist names: the normal one with
electrical terminals only (e.g., "bjt") and an electrothermal model
that has an additional pair of thermal terminals. The voltage in this
thermal port is the difference between the device temperature and the
ambient temperature. The current is proportional to the power
dissipated in the device. The netlist name for the electrothermal
model is formed by adding "_t" to the original name (e.g., "bjt_t").

To implement an automatic electrothermal model, set the following
attribute::

    makeAutoThermal = True

The ``process_params()`` function must be modified to accept an
additional argument as follows::

    def process_params(self, thermal = False):
        # Set flag to re-add thermal port
        self.__addThermalPorts = True
        ...
	self.csOutPorts = [(tBi, 2), (tBi, 0), (0, 2), (tref, tib)]
        self.controlPorts = [(tBi, 2), (tBi, 0), (tib, tref)]
        ...
        if not thermal:
            # Calculate temperature-dependent variables
            self.set_temp_vars(self.temp)

The ``thermal`` flag is set to ``True`` for electrothermal devices. In
this example the temperature-dependent variables are not calculated
during parameter processing if ``thermal == True`` since this
calculation would be redundant. Set the ``__addThermalPorts`` flag to
``True`` in this function if one of ``csOutPorts`` or ``controlPorts``
is changed/reassigned.

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

Must provide the following arguments/functions: 

1. At least one (perhaps more) of the source flags set to ``True``::

        # isDCSource = True
        # isTDSource = True
        # isFDSource = True

2. The ``sourceOutput`` argument that contains tuple with output
   port. Voltage sources are implemented using a gyrator and a current
   source. Example::

     sourceOutput = (0, 1) # for a current source

3. Implement at least one of the following source-related functions::

       def get_DCsource(self):
           # used if isDCSource = True
           # return current value
	   pass
    
       def get_TDsource(self, ctime):
           """
           ctime is the current time
           """
           # used if isTDSource = True
           # return current at ctime
	   pass
      
       def get_FDsource(self):
           """
	   Returns a tuple with a frequency and a current phasor vectors

	   	   (fvec, currentVec)

           """
           # used if isFDSource = True. fvec is defined by the source
           # parameters. 

	   # Example for cos wave:
	   fvec = np.array([self.freq])
	   currentVec = np.array([self.magnitude], dtype=complex)
	   return (fvec, currentVec)

   These functions are used with the following conventions:

     * The DC component is the only one that is active for OP or DC
       analyses. 

     * The DC component is always added to the contribution of the
       other sources. Do not include DC components in the other
       functions.

     * Some analyses (such as some forms of envelope-following) may
       require combined time/frequency or multiple time
       dimensions. The interface may have to be extended to handle
       that. The safest approach seems to be to define a new function
       for each case.

   Optionally, some time-domain sources may implement the following
   function to help controlling time-step size::
   
       def get_next_event(self, ctime):
           """
           Returns time of next discontinuity in function/derivative
           """
           pass

   Also optionally, frequency-domain sources may implement the
   following function to be used for AC analysis::

       def get_AC(self):
           """ 
           Returns AC magnitude and phase
           """
           return cm.rect(self._acmag, self._phase)


Linear frequency-defined 
------------------------

If the attribute ``isFreqDefined = True``, then the model must
also include the following attribute with the port definitions for the
frequency-domain part of the device::

  fPortsDefinition = [(0, 1), (2, 3)]

The format of this list is one tuple per port. In the example above,
there are two ports. The positive terminals are 0 and 2. The other
terminals, 1 and 3 are (local) port references.

The Y/G parameters are calculated in the following functions::

    def get_Y_matrix(self, fvec):
        """
        Documentation 
	fvec is a frequency vector/scalar, but frequency can not be zero
        """
        # For scalar fvec returns Y matrix
	# For vector should return 3-D np.array. The frequency
        # index is the last.
        return ymatrix
    

    def get_G_matrix(self):
        """
        Returns a matrix with the DC G parameters
	"""	
	return ymatrix

``get_ymatrix()`` should work when given for both scalar and vector
frequencies and should take advantage of the vectorization facilities
in numpy. It may not work at DC, that is why ``get_gmatrix()`` is also
needed.

Devices Package Documentation
-----------------------------
    
.. automodule:: devices
       :members:


Re-Generating Catalogs
----------------------

Erase one catalog file (for example, ``device_library.rst``) in the
documentation directory (``doc/``) and re-make the documentation. All
catalogs should be automatically re-generated.
