
Simulator Design Considerations
===============================

Design History
--------------

The design was originally inspired from the experience with fREEDA
[freeda]_ and carrot [carrot]_ plus some ideas from other
simulators and improvements that take advantage of the flexibility in
Python. 


Agnostic Simulation-Approach Circuit Representation
---------------------------------------------------

The internal circuit representation and the device library attempt to
be agnostic regarding to simulation methods. The ``Circuit`` class has
no methods to obtain nodal voltages or calculate the Jacobian of the
circuit. This is delegated to other classes that handle particular
circuit analysis approaches such as nodal, sparse tableau, etc.

Element Interfaces
++++++++++++++++++

All elements but frequency-defined ones are modeled using the
current-source approach [aplac2]_. The interface essentially describes
a subcircuit composed of independent and voltage-controlled
current/charge sources. If the nodal approach is used this has the
advantage that an analysis can be implemented by just considering a
few fundamental components. For example there is no need for MNAM
stamps as the NAM can handle all elements.  This approach is also very
flexible. Both fREEDA-style state-variable and spice-style nonlinear
models have been implemented.

Frequency-defined elements are modeled by their Y matrix (or
equivalent, see S-parameter-based transmission line model). The
interface returns a 3-D matrix with the parameters at all frequencies.
It is possible to conceive a more compact hybrid representation for
some devices that combines the current-source approach with smaller Y
matrices. However this may unnecessarily complicate the implementation
of simpler devices. This interface is being reviewed and may change in
the future. 

Internal Terminal Handling
++++++++++++++++++++++++++

Internal terminals are not tracked directly by the circuit. One of the
advantages of this is that a device can process parameters
independently of the containing circuit (a reference to the circuit is
no longer needed in ``process_params()``). Another advantage is that
the terminal name is just the internal connection number and does not
need to be unique.

By default the last external terminal in a device is taken as the
local reference. Internal voltages are always referred to that local
reference and the corresponding variable unit is taken from the
internal terminal. 

The ``Circuit`` class has now a function to retrieve all internal
terminals, which (as explained above) are not present in the
``termDict`` dictionary.


Internal Terminal Indexing
++++++++++++++++++++++++++

Internal terminals could be internally indexed by its position in the
terminal list, but for devices that are derived from a base device
(such as autothermal devices) the internal terminal indexes change if
the number of external terminals change. This problem can be avoided
by indexing internal terminals separately and refer to internal
terminals as ``(self.numTerms + number)`` instead of using fixed
numbers. A more convenient solution currently implemented is to have
the ``add_internal_terminal()`` and ``add_reference_term()`` functions
return the internal terminal index.

This solution however would not work for a 'hardwired'
subcircuit. Hardwired subcircuits would also present other problems
with parameter lists and can easily be replacede with regular
subcircuits. For this reason no support for them is planned. Just for
the record a possible solution is presented here: have a function in
the Element class that returns terminal indexes. External terminal
indexes would also change in this case::

     class Element:
     	   ...
	   def term_idx(n):
	       return self.baseIndex + n

	   def int_term_idx(n):
	       return self.baseIndex + self.numTerms + n

``baseIndex`` is set to the starting terminal number for a particular
element. 


Separation of current and charge return vectors in eval_cqs()
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The return of ``eval_cqs()`` for device models is normally a tuple of
the form ``(iVec, qVec)``, where ``iVec`` is a vector of currents and
``qVec`` a vector of charges. Any of the vectors may be empty in cases
where there are no currents or charges to calculate. For this reason
sometimes this approach introduces some overhead compared to having
both currents and charges combined in a single vector. However the
current output format is easier to handle at a higher level and thus
it is preferred.

Still, the ``eval()`` and ``eval_and_deriv()`` functions lump currents
and charges in a single vector because this is a lot more efficient
when implemented using AD tapes (at least with the current
interface). However, the analysis code could bypass those completely
and generate custom AD tapes for greater efficiency (trying this is
one of the near-term goals).


Model, netlist variables and sensitivities considerations
---------------------------------------------------------

Currently device parameters can be initially set by three different
means, in order of priority:

1. Explicitly specified in the device line

2. Specified in a global model

3. The default value 

The explicitly-specified parameters are kept in a separated dictionary
in each element. Global netlist variables can be used to set the
values of parameters specified in the device line or the model
line. 

1. If a string is specified as the value of a numeric parameter value,
   then it is marked as a potential variable.

2. Variables are specified in a ``.vars`` statement in the netlist and
   are assumed to be numeric/vectors

3. When the circuit/analysis is initialized, elements/models/analysis
   check the global netlist variable dictionary to find and set the
   variable value.  If the variable is not found raise an
   exception. One problem with this is that by that time the netlist
   line number is lost and the diagnostic message is not as good.

Another difficulty is how to update dependent parameters when a
variable value is changed. This would require to repeat the whole
process for all models/elements as there is no way to know which ones
are affected.  A change in variables/model/element parameters is
likely to happen in sweeps, sensitivity and optimization calculations.
From the above considerations the current solution requires
re-visiting all elements and re-generating all equations.  One work
around is to create a list of elements to be updated when needed in
the analysis code.


.. include:: ../TODO


