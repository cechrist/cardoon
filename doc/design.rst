
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
advantages that a device can process parameters independently of the
containing circuit (a reference to the circuit is no longer needed in
``process_params()``). Another advantage is that the terminal name is
just the internal connection number and does not need to be unique.

By default the last external terminal in a device is taken as the
local reference. Internal voltages are always referred to that local
reference and the corresponding variable unit is taken from the
internal terminal. 

The ``Circuit`` class has now a function to retrieve all internal
terminals, which (as explained above) are not present in the
``termDict`` attribute.



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




.. include:: ../TODO


