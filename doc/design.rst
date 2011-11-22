
Simulator Design Considerations
===============================


Internal Terminal Handling
--------------------------

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
-------------------------------------------------------------

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


Noise/Frequency-Defined Functions
---------------------------------

These functions should work when given for both scalar and vector
frequencies. They should take advantage of the vectorization
facilities in numpy. 

At this time the most convenient type for returning Y matrices for
multiple frequencies seem to be dense 3-D matrices. Sparse formats may
be good for some models but are inefficient for one of the most common
cases: a N-port device with parameters read from a file.

This interface will be reviewed when the AC analysis is implemented.


.. include:: ../TODO


