
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

Separation of current and charge return vectors
-----------------------------------------------

Currently the return of ``eval_cqs()`` is one vector combining
currents and charges. It would be better if the two could be separated
but we have to decide what to do when one of the two vectors is empty.


Noise/Frequency-Defined Functions
---------------------------------

These functions should work when given for both scalar and vector
frequencies. They should take advantage of the vectorization
facilities in numpy. 

At this time the most convenient type for returning Y matrices for
multiple frequencies seem to be dense 3-D matrices.

.. include:: ../TODO


