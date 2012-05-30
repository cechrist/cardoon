
Simulator Design Considerations
===============================

Design History
--------------

Work on the project started in September 2011. The design was
originally inspired from the experience with fREEDA [freeda]_ and
carrot [carrot]_ plus some ideas from other simulators and
improvements that take advantage of the flexibility in Python.


Formulation-Independent Circuit Representation
----------------------------------------------

The internal circuit representation and the device library attempt to
be independent of the formulation used for simulation. The ``Circuit``
class has no methods to obtain nodal voltages or calculate the
Jacobian of the circuit. This is delegated to other classes that
handle particular circuit analysis approaches such as nodal,
port-based, linear/nonlinear separation, *etc.*

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
^^^^^^^^^^^^^^^^^^^^^^^^^^

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
with parameter lists. As hardwired subcircuits can be replaced with
regular subcircuits, no support for them is planned. Just for the
record a possible solution is presented here: have a function in the
Element class that returns terminal indexes. External terminal indexes
would also change in this case::

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

Efficiency notes for nodal analysis
-----------------------------------

For efficiency indexing individual elements in arrays from Python
should be avoided as much as possible.  Advanced numpy indexing to
avoid loops for each element of the matrices was tried but
unfortunately ``G[obj] += Jac`` does not work when some of the
elements selected by ``obj`` are repeated (see tentative numpy
tutorial at http://www.scipy.org/NumPy_Tutorial). Possible approaches
to overcome this are the following:

1. Use a few optimized functions doing the inner loops (perhaps using
   cython) or implementing fancy indexing with sparse matrices.

2. Create a giant AD tape for the whole circuit. The nodalAD module
   implements this. This approach is simpler but in practice test
   results (see nodalAD module) show that it does not improve
   efficiency. Also tape generation seems to require a dense matrix
   multiplication (G * x). This rules out the approach for large
   circuits.

Currently the first approach is being used for both dense and sparse
matrices and is described here: nonlinear (and frequency-defined)
elements are added new attribute vectors ``nD_?pos`` and ``nD_?neg``
that contain the non-reference terminal RC numbers where the internal
current sources are connected. This requires some pre-processing but
in this way we can avoid ``if`` statements in inner loop
functions. For regular linear transconductances this does not seem
necessary as we only have to fill the matrix once.

A profile transient analysis of the soliton line with a matrix size of
3022 (using pysparse) seems to indicate that about half of the time is
spent building and half factoring the matrix. At this time the main
Jacobian is (almost) created and factored from scratch at every
iteration.  The results shown in the table were obtained in a netbook
with an Intel Atom processor. Note the time to evaluate nonlinear
models (``eval_and_deriv``) is only about 25% of the time to build the
matrix.

   ======  =======  =======  =======  ======= ===========================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   ======  =======  =======  =======  ======= ===========================
        1    0.000    0.000   23.301   23.301 <string>:1(<module>)
        1    0.000    0.000   23.301   23.301 profile:0(run_analyses(analysisQueue))
        1    0.000    0.000   23.301   23.301 cardoon:26(run_analyses)
        1    0.036    0.036   23.301   23.301 tran.py:66(run)
      100    0.000    0.000   19.181    0.192 fsolve.py:23(solve)
      100    0.012    0.000   19.181    0.192 nodalSP.py:257(solve_simple)
      100    0.132    0.001   19.169    0.192 fsolve.py:59(fsolve_Newton)
      257    0.056    0.000   18.197    0.071 nodalSP.py:260(get_deltax)
      253    3.516    0.014    8.993    0.036 nodalSP.py:729(get_i_Jac)
      257    0.016    0.000    8.921    0.035 nodalSP.py:246(_get_deltax)
      257    8.453    0.033    8.453    0.033 :0(factorize)
    12126    0.924    0.000    2.096    0.000 cppaddev.py:143(eval_and_deriv)
   ======  =======  =======  =======  ======= ===========================

Of course things may be different with other circuits but for this
case it may be expected to speed the simulator at least by a factor of
two if matrix creation and factorization are optimized.

Currently voltages are stored in terminals only after the final
solution is found. The main reason for this is efficiency as it is
less work to operate directly from the vector of unknowns in the
equation-solving routine.

Move sparse matrix implementation to scipy
++++++++++++++++++++++++++++++++++++++++++

The scipy library offers more flexibility than pysparse. It can handle
complex matrices and perform separated symbolic and numeric
factorizations. 


Temperature handling in electrothermal simulations
--------------------------------------------------

The nodal voltage in the thermal port of electrothermal models
represents the difference between the device temperature and the
ambient temperature. In this way a zero difference is usually a good
guess for the nonlinear solver and the numerical solution is more
robust.

Some issues that should be addressed in the future include:

  * Automatically convert the temperature difference to the actual
    temperature for plotting and saving (the same scheme should also
    work for currents in voltage sources for example).

  * Propagate units across terminals in thermal circuits using an
    algorithm similar to freeda's local reference group checking.


.. include:: ../TODO


