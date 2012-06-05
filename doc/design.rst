
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

Currently the first approach is being used for dense matrices and is
described here: nonlinear (and frequency-defined) elements are added
new attribute vectors ``nD_?pos`` and ``nD_?neg`` that contain the
non-reference terminal RC numbers where the internal current sources
are connected. This requires some pre-processing but in this way we
can avoid ``if`` statements in inner loop functions. For regular
linear transconductances this does not seem necessary as we only have
to fill the matrix once.

For the sparse-matrix nodal implementation, further pre-processing is
needed to create the index vectors to directly fill the main
matrix. The main matrix is built in triplet format (Scipy coo_matrix
format). Values from nonlinear elements are directly inserted in the
coo_matrix ``data`` array::

    M.data[mpidx + self._mbase] = Jac.flat[jacpidx]
    self._mbase += len(mpidx)
    M.data[mnidx + self._mbase] = -Jac.flat[jacnidx]
    self._mbase += len(mnidx)
 
It seems that this method is the fastest that can be achieved without
resorting to compilation.  According to my tests fancy indexing is a
little faster (very little) than using ``nimpy.put()`` and
``numpy.take()`` to fill the matrix. This does not agree with the
comments in http://wesmckinney.com/blog/?p=215 . Perhaps this is
different for other architectures or versions of the program?

Further optimization would possibly require access to the low-level
interface to the C libraries plus compilation of some functions used
in inner loops. 


Profiler results for soliton.net
++++++++++++++++++++++++++++++++

A profile transient analysis of the soliton line with a matrix size of
3022 was performed using scipy matrices, optimized matrix filling and
SuperLU for matrix decomposition. The number of time steps is 501. The
results shown in the table were obtained in a netbook with an Intel
Atom processor:

   ======  =======  =======  =======  ======= ===========================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   ======  =======  =======  =======  ======= ===========================
        1    0.000    0.000   59.380   59.380 cardoon:27(run_analyses)
        1    0.152    0.152   59.380   59.380 tran.py:73(run)
      501    0.040    0.000   44.099    0.088 fsolve.py:23(solve)
      501    0.016    0.000   44.059    0.088 nodal.py:420(solve_simple)
      501    0.260    0.001   44.043    0.088 fsolve.py:59(fsolve_Newton)
      505    0.068    0.000   40.807    0.081 nodal.py:423(get_deltax)
      501    7.340    0.015   25.678    0.051 nodalSP.py:731(get_i_Jac)
      505    0.132    0.000   14.793    0.029 nodalSP.py:204(_get_deltax)
    47282    7.824    0.000    9.421    0.000 nodalSP.py:165(_set_Jac)
      505    0.036    0.000    9.129    0.018 linsolve.py:131(factorized)
      505    0.096    0.000    9.093    0.018 linsolve.py:103(splu)
      505    8.597    0.017    8.597    0.017 :0(dgstrf)
      501    2.756    0.006    7.720    0.015 nodalSP.py:651(update_q)
     1004    4.904    0.005    4.904    0.005 :0(solve)
    23782    1.996    0.000    4.236    0.000 cppaddev.py:143(eval_and_deriv)
   230703    3.812    0.000    3.812    0.000 :0(len)
      505    0.136    0.000    3.084    0.006 coo.py:238(tocsc)
     1020    2.680    0.003    2.680    0.003 :0(max)
      499    0.036    0.000    2.492    0.005 nodalSP.py:223(get_chord_deltax)
     8013    0.408    0.000    2.476    0.000 nodalSP.py:43(set_quad)
    70829    2.388    0.000    2.388    0.000 nodal.py:52(set_i)
    23782    1.604    0.000    2.100    0.000 adfun.py:208(jacobian)
    18224    1.252    0.000    2.076    0.000 nodalSP.py:34(triplet_append)
   ======  =======  =======  =======  ======= ===========================

There seem to be no obvious major bottlenecks. Most of the simulation
time is spent building the matrix and evaluating nonlinear devices
(``get_i_Jac``), followed by linear system solving
(``_get_deltax``). Note than nonlinear device evaluation
``eval_and_deriv`` takes a small percentage of the total simulation
time, in part because there are few nonlinear devices in this circuit.

Profiler results for sum741_profile.net
+++++++++++++++++++++++++++++++++++++++

A profile transient analysis of the soliton line with a matrix size of
189 was performed using scipy matrices, optimized matrix filling and
SuperLU for matrix decomposition. The number of time steps is 101. The
results shown in the table were obtained in a netbook with an Intel
Atom processor:

   ======  =======  =======  =======  ======= ===========================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   ======  =======  =======  =======  ======= ===========================
        1    0.000    0.000   13.149   13.149 cardoon:27(run_analyses)
        1    0.020    0.020   13.149   13.149 tran.py:73(run)
      101    0.012    0.000   11.137    0.110 fsolve.py:23(solve)
      101    0.012    0.000   11.125    0.110 nodal.py:420(solve_simple)
      101    0.084    0.001   11.113    0.110 fsolve.py:59(fsolve_Newton)
      257    0.064    0.000   10.849    0.042 nodal.py:423(get_deltax)
      203    1.844    0.009    6.748    0.033 nodalSP.py:731(get_i_Jac)
      257    0.032    0.000    2.716    0.011 nodalSP.py:204(_get_deltax)
    11960    2.076    0.000    2.424    0.000 nodalSP.py:165(_set_Jac)
     6708    0.604    0.000    1.660    0.000 cppaddev.py:143(eval_and_deriv)
       54    0.192    0.004    1.320    0.024 nodalSP.py:438(get_i_Jac)
      257    0.088    0.000    1.320    0.005 coo.py:238(tocsc)
      257    0.008    0.000    1.284    0.005 linsolve.py:131(factorized)
      257    0.040    0.000    1.276    0.005 linsolve.py:103(splu)
      257    1.008    0.004    1.008    0.004 :0(dgstrf)
      101    0.284    0.003    0.992    0.010 nodalSP.py:651(update_q)
     6708    0.764    0.000    0.888    0.000 adfun.py:208(jacobian)
    54474    0.800    0.000    0.800    0.000 :0(len)
    14586    0.768    0.000    0.768    0.000 nodal.py:52(set_i)
      358    0.088    0.000    0.740    0.002 base.py:270(__mul__)
      260    0.112    0.000    0.528    0.002 compressed.py:20(__init__)
      358    0.064    0.000    0.384    0.001 compressed.py:259(_mul_vector)
     9334    0.376    0.000    0.376    0.000 nodal.py:41(set_xin)
   ======  =======  =======  =======  ======= ===========================

Similar observations can be made for this case, except that in this
circuit most of the devices are nonlinear.


Profiler results using pysparse library (old implementation)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Note: the code used in this profile is not the default currently used
in the program.  A profile transient analysis of the soliton line with
a matrix size of 3022 (using pysparse) seems to indicate that about
half of the time is spent building and half factoring the matrix. At
this time the main Jacobian is (almost) created and factored from
scratch at every iteration.  The results shown in the table were
obtained in a netbook with an Intel Atom processor. Note the time to
evaluate nonlinear models (``eval_and_deriv``) is only about 25% of
the time to build the matrix.

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


