
Nodal Analyses Reference
========================

Introduction
------------

In the following discussion we assume that :math:`x` is the vector of
nodal variables. We denote matrices and frequency-domain quantities
with uppercase letters. 

This section documents the formulation used for nodal-based
analyses. For all analyses there are four types of devices to
consider:

Linear VCCS/QS

    Contribute the :math:`G` and :math:`C` matrices,
    respectively. Elements provide the following attributes to build
    these matrices: ``linearVCCS`` and ``linearVCQS``.
  
Nonlinear VCCS/QS

    Contribute :math:`i(v, v_T)` and :math:`q(v, v_T)` vector
    functions. Here :math:`v` is a vector with controlling port
    voltages and :math:`v_T` is a vector with time-delayed control
    voltages. For the :math:`j^{th}` nonlinear device, vector
    :math:`[i_j , q_j]^T` and its Jacobian are returned by
    ``eval_and_deriv()``. The Jacobian is returned in a single matrix
    with the following format:

    .. math::
 
        \left[ \begin{array}{cc} di_j/dv & di_j/dv_T \\
     	                         dq_j/dv & dq_j/dv_T 
				 \end{array} \right]

    In nodal analyses, controlling port voltages (:math:`v`) are
    obtained as differences between elements of the nodal variables
    vector (:math:`x`). In compact form this can be denoted as:

    .. math::
 
	v = T_{cv} x \; ,

    where :math:`T_{cv}` is an incidence matrix. All nonzero entries
    in :math:`T_{dv}` are either equal to 1 or -1 and there are at
    most two nonzero entries per row. This matrix is implicitly stored
    in the simulator with the control port definitions in each
    element. Samples of time-delayed control voltages can be
    calculated in a similar way using a :math:`T_{dv}` incident matrix
    but the effect of time delays depends on the analysis type.

    Vector functions :math:`i(v, v_T)` and :math:`q(v, v_T)` are
    mapped into the nodal equations using incidence matrices
    :math:`T_i` and :math:`T_q`, which have the same properties as
    the transpose of :math:`T_{cv}`.
  
Frequency-defined devices

    Contribute the complex, frequency-dependent :math:`Y(f)`
    matrix. The corresponding impulse-response matrix is denoted
    :math:`Y(t)`. The following methods are used to build this
    matrix: ``get_Y_matrix()`` for :math:`f \neq 0` and
    ``get_G_matrix()`` for :math:`f = 0`. More details will be
    provided when time-domain support for this kind of device is
    implemented.
  
Sources

    Contribute a source vector. There are 3 types of source: DC
    (:math:`s_{DC}`), time-domain (:math:`s(t)`) and frequency-domain
    (:math:`S(f)`). The following element methods are used:
    ``get_DCsource()``, ``get_TDsource()``, ``get_FDsource()`` and
    ``get_AC()``.


OP/DC Equations
---------------

For DC analysis, time-delayed control voltages are calculated as
regular control voltages; :math:`C` and :math:`q(v, v_T)` are ignored
and only the DC component of sources (:math:`s_{DC}`) and the DC
conductance matrix of frequency-defined elements (:math:`G_0`) are
considered.

The analysis solves the following nonlinear equation iteratively
using Newton's method:

.. math::

    (G + G_0) x + T_i i(T_{cv} x, T_{dv} x) = s_{DC} \; .

Let :math:`G_1 = G + G_0`. The iteration is defined by linearizing
:math:`i(T_{cv} x, T_{dv} x)` as follows:

.. math::

    G_1 x^{k+1} + T_i i^k + J^k_i \, (x^{k+1} - x^k) = s_{DC} \; ,

where the :math:`k` superscript indicates the iteration number and 

.. math::

     J^k_i = T_i \left(
                 \frac{di^k}{dv} T_{cv} + \frac{di^k}{dv_T} T_{dv} 
                 \right) \; .

The :math:`x^k` vector is assumed to be known for each iteration and
:math:`x^{k+1}` is the unknown to solve for. The initial guess
(:math:`x_0`) is set to the values suggested by the nonlinear devices,
if any, or otherwise to zero. The previous equation can be re-arranged
as as the following system of linear equations:

.. math::

     (G_1 + J^k_i) \, \Delta x^{k+1} = 
            s_{DC} - G_1 x^k - T_i i^k \; ,

with :math:`\Delta x^{k+1} = x^{k+1} - x^k`.  This equation can be
seen as the nodal equation of a linearized circuit obtained by
substituting all devices by transconductances in parallel with current
sources that are dependent of the current approximation for the nodal
voltages (:math:`x^k`). 

Iterations stop when

.. math::

   | \Delta x^{k+1}_j | < 
       \mbox{reltol} \; \max(|x^{k+1}_j|,  |x^k_j|) + \mbox{abstol}

for each nodal variable (:math:`j`). In addition, if ``errfunc`` is
set to ``True``, the following check is made:

.. math::

    |(G + G_0) x + T_i i(T_{cv} x, T_{dv} x) - s_{DC}| < \mbox{abstol}


AC Formulation
--------------

Here :math:`X(f)` is the nodal vector in frequency domain. A DC
operating point analysis is performed first to obtain the incremental
conductances and capacitances given by :math:`di/dx` and
:math:`dq/dx`, respectively.  The analysis solves for :math:`X(f)` for
all requested frequencies using the following equation:

.. math::

    \left[ (G + J_i + j 2 \pi f (C + J_q) 
           + Y(f) \right] \, X(f) = S(f)

The :math:`J_i` and :math:`J_q` Jacobian matrices are calculated as
follows:

.. math::

    J_i = T_i \left(
           \frac{di}{dv} T_{cv} + \frac{di}{dv_T} E_{\tau} T_{dv} 
    	   \right) \; ,

    J_q = T_q \left(
           \frac{dq}{dv} T_{cv} + \frac{dq}{dv_T} E_{\tau} T_{dv} 
    	   \right) \; ,

where :math:`E_{\tau}` is a diagonal matrix that includes the effect
of time delays in each control port:

.. math::

     E_{\tau} = \mbox{diag}\left( [\exp(-j \omega \tau_1), 
     	      	\exp(-j \omega \tau_2), \dots , \exp(-j \omega \tau_m)]
		\right)
     	      

Transient Analysis Equations
----------------------------

Transient analysis solves the following nonlinear
algebraic-integral-differential equation given an initial condition,
:math:`x(0) = x_0`:

.. math::

    G x + C \dot{x} + T_i i(T_{cv} x, d(T_{dv} x)) + 
      T_q \dot{q}(T_{cv} x, d(T_{dv} x)) + 
      \int_{0}^\infty Y(\tau) x(t - \tau) d\tau
      = s_{DC} + s(t)  \; ,

where dotted quantities indicate derivative with respect to time and
:math:`d()` is a vector function that applies a (possibly different)
time delay to each control voltage.  The delay function (:math:`d()`)
is implemented by storing all time-delayed control port voltages and
using an interpolation function to find the voltage at the desired
time in the past. 

The initial condition (:math:`x(0)`) is usually obtained by
calculating the operating point of the circuit using the OP analysis.

An integration method (such as Backward Euler (BE) or Trapezoidal
Integration) is applied to transform the differential equation into a
difference equation by discretizing time and approximating derivatives
with respect to time. Here we assume the time step (:math:`h`) is
constant.  For example, using the BE rule:

.. math::

    \dot{q}(t_n) = \dot{q}_n \approx \frac{q_n - q_{n-1}}{h} \; ,

here, the subscript :math:`n` denotes the time sample number. For
implicit methods in general,

.. math::

    \dot{q_n} \approx a_0 q_n - f_{n-1}(q) \; ,

with :math:`f_{n-1}(x)` being a function that depends on the previous
samples of :math:`q`:

.. math::

    f_{n-1}(q) = a_1 q_{n-1} + a_2 q_{n-2} + \dots \; .

The :math:`a_i; i=0,1,\dots` coefficients depend on the time step size
and the integration method. Substituting dotted variables and
discretizing the convolution operation the resulting circuit equation
is the following:

.. math::

    G' x_n + i'(x_n) = s' \; ,

with

.. math::

   G' = G + Y_0 + a_0 C

   i'(x_n) = T_i i(T_{cv} x_n, d(T_{dv} x_n)) + 
             a_0 T_q q(T_{cv} x_n, d(T_{dv} x_n))

   s' = s_n - \sum_{m=1}^\infty \textbf{Y}_m x_{n-m} 
             + f_{n-1}(C x + T_q q) \; ,

where :math:`Y_m = Y(t_m)` and :math:`Y_0 = Y(0)`. Note that
:math:`s'` is a known vector at the :math:`n^{th}` time step. This is
the equation of a DC circuit with a conductance matrix equal to
:math:`G'`, a set of nonlinear currents given by the :math:`i'(x_n)`
function and a source vector given by :math:`s'`. The unknown
(:math:`x_n`) is iteratively solved using Newton's Method (similarly
as in OP/DC analysis). Iterations are defined by linearizing
:math:`i'(x)` as follows:

.. math::

    G' x^{k+1}_n + i'(x^k_n) + J^k_n \Delta x^{k+1}_n
        = s' \; ,

where the :math:`k` subscript denotes the Newton iteration number,
:math:`\Delta x^{k+1}_n = x^{k+1}_n - x^k_n` and :math:`J^k_n =
di'(x^k_n)/dx`.  This equation is re-arranged as follows:

.. math::

    ( G' + J^k_n ) \Delta x^{k+1}_n =
      s' - G' x^k_n - i'(x^k_n) \; ,

as the right-hand side of this equation is known at the :math:`k^{th}`
iteration, :math:`x^{k+1}_n` can be found by solving a linear system
of equations. 

The criterion to stop iterations is the same as in the DC analysis.
