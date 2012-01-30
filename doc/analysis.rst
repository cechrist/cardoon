
Nodal Analyses Reference
========================

In the following discussion we assume that :math:`x` is the vector of
nodal variables. We denote matrices and frequency-domain quantities
with uppercase letters. 

This section documents the formulation used for nodal-based
analyses. For all analyses there are four types of devices to
consider:

  1. Linear VCCS/QS: contribute the :math:`G` and :math:`C` matrices,
     respectively. Elements provide the following attributes to build
     these matrices: ``linearVCCS`` and ``linearVCQS``.
  
  2. Nonlinear VCCS/QS: contribute :math:`i(x)` and :math:`q(x)`
     vector functions. These and their Jacobian (:math:`di/dx` and
     :math:`dq/dx`) are returned by ``eval_and_deriv()`` in each
     element.
  
  3. Frequency-defined devices: contribute the complex,
     frequency-dependent :math:`Y(f)` matrix. The corresponding
     impulse-response matrix is denoted :math:`Y(t)`. The following
     methods are used to build this matrix: ``get_Y_matrix()`` for
     :math:`f \neq 0` and ``get_G_matrix()`` for :math:`f = 0`.
  
  4. Sources: contribute a source vector. There are 3 types of source:
     DC (:math:`s_{DC}`), time-domain (:math:`s(t)`) and
     frequency-domain (:math:`S(f)`). The following element methods
     are used: ``get_DCsource()``, ``get_TDsource()``,
     ``get_FDsource()`` and ``get_AC()``.


OP/DC Equations
---------------

For DC analysis, :math:`C` and :math:`q(x)` are ignored and only the
DC component of sources (:math:`s_{DC}`) and the DC conductance matrix
of frequency-defined elements (:math:`G_0`) are considered.

The analysis solves the following nonlinear equation iteratively
using Newton's method:

.. math::

    (G + G_0) x + i(x) = s_{DC}

Let :math:`G_1 = G + G_0`. The iteration is defined by
linearizing :math:`i(x)` as follows:

.. math::

    G_1 x_{k+1} + i_k + \frac{di_k}{dx} \, (x_{k+1} - x_k) = s_{DC} \; ,

where the :math:`k` suffix indicates the iteration number. :math:`x_k`
is assumed to be known and :math:`x_{k+1}` is the unknown to solve
for. The initial guess (:math:`x_0`) is set to the values suggested by
the nonlinear devices, if any, or otherwise to zero. The previous
equation can be re-arranged as as the following system of linear
equations:

.. math::

     (G_1 + \frac{di_k}{dx}) \, x_{k+1} = 
            s_{DC} - i_k + \frac{di_k}{dx} \, x_k \; ,

This equation can be seen as the nodal equation of a circuit obtained
by substituting the nonlinear devices by current sources and
transcunductances that are dependent of the current approximation for
the nodal voltages (:math:`x_k`). Iterations stop when

.. math::

   | x^{k+1}_n - x^k_n | < \epsilon


AC Formulation
--------------

Here :math:`X(f)` is the nodal vector in frequency domain. A DC
operating point analysis is performed first to obtain the incremental
conductances and capacitances given by :math:`di/dx` and
:math:`dq/dx`, respectively.  The analysis solves for :math:`X(f)` for
all requested frequencies using the following equation:

.. math::

    \left[ (G + \frac{di}{dx}) + j 2 \pi f (C + \frac{dq}{dx}) 
           + Y(f) \right] \, X(f) = S(f)



Transient Analysis Equations
----------------------------

Transient analysis solves the following nonlinear
algebraic-integral-differential equation:

.. math::

    G x + C \dot{x} + i(x) + \dot{q}(x) + 
      \int_{0}^\infty Y(\tau) x(t - \tau) d\tau
      = s_{DC} + s(t)  \; ,

given an initial condition, :math:`x(0) = x_0`, usually obtained by
calculating the operating point of the circuit using the OP
analysis. The dotted quantities indicate derivative with respect to
time. 

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

   i'(x_n) = i(x_n) + a_0 q(x_n)

   s' = s_n - \sum_{m=1}^\infty \textbf{Y}_m x_{n-m} 
             + f_{n-1}(C x + q) \; ,

where :math:`Y_m = Y(t_m)` and :math:`Y_0 = Y(0)`. Note that
:math:`s'` is a known vector at the :math:`n^{th}` time step. This is
the equation of a DC circuit with a conductance matrix equal to
:math:`G'`, a set of nonlinear currents given by the :math:`i'(x_n)`
function and a source vector given by :math:`s'`. The unknown
(:math:`x_n`) is iteratively solved using Newton's Method (similarly
as in OP/DC analysis). Iterations are defined by linearizing
:math:`i'(x)` as follows:

.. math::

    G' x^{k+1}_n + i'(x^k_n) + \frac{di'}{dx} (x^{k+1}_n - x^k_n)
        = s' \; ,

where the :math:`k` subscript denotes the Newton iteration number.
This equation is re-arranged as follows:

.. math::

    \left( G' + \frac{di'}{dx} \right) x^{k+1}_n =
      s' - i'(x^k_n) + \frac{di'}{dx} x^k_n \; ,

as the right-hand side of this equation is known at the :math:`k^{th}`
iteration, :math:`x^{k+1}_n` can be found by solving a linear system
of equations. Iterations stop when

.. math::

   | x^{k+1}_n - x^k_n | < \epsilon

