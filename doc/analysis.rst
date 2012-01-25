
Nodal Analyses Reference
========================

In the following discussion we assume that :math:`v` is the vector of
nodal variables. We denote matrices and frequency-domain quantities
with uppercase letters. 

This section documents the formulation used for nodal-based
analyses. For all analyses there are 4 types of devices to consider:

  1. Linear VCCS/QS: contribute the :math:`G` and :math:`C` matrices,
     respectively.
  
  2. Nonlinear VCCS/QS: contribute :math:`i(v)` and :math:`q(v)`
     vector functions. These and their Jacobian (:math:`di/dv` and
     :math:`dq/dv`) are returned by ``eval_and_deriv()``.
  
  3. Frequency-defined devices: contribute the complex,
     frequency-dependent :math:`Y(f)` matrix. The corresponding
     impulse-response matrix is denoted :math:`Y(t)`
  
  4. Sources: contribute a source vector. There are 3 types of source:
     DC (:math:`s_{DC}`), time-domain (:math:`s(t)`) and
     frequency-domain (:math:`S(f)`).


OP/DC Equations
---------------

For DC analysis, :math:`C` and :math:`q(v)` are ignored and only the
DC component of sources (:math:`s_{DC}`) are considered.

The analysis solves the following nonlinear equation iteratively
using Newton's method:

.. math::

    G v + i(v) = s_{DC}

The iteration is defined by linearizing :math:`i(v)` as follows:

.. math::

    G v_{k+1} + i_k + \frac{di_k}{dv} \, (v_{k+1} - v_k) = s_{DC} \; ,

where the :math:`k` suffix indicates the iteration number. :math:`v_k`
is assumed to be known and :math:`v_{k+1}` is the unknown to solve
for. The initial guess (:math:`v_0`) is set to the values suggested by
the nonlinear devices, if any, or otherwise to zero. The previous
equation can be re-arranged as as the following system of linear
equations:

.. math::

     (G + \frac{di_k}{dv}) \, v_{k+1} = 
            s_{DC} - i_k + \frac{di_k}{dv} \, v_k \; ,

This equation can be seen as the nodal equation of a circuit obtained
by substituting the nonlinear devices by current sources and
transcunductances that are dependent of the current approximation for
the nodal voltages (:math:`v_k`). Iterations stop when

.. math::

   | v^{k+1}_n - v^k_n | < \epsilon


AC Formulation
--------------

Here :math:`V(f)` is the nodal vector in frequency domain. The
incremental conductances and capacitances are given by
:math:`\frac{di}{dv}` and :math:`\frac{dq}{dv}`, respectively.  The
analysis solves for :math:`V(f)` for all requested frequencies using
the following equation:

.. math::

    \left[ (G + \frac{di}{dv}) + j 2 \pi f (C + \frac{dq}{dv}) 
           + Y(f) \right] \, V(f) = S(f)



Planned Transient Analysis Equations
------------------------------------

Transient analysis solves the following nonlinear
algebraic-integral-differential equation:

.. math::

    G v + C \dot{v} + i(v) + \dot{q}(v) + 
      \int_{0}^\infty Y(\tau) v(t - \tau) d\tau
      = s_{DC} + s(t)  \; ,

given an initial condition, :math:`v(0) = v_0`, usually obtained by
calculating the operating point of the circuit using the OP
analysis. The dotted quantities indicate derivative with respect to
time.  An integration method (such as Backward Euler (BE) or
Trapezoidal Integration) is applied to transform the differential
equation into a difference equation by discretizing time and
approximating derivatives with respect to time. For example, using the
BE rule:

.. math::

    \dot{v}(t_n) = \dot{v}_n \approx \frac{v_n - v_{n-1}}{h} \; ,

here, the subscript :math:`n` denotes the time sample number and
:math:`h` is the time step size. For implicit methods in general,

.. math::

    \dot{v_n} \approx a_0 v_n + f_{n-1}(v) \; ,

with :math:`f_{n-1}(v)` being a function that depends on the previous
samples of :math:`v`:

.. math::

    f_{n-1}(v) = a_1 v_{n-1} + a_2 v_{n-2} + \dots \; .

The :math:`a_i; i=0,1,\dots` coefficients depend on the time step size
and the integration method. Substituting dotted variables and
discretizing the convolution operation the resulting circuit equation
is the following:

.. math::

    G' v_n + i'(v_n) = s' \; ,

with

.. math::

   G' = G + Y_0 + a_0 C

   i'(v_n) = i(v_n) + a_0 q(v_n)

   s' = s_n - \sum_{m=1}^\infty \textbf{Y}_m v_{n-m} 
             - C f_{n-1}(v) + f_{n-1}(q) \; ,

where :math:`Y_m = Y(t_m)` and :math:`Y_0 = Y(0)`. Note :math:`s'` is
a known vector at the :math:`n^{th}` time step. This is the equation
of a DC circuit with a conductance matrix equal to :math:`G'`, a set
of nonlinear currents given by the :math:`i'(v_n)` function and a
source vector given by :math:`s'`. The unknown (:math:`v_n`) is
iteratively solved using Newton's Method (similarly as in OP/DC
analysis). Iterations are defined by linearizing :math:`i'(v)` as
follows:

.. math::

    G' v^{k+1}_n + i'(v^k_n) + \frac{di'}{dv} (v^{k+1}_n - v^k_n)
        = s' \; ,

where the :math:`k` subscript denotes the Newton iteration number.
This equation is re-arranged as follows:

.. math::

    \left( G' + \frac{di'}{dv} \right) v^{k+1}_n =
      s' - i'(v^k_n) + \frac{di'}{dv} v^k_n \; ,

as the right-hand side of this equation is known at the :math:`k^{th}`
iteration, :math:`v^{k+1}_n` can be found by solving a linear system
of equations. Iterations stop when

.. math::

   | v^{k+1}_n - v^k_n | < \epsilon

