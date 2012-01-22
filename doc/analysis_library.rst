========================
Analysis Library Catalog
========================
 
ac: AC Sweep
------------

Calculates a AC sweep of a circuit using the nodal approach. After
the analysis is complete, nodal voltages are saved in circuit and
terminals with the ``aC_`` prefix.  After this the analysis drops
to an interactive shell if the ``shell`` global variable is set to
``True``.

An OP analysis is performed first to obtain the Jacobian from
nonlinear devices. Convergence parameters for the Newton method
are controlled using the global variables in ``.options``.

One plot window is generated for each ``.plot`` statement. Use the
following request types for this analysis: ``ac_mag``,
``ac_phase`` or ``ac_dB``.

Example::

    .analysis ac start=100. stop=1MEG num=100 log=True
    .plot ac_dB 153 151 23
    .plot ac_phase 23

Formulation
+++++++++++

In the following discussion we assume that :math:`v` is the vector of
nodal variables in time domain and :math:`V(f)` is the same vector in
frequency domain.  There are 4 types of devices to consider for AC
analysis (so far, must add later time-delayed nonlinear CS):

  1. Linear VCCS/QS: considered in the real, frequency-independent
  :math:`G` and :math:`C` matrices, respectively.
  
  2. Nonlinear VCCS/QS: an OP analysis is performed to obtain the
  operating point. The Jacobian returned by ``eval_and_deriv()`` is used
  to generate the real, frequency-independent :math:`dI/dv` and
  :math:`dQ/dv` matrices, respectively.
  
  3. Frequency-defined devices: contribute the complex,
  frequency-dependent :math:`Y(f)` matrix.
  
  4. Sources: contribute a complex frequency-dependent vector,
  :math:`S(f)`. 

The analysis solves for :math:`V(f)` for all requested frequencies
using the following equation:

.. math::

    \left[ (G + \frac{dI}{dv}) + 
         j 2 \pi f (C + \frac{dQ}{dv}) + Y(f) \right] \, V(f) = S(f)



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 log          0                         Use logarithmic scale                                
 num          50                        Number of points in sweep                            
 shell        0                         Drop to ipython shell after calculation              
 start        1.0          Hz           Frequency sweep start value                          
 stop         10.0         Hz           Frequency sweep stop value                           
 =========== ============ ============ ===================================================== 

dc: DC Sweep
------------

Calculates a DC sweep of a circuit using the nodal approach. After
the analysis is complete, nodal voltages are saved in circuit and
terminals with the ``dC_`` prefix.  After this the analysis drops
to an interactive shell if the ``shell`` global variable is set to
``True``.

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``.

One plot window is generated for each ``.plot`` statement. Use
``dc`` request type for this analysis.

Example::

    .analysis dc device=vsin:v1 param=vdc start=-2. stop=2. num=50 

    # Some options that affect convergence properties
    .options maxiter=300 gyr=1e-5 maxdelta=5.
    
    .plot dc 153 151 23

Formulation
+++++++++++

The formulation is the same as for the OP analysis.



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 device                                 Instance name of device to sweep variable            
 num          50                        Number of points in sweep                            
 param                                  Device parameter to sweep                            
 shell        0                         Drop to ipython shell after calculation              
 start        0.0          V            Sweep start value                                    
 stop         0.0          V            Sweep stop value                                     
 verbose      0                         Show iterations for each point                       
 =========== ============ ============ ===================================================== 

op: DC Operating Point
----------------------

Calculates the DC operating point of a circuit using the nodal
approach. After the analysis is complete, nodal voltages are saved
in circuit and terminals with the ``nD_`` prefix.  After this the
analysis drops to an interactive shell if the ``shell`` global
variable is set to ``True``.

By default the voltage at all external voltages is printed after
the analysis is complete. Optionally the operating points of
nonlinear elements can be printed. 

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``.

Example::

    .analysis op intvars=1 shell=1

OP/DC formulation
+++++++++++++++++

In the following discussion we assume that :math:`v` is the vector
of nodal variables. There are 4 types of devices to consider for
DC operating point:

  1. Linear VCCS: considered in the :math:`G` matrix.
  
  2. Nonlinear VCCS: considered in the :math:`i(v)` vector. This and
     its Jacobian (:math:`di/dv`) are returned by
     ``eval_and_deriv()``.
  
  3. Frequency-defined devices: their DC contribution is added to
     the :math:`G` matrix.
  
  4. Sources: contribute the source vector, :math:`s`. 

The analysis solves the following nonlinear equation iteratively
using Newton's method:

.. math::

    G v + i(v) - s = 0

The iteration is defined by linearizing :math:`i(v)` as follows:

.. math::

    G v_{k+1} + i_k + \frac{di_k}{dv} \, (v_{k+1} - v_k) - s = 0 \; ,

where the :math:`k` suffix indicates the iteration
number. :math:`v_k` is assumed to be known and :math:`v_{k+1}` is
the unknown to solve for. The initial guess (:math:`v_0`) is set
to the values suggested by the nonlinear devices, if any, or
otherwise to zero. The previous equation can be re-arranged as as
the following system of linear equations:

.. math::

     (G + \frac{di_k}{dv}) \, v_{k+1} = 
            s - i_k + \frac{di_k}{dv} \, v_k \; ,

This equation can be seen as the nodal equation of a circuit
obtained by substituting the nonlinear devices by current sources
and transcunductances that are dependent of the current
approximation for the nodal voltages (:math:`v_k`).



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 elemop       0                         Print element operating points                       
 fullAD       0                         Use CPPAD for entire nonlinear part                  
 intvars      0                         Print internal element nodal variables               
 shell        0                         Drop to ipython shell after calculation              
 =========== ============ ============ ===================================================== 

testdev: Test Equations Of a Nonlinear Device
---------------------------------------------

One advantage of using this method over a DC sweep is that no
Newton iterations are needed. The following internal functions are
tested here:

* process_params()
* set_temp_vars()
* eval_cqs()
* eval()
* get_OP()
* power() (for electrothermal models)

After completion the analysis drops to an interactive shell if the
``shell`` global variable is set to ``True``

Example::

    .analysis testdev plot=1 ports_bias = [3V, 3.V, 0V] sweep_port=1 \ 
    	  start = 0V stop= 3V sweep_num=1000 device = mosekv:m1 \ 
    	  param = temp param_val = [-10, 27, 50]



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 device                                 Instance name of device to test                      
 param                                  Parameter for outer sweep                            
 param_val    []                        Vector with parameter values to sweep                
 plot         1                         Auto-plot currents and charges                       
 ports_bias   []           V            Vector with default values of port voltages          
 shell        0                         Drop to ipython shell after calculation              
 start        0.0          V            Sweep start value                                    
 stop         0.0          V            Sweep stop value                                     
 sweep_num    0                         Number of points in sweep                            
 sweep_port   0                         Port number to be swept, starting from zero          
 useAD        1                         Use automatic differentiation                        
 =========== ============ ============ ===================================================== 

