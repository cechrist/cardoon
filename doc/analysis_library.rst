:tocdepth: 2

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

AC formulation documented in :doc:`analysis`

Example::

    .analysis ac start=100. stop=1MEG num=100 log=True
    .plot ac_dB 153 151 23
    .plot ac_phase 23



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

The following parameters can be swept: 

  * Any device parameter of ``float`` type (device name must be
    specified in this case)

  * Global temperature (no device specified): sweep temperature of
    all devices that do not explicitly have ``temp`` set.

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``. The type of matrix used in
this analysis is controlled by the ``sparse`` option. Global
options are documented in :doc:`global_vars`. 

One plot window is generated for each ``.plot`` statement. Use
``dc`` request type for this analysis.

DC analysis formulation is documented in :doc:`analysis`, and
internal classes and functions used in this analysis are
documented in :doc:`analyses_classes`.

Examples::

    # Device parameter sweep
    .analysis dc device=vsin:v1 param=vdc start=-2. stop=2. num=50 

    # Global temperature sweep
    .analysis dc param=temp start=-20C stop=80C 

    # Some options that affect convergence properties
    .options maxiter=300 gyr=1e-5 maxdelta=5.
    
    .plot dc 153 151 23



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 device                                 Instance name of device to sweep variable            
 num          50                        Number of points in sweep                            
 param                                  Parameter to sweep                                   
 shell        0                         Drop to ipython shell after calculation              
 start        0.0          (variable)   Sweep start value                                    
 stop         0.0          (variable)   Sweep stop value                                     
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
the global variables in ``.options``. The type of matrix used in
this analysis is controlled by the ``sparse`` option. Global
options are documented in :doc:`global_vars`. 

OP analysis formulation is documented in :doc:`analysis`, and
internal classes and functions used in this analysis are
documented in :doc:`analyses_classes`.

Example::

    .analysis op intvars=1 shell=1



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 elemop       0                         Print element operating points                       
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

tran: Transient Analysis
------------------------

Solves nodal equations starting from ``t=0`` to ``tstop`` with a
fixed time step equal to ``tstep``. Two integration methods are
supported: Backwards Euler (``im = BE``) and trapezoidal
(``im=trap``). Support for frequency-defined elements and time
delays is not yet included.

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``. The type of matrix used in
this analysis is controlled by the ``sparse`` option. Global
options are documented in :doc:`global_vars`. 

One plot window is generated for each ``.plot`` statement. Use
``tran`` request type for this analysis. By default, only results
for nodes listed in ``.plot`` statements are saved. To save all
nodal variables set ``saveall`` to 1.

Transient analysis formulation is documented in :doc:`analysis`,
and internal classes and functions used in this analysis are
documented in :doc:`analyses_classes`.

Example::

    .analysis tran tstop=1ms tstep=.01ms im=BE

    .plot tran vin vout



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 im           BE                        Integration method                                   
 saveall      0                         Save all nodal voltages                              
 shell        0                         Drop to ipython shell after calculation              
 tstep        1.0e-05      s            Time step size                                       
 tstop        0.001        s            Simulation stop time                                 
 verbose      0                         Show iterations for each point                       
 =========== ============ ============ ===================================================== 

