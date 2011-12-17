========================
Analysis Library Catalog
========================
 
dc
--


DC Sweep Calculation

Calculates a DC sweep of a circuit using the nodal approach. Nodal
voltages are saved after the analysis is complete.

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``.

After completion the analysis drops to an interactive shell if the
``shell`` global variable is set to ``True``


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 device                                 Instance name of device to sweep variable            
 fullAD       0                         Use CPPAD for entire nonlinear part                  
 param                                  Device parameter to sweep                            
 start        0.0          V            Sweep start value                                    
 stop         0.0          V            Sweep stop value                                     
 sweep_num    50                        Number of points in sweep                            
 verbose      0                         Show iterations for each point                       
 =========== ============ ============ ===================================================== 

op
--


DC Operating Point Calculation

Calculates the DC operating point of a circuit using the nodal
approach. Nodal voltages and nonlinear device operating points are
saved after the analysis is complete.

By default the voltage at all external voltages is printed after
the analysis is complete. Optionally the operating points of
nonlinear elements can be printed. 

Convergence parameters for the Newton method are controlled using
the global variables in ``.options``.

After completion the analysis drops to an interactive shell if the
``shell`` global variable is set to ``True``


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 elemop       0                         Print element operating points                       
 fullAD       0                         Use CPPAD for entire nonlinear part                  
 intvars      0                         Print internal element nodal variables               
 =========== ============ ============ ===================================================== 

testdev
-------


Test equations of one nonlinear device

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
 start        0.0          V            Sweep start value                                    
 stop         0.0          V            Sweep stop value                                     
 sweep_num    0                         Number of points in sweep                            
 sweep_port   0                         Port number to be swept, starting from zero          
 useAD        1                         Use automatic differentiation                        
 =========== ============ ============ ===================================================== 

