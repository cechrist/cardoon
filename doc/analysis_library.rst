========================
Analysis Library Catalog
========================
 
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

