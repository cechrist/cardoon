==============
Device Library
==============
 
bjt
---


    Gummel-Poon intrinsic BJT model

    Terminal order: 0 Collector, 1 Base, 2 Emitter, (3 Bulk, not included)::

                      
      C (0) o----,         4----o  E (2)
                  \       /
                   \     /
                  ---------
                      |
                      o 
   
                      B (1)

    Can be used for NPN or PNP transistors.

    Internally may add up to 2 additional nodes (plus gnd) if rb is
    not zero: Bi(3) for the internal base node and, if rbm is
    specified, ib(4) to measure the internal base current and
    calculate Rb(ib)
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 xcjc       1.0                       Fraction of cbc connected internal to rb             
 itf        0.0          A            Transit time dependency on ic                        
 tr         0.0          S            Ideal reverse transit time                           
 eg         1.11         eV           Badgap voltage                                       
 tnom       27.0         C            Nominal temperature                                  
 area       1.0                       Current multiplier                                   
 cjc        0.0          F            Base collector zero bias p-n capacitance             
 nc         2.0                       Base-collector leakage emission coefficient          
 vaf        0.0          V            Forward early voltage                                
 ne         1.5                       Base-emitter leakage emission coefficient            
 nf         1.0                       Forward current emission coefficient                 
 cje        0.0          F            Base emitter zero bias p-n capacitance               
 xtf        0.0                       Transit time bias dependence coefficient             
 xtb        0.0                       Forward and reverse beta temperature coefficient     
 rb         0.0          W            Zero bias base resistance                            
 var        0.0          V            Reverse early voltage                                
 irb        0.0          A            Current at which rb falls to half of rbm             
 type       npn                       Type (npn or pnp)                                    
 xti        3.0                       IS temperature effect exponent                       
 ikr        0.0          A            Corner for reverse-beta high current roll off        
 bf         100.0                     Ideal maximum forward beta                           
 tf         0.0          S            Ideal forward transit time                           
 fc         0.5                       Forward bias depletion capacitor coefficient         
 ikf        0.0          A            Forward-beta high current roll-off knee current      
 br         1.0                       Ideal maximum reverse beta                           
 isat       1.0e-16      A            Transport saturation current                         
 nr         1.0                       Reverse current emission coefficient                 
 mje        0.33                      Base emitter p-n grading factor                      
 temp       None         C            Device temperature                                   
 mjc        0.33                      Base collector p-n grading factor                    
 vtf        0.0          V            Transit time dependency on vbc                       
 vjc        0.75         V            Base collector built in potential                    
 rbm        0.0          W            Minimum base resistance                              
 vje        0.75         V            Base emitter built in potential                      
 isc        0.0          A            Base collector leakage saturation current            
 ise        0.0          A            Base-emitter leakage saturation current              
 ========= ============ ============ ===================================================== 

bjt_t
-----

Electro-thermal version of bjt (extra thermal port)

cap
---


    Linear Capacitor::

                   || C
      0 o----------||---------o 1
                   ||

    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 c          0.0          F            Capacitance                                          
 ========= ============ ============ ===================================================== 

diode
-----


    Diode device (based on spice model)::
    
               o  1                           
               |                            
             --+--
              / \     
             '-+-' 
               |                          
               o  0 

    Includes depletion and diffusion charges.
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 tnom       27.0         C            Nominal temperature                                  
 af         1.0                       Flicker noise exponent                               
 vj         1.0          V            Built-in junction potential                          
 bv         0.0          V            Breakdown voltage                                    
 fc         0.5                       Coefficient for forward-bias depletion capacitance   
 isat       1.0e-14      A            Saturation current                                   
 kf         0.0                       Flicker noise coefficient                            
 temp       None         C            Device temperature                                   
 area       1.0                       Area multiplier                                      
 tt         0.0          s            Transit time                                         
 eg0        1.11         eV           Energy bandgap                                       
 m          0.5                       PN junction grading coefficient                      
 rs         0.0          Ohms         Series resistance                                    
 n          1.0                       Emission coefficient                                 
 ibv        1.0e-10      A            Current at reverse breakdown voltage                 
 cj0        0.0          F            Zero-bias depletion capacitance                      
 xti        3.0                       Is temperature exponent                              
 ========= ============ ============ ===================================================== 

diode_t
-------

Electro-thermal version of diode (extra thermal port)

idc
---


    DC current source. 

    Includes temperature dependence::

                    ______ 
                   /      \ idc
        0 o-------+  --->  +---------o 1
                   \______/  

    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 tnom       27.0         C            Nominal temperature                                  
 tc2        0.0          1/C^2        Current temperature coefficient 2                    
 idc        0.0          A            DC current                                           
 temp       None         C            Device temperature                                   
 tc1        0.0          1/C          Current temperature coefficient 1                    
 ========= ============ ============ ===================================================== 

ind
---


    Linear inductor::

                 __  __  __  _ 
        0       /  \/  \/  \/ \          1
          o----+   /\  /\  /\  +-------o    External view
                  (_/ (_/ (_/  

    Internal implementation uses a gyrator (adds one internal node
    plus uses gnd)::

                                          2
        0  o---------+            +----------------+
                     | gyr V2     |                |
          +         /|\          /^\               |
        Vin        | | |        | | | gyr Vin    ----- gyr^2 * L
          -         \V/          \|/             -----
                     |            |                |
        1  o---------+            +------+---------+
                                         |
                                        --- (terminal 3 here)
                                         V
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 l          0.0          H            Inductance                                           
 ========= ============ ============ ===================================================== 

mosacm
------


    Implements a simplified ACM MOSFET model. 

    Only (some) DC equations are considered for now.
    Terminal order: 0 Drain, 1 Gate, 2 Source, 3 Bulk::

               Drain 0
                       o
                       |
                       |
                   |---+
                   |
      Gate 1 o-----|<-----o 3 Bulk
                   |
                   |---+
                       |
                       |
                       o
              Source 2
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 phi        0.55         V            Surface Potential                                    
 vsat       80000.0      m/s          Saturation Velocity                                  
 temp       None         C            Device temperature                                   
 tox        7.5e-09      m            Oxide Thickness                                      
 l          1.0e-05      m            Channel length                                       
 vt0        0.532        V            Threshold Voltage                                    
 kp         0.0005106    A/V^2        Transconductance Parameter                           
 w          1.0e-05      m            Channel width                                        
 theta      0.814        1/V          Mobility Saturation Parameter                        
 gamma      0.631        V^(1/2)      Bulk Threshold Parameter                             
 ========= ============ ============ ===================================================== 

mosacm_t
--------

Electro-thermal version of mosacm (extra thermal port)

mosekv
------


    Intrinsic EPFL EKV 2.6 MOSFET::

        Terminal order: 0 Drain, 1 Gate, 2 Source, 3 Bulk
        
                 Drain 0
                         o
                         |
                         |
                     |---+
                     |
        Gate 1 o-----|<-----o 3 Bulk
                     |
                     |---+
                         |
                         |
                         o
                Source 2

    Mostly based on [1], but some updates from a later revision (dated
    1999) are also included.
    
    [1] The EPFL-EKV MOSFET Model Equations for Simulation, Technical
    Report, Model Version 2.6, June, 1997, Revision I, September,
    1997, Revision II, July, 1998, Bucher, Christophe Lallement,
    Christian Enz, Fabien Theodoloz, Francois Krummenacher,
    Electronics Laboratories, Swiss Federal Institute of Technology
    (EPFL), Lausanne, Switzerland
    
    This implementation includes accurate current interpolation
    function (optional), works for negative VDS and includes
    electrothermal model, DC operating point paramenters and noise
    equations.
    
    Code originally based on freeda 1.4 implementation
    <http://www.freeda.org>::
    
        // Element information
        ItemInfo Mosnekv::einfo =
        {
          "mosnekv",
          "EPFL EKV MOSFET model",
          "Wonhoon Jang",
          DEFAULT_ADDRESS"transistor>mosfet",
          "2003_05_15"
        };
    
    Parameter limit checking, simple capacitance calculations for
    operating point are not yet implemented.
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 xj         1.0e-07      m            Junction depth                                       
 ekvint     0                         Interpolation function (0: accurate, 1: simple)      
 af         1.0                       Flicker noise exponent                               
 tcv        0.001        V/K          Threshold voltage temperature coefficient            
 avto       0.0          Vm           Area related threshold voltage mismatch parameter    
 ibbt       0.0009       1/K          Temperature coefficient for IBB                      
 tnom       27.0         C            Nominal temperature of model parameters              
 ucex       0.8                       Longitudinal critical field temperature exponent     
 lk         2.9e-07      m            Reverse short channel effect characteristic length   
 leta       0.1                       Short-channel effect coefficient                     
 q0         0.0          A.s/m^2      Reverse short channel effect peak charge density     
 tox        None         m            Oxide thickness                                      
 u0         None         cm^2/(V.s)   Low-field mobility                                   
 np         1.0                       Parallel multiple device number                      
 theta      0.0          1/V          Mobility recuction coefficient                       
 ns         1.0                       Serial multiple device number                        
 type       n                         N- or P-channel MOS (n or p)                         
 ucrit      2.0e+06      V/m          Longitudinal critical field                          
 phi        0.7          V            Bulk Fermi potential                                 
 ibn        1.0                       Saturation voltage factor for impact ionization      
 vmax       None         m/s          Saturation velocity                                  
 dw         0.0          m            Channel width correction                             
 vfb        None         V            Flat-band voltage                                    
 e0         1.0e+12      V/m          Mobility reduction coefficient                       
 agamma     0.0          V^(1/2)m     Area related body effect mismatch parameter          
 Lambda     0.5                       Channel-length modulation                            
 dl         0.0          m            Channel length correction                            
 kf         0.0                       Flicker noise coefficient                            
 temp       None         C            Device temperature                                   
 satlim     54.5982                   Ratio defining the saturation limit if/ir            
 nsub       None         1/cm^3       Channel doping                                       
 ibb        3.0e+08      V/m          Second impact ionization coefficient                 
 akp        0.0          m            Area related gain mismatch parameter                 
 l          1.0e-06      m            Gate length                                          
 vt0        0.5          V            Long_channel threshold voltage                       
 bex        -1.5                      Mobility temperature exponent                        
 kp         5.0e-05      A/V^2        Transconductance parameter                           
 w          1.0e-06      m            Gate width                                           
 iba        0.0          1/m          First impact ionization coefficient                  
 weta       0.25                      Narrow-channel effect coefficient                    
 cox        0.0007       F/m^2        Gate oxide capacitance per area                      
 gamma      1.0          V^1/2        Body effect parameter                                
 ========= ============ ============ ===================================================== 

mosekv_t
--------

Electro-thermal version of mosekv (extra thermal port)

res
---


    Linear Resistor::

                    R
      0 o--------/\/\/\/---------o 1

    If the electro-thermal version is used (res_t), the device is
    nonlinear.
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 tnom       27.0         C            Nominal temperature                                  
 rsh        0.0          Ohms         Sheet resistance                                     
 temp       None         C            Device temperature                                   
 tc2        0.0          1/C^2        Temperature coefficient 2                            
 l          0.0          m            Lenght                                               
 tc1        0.0          1/C          Temperature coefficient 1                            
 r          0.0          Ohms         Resistance                                           
 w          0.0          m            Width                                                
 narrow     0.0          m            Narrowing due to side etching                        
 ========= ============ ============ ===================================================== 

res_t
-----

Electro-thermal version of res (extra thermal port)

svdiode
-------


    State-Variable-Based Diode device (based on Spice model)::

            o  1                           
            |                            
          --+--
           / \     
          '-+-'
            |                          
            o  0    	                  

    Internally represented as::

        0  o
           |
           \ 
           / Rs
           \ 
           / 
           |                                     2
        4  o---------+                  +----------------+
                     | i(x)+dq/dt       |                |
          +         /|\                /|\ gyr vin      /^\ 
        vin        | | |              | | |            | | | gyr v(x)
          -         \V/                \V/              \|/  
                     |                  |                |
        1  o---------+                  +------+---------+
                                               |
                                              --- (terminal 3 is gnd)
                                               V

    Terminal 4 not present if Rs = 0

    Implementation includes depletion and diffusion charges.
    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 tnom       27.0         C            Nominal temperature                                  
 af         1.0                       Flicker noise exponent                               
 vj         1.0          V            Built-in junction potential                          
 bv         0.0          V            Breakdown voltage                                    
 fc         0.5                       Coefficient for forward-bias depletion capacitance   
 isat       1.0e-14      A            Saturation current                                   
 kf         0.0                       Flicker noise coefficient                            
 temp       None         C            Device temperature                                   
 area       1.0                       Area multiplier                                      
 tt         0.0          s            Transit time                                         
 eg0        1.11         eV           Energy bandgap                                       
 m          0.5                       PN junction grading coefficient                      
 rs         0.0          Ohms         Series resistance                                    
 n          1.0                       Emission coefficient                                 
 ibv        1.0e-10      A            Current at reverse breakdown voltage                 
 cj0        0.0          F            Zero-bias depletion capacitance                      
 xti        3.0                       Is temperature exponent                              
 ========= ============ ============ ===================================================== 

svdiode_t
---------

Electro-thermal version of svdiode (extra thermal port)

vdc
---


    DC voltage source. 

    Includes temperature dependence in vdc only::
   
                   ______ 
                  /      \ vdc       Rint
       0 o-------(  -  +  )--------/\/\/\/\--------o 1
                  \______/ 
   
    Implemented using a gyrator if Rint is zero::

                                  2       V2
        0  o---------+            +----------------+
                     | gyr V2     |                |
          +         /|\          /|\              /^\ 
        vin        | | |        | | | gyr vin    | | | gyr vdc
          -         \V/          \V/              \|/  
                     |            |                |
        1  o---------+            +------+---------+
                                  3      |
                                        --- (terminal 3 here)
                                         V  

    

Parameters
++++++++++

 ========= ============ ============ ===================================================== 
 Name       Default      Unit         Description                                          
 ========= ============ ============ ===================================================== 
 tnom       27.0         C            Nominal temperature                                  
 temp       None         C            Device temperature                                   
 rint       0.0          Ohms         Internal resistance                                  
 tc2        0.0          1/C^2        Voltage temperature coefficient 2                    
 vdc        0.0          V            DC current                                           
 tc1        0.0          1/C          Voltage temperature coefficient 1                    
 ========= ============ ============ ===================================================== 

