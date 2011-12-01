======================
Device Library Catalog
======================
 
bjt
---


Gummel-Poon intrinsic BJT model

This implementation based mainly on previous implementation in
carrot and some equations from Pspice manual.

Terminal order: 0 Collector, 1 Base, 2 Emitter::

                  
      C (0) o----,         4----o  E (2)
                  \       /
                   \     /
                  ---------
                      |
                      o 
    
                      B (1)

Can be used for NPN or PNP transistors.

Bulk connection, RC, RE are not included here.

Netlist examples::

    bjt:q1 2 3 4 model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m

    # Electro-thermal version
    bjt_t:q2 2 3 5 pout gnd model = mypnp

    # Model statement
    .model mypnp bjt_t (type=pnp isat=5e-17 cje=60fF vje=0.83 mje=0.35)

Internal Topology
+++++++++++++++++

Internally may add 2 additional nodes (plus gnd) if rb is not
zero: Bi(3) for the internal base node and ib(4) to measure the
internal base current and calculate Rb(ib). The possible
configurations are described here.

1. If RB == 0::

                     +----------------+--o 0 (C)
                     |                |
                    /^\               |
                   ( | ) ibc(vbc)     |
                    \|/               |       
                     |               /|\       
     (B) 1 o---------+              ( | ) ice    
                     |               \V/      
                    /|\               |       
                   ( | ) ibe(vbe)     |
                    \V/               |
                     |                |
                     +----------------+--o 2 (E)

2. If RB != 0::

                                 +----------------+--o 0 (C)
                                 |                |
                                /^\               |
                   ib          ( | ) ibc(vbc)     |
                                \|/               |       
                 ,---,           |               /|\       
     (B) 1 o----( --> )----------+ 3 (Bi)       ( | ) ice    
                 `---`           |               \V/      
                                /|\               |       
                               ( | ) ibe(vbe)     |
                                \V/               |
                                 |                |
                 gyr v13         +----------------+--o 2 (E)
                              
                  ,---,       
             +---( <-- ) -----+
             |    `---`       |
     lref    |                | ib/gyr
     (5) ,---+                |
         |   |    ,---,       | 4 (ib)
         |   +---( --> )------+
         |        `---`       
        ---
         V      gyr ib Rb(ib)
                                       
Charge sources are connected between internal nodes defined
above. If xcjc is not 1 but RB is zero, xcjc is ignored.



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 area         1.0                       Current multiplier                                   
 bf           100.0                     Ideal maximum forward beta                           
 br           1.0                       Ideal maximum reverse beta                           
 cjc          0.0          F            Base collector zero bias p-n capacitance             
 cje          0.0          F            Base emitter zero bias p-n capacitance               
 eg           1.11         eV           Badgap voltage                                       
 fc           0.5                       Forward bias depletion capacitor coefficient         
 ikf          0.0          A            Forward-beta high current roll-off knee current      
 ikr          0.0          A            Corner for reverse-beta high current roll off        
 irb          0.0          A            Current at which rb falls to half of rbm             
 isat         1.0e-16      A            Transport saturation current                         
 isc          0.0          A            Base collector leakage saturation current            
 ise          0.0          A            Base-emitter leakage saturation current              
 itf          0.0          A            Transit time dependency on ic                        
 mjc          0.33                      Base collector p-n grading factor                    
 mje          0.33                      Base emitter p-n grading factor                      
 nc           2.0                       Base-collector leakage emission coefficient          
 ne           1.5                       Base-emitter leakage emission coefficient            
 nf           1.0                       Forward current emission coefficient                 
 nr           1.0                       Reverse current emission coefficient                 
 rb           0.0          W            Zero bias base resistance                            
 rbm          0.0          W            Minimum base resistance                              
 temp         None         C            Device temperature                                   
 tf           0.0          S            Ideal forward transit time                           
 tnom         27.0         C            Nominal temperature                                  
 tr           0.0          S            Ideal reverse transit time                           
 type         npn                       Type (npn or pnp)                                    
 vaf          0.0          V            Forward early voltage                                
 var          0.0          V            Reverse early voltage                                
 vjc          0.75         V            Base collector built in potential                    
 vje          0.75         V            Base emitter built in potential                      
 vtf          0.0          V            Transit time dependency on vbc                       
 xcjc         1.0                       Fraction of cbc connected internal to rb             
 xtb          0.0                       Forward and reverse beta temperature coefficient     
 xtf          0.0                       Transit time bias dependence coefficient             
 xti          3.0                       IS temperature effect exponent                       
 =========== ============ ============ ===================================================== 

bjt_t
-----

Electro-thermal version of bjt (extra thermal port)

cap
---


Linear Capacitor::

               || C
  0 o----------||---------o 1
               ||

Netlist example::

    cap:c1 1 2 c=10uF



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 c            0.0          F            Capacitance                                          
 =========== ============ ============ ===================================================== 

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

Netlist examples::

    diode:d1 1 0 isat=10fA cj0=20fF

    # Electrothermal device
    diode_t:d2 2 3 1000 gnd cj0=10pF tt=1e-12 rs=100 bv = 4.

    # Model statement
    .model dmodel1 diode (cj0 = 10pF tt=1ps)

Internal Topology
+++++++++++++++++

The internal representation is the following::

    0  o
       |
       \ 
       / Rs
       \ 
       / 
       |                      
    2  o---------+            
                 | i(vin)+dq/dt 
      +         /|\           
    vin        | | |          
      -         \V/           
                 |            
    1  o---------+            
                              
                              
                              

Terminal 2 not present if Rs = 0


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 af           1.0                       Flicker noise exponent                               
 area         1.0                       Area multiplier                                      
 bv           0.0          V            Breakdown voltage                                    
 cj0          0.0          F            Zero-bias depletion capacitance                      
 eg0          1.11         eV           Energy bandgap                                       
 fc           0.5                       Coefficient for forward-bias depletion capacitance   
 ibv          1.0e-10      A            Current at reverse breakdown voltage                 
 isat         1.0e-14      A            Saturation current                                   
 kf           0.0                       Flicker noise coefficient                            
 m            0.5                       PN junction grading coefficient                      
 n            1.0                       Emission coefficient                                 
 rs           0.0          Ohms         Series resistance                                    
 temp         None         C            Device temperature                                   
 tnom         27.0         C            Nominal temperature                                  
 tt           0.0          s            Transit time                                         
 vj           1.0          V            Built-in junction potential                          
 xti          3.0                       Is temperature exponent                              
 =========== ============ ============ ===================================================== 

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

Netlist example::

    idc:vdd gnd 4 idc=2mA



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 idc          0.0          A            DC current                                           
 tc1          0.0          1/C          Current temperature coefficient 1                    
 tc2          0.0          1/C^2        Current temperature coefficient 2                    
 temp         None         C            Device temperature                                   
 tnom         27.0         C            Nominal temperature                                  
 =========== ============ ============ ===================================================== 

ind
---


Linear inductor::

             __  __  __  _ 
    0       /  \/  \/  \/ \          1
      o----+   /\  /\  /\  +-------o    External view
              (_/ (_/ (_/  

Netlist example::

    ind:l1 1 0 l=3uH


Internal Topology
+++++++++++++++++

Internal implementation uses a gyrator (adds one internal node)::

                                    il/gyr  (2)
    0  o---------+            +----------------+
                 | gyr V23    |                |
      +         /|\          /^\               |
    Vin        | | |        | | | gyr Vin    ----- gyr^2 * L
      -         \V/          \|/             -----
                 |            |                |
    1  o---------+            +----------------+
                                      |
                                     --- lref (3)
                                      V




Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 l            0.0          H            Inductance                                           
 =========== ============ ============ ===================================================== 

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

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 gamma        0.631        V^(1/2)      Bulk Threshold Parameter                             
 kp           0.0005106    A/V^2        Transconductance Parameter                           
 l            1.0e-05      m            Channel length                                       
 phi          0.55         V            Surface Potential                                    
 temp         None         C            Device temperature                                   
 theta        0.814        1/V          Mobility Saturation Parameter                        
 tox          7.5e-09      m            Oxide Thickness                                      
 vsat         80000.0      m/s          Saturation Velocity                                  
 vt0          0.532        V            Threshold Voltage                                    
 w            1.0e-05      m            Channel width                                        
 =========== ============ ============ ===================================================== 

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

Code originally based on fREEDA 1.4 implementation
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

Netlist examples::

    mosekv:m1 2 3 4 gnd w=30e-6 l=1e-6 type = n ekvint=0

    # Electro-thermal version
    mosekv_t:m1 2 3 4 gnd 1000 gnd w=30e-6 l=1e-6 type = n

    # Model statement
    .model ekvn mosekv (type = n kp = 200u theta = 0.6)

Internal Topology
+++++++++++++++++

The internal topology is the following::

                                  +-------------+--o 0 (D)
                                  |             |
                                  |             |
                                -----           |
                                ----- qd        |       
                                  |            /|\       
     (G) 1 o---------+            |           | | | ids    
                     |            |            \V/      
                     |            |             |       
                   -----          |             |
                   ----- qg       |      qs     |
                     |            |      ||     |
     (B) 4 o---------+------------+------||-----+--o 2 (S)
                                         ||

The impact ionization current is normally added to the drain
current, but if the device is in reverse (Vds < 0 for N-channel)
mode, it is added to the source current.


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 Lambda       0.5                       Channel-length modulation                            
 af           1.0                       Flicker noise exponent                               
 agamma       0.0          V^(1/2)m     Area related body effect mismatch parameter          
 akp          0.0          m            Area related gain mismatch parameter                 
 avto         0.0          Vm           Area related threshold voltage mismatch parameter    
 bex          -1.5                      Mobility temperature exponent                        
 cox          0.0007       F/m^2        Gate oxide capacitance per area                      
 dl           0.0          m            Channel length correction                            
 dw           0.0          m            Channel width correction                             
 e0           1.0e+12      V/m          Mobility reduction coefficient                       
 ekvint       0                         Interpolation function (0: accurate, 1: simple)      
 gamma        1.0          V^1/2        Body effect parameter                                
 iba          0.0          1/m          First impact ionization coefficient                  
 ibb          3.0e+08      V/m          Second impact ionization coefficient                 
 ibbt         0.0009       1/K          Temperature coefficient for IBB                      
 ibn          1.0                       Saturation voltage factor for impact ionization      
 kf           0.0                       Flicker noise coefficient                            
 kp           5.0e-05      A/V^2        Transconductance parameter                           
 l            1.0e-06      m            Gate length                                          
 leta         0.1                       Short-channel effect coefficient                     
 lk           2.9e-07      m            Reverse short channel effect characteristic length   
 np           1.0                       Parallel multiple device number                      
 ns           1.0                       Serial multiple device number                        
 nsub         None         1/cm^3       Channel doping                                       
 phi          0.7          V            Bulk Fermi potential                                 
 q0           0.0          A.s/m^2      Reverse short channel effect peak charge density     
 satlim       54.5982                   Ratio defining the saturation limit if/ir            
 tcv          0.001        V/K          Threshold voltage temperature coefficient            
 temp         None         C            Device temperature                                   
 theta        0.0          1/V          Mobility recuction coefficient                       
 tnom         27.0         C            Nominal temperature of model parameters              
 tox          None         m            Oxide thickness                                      
 type         n                         N- or P-channel MOS (n or p)                         
 u0           None         cm^2/(V.s)   Low-field mobility                                   
 ucex         0.8                       Longitudinal critical field temperature exponent     
 ucrit        2.0e+06      V/m          Longitudinal critical field                          
 vfb          None         V            Flat-band voltage                                    
 vmax         None         m/s          Saturation velocity                                  
 vt0          0.5          V            Long_channel threshold voltage                       
 w            1.0e-06      m            Gate width                                           
 weta         0.25                      Narrow-channel effect coefficient                    
 xj           1.0e-07      m            Junction depth                                       
 =========== ============ ============ ===================================================== 

mosekv_t
--------

Electro-thermal version of mosekv (extra thermal port)

res
---


Resistor::

                R
  0 o--------/\/\/\/---------o 1

Normally a linear device. If the electro-thermal version is used
(res_t), the device is nonlinear.

Netlist examples::

    # Linear resistor (2 terminals)
    res:r1 1 2 r=1e3 tc1=10e-3

    # Electro-thermal resistor (nonlinear, 4 terminals)
    res_t:r1 1 2 3 4 r=1e3 tc1=10e-3



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 l            0.0          m            Lenght                                               
 narrow       0.0          m            Narrowing due to side etching                        
 r            0.0          Ohms         Resistance                                           
 rsh          0.0          Ohms         Sheet resistance                                     
 tc1          0.0          1/C          Temperature coefficient 1                            
 tc2          0.0          1/C^2        Temperature coefficient 2                            
 temp         None         C            Device temperature                                   
 tnom         27.0         C            Nominal temperature                                  
 w            0.0          m            Width                                                
 =========== ============ ============ ===================================================== 

res_t
-----

Electro-thermal version of res (extra thermal port)

svbjt
-----


State-variable-based Gummel-Poon intrinsic BJT model based

This implementation based mainly on previous implementation in
carrot and some equations from Pspice manual, with the addition of
the state-variable definitions.

Terminal order: 0 Collector, 1 Base, 2 Emitter, (3 Bulk, not included)::

                  
  C (0) o----,         4----o  E (2)
              \       /
               \     /
              ---------
                  |
                  o 

                  B (1)

Can be used for NPN or PNP transistors.

Bulk connection, RC, RE are not included for now.

Netlist examples::

    bjt:q1 2 3 4 model = mypnp isat=4e-17 bf=147 vaf=80 ikf=4m

    # Electro-thermal version
    bjt_t:q2 2 3 5 pout gnd model = mypnp

    # Model statement
    .model mypnp bjt_t (type=pnp isat=5e-17 cje=60fF vje=0.83 mje=0.35)

Internal Topology
+++++++++++++++++

The state variable formulation is achieved by replacing the BE and
BC diodes (Ibf, Ibr) with state-variable based diodes. This
requires two additional variables (nodes) but eliminates large
positive exponentials from the model::

                              3 (x2)
                  +--------------------------+
                  |                          |
                 /|\                        /^\ 
                ( | ) gyr v2               ( | ) gyr vbc(x)
                 \V/                        \|/  
         lref     |                          |
         (5) ,----+--------------------------+ 
             |    |                          |               
             |   /^\                        /|\              
             |  ( | ) gyr v1               ( | ) gyr vbe(x)  
            ---  \|/                        \V/  
             V    |                          |
                  +--------------------------+
                               4 (x1)               
                                              
All currents/charges in the model are functions of voltages v3
(x2) and v4 (x1). Note that vbc and vbe are now also functions of
x1, x2.

In addition we may need 2 additional nodes (plus gnd) if rb is not
zero: Bi(3) for the internal base node and ib(4) to measure the
internal base current and calculate Rb(ib).

1. If RB == 0::

                       +----------------+--o 0 (C)
                -      |                |
                      /^\               |
               v2    ( | ) ibc(x2)      |
                      \|/               |       
                +      |               /|\       
       (B) 1 o---------+              ( | ) ice(x1,x2)
                +      |               \V/      
                      /|\               |       
               v1    ( | ) ibe(x1)      |
                      \V/               |
                -      |                |
                       +----------------+--o 2 (E)

2. If RB != 0 and IRB != 0::

                                 +----------------+--o 0 (C)
                            -    |                |
                                /^\               |
                gyr v75    v2  ( | ) ibc(x2)      |
                                \|/               |       
                 ,---,      +    |               /|\       
     (B) 1 o----( --> )----------+ 6 (Bi)       ( | ) ice(x1,x2)
                 `---`      +    |               \V/      
                                /|\               |       
                           v1  ( | ) ibe(x1)      |
                                \V/               |
                            -    |                |
                 gyr v16         +----------------+--o 2 (E)
                              
                  ,---,       
             +---( <-- ) -----+
             |    `---`       |
      lref   |                | ib/gyr
      (5) ,--+                |
          |  |    ,---,       | 7 (ib)
          |  +---( --> )------+
          |       `---`       
         --- 
          V     gyr ib Rb(ib)
                                       
Charge sources are connected between internal nodes defined
above. If xcjc is not 1 but RB is zero, xcjc is ignored.


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 area         1.0                       Current multiplier                                   
 bf           100.0                     Ideal maximum forward beta                           
 br           1.0                       Ideal maximum reverse beta                           
 cjc          0.0          F            Base collector zero bias p-n capacitance             
 cje          0.0          F            Base emitter zero bias p-n capacitance               
 eg           1.11         eV           Badgap voltage                                       
 fc           0.5                       Forward bias depletion capacitor coefficient         
 ikf          0.0          A            Forward-beta high current roll-off knee current      
 ikr          0.0          A            Corner for reverse-beta high current roll off        
 irb          0.0          A            Current at which rb falls to half of rbm             
 isat         1.0e-16      A            Transport saturation current                         
 isc          0.0          A            Base collector leakage saturation current            
 ise          0.0          A            Base-emitter leakage saturation current              
 itf          0.0          A            Transit time dependency on ic                        
 mjc          0.33                      Base collector p-n grading factor                    
 mje          0.33                      Base emitter p-n grading factor                      
 nc           2.0                       Base-collector leakage emission coefficient          
 ne           1.5                       Base-emitter leakage emission coefficient            
 nf           1.0                       Forward current emission coefficient                 
 nr           1.0                       Reverse current emission coefficient                 
 rb           0.0          W            Zero bias base resistance                            
 rbm          0.0          W            Minimum base resistance                              
 temp         None         C            Device temperature                                   
 tf           0.0          S            Ideal forward transit time                           
 tnom         27.0         C            Nominal temperature                                  
 tr           0.0          S            Ideal reverse transit time                           
 type         npn                       Type (npn or pnp)                                    
 vaf          0.0          V            Forward early voltage                                
 var          0.0          V            Reverse early voltage                                
 vjc          0.75         V            Base collector built in potential                    
 vje          0.75         V            Base emitter built in potential                      
 vtf          0.0          V            Transit time dependency on vbc                       
 xcjc         1.0                       Fraction of cbc connected internal to rb             
 xtb          0.0                       Forward and reverse beta temperature coefficient     
 xtf          0.0                       Transit time bias dependence coefficient             
 xti          3.0                       IS temperature effect exponent                       
 =========== ============ ============ ===================================================== 

svbjt_t
-------

Electro-thermal version of svbjt (extra thermal port)

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

This model has better convergence properties. Externally it
behaves exactly like the regular diode device. 

Implementation includes depletion and diffusion charges. 

Netlist examples::

    svdiode:d1 1 0 isat=10fA cj0=20fF

    # Electrothermal device
    svdiode_t:d2 2 3 1000 gnd cj0=10pF tt=1e-12 rs=100 bv = 4.

    # Model statement
    .model dmodel1 svdiode (cj0 = 10pF tt=1ps)

Internal Topology
+++++++++++++++++

The internal representation is the following::

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
    1  o---------+                  +--------+-------+
                                             |
                                            --- lref (3)
                                             V

Terminal 4 not present if Rs = 0


Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 af           1.0                       Flicker noise exponent                               
 area         1.0                       Area multiplier                                      
 bv           0.0          V            Breakdown voltage                                    
 cj0          0.0          F            Zero-bias depletion capacitance                      
 eg0          1.11         eV           Energy bandgap                                       
 fc           0.5                       Coefficient for forward-bias depletion capacitance   
 ibv          1.0e-10      A            Current at reverse breakdown voltage                 
 isat         1.0e-14      A            Saturation current                                   
 kf           0.0                       Flicker noise coefficient                            
 m            0.5                       PN junction grading coefficient                      
 n            1.0                       Emission coefficient                                 
 rs           0.0          Ohms         Series resistance                                    
 temp         None         C            Device temperature                                   
 tnom         27.0         C            Nominal temperature                                  
 tt           0.0          s            Transit time                                         
 vj           1.0          V            Built-in junction potential                          
 xti          3.0                       Is temperature exponent                              
 =========== ============ ============ ===================================================== 

svdiode_t
---------

Electro-thermal version of svdiode (extra thermal port)

tlinps4
-------


4-terminal physical transmission line model using scattering
parameters::

         0 o===================================o 2
                           Z0
         1 o===================================o 3


This model is similar to tlinpy4, but it is more robust and can
handle lossless lines, even at DC, but internally requires 2
additional ports to keep track of v1+ and v2+. This model is more
suitable for convolution as the S parameters are better behaved
than the Y parameters.

Netlist Examples::

  tlinps4:tl1 in gnd out gnd z0mag=100. length=0.3m
  .model c_line tlins4 (z0mag=75.00 k=7 fscale=1.e10 alpha = 59.9)

Internal Topology
+++++++++++++++++

The model is symmetric. The schematic for Port 1 is shown here::

           I1                              v1+ + v1-          v1-
          --->                               ---->     v1+   ---->
      0 o--------,                          ,------------+----------,  4
   +             |                          |            |          |  
                 |                          |           ,-,  s12 v2+|  
  V1            /|\ (v1+ - s12 v2+)/Z0     /^\          | |        /|\ 
               ( | )                      ( | )       1 | |       ( | )
   -            \V/                    V1  \|/          '-'        \V/ 
                 |                          |            |          |  
      1 o--------+                          +---------+--+----------'   
                                                      |
                                                     --- lref (6)
                                                      V


Note: for a matched transmission line, s11 = s22 = 0 and s12 =
s21. The equivalent 'Y' matrix is::

           |              1/Z0    -s12/Z0 |
           |                              |
           |             -s21/Z0    1/Z0  |           
       Y = |                              |
           | -1            1        s12   |
           |                              |
           |        -1    s21        1    |



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 alpha        0.1          dB/m         Attenuation                                          
 fscale       0.0          Hz           Scaling frequency for attenuation                    
 k            1.0                       Effective relative dielectric constant               
 length       0.1          m            Line length                                          
 tand         0.0                       Loss tangent                                         
 z0mag        50.0         Ohms         Magnitude of characteristic impedance                
 =========== ============ ============ ===================================================== 

tlinpy4
-------


4-terminal physical transmission line model using Y parameters::


         0 o===================================o 2
                           Z0
         1 o===================================o 3


Code derived from fREEDA tlinp4 element. fREEDA implementation by
Carlos E. Christoffersen, Mete Ozkar, Michael Steer

Two models are supported dependent on the secting of nsect: When
``nsect = 0`` (not set) the frequency-domain model is enabled.
When ``nsect > 0`` the transmission line is expanded in 
``nsect`` RLCG subsections.

Netlist Examples::

  tlinpy4:tl1 in gnd out gnd z0mag=100. length=0.3m
  .model c_line tlinpy4 (z0mag=75.00 k=7 fscale=1.e10 alpha = 59.9)


Internal Topology
+++++++++++++++++

The internal schematic is the following::
             
      0 o----+------,               ,-----+-------o 2
   +         |      |               |     |              +
            ,-,     |               |    ,-, 
  v1        | |    /|\ y12 v2      /|\   | |             v2
        y11 | |   ( | )           ( | )  | | y22
   -        '-'    \V/      y21 v1 \V/   '-'             -
             |      |               |     |  
      1 o----+------'               '-----+-------o 3

                   y11 = y22 , y12 = y21



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 alpha        0.1          dB/m         Attenuation                                          
 fopt         0            Hz           Optimum frequency for discrete approximation         
 fscale       0.0          Hz           Scaling frequency for attenuation                    
 k            1.0                       Effective relative dielectric constant               
 length       0.1          m            Line length                                          
 nsect        0                         Enable discrete approximation with n sections        
 tand         0.0                       Loss tangent                                         
 z0mag        50.0         Ohms         Magnitude of characteristic impedance                
 =========== ============ ============ ===================================================== 

vdc
---


DC voltage source. 

Includes temperature dependence in vdc only::
                      
               ,---,  vdc       Rint
   0 o--------( - + )---------/\/\/\/\--------o 1
               '---'  

Netlist example::

    vdc:vdd 1 0 vdc=3V


Internal Topology
+++++++++++++++++

Implemented using a gyrator if Rint is zero::

                                   i/gyr (2)
    0  o---------+            +----------------+
                 | gyr V23    |                |
      +         /|\          /|\              /^\ 
    vin        | | |        | | | gyr vin    | | | gyr vdc
      -         \V/          \V/              \|/  
                 |            |                |
    1  o---------+            +----------------+
                                      |
                                     --- lref (3)
                                      V



Parameters
++++++++++

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 rint         0.0          Ohms         Internal resistance                                  
 tc1          0.0          1/C          Voltage temperature coefficient 1                    
 tc2          0.0          1/C^2        Voltage temperature coefficient 2                    
 temp         None         C            Device temperature                                   
 tnom         27.0         C            Nominal temperature                                  
 vdc          0.0          V            DC current                                           
 =========== ============ ============ ===================================================== 

