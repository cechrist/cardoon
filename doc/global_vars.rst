===================================
Global Options and Constants Tables
===================================
 
Options Variables
-----------------

The following variables can be set from the netlist using the
``.options`` keyword. Available options are listed below:

.. include:: options_table

Global options can be accessed within an interactive shell in the
program. For example::

  cechrist@moon:~/$ cardoon -i
  
  Cardoon Circuit Simulator 0.4 release 0.4.1.dev
  Type CTR-D to exit
  In <1>: glVar.gyr
  Out<1>: 0.001
  In <2>: glVar.abstol
  Out<2>: 1e-08
  In <3>: 

Physical and Mathematical Constants
-----------------------------------

Physical constants are used within models and they are also available
from interactive shells::

    In <4>: from globalVars import const
    In <5>: const.epsilon0
    Out<5>: 8.8541878170000005e-12

Table of physical and mathematical constants:

.. include:: constants_table

