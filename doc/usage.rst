

Usage
=====

From the command line
---------------------

In the most common usage mode, the main program reads a netlist file,
builds the circuit described there and runs any specified
analyses. Optionally some analyses can drop to an ipython shell with
access to internal variables after calculations are finished for
interactive work.

Sample session
++++++++++++++

After the source is installed, change into the cardoon directory and
run the program with no arguments for a brief help message::

    cechrist@venus:~$ python -m cardoon
    
    Cardoon Circuit Simulator 0.6 release 0.6.0.dev
    Usage:
            cardoon <netlistname>        : Process netlist file
            cardoon -c                   : Generate catalogs
            cardoon -x <script> <args>   : execute Python script

In this document the first form of execution listed above is
described. The last two forms allow execution of arbitrary Python code
with access to the internal simulator classes and are not documented
here (for now).

To run one netlist, change into the ``cardoon/workspace`` directory
and execute the program as follows::

    cechrist@venus:~/src/cardoon/workspace$ python -m cardoon bias_npn.net 
    
    Cardoon Circuit Simulator 0.6 release 0.6.0.dev
    ******************************************************
                 Operating point analysis
    ******************************************************
    
     # Test of a transistor device 
    
    Using dense matrices
    
    Number of iterations =  17
    Residual =  5.72923228467e-06
    
     Node      |  Value               | Unit 
    ----------------------------------------
    1          |        3.49950994522 | V
    10         |        11.9158367321 | V
    2          |       0.800179840214 | V
    3          |                 10.0 | V
    gnd        |                  0.0 | V
    
    Element:  svbjt:q1
    
        Internal nodal variables:
    
        et         :     0.00425607912346 V
        ct         :        3.47695418943 V
        x2         :       -2.69945983799 s.v.
        x1         :        290.954259747 s.v.
        Bi         :       0.800050107221 V
        ib         :      0.0958314599978 0.001 A
    
        Operating point info:
    
        IB         : 9.58314599978e-05
        IC         : 0.00841632678691
        IE         : -0.00851215824691
        Power      : 0.0294934732763
        Temp       : 29.1439
        VBE        : 0.795794028098
        VCE        : 3.49525386609
    
    
    op analysis time: 0.01 s
    
    Press [Enter] to exit ...


Interactive Session
-------------------

Interactive use is in experimental mode and likely to change. The main
functions are provided in the ``simulator`` module. For example the
``ring_osc_ahkab.net`` from the ``/examples`` directory is run as
follows::

    In [1]: import cardoon.simulator as cs
    
    In [2]: cs.run
    cs.run_analyses  cs.run_netlist   
    
    In [2]: cs.run_netlist('ring_osc_ahkab.net')
    ******************************************************
                     Transient analysis
    ******************************************************
    
     ### RING OSCILLATOR with 3 inverters ### 
    
    Using dense matrices
    
    Calculating DC operating point ... Succeded.
    
    System dimension: 7
    Number of samples: 200
    Integration method: trap
    Printing one dot every 50 samples:
    ...
    Average iterations: 2.345
    Average residual: 1.31511609988e-06
    
    tran analysis time: 1.24 s

    In [3]: 

(a plot should pop-up).

.. automodule:: cardoon.simulator
       :members:


Netlist Format
--------------

A very brief description is provided here. The netlist syntax
resembles somewhat the syntax used in other simulators such as spice,
fREEDA and QUCS, but at least for now it has some simplifications. The
netlist is case-sensitive. Each line specifies one circuit element, an
analysis to perform or another command. In general lines can be
entered in any order. Order is important only to define subcircuit
blocks.

**Title line**

  The first line in the netlist is used to set the title and does not
  follow any particular syntax.

**End of netlist line**

  It a ``.end`` line is read, the parser stops reading the file and any
  additional lines are ignored.

**Line continuation**

  The backslash ("\\") at the end of a line means that the line must
  be joined with the next one. The following is taken as single line::

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40]

  This is different from spice syntax but it is easier to read from
  the parser.

**Parameters**

  Parameters can be ``float`` or ``int`` numbers, strings (``str``) or
  numerical vectors. Spice suffixes (``uF``, ``mA``, ``kHz``, ``GHz``,
  *etc.*) can be used to specify multipliers::

      model= mynpn v1 = 1kOhm r2 = 1e2MEG

  Some devices (such as the memristor) accept an expression as a
  parameter. Expressions must be enclosed in single quotes (') and can
  contain parenthesis and white spaces. As expressions are evaluated
  directly by the Python parser, Python syntax must be used. Constants
  must be written as numbers and standard Spice suffixes *can not* be
  used inside expressions. Mathematical functions are available but
  must be preceded by the ``np.`` prefix. These restrictions may be
  relaxed in the future::

      mem:m1 2 0 m = '1e3 * (np.cosh(1e7 * q)-1.)' 

**Element lines**

  General format::

      <element type>:<name> <node list> [<model>] <parameter list>

  Node names can be strings or numbers. A terminal named ``gnd`` (or
  ``0``) is assumed to be the global reference node for all
  circuits/subcircuits.

  <model> is optional. Parameters specified in the element line
  override parameters in model. In the following example, ``tc1`` is
  set to 1e-5::

      res:r1 1 gnd model = mymodel r=50. tc1=1e-5
      .model mymodel res (tc1=1e-4)

  Elements are documented in the :doc:`device_library`.

**Analysis lines**

  General format::

     .analysis <analysis type> <parameter list>

  Available analyses are documented in the :doc:`analysis_library`.

  Examples::

      .analysis ac start=.1GHz stop=10GHz sweep_num=200 log=True shell=0

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40] 

**Global options**

  General format (similar to spice's options):: 

      .options <parameter list>
   
  Example::
   
       .options temp=29.1439 gyr=1e-3

  Global options are documented in the :doc:`global_vars`.   
   
**Subcircuits**

  Subcircuits use a syntax similar to spice. general form for
  subcircuit definition::

    .subckt <name> <list of external nodes> 
    
    .ends

  The global reference node (``gnd`` or ``0``) can not be included as
  an external node, but if present in the subckt definition it is
  **assumed to be connected to the ``gnd`` node of other
  circuits/subcircuits**.  Example::

      res:r1 2 gnd r=40.
      x1 2 3 parasitic1
      x2 3 4 parasitic1

      .subckt parasitic1 in out
      res:r1 in out r=1kOhm
      cap:c2 out gnd c=1nH
      .ends

  Here ``gnd`` in the ``parasitic1`` definition is the same node as
  ``gnd`` in the main circuit.

**Include files**

  General format::

       .include <filename>

  The file is inserted as a part of the netlist in the position of the
  ``.include`` statement.

**Netlist variables**
	  
  Examples::

       .vars freq = 1GHz iin = .5mA
       .vars portVolt1 = [1, 2, 0.]
       idc:i1 gnd 20 idc=iin

  Numeric/vector netlist variables are defined with the ``.vars``
  keyword. Many occurences of this keyword may appear in the
  netlist. No checking is made for repeated definitions. The last
  definition overwrites any previous one.
  
  Netlist variables can be used as parameter values for element, model
  and analysis lines. ``.var`` definitions can be placed anywhere in the
  netlist.

**Output commands**

  There are two output commands: ``.plot`` and ``.save``. Both of them
  use the same syntax. Examples::

    .plot dc in out
    .plot tran 5 out3
    .plot tran vdc:amp1:i
    # In general:
    .plot <type> <list of terminals>

  In the examples, ``dc`` and ``tran`` are the type of output to
  plot. Some possible types are the following: ``dc``, ``ac_mag``,
  ``ac_phase``, ``tran``. Check the :doc:`analysis_library` to see what
  types of requests are accepted by each analysis.  

  Terminals can be external or internal. For external terminals just
  specify the terminal name.  Internal terminals are specified as
  follows::

    <element type>:<name>:<internal terminal name>
    # Example: 'x1' internal terminal from 'svbjt:q1'
    svbjt:q1:x1

  Check the internal topology of each device in the
  :doc:`device_library` to find the internal terminal names for each
  device. In the documentation external terminals are numbered,
  starting with ``0`` and internal terminals have alphanumeric
  labels. Internal reference terminals (i.e., ``tref``) are not
  accessible.  In the following example the internal terminal name is
  'i'::

        0                              i/gyr      Term: i
           o---------+            +----------------+
                     | gyr V(i)   |                |
          +         /|\          /|\              /^\ 
        vin        ( | )        ( | ) gyr vin    ( | ) gyr vdc
          -         \V/          \V/              \|/  
                     |            |                |
           o---------+            +----------------+
        1                                 |
                                         --- tref
                                          V

  Each recognized plot line generates a new figure. Results stored in
  terminals listed in a single plot line are grouped in a single
  figure. If an analysis does not recognize a request type, the
  request is ignored.

  ``.save`` statements save the requested information in a numpy
  ``.npz`` file. The file name is formed as follows by taking the main
  netlist file name minus `.net` plus ``_<request name>.npz``. For
  example, if the netlist file name is ``vsin.net``, the file created
  for an ``ac`` request is ``vsin_ac.npz``. Data saved in this file
  can be loaded in a python session using the numpy ``load`` function
  as follows::

    >>> import numpy as np
    >>> l=np.load('vsin_ac.npz')
    >>> l.files
    ['1', '2', 'xaxis']
    >>> l['1']
    array([ 1.00000000 -6.28318278e-06j,  0.99999999 -7.22413241e-06j,
            0.99999999 -8.30599536e-06j,  0.99999999 -9.54987410e-06j,
	    ...
      

**Electrothermal devices**

  Refer to the :doc:`device_library` to find which devices support
  electrothermal models. The netlist name for an electrothermal model
  is formed by adding "_t" to the original name (e.g., ``bjt_t``).  An
  electrothermal model has an additional pair of thermal
  terminals. The voltage in this thermal port is the difference
  between the device temperature and the ambient temperature. The
  current is proportional to the power dissipated in the device.



Generating this documentation
-----------------------------

The main documentation files are kept in the ``doc``
directory. Documentation can be generated in html or LaTeX formats
(other formats are possible but not tested).  The documentation can be
generated as follows::

    cd doc
    make html

The device or analysis catalogs are not checked for dependencies. To
force re-generation of those, you can just remove
``device_library.rst`` (or run ``cardoon -c`` in the doc directory)
and re-make the documentation. The ``latex`` targets can be used to
generate the documentation in latex format.
