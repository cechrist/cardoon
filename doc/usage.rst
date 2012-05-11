

Installation and Usage
======================

Installation
------------

Cardoon is being developed in Debian GNU/Linux
<http://www.debian.org>, but it should work in other environments if
the Python interpreter and the libraries are available.  In addition
to a Python interpreter (version 2.6 or possibly 2.7), several
free libraries are required to run cardoon. In Debian Linux or
Debian-based distributions (Ubuntu, Linux Mint and others) all of
these libraries except pycppad can be installed as packages::

  apt-get install python-numpy python-sparse python-matplotlib \
  python-pyparsing ipython python-sphinx pylint 

To install pycppad install cppad first by following the instructions
on their websites. The corresponding websites for each library are
listed below:

* pycppad: Automatic differentiation

  - cppad:  http://www.coin-or.org/CppAD/Doc/cppad.xml

  - pycppad:  http://www.seanet.com/~bradbell/pycppad/pycppad.xml 

* numpy:  http://numpy.scipy.org/ (matrix and vector support)

* pysparse:  http://pysparse.sourceforge.net/ (sparse matrix support)

* matplotlib:  http://matplotlib.sourceforge.net/ (plotting)

* pyparsing:  http://pyparsing.wikispaces.com/ (parser)

* ipython:  http://ipython.org/ (iterative shells)

* Sphinx: http://sphinx.pocoo.org/ (documentation)

* pyreverse from the pylint package: http://www.logilab.org/2560 (to
  generate UML diagrams)

You may also need to install *git* to fetch the source code from
the github repository::

    git clone git://github.com/cechrist/cardoon.git

Alternatively the source code can be downloaded as a zip file from:

https://github.com/cechrist/cardoon

Just unpack the zip file in some directory.

Usage
-----

At this moment the main program reads a netlist file, builds the
circuit described there and runs any specified analyses. Optionally
some analyses can drop to an ipython shell with access to internal
variables after calculations are finished for interactive work.

After the source is installed, change into the cardoon directory and
run the program with no arguments for a brief help message::

    cechrist@moon:~/wd/cardoon$ python cardoon.py 
    Usage:
            cardoon <netlistname>  : Process netlist file
            cardoon -c             : Generate catalogs
            cardoon -i             : immediatly drop to Ipython shell

Then you can change into the test/devices directory and try some of
the netlists there::

    cechrist@moon:~/wd/cardoon$ cd test
    cechrist@moon:~/wd/cardoon/test$ python ../cardoon.py simple.net
    Trying Simple Newton's method
    ******************************************************
                 Operating point analysis
    ******************************************************
    Number of iterations =  3
    Residual =  2.20086083382e-08
    
     Node      |  Value               | Unit 
    ----------------------------------------
    1          |             0.672162 | V
    2          |              2.67216 | V
    gnd        |                  0.0 | V
    
    Element:  svdiode:d1
     Variable  |  Value 
    -------------------------
        Cd     | 6.69110617079e-05
        ID     | 0.00193278379487
      Sshot    | 6.19332180271e-22
     Sthermal  | 0.0
        VD     | 0.672162035866
        gd     | 0.0747260598229
     kFliker   | 0.0
        x      | 4.38235965115


Netlist Format
--------------

A very brief description is provided here. The netlist syntax
resembles somewhat the syntax used in other simulators such as spice,
fREEDA and QUCS, but at least for now it has some simplifications. The
netlist is case-sensitive. Each line specifies one circuit element, an
analysis to perform or another command.

* The backslash ("\\") at the end of a line means that the line must
  be joined with the next one. The following is taken as single line::

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40]

  This is different from spice syntax but it is easier to read from
  the parser.

* Parameters can be ``float`` or ``int`` numbers, strings (``str``) or
  numerical vectors. All spice suffixes can be used to specify
  multipliers::

      model= mynpn v1 = 1kOhm r2 = 1e2MEG

* Element lines::

      <element type>:<name> <node list> [<model>] <parameter list>

  <model> is optional. Parameters specified in the element line
  override parameters in model. In the following example, ``tc1`` is
  set to 1e-5::

      res:r1 1 gnd model = mymodel r=50. tc1=1e-5
      .model mymodel res (tc1=1e-4)

  Elements are documented in the :doc:`device_library`.

* Analysis lines::

     .analysis <analysis type> <parameter list>

  Available analyses are documented in the :doc:`analysis_library`.

  Examples::

      .analysis ac start=.1GHz stop=10GHz sweep_num=200 log=True shell=0

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40] 

* Global options (similar to spice's options):: 

      .options <parameter list>
   
  Example::
   
       .options temp=29.1439 gyr=1e-3

  Global options are documented in the :doc:`global_vars`.   
   
* Subcircuits use a syntax similar to spice::

      x1 2 3 4 X1
      x2 2 gnd 3 X1

      .subckt X1 in out gnd
      res:r1 in out r=1kOhm
      cap:c2 out gnd c=1nH
      .ends

* Include files::

       .include <filename>


* Netlist variables::

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

* Output commands: the only command to control output defined so far
  is ``.plot``. Examples::

    .plot dc in out
    .plot tran 5 out3
    # In general:
    .plot <type> <list of terminals>

  In the examples, ``dc`` and ``tran`` are the type of output to
  plot. Some possible types are the following: ``dc``, ``ac_mag``,
  ``ac_phase``, ``tran``. Check the :doc:`analysis_library` to see what
  types of requests are accepted by each analysis.  Each recognized
  plot line generates a new figure. Results stored in terminals listed
  in a single plot line are grouped in a single figure. If an analysis
  does not recognize a request type, the request is ignored.


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
