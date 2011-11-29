

Installation and Usage
======================

Installation
------------

In addition to a Python interpreter (version 2.6 or possibly 2.7), the
following free libraries are required to run cardoon. All, except
pycppad can be installed directly from Debian/Ubuntu packages. To
install pycppad, just follow the instructions on its website.

* numpy:  http://numpy.scipy.org/

* pycppad:  http://www.seanet.com/~bradbell/pycppad/pycppad.xml

* matplotlib:  http://matplotlib.sourceforge.net/

* pyparsing:  http://pyparsing.wikispaces.com/

* ipython:  http://ipython.org/

You may also need to install *git* to fetch the source code from
the github repository::

    git clone git://github.com/cechrist/cardoon.git

Alternatively the source code can be downloaded as a zip file from:

https://github.com/cechrist/cardoon

Just unpack the zip file in some directory.

Usage
-----

At this moment the main program is very simple. It just reads a
netlist file, builds the circuit described there and runs any
specified analyses.

After the source is installed, change into the cardoon directory and
run the program with no arguments for a brief help message::

    cechrist@moon:~/wd/cardoon$ python cardoon.py 
    Usage:
            cardoon <netlistname>  : Process netlist file
            cardoon -c             : Generate catalogs
            cardoon -i             : drop to Ipython shell

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

A very brief and incomplete description is provided here. The netlist
syntax resembles somewhat the syntax used in other simulators such as
spice, fREEDA and QUCS, but at least for now it has some
simplifications. The netlist is case-sensitive. Each line specifies
one circuit element, an analysis to perform or another command.

#. The backslash ("\\") at the end of a line means that the line must
   be joined with the next one. The following is taken as single
   line::

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40]

   This is different from spice syntax but it is easier to read from
   the parser.

#. Parameters can be ``float`` or ``int`` numbers, strings (``str``)
   or numerical vectors. All spice suffixes can be used to specify
   multipliers::

      model= mynpn v1 = 1kOhm r2 = 1e2MEG

#. Element lines::

      <element type>:<name> <node list> [<model>] <parameter list>

   <model> is optional. Parameters specified in the element line
   override parameters in model. In this example, ``tc1`` is set to
   1e-5::

      res:r1 1 gnd model = mymodel r=50. tc1=1e-5
      .model mymodel res (tc1=1e-4)

#. Analysis lines::

     .analysis <analysis type> <parameter list>

  Example::

      .analysis testdev plot=1 ports_bias = [.7V] sweep_port=0 \
      start = .1V stop= .8V sweep_num=1100 device = diode:d2 \
      param = temp param_val = [0., 27, 40]

#. Global variables:: 

      .options <parameter list>
   
   Example::
   
       .options temp=29.1439 gyr=1e-3
   
   List of global variables (check globalVars.py for an updated list)

 =========== ============ ============ ===================================================== 
 Name         Default      Unit         Description                                          
 =========== ============ ============ ===================================================== 
 abstol       1.0e-08                   Absolute tolerance                                   
 gyr          0.01         S            Default gain in internal gyrators                    
 maxiter      100                       Maximum number of Newton iterations                  
 reltol       1.0e-08                   Relative tolerance                                   
 shell        0                         Drop to ipython shell after calculation              
 temp         27.0         C            Ambient temperature                                  
 =========== ============ ============ =====================================================  

#. Subcircuits use a syntax similar to spice::

      x1 2 3 4 X1
      x2 2 gnd 3 X1

      .subckt X1 in out gnd
      res:r1 in out r=1kOhm
      cap:c2 out gnd c=1nH
      .ends

#. Include files::

       .include <filename>

For now there are no output commands defined.


Generating this documentation
-----------------------------

The main documentation files are kept in the ``doc`` directory. To
generate the documentation in html or LaTeX formats (other formats are
possible but not tested) the following packages are needed:

* Sphinx (http://sphinx.pocoo.org/)

* pyreverse from the pylint package yo generate UML diagrams
  (http://www.logilab.org/2560)

The documentation can be generated as follows::

    cd doc
    make html

The device or analysis catalogs are not checked for dependencies. To
force re-generation of those, you can just remove
``device_library.rst`` and re-make the documentation. The ``latex``
targets can be used to generate the documentation in latex format.
