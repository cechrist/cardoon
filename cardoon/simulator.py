"""
Simulator High-Level Functions
------------------------------

.. module:: simulator
.. moduleauthor:: Carlos Christoffersen

High-level functions for interactive use. Example::

  from cardoon import simulator as cs
  help(cs)
  cs.run_netlist('oscillator.net')


"""
from __future__ import print_function
import time
from globalVars import const, glVar
import circuit as cir
from netlistparser import parse_file
import analyses
import devices

copyright = u'2011, 2012, 2013, Carlos Christoffersen and others'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.6'
# The full version, including alpha/beta/rc tags.
release = '0.6.0.dev'

license = """
Copyright (C) {0}
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
""".format(copyright)


def reset_all():
    """
    Clean all circuits and reset global options
    """
    # Set global variables to default values
    glVar.reset()
    glVar.set_attributes()
    # Erase circuits
    cir.reset_allckt()


def parse_net(filename, ckt = None):
    """
    Parse netlist file given in filename. 

    If ckt is not given parses the circuit named 'main'. Returns
    analysisQueue with analyses to be performed.
    """
    if not ckt:
        ckt = cir.get_mainckt()
    return parse_file(filename, ckt)


def run_analyses(analysisQueue, ckt = None):
    """
    Run all analyses in analysisQueue applied to the provided circuit
    (or 'main' if no circuit provided). If queue is empty just print
    circuit.
    """
    if not ckt:
        ckt = cir.get_mainckt()
    # Perform requested analyses
    if analysisQueue:
        for an in analysisQueue:
            try:
                start = time.clock()
                an.run(ckt)
                elapsed = time.clock() - start
                print('{0} analysis time: {1} s\n'.format(an.anType, elapsed))
            except analyses.AnalysisError as ae:
                print(ae)
    else:
        print('Nothing to do. Printing circuit:\n')
        print(ckt.netlist_string())
        print(ckt.globals_to_str())


def new_elem(elemType, name, **kwargs):
    """
    Returns a new element instance

    elemType: type of device (diode, res, cap, ind, etc.)
    name: instance name (without type)

    The remaining arguments should have the following format:
    <param name> = <param value>. Type of <param value> should
    exactly match the expected parameter type.
    
    Sample usage: newElem('diode', 'd4', isat=2.1e-15, cj0=1e-12)

    """
    elem = devices.devClass[elemType](name)
    elem.set_params(**kwargs)
    return elem


def run_netlist(fileName, circuitName=None, reset=True):
    """
    Read circuit and run analyses in netlist given in fileName 

    circuitName: name of circuit instance to be created, if None, use
    main circuit

    If reset == True, clear simulator state before reading netlist,
    otherwise append to any existing data.
    """
    if reset:
        reset_all()

    if circuitName:
        ckt = cir.Circuit(circuitName)
    else:
        ckt = cir.get_mainckt()        

    analysisQueue = parse_net(fileName, ckt)
    run_analyses(analysisQueue, ckt)



