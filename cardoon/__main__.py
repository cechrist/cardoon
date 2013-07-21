#!/usr/bin/python
"""
Main simulator program

For now this is quite rudimentary. 
"""
from __future__ import print_function
import sys
import simulator
import catalogs
from netlistparser import ParseError
from circuit import CircuitError
from paramset import ParamError
from simulator import version, release

print('\nCardoon Circuit Simulator {0} release {1}'.format(version, release))
if len(sys.argv) < 2:
    print('Usage:')
    print('        cardoon <netlistname>        : Process netlist file')
    print('        cardoon -c                   : Generate catalogs')
    print('        cardoon -x <script> <args>   : execute Python script')
    exit(1)

if sys.argv[1] == '-c':
    catalogs.make_catalogs()
elif sys.argv[1] == '-x':
    # Execute python script
    with open(sys.argv[2], 'r') as f:
        sys.argv = sys.argv[2:]
        exec(f)
else:
    # Use 'main' circuit
    try:
        # No need to reset anything since we are just starting
        simulator.run_netlist(sys.argv[1], reset=False)
        print('Press [Enter] to exit ...')
        c = raw_input()
    except (ParseError, CircuitError, ParamError) as ex:
        print(ex)
        exit(1)


