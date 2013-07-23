#!/usr/bin/python
"""
Main simulator program

For now this is quite rudimentary. 
"""
from __future__ import print_function
import sys
import warnings
import simulator
import catalogs
from netlistparser import ParseError
from circuit import CircuitError
from paramset import ParamError
from simulator import version, release, license

# Comment this out to see all warnings 
warnings.filterwarnings('ignore', category=RuntimeWarning)
# Show user warnings
warnings.filterwarnings('always', category=UserWarning)


print('\nCardoon Circuit Simulator {0} release {1}'.format(version, release))
if (len(sys.argv)) < 2 or (sys.argv[1] == '-h'):
    print('Usage:')
    print('        cardoon <netlistname>        : Process netlist file')
    print('        cardoon -c                   : Generate catalogs')
    print('        cardoon -h                   : Print this message')
    print('        cardoon -l                   : Print license\n')
    exit(1)

if sys.argv[1] == '-c':
    catalogs.make_catalogs()
elif sys.argv[1] == '-l':
    print(license)
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


