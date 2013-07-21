#!/usr/bin/python
#
# simple profiling script
#
# Important: set PYTHONPATH to point to the main directory before running

import pstats
import profile
import sys
import cardoon.simulator as cs

if len(sys.argv) != 2:
    print('Usage: run -i profile_netlist.py <netlist file>')
else:
    analysisQueue = cs.parse_net(sys.argv[1])
    profile.run('cs.run_analyses(analysisQueue)','foo.out')
    p = pstats.Stats('foo.out')
    p.strip_dirs().sort_stats(-1)
    p.sort_stats('cumulative').print_stats(25)
