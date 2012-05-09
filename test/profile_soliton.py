#!/usr/bin/python

# This script is intended to be run within the cardoon shell:
#
#  cardoon -i
#  run -i profile_soliton.py

import pstats
import profile
analysisQueue = parse_net('soliton.net')
profile.run('run_analyses(analysisQueue)','foo.out')
p = pstats.Stats('foo.out')
p.strip_dirs().sort_stats(-1)
p.sort_stats('cumulative').print_stats(25)
