"""
Cardoon electronic simulator library
------------------------------------

.. module:: cardoon
.. moduleauthor:: Carlos Christoffersen

The following modules and packages are available:

  simulator: high-level functions 
  circuit: circuit representation
  globalVars: global variables and constants
  paramset: parameter handling for netlist items
  netlistparser: functions to parse netlists
  catalogs: funtions to generate different catalogs
  
  devices: device models
  analyses: nodal analysis, DC, tran, sensitivities, etc.

"""

from simulator import version, release, copyright, license, run_netlist

__all__ = ['simulator', 'circuit', 'analyses', 'devices']



