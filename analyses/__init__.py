"""
This package contains all analyses modules.  Analysis classes are
imported into a dictionary (similarly as with the 'devices'
package). Keys are the device types. So to create a new device just
use:

analyses.anClass['analysisType']

----------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

# Regular 'analysis' modules listed here
netElemList = ['testdev']

# Add here any modules to be imported in addition to netElemList
__all__ = netElemList 

from analysis import AnalysisError

anClass = {}

for modname in netElemList:
    module = __import__(modname, globals(), locals())
    anClass[module.Analysis.anType] = module.Analysis




