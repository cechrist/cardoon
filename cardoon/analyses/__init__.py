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
analysisList = ['testdev', 'op', 'dc', 'ac', 'tran']

# Add here any modules to be imported in addition to analysisList
__all__ = analysisList 

from analysis import AnalysisError

# Dictionary with existing analysis types. Key is netlist name.
anClass = {}
# This is a list of valid request types. Each analysis module must
# define a list with allowed types (reqTypes).
validReqTypes = list()

for modname in analysisList:
    module = __import__(modname, globals(), locals())
    anClass[module.aClass.anType] = module.aClass
    try: 
        validReqTypes += module.reqTypes
    except AttributeError:
        # do nothing if attribute does not exist
        pass

#---------------------------------------------------------------------
class OutRequest:
    """
    Holds a set of output variables requests

    Output request consist in: 

      1. Type of of request (``type``): dc, ac_*, tran, etc. These are
         defined by each analysis.

      2. List of variables (``varlist``): for external terminals,
         these are strings with terminal name. For internal terminals,
         a list with device and terminal names.

    After initialization the circuit adds a list of terminals in the
    ``termlist`` attribute.
    """

    def __init__(self, reqtype, varlist):
        if reqtype not in validReqTypes:
            raise CircuitError(
                'Not a valid output request type: {0}'.format(reqtype)
                + '\nValid types: {0}'.format(validReqTypes))
        self.type = reqtype
        self.varlist = varlist



