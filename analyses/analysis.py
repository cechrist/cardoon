"""
This module defines some handy classes/functions for analyses

----------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

from globalVars import glVar
from IPython.Shell import IPShellEmbed

class AnalysisError(Exception):
    """
    Exception class to be used for analyses. 
    """
    pass


def ipython_drop(glo, loc):
    """
    Add this function in your code to drop to an ipython shell

    (provided that glVar.shell == True)

    glo: globals()
    loc: locals()
    """
    if glVar.shell: 
        args = ['-pi1','In <\\#>: ','-pi2','   .\\D.: ',
                '-po','Out<\\#>: ','-nosep']
        ipshell = IPShellEmbed(
            args, 
            banner = 'Dropping into IPython, type CTR-D to exit',
            exit_msg = 'Leaving Interpreter, back to program.')
        ipshell(global_ns = glo, local_ns = loc)

