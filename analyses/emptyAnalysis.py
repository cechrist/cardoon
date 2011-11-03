"""
Generic analysis module code to use as a model to create new ones

-------------------------------------------------------------------

Mandatory attributes
====================

Argument: anType = 'string' 

Specifies the netlist name of the analysis, for example 'tran'
transient analysis.

Mandatory Functions
===================

def run(self, circuit): self explanatory

-------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
import circuit as cir
from globalVars import glVar
from analysis import AnalysisError

class Analysis:
    """
    Class name must be "Analysis"

    Document here what the analysis does and explain parameters, if
    needed.
    """

    # antype is the netlist name of the analysis: .analysis tran tstart=0 ...
    anType = "empty_analysis"

    # Define parameters as follows
    paramList = (
        ('tstart', 'Initial time', 's', 0.),
        ('tstop', 'Stop time', 's', 0.),
        ('tstep', 'Initial time step', 's', 1e-3)
        )

    def run(self, circuit):
        """
        Calling this functions performs the analysis on the given
        circuit. Throw AnalysisError if something goes wrong.
        """
        # Parameters are available in self.pset.attribute
        pass




# Here you can add additional functions and classes that only are
# visible withing this module.

