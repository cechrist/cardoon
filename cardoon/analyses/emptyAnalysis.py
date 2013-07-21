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

"""

import numpy as np
import circuit as cir
from cardoon.globalVars import glVar
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
    paramDict = dict(
        intvars = ('Print internal element nodal variables', '', bool, False),
        elemop = ('Print element operating points', '', bool, False)
        )

    def __init__(self):
        # Just init the base class
        ParamSet.__init__(self, self.paramDict)

    def run(self, circuit):
        """
        Calling this functions performs the analysis on the given
        circuit. Throw AnalysisError if something goes wrong.
        """
        # Parameters are available in self.pset.attribute
        pass




# Here you can add additional functions and classes that only are
# visible withing this module.

