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
import numpy as np
import matplotlib.pyplot as plt

class AnalysisError(Exception):
    """
    Exception class to be used for analyses. 
    """
    pass


def ipython_drop(msg, glo, loc):
    """
    Add this function in your code to drop to an ipython shell

    (provided that glVar.shell == True)

    msg: informative message to print when shell is invoked
    glo: globals()
    loc: locals()
    """
    args = ['-pi1','In <\\#>: ','-pi2','   .\\D.: ',
            '-po','Out<\\#>: ','-nosep']
    ipshell = IPShellEmbed(
        args, 
        banner = 'Dropping into IPython, type CTR-D to exit' + msg,
        exit_msg = 'Leaving Interpreter, back to program.')
    ipshell(global_ns = glo, local_ns = loc)


def process_requests(circuit, reqtype, xaxis, xlabel, attribute, 
                     f1 = lambda x:x, unitOverride = '', log = False):
    """
    Process plot and save output requests

    Arguments:

    circuit: circuit instance with requests
    reqtype: tran, dc, etc.
    xaxis: numpy array to use for x axis
    xlabel: Label for x-axis
    attribute: attribute to extract from terminals (numpy array)
    f1: function to apply to attribute
    unitOverride: override terminal units
    log: use logarithmic scale
    """
    if log:
        pltfunc = plt.semilogx
    else:
        pltfunc = plt.plot

    flag = False
    for outreq in circuit.plotReqList:
        if outreq.type == reqtype:
            flag = True
            plt.figure()
            plt.grid(True)
            plt.xlabel(xlabel)
            for term in outreq.termlist:
                if unitOverride:
                    unit = unitOverride
                else:
                    unit = term.unit
                pltfunc(xaxis, f1(getattr(term, attribute)), 
                        label = 'Term: {0} [{1}]'.format(term.nodeName, unit)) 
            if len(outreq.varlist) == 1:
                plt.ylabel(
                    'Term: {0} [{1}]'.format(term.nodeName, unit))
            else:
                plt.legend()
    if flag:
        plt.show()
    flag = False
    if circuit.saveReqList:
        saveDict = {'xaxis': xaxis}
        for outreq in circuit.saveReqList:
            if outreq.type == reqtype:
                flag = True
                for term in outreq.termlist:
                    saveDict[term.nodeName] = f1(getattr(term, attribute))
    if flag:
        # Take filename minus '.net' and append request name
        outfile = circuit.filename.split('.net')[0] + '_' + reqtype
        # One or more variables should be saved
        np.savez(outfile, **saveDict)
