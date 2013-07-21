"""
Catalog Generation Functions
----------------------------

.. module:: catalogs
.. moduleauthor:: Carlos Christoffersen

"""
from __future__ import print_function
import os
from globalVars import const, glVar
import paramset as ps
import analyses
import devices

def device_catalog():
    """
    Generates device library catalog

    file: device_library.rst 
    """
    with open('device_library.rst', 'w') as f:
        print(":tocdepth: 3\n", file=f)
        print("======================", file=f)
        print("Device Library Catalog", file=f)
        print("======================", file=f)
        print(" ", file=f)

        # sort devices according to category
        catDict = dict()
        for key in devices.devClass:
            dev = devices.devClass[key]
            try:
                catDict[dev.category].append(dev)
            except KeyError:
                catDict[dev.category] = [dev]
                
        # Generate one Section per category
        for category in sorted(catDict.iterkeys()):
            print(category, file=f)
            print('=' * len(category) + '\n', file=f)
            devlist = catDict[category]
            devlist.sort(key = lambda x: x.devType)
            for dev in devlist:
                key = dev.devType
                if key[-2:] == '_t':
                    #print('-' * len(key) + '\n', file=f)
                    ts = """
Electro-thermal version
+++++++++++++++++++++++

Electro-thermal version with extra thermal port: **{0}**
"""
                    print(ts.format(key), file=f)
                    continue
                # Print doc string
                doc = dev.__doc__
                if hasattr(dev,'extraDoc'):
                    doc += dev.extraDoc
                doc = doc.split(os.linesep)
                # Remove empty first line 
                doc.pop(0)
                # Add header with netlist name to title
                print(key + ': ' + doc.pop(0)[4:], file=f)
                print('-' * (len(key)+2) + doc.pop(0)[4:], file=f)
                for line in doc:
                  print(line[4:], file=f)
    
                print('\nParameters', file=f)
                print('++++++++++\n', file=f)
                # create parameter set from device dictionary
                pset = ps.ParamSet(dev.paramDict)
                print(pset.describe_parameters(), file=f)
            

def analysis_catalog():
    """
    Generates analysis library catalog

    file: analysis_library.rst 
    """
    with open('analysis_library.rst', 'w') as f:
        print(":tocdepth: 2\n", file=f)
        print("========================", file=f)
        print("Analysis Library Catalog", file=f)
        print("========================", file=f)
        print(" ", file=f)
        # loop through all analyses
        for key in sorted(analyses.anClass.iterkeys()):
            val = analyses.anClass[key]
            # Print doc string
            doc = val.__doc__.split(os.linesep)
            # Remove empty first line 
            doc.pop(0)
            # Add header with netlist name to title
            print(key + ': ' + doc.pop(0)[4:], file=f)
            print('-' * (len(key)+2) + doc.pop(0)[4:], file=f)
            for line in doc:
              print(line[4:], file=f)
            print('\nParameters', file=f)
            print('++++++++++\n', file=f)
            # create parameter set from device dictionary
            pset = ps.ParamSet(val.paramDict)
            print(pset.describe_parameters(), file=f)


def options_const_tables():
    """
    Generates tables for global variables and constants

    files: constants_table, options_table
    """
    with open('constants_table', 'w') as f:
        print(const.describe_parameters(), file=f)
    with open('options_table', 'w') as f:
        print(glVar.describe_parameters(), file=f)


def make_catalogs():
    """
    Generates all catalogs and tables
    """
    print('Generating catalogs ...')
    device_catalog()
    analysis_catalog()
    options_const_tables()

