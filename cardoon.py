#!/usr/bin/python
"""
Main simulator program

For now this is quite rudimentary. 
"""
from __future__ import print_function
import circuit as cir
from netlistparser import parse_file, ParseError, analysisQueue
from paramset import ParamError
import analyses
import os

def parse_net(filename, ckt = None):
    """
    Parse netlist file given in filename. If ckt is not given it
    parses in the circuit named 'main'. Returns analysisQueue with
    analyses to be performed.
    """
    if not ckt:
        ckt = cir.get_mainckt()
    parse_file(filename, ckt)
    return analysisQueue


def run_analyses(analysisQueue, ckt = None):
    """
    Run all analyses in analysisQueue applied to the provided circuit
    (or 'main' if no circuit provided). If queue is empty just print
    circuit.
    """
    if not ckt:
        ckt = cir.get_mainckt()
    # Perform requested analyses
    if analysisQueue:
        for an in analysisQueue:
            try:
                an.run(ckt)
            except analyses.AnalysisError as ae:
                print(ae)
    else:
        print('Nothing to do. Printing circuit:\n')
        print(ckt.netlist_string())
        print(ckt.globals_to_str())

def device_catalog():
    """
    Generates device library catalog: device_library.rst 
    """
    import devices
    import paramset as ps

    with open('device_library.rst', 'w') as f:
        print("======================", file=f)
        print("Device Library Catalog", file=f)
        print("======================", file=f)
        print(" ", file=f)
        # loop through all devices
        for key in sorted(devices.devClass.iterkeys()):
            if key[-2:] == '_t':
                #print('-' * len(key) + '\n', file=f)
                print('\nElectro-thermal version with extra thermal port:', 
                      key, '\n', file=f)
                continue
            val = devices.devClass[key]
            # Print doc string
            doc = val.__doc__
            if hasattr(val,'extraDoc'):
                doc += val.extraDoc
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
            pset = ps.ParamSet(val.paramDict)
            print(pset.describe_parameters(), file=f)
            

def analysis_catalog():
    """
    Generates analysis library catalog: analysis_library.rst 
    """
    import analyses
    import paramset as ps

    with open('analysis_library.rst', 'w') as f:
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


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print('\nUsage:')
        print('        cardoon <netlistname>  : Process netlist file')
        print('        cardoon -c             : Generate catalogs')
        print('        cardoon -i             : drop to Ipython shell')
        exit(1)

    if sys.argv[1] == '-c':
        # generate catalogs
        device_catalog()
        analysis_catalog()
    elif sys.argv[1] == '-i':
        # drop to ipython shell
        import devices
        from globalVars import const, glVar
        from IPython.Shell import IPShellEmbed
        args = ['-pi1','In <\\#>: ','-pi2','   .\\D.: ',
                '-po','Out<\\#>: ','-nosep']
        ipshell = IPShellEmbed(args, 
                               banner = 'Type CTR-D to exit',
                               exit_msg = 'Leaving Interpreter.')
        ipshell()
    else:
        # Use 'main' circuit
        try:
            analysisQueue = parse_net(sys.argv[1])
            run_analyses(analysisQueue)
        except (ParseError, cir.CircuitError, ParamError) as ex:
            print(ex)
            exit(1)


