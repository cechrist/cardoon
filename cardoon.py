#!/usr/bin/python
"""
Main simulator program

For now this is quite rudimentary. 
"""
from __future__ import print_function
from globalVars import const, glVar
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
                    print('\nElectro-thermal version with extra thermal port:', 
                          key, '\n', file=f)
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
    Generates analysis library catalog: analysis_library.rst 
    """
    import analyses
    import paramset as ps

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


if __name__ == "__main__":
    import sys
    import doc
    print('\nCardoon Circuit Simulator {0} release {1}'.format(doc.version,
                                                               doc.release))
    if len(sys.argv) < 2:
        print('Usage:')
        print('        cardoon <netlistname>        : Process netlist file')
        print('        cardoon -c                   : Generate catalogs')
        print('        cardoon -x <script> <args>   : execute Python script')
        print('        cardoon -i                   : drop to Ipython shell')
        exit(1)

    if sys.argv[1] == '-c':
        # generate catalogs
        print('Generating catalogs ...')
        device_catalog()
        analysis_catalog()
        options_const_tables()
    elif sys.argv[1] == '-x':
        # Execute python script
        with open(sys.argv[2], 'r') as f:
            sys.argv = sys.argv[2:]
            exec(f)
    elif sys.argv[1] == '-i':
        # drop to ipython shell
        import devices
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


