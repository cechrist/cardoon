#!/usr/bin/python
"""
Main simulator program
"""

import circuit as cir
from netlistparser import parse_file, ParseError, analysisQueue
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
    #ckt.flatten()
    ckt.init()
    # The following is slow, so we may want to remove/change it later
    #ckt.check_sanity()
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
                print ae
    else:
        print 'Nothing to do. Printing circuit:\n'
        print ckt
        print ckt.models_to_str()

def device_catalog():
    """
    Generates device library catalog: device_library.rst (to stdout)
    """
    import devices
    import paramset as ps
    print "======================"
    print "Device Library Catalog"
    print "======================"
    print " "
    # loop through all devices
    for key in sorted(devices.devClass.iterkeys()):
        print key
        print '-' * len(key) + '\n'
        if key[-2:] == '_t':
            print 'Electro-thermal version of', key[:-2], '(extra thermal port)\n'
            continue
        val = devices.devClass[key]
        # Print doc string
        doc = val.__doc__.split(os.linesep)
        for line in doc:
          print line[4:]
        print '\nParameters'
        print '++++++++++\n'
        # create parameter set from device dictionary
        pset = ps.ParamSet(val.paramDict)
        print pset.describe_parameters()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print '\nUsage:'
        print '        cardoon <netlistname>'
        print '        cardoon -c'
        print '        cardoon -i'
        exit(1)

    if sys.argv[1] == '-c':
        # generate device library catalog
        device_catalog()
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
        except (ParseError, cir.CircuitError) as ex:
            print ex
            exit(1)
        else:
            run_analyses(analysisQueue)

