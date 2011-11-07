#!/usr/bin/python
"""
Main simulator program
"""

import circuit as cir
from netlistparser import parse_file, ParseError, analysisQueue
import analyses

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


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print "Usage: cardoon <netlistname> \n"
        exit(1)
    # Use 'main' circuit
    try:
        analysisQueue = parse_net(sys.argv[1])
    except (ParseError, cir.CircuitError) as ex:
        print ex
        exit(1)
    else:
        run_analyses(analysisQueue)

