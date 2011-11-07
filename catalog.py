#!/usr/bin/python
"""
Generates device library catalog: device_library.rst (to stdout)
"""

#import circuit as cir
#from netlistparser import parse_file, ParseError, analysisQueue
#import analyses
import devices
import paramset as ps

if __name__ == "__main__":

    print "=============="
    print "Device Library"
    print "=============="
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
        print val.__doc__
        print '\nParameters'
        print '++++++++++\n'
        # create parameter set from device dictionary
        pset = ps.ParamSet(val.paramDict)
        print pset.describe_parameters()

