#!/usr/bin/python

# This script can be run within cardoon environment:
#
#  cardoon -x run_tests.py [-g]
#

import sys
import os
import numpy as np

netlists = ['sum741_profile.net', 'soliton.net']
usage = '\nUsage: cardoon -x run_tests.py [-g]'
flag = False

# Check arguments
if len(sys.argv) == 2:
    if sys.argv[1] != '-g':
        print(usage)
    else:
        # Re-generate output data
        flag = True
        print('\n=======================================================')
        print(' Re-generating reference results.')
        print(' Existing results *will be erased*')
        print('\n=======================================================')
elif len(sys.argv) > 1:
    print(usage)



for net in netlists:
    circuit = cir.get_mainckt()
    analysisQueue = parse_net(net, circuit)
    run_analyses(analysisQueue)
    basename = net.split('.net')[0]
    if flag:
        # Move output data to reference file names
        if circuit.saveReqList:
            for outreq in circuit.saveReqList:
                oldname = basename + '_' + outreq.type + '.npz'
                newname = basename + '_' + outreq.type + '_ref.npz'
                try:
                    os.rename(oldname, newname)
                except OSError:
                    pass
    else:
        # Look for *_ref.npz files and compare with corrensponding *.npz
        files = [name for name in os.listdir('.') 
                 if (name.find(basename) == 0) and (name.find('_ref.npz') > 0)]
        residual = 0.
        for reffile in files:
            delta = 0.
            # Look for corresponding output file
            outfile = reffile.replace('_ref','')
            refResult = np.load(reffile)
            result = np.load(outfile)
            for var in refResult.files:
                delta = refResult[var] - result[var]
                res = np.max(abs(delta))
                residual += res
                if not np.all(abs(delta) < (abs(glVar.reltol * refResult[var])
                                            + glVar.abstol)):
                    message = """
Test failed:

Reference file: {0}, Variable: {1}, Residual: {2}

""".format(reffile, var, res)
                    raise Exception(message)
        print('\n=======================================================')
        print(' Success: {0}'.format(net))
        print(' Sum of residuals: {0}'.format(res))
        print('=======================================================\n')

    # Reset everything
    # Set global variables to default values
    glVar.reset()
    glVar.set_attributes()
    # Erase circuits
    cir.reset_allckt()

if not flag:
    print('\n============================================================')
    print('              All tests succesfully completed!')
    print('============================================================\n')
