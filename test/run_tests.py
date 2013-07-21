#!/usr/bin/python
#
# Use the script after program modifications to check for regressions

import sys
import os
import numpy as np

# Make sure we import the right code 
pwd = os.getcwd()
cardoonPath = pwd[:pwd.rfind('test')]
sys.path.insert(0, cardoonPath)
import cardoon.simulator as cs
import cardoon.circuit as cir

flag = False
# Check arguments
if len(sys.argv) >= 2:
    if sys.argv[1] == '-g':
        if len(sys.argv) == 2:
            print('\nUsage: run_tests.py [-g] <netlists>')
        else:
            # Re-generate output data
            netlists = sys.argv[2:]
            flag = True
            print('\n=======================================================')
            print(' Re-generating reference results.')
            print(' *Existing results are being erased*')
            print('\n=======================================================')
    else:
        if sys.argv[1] == '*.net':
            # argument was not expanded (we are probably on Windows)
            files = os.listdir(os.curdir)
            netlists = [f for f in files if f[-4:]=='.net']
        else:
            netlists = sys.argv[1:]
else:
    print('\nNo netlist specified, exiting ...')
    exit(-1)

for net in netlists:
    circuit = cir.get_mainckt()
    analysisQueue = cs.parse_net(net, circuit)
    cs.run_analyses(analysisQueue)
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
        warningflag = False
        # Look for *_ref.npz files and compare with corrensponding *.npz
        files = [name for name in os.listdir('.') 
                 if (name.find(basename) == 0) and (name.find('_ref.npz') > 0)]
        residual = 0.
        if not files:
            print('\nFatal: no reference files for {0}'.format(net))
            print('\n*** Test Failed! ***')
            exit(-1)
        for reffile in files:
            delta = 0.
            # Look for corresponding output file
            outfile = reffile.replace('_ref','')
            try:
                refResult = np.load(reffile)
            except IOError:
                print('\nProblem reading reference file: {0}'.format(reffile))
                print('\n*** Test Failed! ***')
                exit(-1)
            try:
                result = np.load(outfile)
            except IOError:
                print('\nProblem reading output file: {0}'.format(outfile))
                print('\n*** Test Failed! ***')
                exit(-1)

            for var in refResult.files:
                try:
                    delta = refResult[var] - result[var]
                except KeyError:
                    print('\nError: "{0}" not found in {1}'.format(
                            var, outfile))
                    print('\n*** Test Failed! ***')
                    exit(-1)
                res = np.average(abs(delta))
                residual += res
                if not np.all(abs(delta) < (abs(cs.glVar.reltol * 
                                                refResult[var])
                                            + cs.glVar.abstol)):
                    message = """
Reference file: {0}
Variable: {1}
Residual: {2}""".format(reffile, var, res)
                    if res < 1e-3:
                        print('\n=============================================')
                        message = "Warning: \n" + message
                        print(message)
                        print('=============================================')
                        warningflag = True
                    else:
                        message = "\nTest failed: \n" + message
                        raise Exception(message)
        print('\n=======================================================')
        print(' Success: {0}'.format(net))
        print(' Sum of residuals: {0}'.format(residual))
        print('=======================================================\n')

    # Reset everything
    cs.reset_all()

if not flag:
    if warningflag:
        print('\n============================================================')
        print('              All tests completed (with Warnings)')
        print('============================================================\n')
    else:
        print('\n============================================================')
        print('              All tests succesfully completed!')
        print('============================================================\n')
