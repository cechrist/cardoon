"""
Test nodal module: to be run inside cardoon shell (cardoon -i)

Run as follows: %run -i testnodal 
"""
import numpy as np
import analyses.nodal as nd

netlist = raw_input('netlist: ')
parse_net(netlist)
ckt=cir.get_mainckt()
ckt.init()
nd.make_nodal_circuit(ckt)
dc = nd.DCNodal(ckt)
x = dc.get_guess() + 1e-3

x = dc.get_guess()
s = dc.get_source()

for k in range(10):
    (iVec, Jac) = dc.get_i_Jac(x)
    try:
        deltax = np.linalg.solve(Jac, iVec - s)
    except:
        print 'oops'
        deltax = np.dot(np.linalg.pinv(Jac), iVec-s)
    x -= deltax
    print x


