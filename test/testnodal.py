"""
Test nodal module: to be run inside cardoon shell (cardoon -i)

Run as follows: %run -i testnodal 
"""
import numpy as np
import analyses.nodal as nd


parse_net('bias_npn.net')
ckt=cir.get_mainckt()
nodalckt = nd.NodalCircuit(ckt)
x = nodalckt.get_guess() + 1e-3
Jac = np.zeros((nodalckt.dimension, nodalckt.dimension))
iVec = np.zeros(nodalckt.dimension)
s = np.zeros(nodalckt.dimension)
nodalckt.get_DC_source(s)

for k in range(3):
    nodalckt.get_DC_i_Jac(x, iVec, Jac)
    try:
        deltax = np.linalg.solve(Jac, iVec - s)
    except:
        print 'oops'
        deltax = np.dot(np.linalg.pinv(Jac), iVec-s)
    x -= deltax
    print x


