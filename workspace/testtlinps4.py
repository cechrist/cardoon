"""
Test tlinpy4 model
"""
import numpy as np

parse_net('test.net')
ckt=cir.get_mainckt()
t1 = ckt.elemDict['tlinps4:tl2']
y1 = t1.get_Y_matrix(1e9)
print y1
fvec=np.array([1e3, 2e5, 3e7, 2e10])
y2 = t1.get_Y_matrix(fvec)
print y2
print t1.get_G_matrix()
v=np.array([1., 1.001, 0.5, 0.5])
d=t1.get_OP(v)
t1.print_vars()
