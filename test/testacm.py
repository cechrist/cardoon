# Script to test ACM model derivatives. Run with cardoon -i

import numpy as np
import pycppad
from devices import devClass
from devices.mosACM import inv_f, inv_f1
m1=devClass['mosacm']('test')
m1.set_attributes()
m1.process_params()
m1.set_temp_vars(m1.temp)
vport = np.array([3.,2.,1.])
print m1.eval_cqs(vport)
print m1.eval(vport)
print m1.eval_and_deriv(vport)
fout = m1._func.forward(0, vport)
jac = m1._func.jacobian(vport)
