import numpy as np


r = np.eye(PSI.shape[0], PSI.shape[1])
u,s,v = np.linalg.svd(PSI, full_matrices = 1)
PSIinv = np.dot(np.transpose(u),Vec)
PSIinv = np.dot(np.transpose(r), PSIinv)
for x in range (PSI.shape[0]):
    PSIinv[x] *= s[x]
PSIinv = np.dot(np.transpose(v), PSIinv)

INV = np.dot(np.transpose(v), np.transpose(r))
for x in range (PSI.shape[0]):
    PSIinv[x] /= s[x]
INV = np.dot(INV, np.transpose(u))
