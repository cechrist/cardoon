"""
:mod:`csUtil` -- Utility functions for compressed sampling
----------------------------------------------------------

This module contains utility functions for compressed sampling:

* Compression matrix creation

* Recover signal from compressed coefficients

* Get the recovery Jacobian for a given set of compressed coefficients

"""
import numpy as np
import spams

# Compression matrix creation

def get_C(m, nsamples):
    """
    Returns a random m x nsamples compression matrix
    """
    C = np.random.normal(0, np.sqrt(nsamples), size = (m, nsamples))
    col_sums = np.sqrt((C*C).sum(axis=0))
    C = np.asfortranarray(C / col_sums[np.newaxis, :])
    return C

# Recover a signal
def recover(C, xc):
    (m, nsamples) = C.shape
    if len(xc.shape) == 1:
        xComp = np.asfortranarray(xc.reshape(m,1))
        nsig = 1
    else:
        nsig = xc.shape[1]
        xComp = np.asfortranarray(xc)

#    # Use Orthogonal Matching Pursuit algorithm
#    L = m
#    eps = 1e-13
#    lambda1 = 0.
#    numThreads = -1
#    xd = spams.omp(xComp, C, L=L, eps= eps, lambda1 = lambda1, 
#                   return_reg_path = False, 
#                   numThreads = numThreads).todense().getA()
    
    # Homotopy-LARS algorithm 
    lambda1 = 1e-10
    lambda2 = 1e-6
    numThreads = -1
    mode = 1
    xd = spams.lasso(xComp, D=C, 
                     lambda1 = lambda1, lambda2 = lambda2,
                     mode = mode,
                     numThreads = numThreads).todense().getA()

    return xd

# Decompression Jacobian calculation
def get_Jc(C, xc, delta = 1e-4):
    (m, nsamples) = C.shape
    incr = np.zeros(m)
    for i in range(0,m):
        if xc[i]:
            incr[i] = delta*xc[i]
        else:
            incr[i] = delta
    M1 = np.zeros((m, m+1))
    M1[:,0] = xc
    M1[:,1:] = np.tile(xc,(m,1)).T + np.diag(incr)
    xComp = np.asfortranarray(M1)
    
    Xrec = recover(C, xComp)
    J_c = Xrec[:,1:] - np.tile(Xrec[:,0],(m,1)).T
    J_c /= incr
    return J_c
