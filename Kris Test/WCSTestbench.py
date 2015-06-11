# Program for building wavelets
# Working with imported data to have a look at support structures.

# First start by importing needed libraries
import numpy as np
import time
import scipy.sparse as sparse
import scipy.sparse.linalg as spy
# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

# Works with powers of 2.  i.e.  2, 4, 8, 16, 32, 64, ...
numcoeff = 1024

# Convergence at difference
convergence = 1e-8

# Number of steps in test waveform (should be higher than numcoeff)
steps = 2*numcoeff

# Percentage of highest wavelet coefficients to keep (sets rest to 0) 0.3 = 30%
rejectpercent = .1

# Percentage of compressed coefficients to corrupt by adding noise  0.2 = 20%
corrupt = 0.0

# Mean and Standard Deviation of the noise
nmean = 0.
nsd = 100

# Use a conditioned random matrix or force isometry A* A = I (Noisy but works
# can be improved with random matrix "decoupling" or "whitening")
ForceIsometry = False

# Search Speed:  Normally, each iteration of the recovery algorithm
# searches for and adds one new column to the support set (SearchSpeed = 1),
# increasing the search speed (SearchSpeed > 1) will increase the number of
# columns added to the support set per iteration.  However, this will introduce
# error into the result if increased too much.  Noise level affects the maximum
# speed you can use and still get an accurate recovery.
SearchSpeed = 1

# Use SVD Decomposition instead of matrix inverse in recovery algorithm
UseSVDDecomp = True

# Generate separate matrices for sensing and compression?
UseDifferentSensingMatrix = False

# Load iteration data from simulator
xMatrix = np.loadtxt('outputs.txt', delimiter=',')

# ------------------------------------------------------------
# End of Inputs
# ------------------------------------------------------------

ttstart = time.time()
# For now we will only use Haar wavelets.

# ------------------------------------------------------------
# Start with the Haar Basis function.
# ------------------------------------------------------------
# For Haar, the basis function is a box function that encompasses the
# entire period we are measuring over.
Basis = np.ones(numcoeff)

# Next we need the wavelet function.

# This function is derived mathematically by performing the scaling calculation
# with the original box function.
print('Building Wavelet Matrix...')
Wavelet = np.hstack((np.ones(numcoeff/2), (-1)*np.ones(numcoeff/2)))


# Now we need the wavelet transform matrix

# The first row of the matrix will be made up of the Basis functions elements
# and the second row will be the mother wavelet functions elements.

WM = np.empty((numcoeff, numcoeff))
for x in range(numcoeff):
    WM[0,x] = Basis[x]
    WM[1,x] = Wavelet[x]

# The rest of the rows are made up of the Wavelet functions translations.

Scale = 1;
SF = 1;

for x in range(2,numcoeff):
    
    # Generate the scaled wavelet function
    if Scale == 1:
        Scale = 0
        for y in range(numcoeff):
            if (y < numcoeff/(2**(SF+1))):
                Wavelet[y] = 1
            elif ((y >= numcoeff/(2**(SF+1))) and (y<2*numcoeff/(2**(SF+1)))):
                Wavelet[y] = -1
            else:
                Wavelet[y] = 0
            
    for y in range(numcoeff):    
        WM[x, y] = (2**(SF*0.5)) * Wavelet[y]
    
    if (Wavelet[numcoeff-1] < 0):
        Scale = 1
        SF += 1
    else:
        Wavelet = np.roll(Wavelet, numcoeff/(2**SF))

WM /= numcoeff**(0.5)

def GetSamples(Meas, rejectpercent, numcoeff):
    # ------------------------------------------------------------
    # Sample and Transform Waveform
    # ------------------------------------------------------------
    compressed = np.dot(WM, Meas)
    
    
    # Reject lower x% of non-zero samples
    maxcoeff = int(rejectpercent*np.count_nonzero(compressed))
    markedcoeff = 0
    floorcoeff = np.zeros(numcoeff)
    
    sortedaxis = np.argsort(abs(compressed))
    startindex = sortedaxis.size-1
    
    for x in range (0, maxcoeff):
        floorcoeff[sortedaxis[startindex-x]] = compressed[sortedaxis[startindex-x]]
    return compressed
    
    



# ------------------------------------------------------------
# Build the compressed sensing matrix.
# ------------------------------------------------------------
#
print('Building CS Matrix...')
# First determine sparsity factor k (number of non-zero elements in x)
k = int(rejectpercent*numcoeff)


# Use random rows from the Wavelet matrix but not the first row.

m = 2*k

A = np.zeros([m, numcoeff])

mean = 0
sd = 1.0/m

# Generate a compression matrix
# Force A to be a good isometry.  A* A = I
if (ForceIsometry == True):
    mean = 0
    sd = 1.0/m
    A = np.random.normal(mean, sd, [m,numcoeff])
    U,S,V = np.linalg.svd(A,0,1)
    A = np.copy(V)
else:
    A = np.random.normal(mean, sd, [m,numcoeff])
    A += m*np.eye(m, numcoeff)
    # Normalize A
    A/=m

def compress(floorcoeff, A, m, corrupt, nmean, nsd):
    # Compress the signal
    Cvec = np.empty(m)
    
    Cvec = np.dot(A, floorcoeff)
    
    # Generate a sensing matrix
    if (UseDifferentSensingMatrix == True):
        # Force A to be a good isometry.  A* A = I
        if (ForceIsometry == True):
            mean = 0
            sd = 1.0/m
            A = np.random.normal(mean, sd, [m,numcoeff])
            U,S,V = np.linalg.svd(A,0,1)
            A = np.copy(V)
        else:
            A = np.random.normal(mean, sd, [m,numcoeff])
            for x in range (0,m):
                for y in range(0,numcoeff):
                    if (x==y):
                        A[x,y] = np.random.normal(m, sd)
        # Normalize A
        A/=m
        
    # ------------------------------------------------------------
    # Corrupt the x% of the measurements
    # ------------------------------------------------------------
    
    randarray = np.zeros(m)
    
    for x in range(0,m):
        randarray[x] = x
    
    np.random.shuffle(randarray)
    
    for x in range(0,int(m*corrupt)):
        Cvec[randarray[x]] += np.random.normal(nmean,nsd)
    
    return Cvec, A
    
def CSRecover_onenode (Coeff, CSMatrix, UseQuickAlgorithm, Supportprev, testcoeff, SearchSpeed = 1):
    
    """
    ----------------------------------------------------------------
    Recover One Node
    ----------------------------------------------------------------
    
    Recovers the sample voltages for one node using a compressed
    sensing recovery algorithm.
    
    Coeff: Compressed coefficients.
    
    CSMatrix: Sensing matrix used in recovery algorithm.  Should have
                number of compressed coefficients rows and number of 
                samples columns.
                
    SearchSpeed: Speeds up the algorithm by selecting multiple columns
                   of the CSMatrix to the known support per iteration
                   of the recovery algorithm.  A higher speed can
                   introduce error in the recovery and should not exceed
                   number of 3*samples/4.  Default value is 1.
                   
    UseQuickAlgorithm: Speeds up each iteration of the recovery
                         algorithm by replacing the pseudoinverse
                         operation by singular value decomposition.  U
                         and V^H are cancelled out by multiplying by
                         the hermetian of each while the singular value
                         matrix is cancelled out by using the transpose
                         this works because all singular values of a
                         properly generated CSMatrix are approximately 1.
    """

    sparsity = CSMatrix.shape[0]
    numcoeff = CSMatrix.shape[1]
    Support = np.copy(Supportprev)
    xvec = np.zeros([numcoeff])
    SupportMatrix = np.zeros([sparsity,numcoeff])

    F = np.zeros([sparsity,numcoeff])
    r = np.zeros([sparsity, 2])
    r[:,1] = np.transpose(Coeff)
    g = np.empty([numcoeff])
    colnorms = np.empty(numcoeff)

    j=0
    if (SearchSpeed > 3*numcoeff/4):
        print ("Search speed set to high.  Setting to maximum (" + str(3*numcoeff/4) + ")")
        SearchSpeed = 3*numcoeff/4
        
    for x in range (0,numcoeff):
        colnorms[x] = np.linalg.norm(CSMatrix[:,x], 2)
    
    # Calculate the value of xvec
    xold = xvec*0
    DeltaX = 1
    while (DeltaX>10**-8):#abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2))>10**-8):
    
        g = abs(np.dot(np.transpose(CSMatrix), r[:,1]))
    
        # Calculate j and amax
    
        g = np.divide(g,colnorms)
    
        if (SearchSpeed==0):
            j = np.argsort(g)[g.size-1]
            # Update support set 
            Support[j] += 1
            # Add new supported column to AI
            SupportMatrix[:,j]=np.transpose(CSMatrix[:,j])
        else:
            for x in range(0, int((g.size-np.linalg.norm(Support, 0))*SearchSpeed)):
                Support[np.argsort(g)[g.size-1-x]]+=1
                # Add new supported column to AI
                SupportMatrix[:,np.argsort(g)[g.size-1-x]]=np.transpose(CSMatrix[:,np.argsort(g)[g.size-1-x]])
    
        # Calculate new adjustment to X and add the adjustment to it
        if (UseQuickAlgorithm == True):
            Su = np.dot(SupportMatrix, WM)
            u,s,v = np.linalg.svd(Su, full_matrices = 1)
            F[:sparsity, :sparsity] = np.diag(s)
            AL = np.dot(np.transpose(u), Coeff)
            AL = np.dot(np.transpose(F), AL)
            xvec = np.dot(np.transpose(v), AL)
        else:
            xvec = np.dot(np.linalg.pinv(SupportMatrix), Coeff)
    
        # Update r
        newrow = testcoeff - np.dot(CSMatrix,np.dot(WM, xvec))
        r[:,0] = np.copy(r[:,1])
        r[:,1] = np.copy(newrow)
        print("DeltaX: " + str(max(abs(xvec-xold))))
        DeltaX = max(abs(xvec-xold))
        xold= xvec
        print(abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2)))
        
    return xvec, Support
    
import matplotlib.pylab as plt
plt.ion()
Supportprev = np.zeros(numcoeff)
SupportMatrix = np.copy(A)
for x in range(0, Supportprev.shape[0]):
    SupportMatrix[:,x]=np.transpose(A[:,x])
sparsity = A.shape[0]
numcoeff = A.shape[1]
F = np.zeros([sparsity,numcoeff])
u,s,v = np.linalg.svd(A, full_matrices = 1)
F[:sparsity, :sparsity] = np.diag(s)
for i in range (35,xMatrix.shape[0],5):
    floorcoeff = GetSamples(xMatrix[i,:], rejectpercent, numcoeff)
    testcoeff = GetSamples(xMatrix[0,:], rejectpercent, numcoeff)
    Cvec, A = compress(floorcoeff, A, m, corrupt, nmean, nsd)
    TCvec, AB = compress(floorcoeff, A, m, corrupt, nmean, nsd)
    #test = np.dot(A, WM)
    #Cvec = np.dot(test, xMatrix[i,:])
    TVec = np.dot(np.dot(AB, WM), testcoeff)   
    xvec, Support = CSRecover_onenode (Cvec, A, UseSVDDecomp, Supportprev, TVec, SearchSpeed = SearchSpeed)    
    #xvec = np.dot(np.linalg.pinv(np.dot(A, WM)), Cvec)    
    tstop = time.time()
    # ------------------------------------------------------------
    # Reverse Wavelet Transform
    # ------------------------------------------------------------
    #fluncompressed = np.empty(numcoeff)
    fluncompressed = xvec#np.dot(np.transpose(WM), np.transpose(xvec))
    
    
    
    # ------------------------------------------------------------
    # Generate Results
    # ------------------------------------------------------------
    
    orig = np.dot(np.transpose(WM), fluncompressed)
    
    
    # Plot original signal and uncompressed samples
    tstep = 1
    tim = np.zeros(numcoeff)
    for x in range(numcoeff):
        tim[x] = x*tstep
    print ("Total Time: " + str(tstop-ttstart) + "s")
    # Plot Final Result "Regular Format"
    plt.figure(1)
    plt.subplot(211)
    plt.plot(tim, fluncompressed)
    plt.subplot(212)    
    plt.plot(tim, xMatrix[i,:])
    #plt.figure(i/5)    
    #plt.bar(tim, Support-Supportprev)
    #Supportprev = np.copy(Support)
