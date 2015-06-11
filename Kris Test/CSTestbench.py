# Program for building wavelets

# First start by importing needed libraries
import numpy as np
import time

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
corrupt = 0.3

# Mean and Standard Deviation of the noise
nmean = 0.
nsd = 0.001

# Use a conditioned random matrix or force isometry A* A = I (Noisy but works
# can be improved with random matrix "decoupling" or "whitening")
ForceIsometry = False

# Search Speed:  Normally, each iteration of the recovery algorithm
# searches for and adds one new column to the support set (SearchSpeed = 1),
# increasing the search speed (SearchSpeed > 1) will increase the number of
# columns added to the support set per iteration.  However, this will introduce
# error into the result if increased too much.  Noise level affects the maximum
# speed you can use and still get an accurate recovery.
SearchSpeed = 100

# Use SVD Decomposition instead of matrix inverse in recovery algorithm
UseSVDDecomp = True

# Generate separate matrices for sensing and compression?
UseDifferentSensingMatrix = False

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

# ------------------------------------------------------------
# Generate Test Waveform
# ------------------------------------------------------------

print('Generating Waveform...')
# Generate a waveform for compression
import math
f = 1000
t = 1/f
waveform = np.zeros(steps/4)
timeVec = np.empty(steps)
for x in range (steps):
    timeVec[x] = x

base = np.zeros(steps/4)

for x in range(steps/4):
    if (x < steps/8):
        waveform[x] = x
    else:
        waveform[x] = steps/4-x

for x in range(steps/4):
    base[x] = steps/4*math.sin(8*3.141592654*x/steps)

waveform = np.hstack((waveform, base))

for x in range(steps/4):
    if (x < steps/8):
        base[x] = x
    else:
        base[x] = steps/2-x

waveform = np.hstack((waveform, base))

for x in range(steps/4):
    if (x < steps/8):
        base[x] = -1
    else:
        base[x] = 1

waveform = np.hstack((waveform, base))

for x in range(steps):
    waveform[x] += math.sin(8*3.141592654*x*10000/steps)

waveform*=0.00001

# ------------------------------------------------------------
# Sample and Transform Waveform
# ------------------------------------------------------------
samples = np.empty(numcoeff)
compressed = np.empty(numcoeff)
for x in range (0, numcoeff):
    samples[x] = waveform[x*steps/numcoeff]

compressed = np.dot(WM, samples)

# Reject lower x% of non-zero samples
maxcoeff = int(rejectpercent*np.count_nonzero(compressed))
markedcoeff = 0
floorcoeff = np.zeros(numcoeff)

sortedaxis = np.argsort(abs(compressed))
startindex = sortedaxis.size-1

#for x in range (0, maxcoeff):
#    floorcoeff[sortedaxis[startindex-x]] = compressed[sortedaxis[startindex-x]]
floorcoeff = compressed



# ------------------------------------------------------------
# Build the compressed sensing matrix.
# ------------------------------------------------------------
#
print('Building CS Matrix...')
# First determine sparsity factor k (number of non-zero elements in x)
k = int(rejectpercent*numcoeff)
#k = np.linalg.norm(floorcoeff, 0)


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
    

# ------------------------------------------------------------
# Inverse Transform
# ------------------------------------------------------------
# Initialization:
tstart= time.time()
g = np.empty([numcoeff])
T = np.zeros(numcoeff)
xvec = np.zeros([numcoeff])
AI = np.zeros([A.shape[0],A.shape[1]])

F = np.zeros([A.shape[0],A.shape[1]])

rprev = np.zeros(m)
rnext = np.zeros(m)
rprev *= 0
rnext = Cvec
j=0
colnorms = np.empty(numcoeff)
for x in range (0,numcoeff):
    colnorms[x] = np.linalg.norm(A[:,x], 2)

# Calculate the value of xvec
i=0

print('Recovering Signal...') 
import scipy.sparse as sp

while (abs(np.linalg.norm(rprev,2)-np.linalg.norm(rnext,2))>convergence):
    i = i+1

    g = abs(np.dot(np.transpose(A), rnext))

    # Calculate j and amax

    g=np.divide(g,colnorms)

    if (SearchSpeed==1):
        j = np.argsort(g)[g.size-1]
        # Update support set 
        T[j] += 1
        # Add new supported column to AI
        AI[:,j]=np.transpose(A[:,j])
    else:
        for x in range(0, SearchSpeed):
            T[np.argsort(g)[g.size-1-x]]+=1
            # Add new supported column to AI
            AI[:,np.argsort(g)[g.size-1-x]]=np.transpose(A[:,np.argsort(g)[g.size-1-x]])

    # Calculate new adjustment to X and add the adjustment to it
    if (UseSVDDecomp == True):
        u,s,v = np.linalg.svd(AI, full_matrices = 1)
        r = np.eye(u.shape[0], v.shape[0])
        r = sp.dia_matrix((s, 0), shape = (u.shape[0], v.shape[0])).todense()    
        PSIinv = np.dot(np.transpose(u), Cvec)
        PSIinv = np.dot(np.transpose(r), PSIinv)
        xvec = np.dot(np.transpose(v),PSIinv.A1)
    else:
        xvec = np.dot(np.linalg.pinv(AI), Cvec)

    # Update r
    newrow = Cvec - np.dot(A, xvec)
    rprev = np.copy(rnext)
    rnext = newrow

    print(abs(np.linalg.norm(rprev,2)-np.linalg.norm(rnext,2)))

print ("Number of Iterations: " + str(i))
tstop = time.time()
# ------------------------------------------------------------
# Reverse Wavelet Transform
# ------------------------------------------------------------
fluncompressed = np.empty(numcoeff)
import scipy.linalg

fluncompressed = np.dot(WM.transpose(),xvec)
#fluncompressed = np.dot(scipy.linalg.pinv2(np.dot(A, WM)), Cvec)
#fluncompressed = np.dot(np.transpose(WM), np.dot(scipy.linalg.pinv2(A),Cvec))



# ------------------------------------------------------------
# Generate Results
# ------------------------------------------------------------

orig = np.dot(np.transpose(WM), fluncompressed)

import matplotlib.pylab as plt
plt.ion()
# Plot original signal and uncompressed samples
tstep = int(steps/numcoeff)
samplot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, tstep):
        samplot[x*tstep + y] = samples[x]
        
# Plot compressed signal coefficients
complot = np.empty(steps)
origplot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, tstep):
        complot[x*steps/numcoeff + y] = compressed[x]
        origplot[x*tstep + y] = orig[x]

# Plot compressed signal with rejected coefficients

plotcompressed = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, tstep):
        plotcompressed[x*tstep + y] = compressed[x]

# Plot uncompressed recovered signal

uncompressed = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, tstep):
        uncompressed[x*tstep + y] = fluncompressed[x]
tim = np.zeros(numcoeff)
for x in range(numcoeff):
    tim[x] = x*tstep
ttstop = time.time()
print ("Sensing Time: " + str(tstop-tstart) + "s")
print ("Total Time: " + str(ttstop-ttstart) + "s")

# Plot Original Signal
plt.plot(timeVec, waveform)

# Plot Final Result "Boxy Format"
#plt.plot(timeVec, uncompressed)

# Plot Wavelet Transform/UnTransform "Boxy Format"
#plt.plot(timeVec, samplot)

# Plot Compressed Coefficients
#plt.plot(timeVec, complot)

# Plot Final Result "Regular Format"
plt.plot(tim, fluncompressed)


