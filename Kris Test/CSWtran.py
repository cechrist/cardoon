# Program for building wavelets

# First start by importing needed libraries
import numpy as np
import scipy as sp
import scipy.linalg as lin

# Works with powers of 2.  i.e.  2, 4, 8, 16, 32, 64, ...
numcoeff = 2**9

# For now we will only use Haar wavelets.

# Start with the Haar Basis function.

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
print('Generating Waveform...')
# Generate a waveform for compression
import math
f = 1000
t = 1/f
steps = 4096
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


samples = np.empty(numcoeff)
compressed = np.empty(numcoeff)
for x in range (0, numcoeff):
    samples[x] = waveform[x*steps/numcoeff]

compressed = np.dot(WM, samples)

# Reject lower x% of non-zero samples
maxcoeff = int(.3*np.count_nonzero(compressed))
markedcoeff = 0
floorcoeff = np.zeros(numcoeff)

sortedaxis = np.argsort(abs(compressed))
startindex = sortedaxis.size-1

for x in range (0, maxcoeff):
    floorcoeff[sortedaxis[startindex-x]] = compressed[sortedaxis[startindex-x]]

# ------------------------------------------------------------
# Build the compressed sensing matrix.
# ------------------------------------------------------------
#
print('Building CS Matrix...')
# First determine sparsity factor k (number of non-zero elements in x)
k = np.linalg.norm(floorcoeff, 0)


# Use random rows from the Wavelet matrix but not the first row.

# Start at 2k
m = 2*k

A = np.zeros([m, numcoeff])
mean = 0
sd = 1
#for x in range(0, m):
     #A[x,:] = np.random.normal(mean, sd, (numcoeff))

A += np.random.normal(mean, sd, (m,numcoeff))
mean = int(numcoeff/2)
sd = int(numcoeff/5)

#for x in range(0,m):
#    np.random.normal(mean, sd)
#    A[x,:] = np.transpose(WM)[x,:]

#for x in range(0,m):
    #row = np.floor(np.random.normal(mean, sd))
    #A[x,:] = WM[x,:]
q,r = lin.qr(np.random.normal(mean, sd, (numcoeff,numcoeff)))
d =sp. diagonal(r)
ph = d/sp.absolute(d)
q = sp.multiply(q,ph,q)

for x in range(0,m):
    np.random.normal(mean, sd)
    A[x,:] = q[np.random.normal(mean, sd),:]

# Normalize A
A/=np.sqrt(numcoeff)
#for x in range (0, numcoeff):
#    A[:,x] *= np.sqrt(m)/np.linalg.norm(A[:,x],2)
#    A[:,x] = A[:,x]/np.linalg.norm(A[:,x])


# Find the coherence
#M=0
#for x in range (0, numcoeff):
#    for y in range (0, numcoeff):
#        if (x!=y):
#            Mtmp = abs(np.dot(np.transpose(A[:, x]), A[:,y]))
#            if (Mtmp > M):
#                M = Mtmp

# Transform the coefficients
Cvec = np.dot(A, floorcoeff)

# Inverse Transform

# Initialization:
g = np.empty([numcoeff])
T = np.ones(numcoeff)
xvec = np.zeros([numcoeff])
AI = np.zeros([A.shape[0],A.shape[1]])

r = np.zeros([m, 2])
r[:,0] = np.transpose(Cvec)
r[:,1] = np.transpose(Cvec)
j=0
colnorms = np.empty(numcoeff)
for x in range (0,numcoeff):
    colnorms[x] = np.linalg.norm(A[:,x], 2)

# Calculate the value of xvec
i=0
U,S,V = np.linalg.svd(A, 0, 1)
print('Recovering Signal...')    
while (np.linalg.norm(r[:,1],2)>1e-8):
    #break    
    i = i+1
    g = abs(np.dot(np.transpose(A), r[:,0]))
    
    # Calculate j and amax

    g=np.divide(g,colnorms)
    j = np.argsort(g)[g.size-1]
            
    # Update support set 
    T[j] += 1

    # Add new supported column to AI

    AI[:,j]=np.transpose(A[:,j])

    # Calculate new adjustment to X and add the adjustment to it
    xvec = np.dot(np.linalg.pinv(AI), Cvec)
    #xvec = np.dot(np.transpose(AI), Cvec)
    #xvec, residuals, rank, s = np.linalg.lstsq(AI, Cvec)
    # Update r
    newrow = Cvec - np.dot(A, xvec)
    r[:,0] = np.copy(r[:,1])
    r[:,1] = np.copy(newrow)
    
    print(np.linalg.norm(r[:,1],2))


#xvec, residuals, rank, s = np.linalg.lstsq(A, Cvec)
# Reverse Wavelet Transform
fluncompressed = np.empty(numcoeff)

#xvec = np.dot(np.linalg.pinv(A), Cvec)
fluncompressed = np.dot(np.transpose(WM), np.transpose(xvec))

orig = np.dot(np.transpose(WM), floorcoeff)

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


plt.plot(timeVec, waveform)
#plt.plot(timeVec, samplot)
#plt.plot(timeVec, uncompressed)
#plt.plot(timeVec, origplot)
plt.plot(tim, fluncompressed)