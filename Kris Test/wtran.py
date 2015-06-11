# Program for building wavelets

# First start by importing needed libraries
import numpy as np

# Works with powers of 2.  i.e.  2, 4, 8, 16, 32, 64, ...
numcoeff = 16

# For now we will only use Haar wavelets.

# Start with the Haar Basis function.

# For Haar, the basis function is a box function that encompasses the
# entire period we are measuring over.
Basis = np.ones(numcoeff)

# Next we need the wavelet function.

# This function is derived mathematically by performing the scaling calculation
# with the original box function.

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

# Generate a waveform for compression
import math
f = 1000
t = 1/f
steps = 2048
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
    waveform[x] += math.sin(8*3.141592654*x*1000000/steps)


samples = np.empty(numcoeff)
compressed = np.empty(numcoeff)
for x in range (0, numcoeff):
    samples[x] = waveform[x*steps/numcoeff]

compressed = np.dot(WM, samples)


import matplotlib.pylab as plt
plt.ion()
# Plot original signal and uncompressed samples

samplot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, steps/numcoeff):
        samplot[x*steps/numcoeff + y] = samples[x]




# Plot compressed signal coefficients
complot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, steps/numcoeff):
        complot[x*steps/numcoeff + y] = compressed[x]



# Reject half of the samples that represent the higher frequency components
temp = np.zeros(numcoeff)
for x in range (numcoeff/2):
    temp[x] = compressed[x]
    
for x in range(numcoeff/2, numcoeff):
    temp[x] = 0

uncompressed = np.empty(numcoeff)
uncompressed = np.dot(np.transpose(WM), temp)

# Plot uncompressed signal with rejected coefficients

uncomplot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, steps/numcoeff):
        uncomplot[x*steps/numcoeff + y] = uncompressed[x]
print(temp)


# Reject half of the samples that have an absolute value greater than .5
for x in range (numcoeff):
    if (math.fabs(compressed[x]) >= 5):
        temp[x] = compressed[x]
    else:
        temp[x] = 0

fluncompressed = np.empty(numcoeff)
fluncompressed = np.dot(np.transpose(WM), temp)

# Plot uncompressed signal with rejected coefficients

floorplot = np.empty(steps)
for x in range(numcoeff):
    for y in range (0, steps/numcoeff):
        floorplot[x*steps/numcoeff + y] = fluncompressed[x]


#plt.plot(timeVec, waveform)
#plt.plot(timeVec, samplot)
#plt.plot(timeVec, complot)
#plt.plot(timeVec, uncomplot)
#plt.plot(timeVec, floorplot)

