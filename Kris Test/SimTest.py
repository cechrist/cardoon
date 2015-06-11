# Program for building wavelets
# First start by importing needed libraries
import numpy as np
import math as mp
import matplotlib.pylab as plt

plt.ion()

# Works with powers of 2.  i.e.  2, 4, 8, 16, 32, 64, ...
numcoeff = 80

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


#Newton Method
dx = 1
f = 1

# Create Jacobian
J = np.empty(numcoeff)

# Set initial guess
x0 = np.ones(numcoeff)*.2

# Create X matrix
x = np.zeros(numcoeff)
Fx = np.zeros(numcoeff)
itcount = 0

while dx > 1e-5:
    itcount += 1
    # Build Jacobian
    for l in range(numcoeff):
        J[l]=-2*pi/numcoeff*mp.cos(2*pi*(x0[l])/numcoeff)
        Fx[l] = l*1.0/numcoeff - mp.sin(2*pi*(x0[l])/numcoeff)
     
    for l in range(0,numcoeff):
        x[l] = x0[l] -  Fx[l]/J[l]

    #x = np.dot(WM, x)
    #for l in range(1,numcoeff):
    #    if x[l] < .5*x[0]:
    #        x[l] = 0
    
    dx = np.max(abs(Fx))
    x0 = x
    #x0 = np.dot(np.transpose(WM), x)




