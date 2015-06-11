# Program for building wavelets

# First start by importing needed libraries
import scipy.sparse as sp
import numpy as np
import sympy as sy

# Variables that are used later on in the code are introduced here.

s = 2 # Dialation factor (unused)

numcoeff = 10 # Number of filter coefficients

Result = np.empty(numcoeff)  # Result Matrix:  See Discussion Below
Equations = []
temp = "" # Temporary variable for building equations
coeff = np.empty(numcoeff) # Array of wavelet scaling coefficients
FC = np.zeros((numcoeff, numcoeff)) # The filter coefficient matrix
M = np.empty((numcoeff, numcoeff))  # The scaling function matrix used to find
                                    # the scaling function F(x).

# First we need to define our symbols.
symlist = []
for i in range (numcoeff):
	symlist.append('a'+str(i))

SFsymslist = []
for i in range (numcoeff):
	SFsymslist.append('Fx'+str(i))
 

# Replace the symbols lists with the actual symbols
syms = sy.symbols(symlist)
SFsyms = sy.symbols(SFsymslist)

# There are two equations that are used for filter coefficients for all wavelets.
# One is that the sum of all coefficients must be equal to a selected dialation factor (s)
# (usually equal to 2 so it will be set to 2 here) and the coeffients must be orthogonal
# so multiplying the coefficient function by itself shifted by a factor (s*l in this case)
# will yield the result 2*delta with the delta function being centered at l.

# Equation 1:
# We need to build a row vector to be used for the first equation (sum of all filter coefficients
# equals s)

# Each equation is built from strings and converted into a set of symbols.
temp = symlist[0]
for i in range (1,numcoeff):
	temp += '+' + symlist[i]

# Convert equation into a set of symbols for use with the sympy solver.
Equations.append(sy.Eq(sy.sympify(temp), 2))


# Equations 2 to numcoeff/2 + 1:
# We will select l = 2 to N/2-1 for the shift for this equation so we can generate a unique
# equation for each value of l.  If a_(k+2*l) is out of or set (i.e. k+2*l > numcoeff) we will
# not use the value of a_(k+2*l).

# Each equation is built from strings and converted into a set of symbols.
for l in range (numcoeff/2 - 1):
	temp = ""
	for k in range (numcoeff):
		if (k+2*l < numcoeff):
			temp +=  '+' + symlist[k] + '*' + symlist[k+2*l]
	if (l == 0):
		delta = 2
	else:
		delta = 0
	# Convert equation into a set of symbols for use with the sympy solver.
	Equations.append(sy.Eq(sy.sympify(temp), delta))


# Equations numcoeff/2 to numcoeff-1

# These equations are defined by our wavelet basis.  For Debuches its 
# SUM(from -inf to inf) { (-1)^k a_k * k^l } = 0 where k is the number of coefficients 
# (numcoeff) and l is the translation of the wavelet basis (each equation is generated 
# for each translation)

for l in range(numcoeff/2):
	temp = ""
	for k in range(numcoeff):
		if ((-1)**k < 0):
			temp += "-"
		else:
			temp += "+"
		temp += str(k**l) + '*' + symlist[k]
	# Convert equation into a set of symbols for use with the sympy solver.
	Equations.append(sy.Eq(sy.sympify(temp), 0))

# Initial guess vector
guess = np.zeros(numcoeff)
guess[0] = 1

# Now we solve the equations and find the falues for our filter coefficients
Coefficients = sy.nsolve(Equations, syms, guess.tolist())

# And the filter coefficients need to be extracted.  Coefficients will have two entries that are made
# up of each of the filter elements (there are two possible conbinations).  Both have the same
# values but in a different order so we can selecte either.  We will use the first since its in the
# desired order.

for i in range(numcoeff):
	coeff[i] = sy.N(Coefficients[i])

# Now to build the wavelet scaling function.  This is generated using the dialation
# equation.  F(x) = SUM(from -inf to inf) { a_k * F(2x - k) } = 0 where k is an integer
# in the range 0 to m and m is the number of coefficients in equation x.

m = 0

for x in range (numcoeff):
    for k in range (m+1):
        if (2*x-k < numcoeff and k < numcoeff):
            FC[x, 2*x-k] = coeff[k]
    m += 2

# And, FC * F(x) = F(x) so (SF-I)F(x) = 0 with (SF-I) = M

M = FC - np.eye(numcoeff)


# To solve for the scaling function the equations for M*F(x) = zeros() are needed.
SFEquations = []
for i in range (numcoeff):
    # Each equation is built from strings and converted into a set of symbols.
    temp = str(M[i,0]) + '*' + SFsymslist[0]
    for j in range (1,numcoeff):
        temp += '+' + str(M[i,j]) + '*' + SFsymslist[j]

    # Convert equation into a set of symbols for use with the sympy solver.
    SFEquations.append(sy.Eq(sy.sympify(temp), 0))

# Initial guess vector
guess = np.ones(numcoeff)*0

# Now we solve the equations and find the scaling function
Fx = sy.nsolve(SFEquations, SFsyms, guess.tolist())
