import scipy.sparse as sparse
import numpy as np
            
def wavelet_matrix(n):
    """
    --------------------------------------------------------------
    Wavelet Matrix
    --------------------------------------------------------------

    Generate a wavelet transform matrix of size 2^n rows by 2^n
    columns.  Matrix is returned in COO format.
    Currently only HARR basis wavelets are supported.

    """
    # ------------------------------------------------------------
    # Haar Basis function.
    # ------------------------------------------------------------
    # For Haar, the basis function is a box function that encompasses the
    # entire period we are measuring over.
    Basis = np.ones(n)
    
    # Next we need the wavelet function.
    
    # This function is derived mathematically by performing the scaling calculation
    # with the original box function.
    
    Wavelet = np.hstack((np.ones(n/2), (-1)*np.ones(n/2)))
    
    
    # Now we need the wavelet transform matrix
    
    # The first row of the matrix will be made up of the Basis functions elements
    # and the second row will be the mother wavelet functions elements.
    
    data = np.array(np.zeros(n**2))
    row = np.array(np.zeros(n**2))
    col = np.array(np.zeros(n**2))
    
    for x in range(n):
        data[x] = Basis[x]/n**(0.5)
        row[x] = 0
        col[x] = x
        data[x+n] = Wavelet[x]/n**(0.5)
        row[x+n] = 1
        col[x+n] = x
    
    # The rest of the rows are made up of the Wavelet functions translations.
    
    Scale = 1;
    SF = 1;
    
    for x in range(2,n):
        
        # Generate the scaled wavelet function
        if Scale == 1:
            Scale = 0
            for y in range(n):
                if (y < n/(2**(SF+1))):
                    Wavelet[y] = 1
                elif ((y >= n/(2**(SF+1))) and (y<2*n/(2**(SF+1)))):
                    Wavelet[y] = -1
                else:
                    Wavelet[y] = 0
        
        # Record Wavelet into the COO matrix data/row/column arrays.        
        for y in range(n):
            if (Wavelet[y] != 0):
                data[n*x+y] = (2.0**(SF*0.5)) * Wavelet[y]/n**(0.5)
                row[n*x+y] = x
                col[n*x+y] = y
    
        # Translate the wavelet function
        if (Wavelet[n-1] < 0):
            Scale = 1
            SF += 1
        else:
            Wavelet = np.roll(Wavelet, n/(2**SF))
    
    # Some quick cleanup
    del Wavelet, Basis, Scale, SF
    
    # Create and return the COO matrix with the data/row/column arrays        

    return sparse.coo_matrix((data, (row, col)), shape=(n, n))


def cs_matrix(m,N):
    
    """
    ----------------------------------------------------------------
    CS Matrix
    ----------------------------------------------------------------
    
    Build a compressed sensing matrix of size m rows by N columns.
    m is equal to the sparsity of the signal (number of
    non-zero elements) to be compressed/sensed and N is equal to the
    number of elements in the uncompressed signal.
    
    Currently only a random gaussian matrix is supported.
    
    """

    # Use random rows from the Wavelet matrix but not the first row.
    m *= 2
    
    A = np.zeros([m, N])
    
    mean = 0
    sd = 1.0/m
    
    # Generate a compression matrix
    A = np.random.normal(mean, sd, [m,N])
    
    # Offset diagonal entries mean by m.
    A += m*np.eye(m, N)

    # Normalize and return A
    A/=m

    return A

def CSMatrix_Inverse (u, s, v, Vect = 0, Type = 'Quick'):
    
    """
    ----------------------------------------------------------------
    CS Matrix Inverse
    ----------------------------------------------------------------
    
    Calculates the inverse of the compressed sensing (CS) matrix.
    Regular Pseudo Inverse (Type = Pinv) and a quick inverse (Quick)
    is supported.  Quick inverse is the default.
    
    u,s,v :  u and v are the u and v matrices from the singular value
            decomposition (SVD) of the CS matrix.  s is the vector of
            singular values of the SVD of the CS matrix.
    
    Type : Type of inverse.  Pinv for regular Pseudo inverse (returns the
          pseudo inverse of the original matrix baised on its u, s, and v)
          and Quick for the quick inverse.  The quick inverse performs the
          inverse step by step and reduces total calculations (returns the
          inverse operation dot(PINV(PSI), Vect) but is less accurate)
    """
    
    
    if (Type == 'Pinv'):
        
        #  This algorithm will use the SVD technique to calculate
        #  pseudo inverse.
        
        r = np.eye(u.shape[0], v.shape[0])
        r[:u.shape[0], :u.shape[0]] = np.diag(s)
        PSIinv = np.dot(np.transpose(r), np.transpose(u))
        PSIinv = np.dot(np.transpose(v), PSIinv)
        
    else:
        
        r = np.eye(u.shape[0], v.shape[0])
        r[:u.shape[0], :u.shape[0]] = np.diag(s)
        PSIinv = np.dot(np.transpose(u), Vect)
        PSIinv = np.dot(np.transpose(r), PSIinv)
        PSIinv = np.dot(np.transpose(v), PSIinv)
    
    del u, s, v, r
    return PSIinv

def CSRecover_onenode (Coeff, CSMatrix, UseQuickAlgorithm, SearchSpeed = 1):
    
    """
    ----------------------------------------------------------------
    Recover One Node
    ----------------------------------------------------------------
    
    Recovers the sample voltages for one node using a compressed
    sensing recovery algorithm.
    
    Coeff : Compressed coefficients.
    
    CSMatrix : Sensing matrix used in recovery algorithm.  Should have
                number of compressed coefficients rows and number of 
                samples columns.
                
    SearchSpeed : Speeds up the algorithm by selecting multiple columns
                   of the CSMatrix to the known support per iteration
                   of the recovery algorithm.  A higher speed can
                   introduce error in the recovery and should not exceed
                   number of 3*samples/4.  Default value is 1.
                   
    UseQuickAlgorithm : Speeds up each iteration of the recovery
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
    Support = np.zeros(numcoeff)
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
    while (abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2))>10**-8):
    
        g = abs(np.dot(np.transpose(CSMatrix), r[:,1]))
    
        # Calculate j and amax
    
        g = np.divide(g,colnorms)
    
        if (SearchSpeed==1):
            j = np.argsort(g)[g.size-1]
            # Update support set 
            Support[j] += 1
            # Add new supported column to AI
            SupportMatrix[:,j]=np.transpose(CSMatrix[:,j])
        else:
            for x in range(0, SearchSpeed):
                Support[np.argsort(g)[g.size-1-x]]+=1
                # Add new supported column to AI
                SupportMatrix[:,np.argsort(g)[g.size-1-x]]=np.transpose(CSMatrix[:,np.argsort(g)[g.size-1-x]])
    
        # Calculate new adjustment to X and add the adjustment to it
        if (UseQuickAlgorithm == True):
            u,s,v = np.linalg.svd(SupportMatrix, full_matrices = 1)
            F[:sparsity, :sparsity] = np.diag(s)
            AL = np.dot(np.transpose(u), Coeff)
            AL = np.dot(np.transpose(F), AL)
            xvec = np.dot(np.transpose(v), AL)
        else:
            xvec = np.dot(np.linalg.pinv(SupportMatrix), Coeff)
    
        # Update r
        newrow = Coeff - np.dot(CSMatrix, xvec)
        r[:,0] = np.copy(r[:,1])
        r[:,1] = np.copy(newrow)
        
    return xvec
    
def CSRecover (Coeff, Numnodes, CSMatrix, UseQuickAlgorithm, SearchSpeed = 1):

    """
    ----------------------------------------------------------------
    Recover All Nodes
    ----------------------------------------------------------------
    
    Recovers the sample voltages for all nodes using a compressed
    sensing recovery algorithm.
    
    Coeff : Compressed coefficients.
    
    Numnodes : Number of nodes
    
    CSMatrix : Sensing matrix kernel used in recovery algorithm.  Should have
                number of compressed coefficients rows and number of 
                samples columns.  Matrix is used on every node.
                
    SearchSpeed : Speeds up the algorithm by selecting multiple columns
                   of the CSMatrix to the known support per iteration
                   of the recovery algorithm.  A higher speed can
                   introduce error in the recovery and should not exceed
                   number of 3*samples/4.  Default value is 1.
                   
    UseQuickAlgorithm : Speeds up each iteration of the recovery
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
    Support = np.zeros(numcoeff*Numnodes)
    xvec = np.zeros([numcoeff*Numnodes])
    SupportMatrix = np.zeros([sparsity*Numnodes,numcoeff*Numnodes])
    RecoveryMatrix = np.kron(np.eye(Numnodes), CSMatrix)
    
    F = np.zeros([sparsity*Numnodes,numcoeff*Numnodes])
    r = np.zeros([sparsity*Numnodes, 2])
    r[:,1] = np.transpose(Coeff)
    g = np.empty([numcoeff*Numnodes])
    colnorms = np.empty(numcoeff*Numnodes)
    j = 0
    
    if (SearchSpeed > 3*numcoeff/4):
        print ("Search speed set to high.  Setting to maximum (" + str(3*numcoeff/4) + ")")
        SearchSpeed = 3*numcoeff/4
        
    for x in range (0,numcoeff*Numnodes):
        colnorms[x] = np.linalg.norm(RecoveryMatrix[:,x], 2)
    
    # Calculate the value of xvec
    while (abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2))>10**-8):
    
        g = abs(np.dot(np.transpose(RecoveryMatrix), r[:,1]))
    
        # Calculate j and amax
    
        g = np.divide(g,colnorms)
    
        if (SearchSpeed==1):
            for x in range (Numnodes):
                j = np.argsort(g[numcoeff*x:numcoeff*x+numcoeff])[numcoeff-1]
                # Update support set 
                if (g[x*numcoeff+j] != 0):
                    Support[numcoeff*x+j] += 1
                    # Add new supported column to AI
                    SupportMatrix[:,numcoeff*x+j]=np.transpose(RecoveryMatrix[:,numcoeff*x+j])

        else:
            for x in range(0, Numnodes):
                sort = np.argsort(g[numcoeff*x:numcoeff*x+numcoeff])
                for y in range(0, SearchSpeed):
                    j = sort[numcoeff-1-y]
                    Support[numcoeff*x+j] += 1
                    # Add new supported column to AI
                    SupportMatrix[:,numcoeff*x+j]=np.transpose(RecoveryMatrix[:,numcoeff*x+j])
    
        # Calculate new adjustment to X and add the adjustment to it
        if (UseQuickAlgorithm == True):
            u,s,v = np.linalg.svd(SupportMatrix, full_matrices = 1)
            F[:sparsity*Numnodes, :sparsity*Numnodes] = np.diag(s)
            AL = np.dot(np.transpose(u), Coeff)
            AL = np.dot(np.transpose(F), AL)
            xvec = np.dot(np.transpose(v), AL)
        else:
            xvec = np.dot(np.linalg.pinv(SupportMatrix), Coeff)
    
        # Update r
        newrow = Coeff - np.dot(RecoveryMatrix, xvec)
        r[:,0] = np.copy(r[:,1])
        r[:,1] = np.copy(newrow)
        print(abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2)))
    
    return xvec

def CSRecover_csc (Coeff, Numnodes, CSMatrix, UseQuickAlgorithm, SearchSpeed = 1):

    """
    ----------------------------------------------------------------
    Recover All Nodes using Sparse Coo format ***UNFINISHED***
    ----------------------------------------------------------------
    
    Recovers the sample voltages for all nodes using a compressed
    sensing recovery algorithm.  Operates in csc format.
    
    Coeff : Compressed coefficients.
    
    Numnodes : Number of nodes
    
    CSMatrix : Sensing matrix kernel used in recovery algorithm.  Should have
                number of compressed coefficients rows and number of 
                samples columns.  Matrix is used on every node.
                
    SearchSpeed : Speeds up the algorithm by selecting multiple columns
                   of the CSMatrix to the known support per iteration
                   of the recovery algorithm.  A higher speed can
                   introduce error in the recovery and should not exceed
                   number of 3*samples/4.  Default value is 1.
                   
    UseQuickAlgorithm : Speeds up each iteration of the recovery
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
    Support = np.zeros(numcoeff*Numnodes)
    xvec = np.zeros([numcoeff*Numnodes])
    SupportMatrix = sparse.csc_matrix(np.zeros([sparsity*Numnodes,numcoeff*Numnodes]))
    RecoveryMatrix = sparse.csc_matrix(np.kron(np.eye(Numnodes), CSMatrix))
    
    F = sparse.csc_matrix(np.zeros([sparsity*Numnodes,numcoeff*Numnodes]))
    r = np.zeros([sparsity*Numnodes, 2])
    r[:,1] = np.transpose(Coeff)
    g = np.empty([numcoeff*Numnodes])
    colnorms = np.empty(numcoeff*Numnodes)
    j = 0
    
    if (SearchSpeed > 3*numcoeff/4):
        print ("Search speed set to high.  Setting to maximum (" + str(3*numcoeff/4) + ")")
        SearchSpeed = 3*numcoeff/4
        
    colnorms = np.copy(RecoveryMatrix.multiply(RecoveryMatrix).sum(0)[0,:])
    
    # Calculate the value of xvec
    while (abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2))>10**-8):
    
        g = abs(RecoveryMatrix.transpose()*r[:,1])
    
        # Calculate j and amax
    
        g = np.divide(g,colnorms)
    
        if (SearchSpeed==1):
            for x in range (Numnodes):
                j = np.argsort(g[numcoeff*x:numcoeff*x+numcoeff])[numcoeff-1]
                # Update support set 
                if (g[x*numcoeff+j] != 0):
                    Support[numcoeff*x+j] += 1
                    # Add new supported column to AI
                    SupportMatrix[:,numcoeff*x+j]=np.transpose(RecoveryMatrix[:,numcoeff*x+j])

        else:
            for x in range(0, Numnodes):
                sort = np.argsort(g[numcoeff*x:numcoeff*x+numcoeff]).flatten()
                for y in range(0, SearchSpeed):
                    j = sort[numcoeff-1-y]
                    Support[numcoeff*x+j] += 1
                    # Add new supported column to AI
                    SupportMatrix[:,numcoeff*x+j]=RecoveryMatrix[:,numcoeff*x+j]
    
        # Calculate new adjustment to X and add the adjustment to it
        if (UseQuickAlgorithm == True):
            u,s,v = np.linalg.svd(SupportMatrix, full_matrices = 1)
            F[:sparsity*Numnodes, :sparsity*Numnodes] = np.diag(s)
            AL = np.dot(np.transpose(u), Coeff)
            AL = np.dot(np.transpose(F), AL)
            xvec = np.dot(np.transpose(v), AL)
        else:
            xvec = np.dot(np.linalg.pinv(SupportMatrix), Coeff)
    
        # Update r
        newrow = Coeff - np.dot(RecoveryMatrix, xvec)
        r[:,0] = np.copy(r[:,1])
        r[:,1] = np.copy(newrow)
        print(abs(np.linalg.norm(r[:,0],2)-np.linalg.norm(r[:,1],2)))
    
    return xvec

Rec = CSRecover_csc(Cvec3, 3, A, UseQuickAlgorithm = True, SearchSpeed = 100)