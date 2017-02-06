import numpy as np

# This library contains code to calculate Leaf Area Density (LAD) profiles from LiDAR-
# derived 3D point cloud data using a model based on radiative tranfer theory.  The 
# method was developed by Detto et al. and is described in full in their publication: 
# Detto, M., G. Asner, Sonnentag, O. and H. C. Muller-Landau. 2015. Using stochastic
# radiative transfer models to estimate leaf area density from multireturn LiDAR in
# complex tropical forests. Journal of Geophysical Research. 

#----------------------------------------------------------------------------------------
# This function calculates the LAD profile based on the radiative transfer model
# Inputs
# pts :: a vector containing the lidar points.  This will have the following structure:
#            - x, y, z, k, c, A (where k=return number, c=classification, A=scan angle
# zi  :: a vector containing the heights at which LAD will be estimated.
# tl  :: leaf inclination model ('planophile'.'erectophile','spherical')
# n   :: a matrix with the dimensions M(number of layers)
#                                    xS(number of scan angles)
#                                    xK(number of points per return number)
#        Note that this is an optional input arguement, but saves reformulating n if
#        already calculated once (i.e. for testing different leaf inclination models).
#        In the first instance - feed it an empty array, i.e. n=np.array([]).  This is 
#        the default setting.
# Outputs
# u   :: a 1D vector of LAD at heights zi
# n   :: matrix (dimensions MxSxK) as described above
# I   :: matrix (dimensions MxSxK) probability for a beam with angle A to intercept fewer
#        than k leaves up to a depth of zi
# U   :: matrix (dimensions MxSxK) probability for leaf at depth zi to be the kth contact
#        along the beam path
def calculate_LAD(pts,zi,max_k,tl,n=np.array([])):
    #print "Calculating LAD using radiative tranfer model"
    # first unpack pts
    #keep = np.all((pts[:,3]<=max_k,pts[:,4]==1),axis=0)
    keep = pts[:,3]<=max_k
    z0 = np.max(zi) - pts[keep,2]
    R  = pts[keep,3]
    A  = np.abs(pts[keep,5])
    # define other variables
    dz = np.abs(zi[0]-zi[1])
    M  = zi.size
    th = np.unique(A)
    S  = th.size
    K  = int(R.max())

    # calculate n(z,s,k), a matrix containing the number of points per depth, scan angle
    # and return number
    if n.size == 0:
        #print "\tCalculating return matrix"
        n = np.zeros((M,S,K),dtype='float')
        for i in range(0,M):
            use1 = np.all((z0>zi[i],z0<=zi[i]+dz),axis=0)
            for j in range(0,S):
                use2 = A[use1]==th[j]
                for k in range(0,K):
                    n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok

    # calculate penetration functions for each scan angle
    #print "\tCalculating penetration functions"
    I = np.zeros((M,S,K),dtype='float')
    U = np.zeros((M,S,K),dtype='float')
    G = np.zeros((M,S),dtype='float')
    control = np.zeros((M,S))
    n0 = np.zeros(S,dtype='float')
    for i in range(0,S):
        n0[i]=np.sum(np.all((R==1, A==th[i]),axis=0))
        n1=np.sum(n[:,i,:],axis=1)
        for j in range(0,K):
            #I[:,i,j]=1-np.cumsum(np.sum(n[:,i,j],axis=2))/n0[i]
            I[:,i,j]=1-np.cumsum(n[:,i,j])/n0[i]
            U[:,i,j]=n[:,i,j]/n1
    
        #Apply correction factor for limited available return
        U[:,i,0]=U[:,i,0]*I[:,i,K-1]
        control[:,i]=n1>0
        G[:,i]=Gfunction(tl,th[i],zi)

    # Compute LAD from ensemble across scan angles
    #print "\tComputing LAD from ensemble across scan angles"
    p = np.sum(n[:,:,0],axis=1)
    p_indices = np.arange(p.size)
    jj = p_indices[p>0][0]-1
    alpha =np.zeros(M,dtype='float')*np.nan
    beta =np.zeros(M,dtype='float')*np.nan
    U0 =np.zeros(M,dtype='float')*np.nan
    for i in range(0,M):
        use = control[i,:]>0
        w=n0[use]/np.sum(n0[use])
        U0[i] = np.inner((G[i,:]/np.abs(np.cos(np.conj(th)/180*np.pi))),(n0/np.sum(n0)))
        if w.size >0:
            # Eq 8a from Detto et al., 2015
            alpha[i]= 1-np.inner(I[i,use,0],w)
            # Eq 8b
            beta[i] = np.inner((U[i,use,0]*G[i,use]/np.abs(np.cos(np.conj(th[use])/180*np.pi))),w)
    alpha[:jj+1]=0
    alpha_indices=np.arange(alpha.size)
    use = alpha_indices[np.isfinite(alpha)]
    alpha = np.interp(alpha_indices,use,alpha[use])
    alpha[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan
    
    beta[:jj+1]=U0[:jj+1]
    beta_indices = np.arange(beta.size)
    use = beta_indices[np.isfinite(beta)]
    beta = np.interp(beta_indices,use,beta[use]) 
    beta[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan
    # numerical solution
    #print "\tNumerical solution"
    u = np.zeros(M,dtype='float')
    for i in range(jj,M): # check indexing here
        # Eq 6
        u[i] = (alpha[i]-np.inner(beta[:i],u[:i]*dz))/(beta[i]*dz)
        if u[i]<0:
            u[i]=0

    u[np.isfinite(u)==False]=0#np.nan
    return u,n,I,U

#----------------------------------------------------------------------------------------
# This function calculates the Ross G function
# ze is in degrees
# LAD = leaf angle distribution in radians
def Gfunction(LAD,ze,zi):
    ze = np.abs(ze)*np.pi/180.
    th = np.linspace(0.00000001,np.pi/2.,101,endpoint=True)
    A = np.cos(ze)*np.cos(th)
    J = 1./np.tan(th)*1./np.tan(ze)
    use = np.abs(J)<=1
    phi = np.arccos(J)
    A[use] = np.cos(th[use])*np.cos(ze)*(1+(2/np.pi)*(np.tan(phi[use])-phi[use]))
    

    if LAD == 'planophile':
        f=2/np.pi*(1+np.cos(2*th))
        G = np.trapz(A*f,th)
    elif LAD == 'erectophile':
        f=2/np.pi*(1-np.cos(2*th))
        G = np.trapz(A*f,th)
    elif LAD == 'spherical':
        G=0.5
    # options here to include more leaf angle distributions (see paper by Detto et al., 2015)
    return G

    
#----------------------------------------------------------------------------------
# This function is an adapted version of Detto et al. (2015)'s radiative transfer
# model.  In this version, however, I account for the variable attenuation of lidar
# returns in the canopy: not every return reaches the ground and therefore there is
# "missing information" from the LiDAR returns, especially within the lower canopy
# that manifests itself as a negative bias in leaf area lower in the canopy.

def calculate_LAD_DTM(pts,zi,max_k,tl):
    #print "Calculating LAD using radiative tranfer model"
    # first unpack pts
    keep = pts[:,3]<=max_k
    #keep = np.all((pts[:,3]<=max_k,pts[:,4]==1),axis=0)
    z0 = np.max(zi) - pts[keep,2]
    R  = pts[keep,3]
    A  = np.abs(pts[keep,5])
    # define other variables
    dz = np.abs(zi[0]-zi[1])
    M  = zi.size
    th = np.unique(A)
    S  = th.size
    K  = int(R.max())

    # calculate n(z,s,k), a matrix containing the number of points per depth, scan angle
    # and return number
    #print "\tCalculating return matrix"
    n = np.zeros((M,S,K))
    for i in range(0,M):
        use1 = np.all((z0>zi[i],z0<=zi[i]+dz),axis=0)
        for j in range(0,S):
            use2 = A[use1]==th[j]
            for k in range(0,K):
                n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok

    #derive correction factor for different return numbers
    CF = np.zeros(K)
    CF[0]=1.
    for i in range(1,K):
        k = i+1
        N_veg_kprev = float(np.all((pts[:,3]==k-1,pts[:,4]==1),axis=0).sum())
        N_k = float((pts[:,3]==k).sum())
        CF[i]=N_veg_kprev/N_k
        n[:,:,i]*=np.product(CF[:k])

    u,n,I,U = calculate_LAD(pts,zi,max_k,tl,n)
    
    return u,n,I,U
