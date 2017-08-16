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
def calculate_LAD_Detto_original(pts,zi,max_k,tl,n=np.array([])):
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
            #use1 = np.all((z0>zi[i],z0<=zi[i]+dz),axis=0)
            use1 = np.all((z0>=zi[i],z0<zi[i]+dz),axis=0)
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
    #jj = p_indices[p>0][0]-1
    jj = p_indices[p>0][0]
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

    #alpha[:jj+1]=0
    alpha[:jj]=0
    alpha_indices=np.arange(alpha.size)
    use = alpha_indices[np.isfinite(alpha)]
    alpha = np.interp(alpha_indices,use,alpha[use])
    alpha[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan
    
    #beta[:jj+1]=U0[:jj+1]
    beta[:jj]=U0[:jj]
    beta_indices = np.arange(beta.size)
    use = beta_indices[np.isfinite(beta)]
    beta = np.interp(beta_indices,use,beta[use]) 
    beta[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan

    # numerical solution
    #print "\tNumerical solution"
    u = np.zeros(M,dtype='float')
    for i in range(jj,M):
        # Eq 6
        u[i] = (alpha[i]-np.inner(beta[:i],u[:i]*dz))/(beta[i]*dz)
        if u[i]<0:
            u[i]=0
    u[~np.isfinite(u)]=0#np.nan
    return u,n,I,U


# An updated version of Detto's radiative transfer scheme - in this instance, I do not interpolate across
# nodata gaps within the canopy - these generally arise where we have no returns at all, and my preference
# is for the simpler treatment that these are levels with zero LAD
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
            use1 = np.all((z0>=zi[i],z0<zi[i]+dz),axis=0)
            for j in range(0,S):
                use2 = A[use1]==th[j]
                for k in range(0,K):
                    n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok

    ##### New test -> find all scan angles for which there are no returns, and remove these scan angles from the inversion, as
    #                 we don't have the required info at this scan angle to make this calculation.  This happens very rarely -
    #                 i.e. cases where there are very few returns at all at this scan angle.  Future efforts could go further
    #                 here by removing scan angles for which the number of 1st returns falls below a threshold, rather than 0.
    n0_test = np.sum(n[:,:,0],axis=0)
    if np.sum(n0_test==0)>0:
        S_old = th.size
        S = S_old-np.sum(n0_test==0)
        th = th[n0_test>0]
        S = th.size
        # rebuild n, but without scan angles missing all first returns
        n = np.zeros((M,S,K),dtype='float')
        for i in range(0,M):
            use1 = np.all((z0>=zi[i],z0<zi[i]+dz),axis=0)
            for j in range(0,S):
                use2 = A[use1]==th[j]
                for k in range(0,K):
                    n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok        

    ##### New test -> let's add 1 to all 1st return bins for which there are also other returns
    for s in range(0,S):
        n1_test=np.sum(n[:,s,:],axis=1)   # this replicates code used in loop below for calculating penetration functions
        mask = np.all((n1_test>0,n[:,s,0]==0),axis=0) # this finds all layers for which there are no first returns, but 2nd and 3rd returns present.  This leads to a beta value of zero resulting in a nan in the ultimate distribution of u.  I add an arbitrary single return to this canopy layer, since this enables the calculation an estimate of LAD in this layer rather than setting it as zero when it is known that there are reflections at this level.
        n[mask,s,0] = 1

    # calculate penetration functions for each scan angle
    #print "\tCalculating penetration functions"
    I = np.zeros((M,S,K),dtype='float')
    U = np.zeros((M,S,K),dtype='float')
    G = np.zeros((M,S),dtype='float')
    control = np.zeros((M,S))
    n0 = np.zeros(S,dtype='float')
    for i in range(0,S):
        n0[i]=np.sum(n[:,i,0])  # n0 defines the number of first returns for a given scan angle
        n1=np.sum(n[:,i,:],axis=1)
        for j in range(0,K):
            I[:,i,j]=1-np.cumsum(n[:,i,j])/n0[i]
            # account for occasions where there are no returns at a given scan angle
            # i.e. where n1==0 set U equal to 0
            U[n1!=0,i,j]=n[n1!=0,i,j]/n1[n1!=0]
            U[n1==0,i,j]=0
            #U[:,i,j]=n[:,i,j]/n1
        ## This is the original transcription of Detto's code
        ###Apply correction factor for limited available return
        ##U[:,i,0]=U[:,i,0]*I[:,i,K-1]
        ##---------------
        ## Update below - applying correction factor across all available returns, as indicated in
        ## equation 5 from their paper
        for k in range(0,K):
            U[:,i,k]=U[:,i,k]*I[:,i,-1]
        control[:,i]=n1>0
        ##---------------
        G[:,i]=Gfunction(tl,th[i],zi)
    # Compute LAD from ensemble across scan angles
    #print "\tComputing LAD from ensemble across scan angles"
    p = np.sum(n[:,:,0],axis=1)
    p_indices = np.arange(p.size)
    #jj = p_indices[p>0][0]-1
    jj = p_indices[p>0][0]
    alpha =np.zeros(M,dtype='float')*np.nan
    beta =np.zeros(M,dtype='float')*np.nan
    U0 =np.zeros(M,dtype='float')*np.nan
    for i in range(0,M):
        use = control[i,:]>0
        w=n0[use]/np.sum(n0[use])
        U0[i] = np.inner((G[i,:]/np.abs(np.cos(np.conj(th)/180*np.pi))),(n0/np.sum(n0)))
        if w.size >0:
            # Eq 8a from Detto et al., 2015
            alpha[i]= 1.-np.inner(I[i,use,0],w)
            # Eq 8b
            beta[i] = np.inner((U[i,use,0]*G[i,use]/np.abs(np.cos(np.conj(th[use])/180*np.pi))),w)
    
    #### In Detto's original code, interpolation is used to traverse nodata gaps.  I am not sure why this 
    ### approach was used as nodata gaps arise when the number of returns within a given canopy layer is 
    ### zero.  The simplest explanation for this scenario is that there is very low leaf area density and 
    ### hence no interception of the LiDAR pulses.  The interpolation function used by Detto seems likely 
    ### to overestimate leaf area - perhaps by a large amount where there are distinct layers in the
    ### canopy.  In contrast, I set ui  = 0 for these layers.  This is simpler, and matches the modelled
    ### scenario more closely.

    alpha[:jj]=0
    """
    alpha_indices=np.arange(alpha.size)
    use = alpha_indices[np.isfinite(alpha)]
    alpha = np.interp(alpha_indices,use,alpha[use])
    alpha[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan
    """
    alpha[~np.isfinite(alpha)]=0

    beta[:jj]=U0[:jj]
    """
    beta_indices = np.arange(beta.size)
    use = beta_indices[np.isfinite(beta)]
    beta = np.interp(beta_indices,use,beta[use]) 
    beta[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range. I convert extrapolated values to np.nan
    """
    beta[~np.isfinite(beta)]=0
    # numerical solution
    #print "\tNumerical solution"
    u = np.zeros(M,dtype='float')
    for i in range(jj,M):
        # Eq 6
        if beta[i]!=0:
            u[i] = (alpha[i]-np.inner(beta[:i],u[:i]*dz))/(beta[i]*dz)
        if u[i]<0:
            u[i]=0
    #print u
    u[~np.isfinite(u)]=0#np.nan
    return u,n,I,U

#----------------------------------------------------------------------------------------
# This function calculates the Ross G function
# ze is in degrees
# LAD = leaf angle distribution in radians
def Gfunction(LAD,ze,zi):
    ze = np.abs(ze)*np.pi/180.
    if ze == 0:
        ze+=0.0000000001 # arbitrary addition here to prevent error messages.  In practice, this doesn't actually have an impact, because 1/tan(0.00000000001) is >> 1, and therefore this gets filtered out a few lines later.
    th = np.linspace(0.0000000001,np.pi/2.,101,endpoint=True)
    A = np.cos(ze)*np.cos(th)
    J = 1./np.tan(th)*1./np.tan(ze)
    use = np.abs(J)<=1
    phi = np.arccos(J[use])
    A[use] = np.cos(th[use])*np.cos(ze)*(1+(2/np.pi)*(np.tan(phi)-phi))
    

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

def calculate_LAD_DTM(pts,zi,max_k,tl,min_returns = 10):

    # First check -  if there are not sufficient later returns, then the results from their
    # inclusion are likely to be poor.  In the simplest case, there are not enought returns
    # at k=kmax, irrespective of scan angle. Therefore, we decrease kmax to compensate.
    while np.all((np.sum(pts[:,3]==max_k)<min_returns, max_k>1)):
        print '\t\t WARNING: not enough returns for kmax = ', max_k, '. Resetting kmax = ', max_k-1
        max_k-=1

    #print "Calculating LAD using radiative tranfer model"
    # first unpack pts
    keep = pts[:,3]<=max_k
    #keep = np.all((pts[:,3]<=max_k,pts[:,4]==1),axis=0)
    z0 = np.max(zi) - pts[keep,2]
    R  = pts[keep,3]
    Class  = pts[keep,4]
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
        use1 = np.all((z0>=zi[i],z0<zi[i]+dz),axis=0)
        for j in range(0,S):
            use2 = A[use1]==th[j]
            for k in range(0,K):
                n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok
    #print n.sum(axis=0)
    # Second check! 
    # Now loop through the scan angles and aggregate scan angles where there are not sufficient returns
    # to reliably invert profile - this is usually determined by the number of returns for k=kmax
    N_k_per_A = np.sum(n,axis=0)
    N_per_A = np.sum(N_k_per_A,axis=1)
    test = np.min(N_k_per_A,axis=1)
    th_fail = th[test<min_returns]
    th_pass = th[test>=min_returns]
    # if all theta bins fail test, then abandon all theta info
    N_fail = th_fail.size
    if N_fail == S:
        print '\t\t WARNING: not sufficient returns for any of ', S, 'scan angles to invert independently therefore aggregating all scan angles'
        S=1
        th=np.array(np.sum(th*N_per_A)/float(np.sum(N_per_A))).reshape(1)
        A[:] = th
        n = np.sum(n,axis=1).reshape(M,S,K)
        pts[keep,5]=th

    # Otherwise aggregate th_fail into nearest th until no more scan angles fail test.
    elif N_fail>0:
        print '\t\t WARNING: not enough returns for ', N_fail ,'/',S,' scan angles.  Aggregating with neighbouring scan angles'
        for tt in range(0,N_fail):
            # find nearest passing bin
            idx = (np.abs(th_pass-th_fail[tt])).argmin()

            n_Afail = (A==th_fail[tt]).sum()
            idx_Afail = A==th_fail[tt]
            n_Apass = (A==th_pass[idx]).sum()
            idx_Apass = A==th_pass[idx]
            
            # weighted average for th_pass           
            th_pass[idx] = (th_pass[idx]*n_Apass + th_fail[tt]*n_Afail)/float(n_Apass + n_Afail)
            A[np.any((idx_Afail,idx_Apass),axis=0)] = th_pass[idx]

        #Rebuild n    
        pts[keep,5]=A.copy()
        th = np.unique(A)
        S  = th.size
        n = np.zeros((M,S,K))
        for i in range(0,M):
            use1 = np.all((z0>=zi[i],z0<zi[i]+dz),axis=0)
            for j in range(0,S):
                use2 = A[use1]==th[j]
                for k in range(0,K):
                    n[i,j,k]=np.sum(R[use1][use2]==k+1) # check conditional indexing - should be ok
    #print n.sum(axis=0)
    #derive correction factor for different return numbers for each scan angle
    #CF = np.zeros(K)
    CF = np.zeros((S,K))
    CF[:,0]=1.
    for s in range(0,S):
        for k in range(1,K):
            this_k = k+1
            """
            N_veg_kprev = float(np.all((pts[:,3]==this_k-1,pts[:,4]==1),axis=0).sum())
            N_k = float((pts[:,3]==this_k).sum())
            """
            N_veg_kprev = float(np.all((R==this_k-1,Class==1,A==th[s]),axis=0).sum())
            N_k = float(np.all((R==this_k,A==th[s]),axis=0).sum())
            # in the case where there are no returns at return number = k for scan angle s
            # we have no information about vegetation and therefore cannot make a correction
            # This is likely to be rare, and possible future fixes could include exclusion
            # of scan angles that do not have sufficient numbers of returns.
            if N_k == 0:
                CF[s,k]=0
            else:
                CF[s,k]=N_veg_kprev/N_k
            n[:,s,k]*=np.product(CF[s,:this_k])
    #print n.sum(axis=0)
    u,n,I,U = calculate_LAD(pts,zi,max_k,tl,n)
    
    return u,n,I,U



# Overall wrapper for radiative transfer model
def calculate_LAD_rad_DTM_full(sample_pts,max_height,layer_thickness,minimum_height,max_return,leaf_angle_dist='spherical'):
    heights = np.arange(0,max_height+1,layer_thickness)
    u,n,I,U = calculate_LAD_DTM(sample_pts,heights,max_return,leaf_angle_dist)
    LAD_rad=u[::-1]
    mask = heights <= minimum_height
    LAD_rad[mask]=0
    return heights, LAD_rad
