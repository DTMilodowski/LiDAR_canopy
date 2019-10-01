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
    alpha[use[-1]+1:]=np.nan # python's interpolation function extends to end of given x range.
                             # I convert extrapolated values to np.nan

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
def calculate_LAD_old(pts,zi,max_k,tl,n=np.array([])):
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

    # find all scan angles for which there are no returns, and remove these scan angles from the inversion
    # This happens very rarely -i.e. cases where there are very few returns at all at this scan angle.
    # Future efforts could go further here by removing scan angles for which the number of 1st returns
    # falls below a threshold, rather than 0.
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
    n_orig = n.copy()
    for s in range(0,S):
        n1_test=np.sum(n[:,s,:],axis=1)
        mask = np.all((n1_test>0,n[:,s,0]==0),axis=0) # this finds all layers for which there are no first returns, but 2nd and 3rd returns present.  This leads to a beta value of zero resulting in a nan in the ultimate distribution of u.  I add an arbitrary single return to this canopy layer, since this enables the calculation an estimate of LAD in this layer rather than setting it as zero when it is known that there are reflections at this level.
        #n[mask,s,0] = 1

    # calculate penetration functions for each scan angle
    #print "\tCalculating penetration functions"
    I = np.zeros((M,S,K),dtype='float')
    U = np.zeros((M,S,K),dtype='float')
    G = np.zeros((M,S),dtype='float')
    control = np.zeros((M,S))
    n0 = np.zeros(S,dtype='float')
    for i in range(0,S):
        n0[i]=np.sum(n[:,i,0])  # n0 defines the number of first returns for a given scan angle
        n1=np.sum(n[:,i,:],axis=1) # n1 defines the total number of returns for a given scan angle
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
        control[:,i]=n1>0 # for each scan angle control masks bins with no returns
        ##---------------
        G[:,i]=Gfunction(tl,th[i],zi)

    #print '-=-=-=-=-=-=-'
    #print control
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

# An updated version of Detto's radiative transfer scheme - in this instance, I do not interpolate across
# nodata gaps within the canopy - these generally arise where we have no returns at all, and my preference
# is for the simpler treatment that these are levels with zero LAD
# Note that I've switched the naming scheme for n0 and n1 so that now n1 = number of first returns and
# nall = total number of returns
def calculate_LAD(pts,zi,max_k,tl,n=np.array([]),test_sensitivity=True):
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

    # find all scan angles for which there are no returns, and remove these scan angles from the inversion
    # This happens rarely -i.e. cases where there are very few returns at all at this scan angle.
    # Future efforts could go further here by removing scan angles for which the number of 1st returns
    # falls below a threshold, rather than 0.
    first_return_count_by_theta = np.zeros(th.size)
    for ii,theta in enumerate(th):
        first_return_count_by_theta[ii] = np.sum(np.all((A==theta,R==1),axis=0))

    S = np.sum(first_return_count_by_theta>0)
    th = th[first_return_count_by_theta>0]

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

    # Build a set of masks for the set of cases for which there needs to be modification:
    # - CASE 1: no returns at all for scan angle for a given canopy layer. We treat this
    #   as simply no interactions with vegetation on that scan angle, therefore the
    #   weighted contribution towards the PAD for that scan angle is zero
    #
    # - CASE 2: there are second or third returns for a given scan angle, but there are
    #   no first returns. In this scenario, beta is zero, therefore the contribution towards
    #   PAD is undefined. However, the presence of returns does indicate that there is
    #   vegetation. We therefore use linear interpolation along this scan angle to
    #   approximate the PAD.
    #
    case1_flag = np.zeros((M,S),dtype='bool')
    case2_flag = np.zeros((M,S),dtype='bool')
    n1 = n[:,:,0] #  number of first returns (per layer, per scan angle)
    nall = np.zeros((M,S)) # hosts number of all returns
    for ss in range(0,S):
        nall[:,ss]=np.sum(n[:,ss,:],axis=1)
        case1_flag[:,ss] = nall[:,ss]==0
        case2_flag[:,ss] = np.all((nall[:,ss]>0,n1[:,ss]==0),axis=0)

    # calculate penetration functions for each scan angle
    #print "\tCalculating penetration functions"
    I = np.zeros((M,S,K),dtype='float')
    U = np.zeros((M,S,K),dtype='float')
    G = np.zeros((M,S),dtype='float')
    for ss in range(0,S):
        case1_mask = case1_flag[:,ss] # no returns in this layer
        case2_mask = case2_flag[:,ss] # no first returns in this layer
        for kk in range(0,K):
            I[:,ss,kk]=1-np.cumsum(n[:,ss,kk])/np.sum(n1[:,ss])
            # account for occasions where there are no returns at a given scan angle
            # i.e. where n1==0 set U equal to 0
            U[case1_mask==False,ss,kk]=n[case1_mask==False,ss,kk]/nall[case1_mask==False,ss]
            U[case1_mask,ss,kk]=0

        ## This is the original transcription of Detto's code
        ###Apply correction factor for limited available return
        ##U[:,i,0]=U[:,i,0]*I[:,i,K-1]
        ##---------------
        ## Update below - applying correction factor across all available returns, as indicated in
        ## equation 5 from their paper
        for k in range(0,K):
            U[:,ss,k]=U[:,ss,k]*I[:,ss,-1]
        #control[:,i]=nall>0 # for each scan angle control masks bins with no returns == case1_flag
        ##---------------
        G[:,ss]=Gfunction(tl,th[ss],zi)

    # Compute LAD from ensemble across scan angles
    #print "\tComputing LAD from ensemble across scan angles"
    p = np.sum(n[:,:,0],axis=1)
    p_indices = np.arange(p.size)
    jj = p_indices[p>0][0] # the first index for which we have returns (top of canopy)
    alpha =np.zeros(M,dtype='float')#*np.nan
    beta =np.zeros(M,dtype='float')#*np.nan
    U0 =np.zeros(M,dtype='float')#*np.nan

    status = np.zeros(M) # 0 = ok, 1 = CASE 1 (across all scan angles), 2 = CASE 2
    penetration_limit = np.zeros(M) # 0 = ok, 1 = no more first returns left

    nall_total = np.sum(nall,axis=0) # sum up 1st return distributions for weighting
    n1_total = np.sum(n1,axis=0) # sum up 1st return distributions for weighting

    for mm in range(0,M):
        use = case1_flag[mm,:]==False # np.all((case1_flag[mm,:]==False,case2_flag[mm,:]==False),axis=0) # control[i,:]>0
        if use.sum() == 0:
            status[mm] = 1 # CASE 1 (no returns, so beta undefined)
        elif np.sum(n1[mm,:]) == 0:
            status[mm]=2 # CASE 2 (no first returns, therefore beta will be zero)
        if np.sum(n1[mm:,:]) == 0:
            penetration_limit[mm] = 1
        U0[mm] = np.inner((G[mm,:]/np.abs(np.cos(np.conj(th)/180*np.pi))),(nall_total/np.sum(nall_total)))
        if use.sum() > 0:#w.size >0:
            w=n1_total[use]/np.sum(n1_total[use]) # w weights the averaging across the ensemble of scan angles
                                                  # angles based on the proportion of first returns
            # Eq 8a from Detto et al., 2015
            alpha[mm]= 1.-np.inner(I[mm,use,0],w)
            # Eq 8b
            beta[mm] = np.inner((U[mm,use,0]*G[mm,use]/np.abs(np.cos(np.conj(th[use])/180*np.pi))),w)

    ### In Detto's original code, interpolation is used to traverse nodata gaps
    ### in alpha and beta. We have modified this slightly.
    ###
    ### Firstly, alpha represents the weighted average fraction of first returns
    ### absorbed up to and including the layer in question. For instances where
    ### there are no first returns (or indeed no returns in a canopy layer,
    ### alpha should be identical to the layer above, indicating no initial
    ### intersections of LiDAR pulses with canopy material.
    ###
    ### Secondly, for beta, the second term in the equation is undefined if
    ### there are no returns (CASE 1) and zero if there are second and
    ### subsequent returns (CASE 2). In both cases this gives rise to an
    ### undefined LAD. The third term is the fraction of returns with a return
    ### number of kmax that are lower in the canopy than the present layer. This
    ### is calculable if there are returns lower in the canopy. Ultimately,
    ### however, beta is zero if n1 is zero, therefore u is undefined. It isn't
    ### clear that beta can be readily interpolated across an arbitrary gap.
    ###
    ### In practice, as the leaf area density is calculated based on alpha and
    ### beta, we avoid interpolating across nodata gaps until after the leaf
    ### area density calculation - i.e. we interpolate across gaps in the
    ### foliage profile, not in alpha or beta, rather than propagate potential
    ### interpolation artefacts.

    alpha[:jj]=0 # no plant area above the top of canopy, so nothing absorbed so far
    beta[:jj]=U0[:jj]

    # numerical solution
    u = np.zeros(M,dtype='float')
    u[:jj]=0
    for mm in range(jj,M):
        # Eq 6
        if status[mm]==0:
            if beta[mm]>0:
                u[mm] = (alpha[mm]-np.inner(beta[:mm],u[:mm]*dz))/(beta[mm]*dz)
            if u[mm]<0:
                u[mm]=0
        # deal with case of no returns, but above penetration limit, so assume no leaf area density (CASE 1)
        if status[mm]==1 and penetration_limit[mm] == 0 :
            u[mm]=0

    u[~np.isfinite(u)]=np.nan # convert np.inf to np.nan here to avoid warning messages

    if test_sensitivity:
        # test sensitivity of u to small change in point distribution by adding one
        # first return iteratively to each canopy layer. Arbitrary threshold set
        # that if sensitivity is an order of magnitde or more, then should be
        # considered unreliable.
        sensitivity = sensitivity_test(pts,n,max_k,tl,zi,u.copy())
        status[sensitivity>2]=3 # set nan-values if sensitivity is greater than factor of two
        #status[sensitivity>10]=3 # set nan-values if sensitivity is greater than factor of ten

    # Now do the interpolation across nodata gaps.
    u[status==1] = 0. # status 1 -> no returns, but above penetration limit
    u[status==2] = np.nan # status 2 -> no first returns, but subsequent returns in
                          #             layer, indicating vegetation is present.
                          #             Set to nan and then interpolate later
    u[status==3] = np.nan # status 3 -> high sensitivity.
                          #             Set to nan and then interpolate later
    u[penetration_limit==1]=np.nan  # Beyond penetration limit of LiDAR in this
                                    # column, so set to nan. This will need gap-
                                    # filling based on a sample over a wider
                                    # neighbourhood.


    # This bit does the interpolation
    u[:jj]=0
    indices = np.arange(u.size)
    use = indices[np.isfinite(u)]
    u = np.interp(indices,use,u[use])
    u[use[-1]+1:]=np.nan    # python's interpolation function by default extends
                            # to the end of the  given x range, which is not
                            # desired here, so reset to nodata.

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

def calculate_LAD_DTM(pts,zi,max_k,tl,min_returns = 10,test_sensitivity=True):

    # First check -  if there are not sufficient later returns, then the results from their
    # inclusion are likely to be poor.  In the simplest case, there are not enought returns
    # at k=kmax, irrespective of scan angle. Therefore, we decrease kmax to compensate.
    while np.all((np.sum(pts[:,3]==max_k)<min_returns, max_k>1)):
        #print('\t\t WARNING: not enough returns for kmax = ', max_k, '. Resetting kmax = ', max_k-1)
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
    pts_z  = pts[keep,2]
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
        #(print '\t\t WARNING: not sufficient returns for any of ', S, 'scan angles to invert independently therefore aggregating all scan angles')
        S=1
        th=np.array(np.sum(th*N_per_A)/float(np.sum(N_per_A))).reshape(1)
        A[:] = th
        n = np.sum(n,axis=1).reshape(M,S,K)
        pts[keep,5]=th

    # Otherwise aggregate th_fail into nearest th until no more scan angles fail test.
    elif N_fail>0:
        #print('\t\t WARNING: not enough returns for ', N_fail ,'/',S,' scan angles.  Aggregating with neighbouring scan angles')
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

    #derive correction factor for different return numbers for each scan angle
    #CF = np.zeros(K)
    CF = np.zeros((S,K))
    CF[:,0]=1.
    for s in range(0,S):
        for k in range(1,K):
            this_k = k+1
            #N_veg_kprev = float(np.all((R==this_k-1,Class==1,A==th[s]),axis=0).sum())
            N_veg_kprev = float(np.all((R==this_k-1,pts_z>=2,A==th[s]),axis=0).sum())
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

    u,n,I,U = calculate_LAD(pts,zi,max_k,tl,n,test_sensitivity=test_sensitivity)

    return u,n,I,U

def calculate_LAD_DTM_old(pts,zi,max_k,tl,min_returns = 10):

    # First check -  if there are not sufficient later returns, then the results from their
    # inclusion are likely to be poor.  In the simplest case, there are not enought returns
    # at k=kmax, irrespective of scan angle. Therefore, we decrease kmax to compensate.
    while np.all((np.sum(pts[:,3]==max_k)<min_returns, max_k>1)):
        #print '\t\t WARNING: not enough returns for kmax = ', max_k, '. Resetting kmax = ', max_k-1
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
    pts_z  = pts[keep,2]

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
        #print '\t\t WARNING: not sufficient returns for any of ', S, 'scan angles to invert independently therefore aggregating all scan angles'
        S=1
        th=np.array(np.sum(th*N_per_A)/float(np.sum(N_per_A))).reshape(1)
        A[:] = th
        n = np.sum(n,axis=1).reshape(M,S,K)
        pts[keep,5]=th

    # Otherwise aggregate th_fail into nearest th until no more scan angles fail test.
    elif N_fail>0:
        #print '\t\t WARNING: not enough returns for ', N_fail ,'/',S,' scan angles.  Aggregating with neighbouring scan angles'
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
            N_veg_kprev_i = float(np.all((R==this_k-1,Class==1,A==th[s]),axis=0).sum())
            N_veg_kprev = float(np.all((R==this_k-1,pts_z>=2,A==th[s]),axis=0).sum())
            N_k = float(np.all((R==this_k,A==th[s]),axis=0).sum())
            #print s, k, n[:,s,k-1].sum(), N_k, N_veg_kprev,N_veg_kprev_i
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
    u,n,I,U = calculate_LAD_old(pts,zi,max_k,tl,n)

    return u,n,I,U

"""
# sensitivity_test
# This function tests the sensitivity of the density estimate at each depth
# level by iterating through each layer, increasing the number of first returns
# in that layer by one, and assessing the change in density. If the new density
# is much lower than the original (usually an indication that there are only one
# or two first returns) then the density estimate is very unreliable, and it is
# recommended that the resolution is coarsened to the point at which it becomes
# stable.
#
# Returns "sensitivity" metric, which is the ratio of the original density to
# the density calculated with modified distribution. High = unreliable
"""
def sensitivity_test(pts,n,max_k,tl,zi,u):
    levels = n.shape[0]
    sensitivity = np.zeros(levels)

    # Full version, which iterates through the canopy layers
    for mm in range(levels):
        n_test=n.copy()
        n_test[mm,:,0]+=1 # add extra first return to this canopy depth
        u_test,n_test,I_test,U_test = calculate_LAD(pts,zi,max_k,tl,n=n_test,test_sensitivity=False)
        sensitivity[mm]=(u/u_test)[mm]

    # alternative version which doesn't require iterating 80 times, but gives
    # similar results
    #n_test=n.copy()
    #n_test[:,:,0]+=1 # add extra first return to this canopy depth
    #u_test,n_test,I_test,U_test = calculate_LAD(pts,zi,max_k,tl,n=n_test,test_sensitivity=False)
    #sensitivity=(u/u_test)

    return sensitivity

# Overall wrapper for radiative transfer model
def calculate_LAD_rad_DTM_full(sample_pts,max_height,layer_thickness,minimum_height,max_return,leaf_angle_dist='spherical'):
    heights = np.arange(0,max_height+1,layer_thickness)
    u,n,I,U = calculate_LAD_DTM(sample_pts,heights,max_return,leaf_angle_dist)
    LAD_rad=u[::-1]
    mask = heights <= minimum_height
    LAD_rad[mask]=0
    return heights, LAD_rad

def calculate_LAD_rad_DTM_test_ensemble(sample_pts,max_height,layer_thickness,minimum_height,max_return,leaf_angle_dist='spherical'):
    heights = np.arange(0,max_height+1,layer_thickness)

    # sample pts in form x,y,z,return,class,scan angle,gps_time
    scan_angles = np.unique(sample_pts[:,5])
    n_A = scan_angles.size
    PAIs = np.zeros(n_A)
    n_returns = np.zeros(n_A)
    u,n,I,U = calculate_LAD_DTM(sample_pts,heights,max_return,leaf_angle_dist,min_returns = 1)
    LAD_rad=u[::-1]
    mask = heights <= minimum_height
    LAD_rad[mask]=0
    PAIs=LAD_rad.sum()
    print("--------")
    print(PAIs)
    return 0#heights, LAD_rad
