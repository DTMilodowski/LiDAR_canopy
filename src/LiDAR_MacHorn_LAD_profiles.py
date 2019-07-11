"""
================================================================================
LiDAR_MacHorn_LAD_profiles.py
--------------------------------------------------------------------------------
Functions that are required to construct PAD profiles using the MacArthur-Horn
approach.
Note that point cloud information currently stored in columns are as follows:
0 = x
1 = y
2 = z
3 = return number
4 = classification
5 = scan angle
6 = GPS time
7 = number of returns recorded for pulse
"""
import numpy as np
import laspy as las
import LiDAR_tools as lidar
import auxilliary_functions as aux

"""
# bin lidar returns into height bins
"""
def bin_returns(pts_in, max_height, layer_thickness):

    mask= pts_in[:,3]==1
    pts=pts_in[mask,:]

    # calculate n ground points
    n_ground_returns = np.sum(pts[:,4]==2)
    # filter to consider only veg returns (class 1 = unclassified (usual
    # scenario); class 3 = low vegetation; class 4 = medium vegetation;
    # class 5 = high vegetation)
    veg_mask = np.any((pts[:,4]==1,pts[:,4]==3,pts[:,4]==4,pts[:,4]==5),axis=0)
    can_pts = pts[veg_mask,:]

    # now set up bins
    lower_lims=np.arange(0,max_height,layer_thickness)
    upper_lims = lower_lims + layer_thickness
    n_bins = lower_lims.size
    profile = np.zeros(n_bins)
    heights = lower_lims+layer_thickness
    n_returns = can_pts.shape[0]

    # bin data
    bin = can_pts[:,2]//layer_thickness
    bin = bin.astype(int)

    for i in range(0,n_returns):
        if bin[i]<n_bins and bin[i]>=0: # in case there are some spuriously high returns outside the permitted range
            profile[bin[i]]+=1

    return heights,profile,n_ground_returns

"""
# bin lidar returns into height bins, but this time weight the returns by the
# number of returns per pulse
"""
def bin_returns_weighted_by_num_returns(pts_in, max_height, layer_thickness):

    pts = pts_in
    weights = 1/pts[:,7]

    # calculate n ground points
    n_ground_returns = np.sum(weights[pts[:,4]==2])
    # filter to consider only veg returns (class 1 = unclassified (usual
    # scenario); class 3 = low vegetation; class 4 = medium vegetation;
    # class 5 = high vegetation)
    veg_mask = np.any((pts[:,4]==1,pts[:,4]==3,pts[:,4]==4,pts[:,4]==5),axis=0)
    can_pts = pts[veg_mask,:]
    can_weights = weights[veg_mask]

    # now set up bins
    lower_lims=np.arange(0,max_height,layer_thickness)
    upper_lims = lower_lims + layer_thickness
    n_bins = lower_lims.size
    profile = np.zeros(n_bins)
    heights = lower_lims+layer_thickness
    n_returns = can_pts.shape[0]

    # bin data
    bin = can_pts[:,2]//layer_thickness
    bin = bin.astype(int)

    for i in range(0,n_returns):
        if bin[i]<n_bins and bin[i]>=0: # in case there are some spuriously high returns outside the permitted range
            #profile[bin[i]]+=1
            profile[bin[i]]+=can_weights[i]

    return heights,profile,n_ground_returns


# Use MacArthur-Horn method to estimate LAD profile from the lidar return profile.  See methods described by Stark et al., Ecology Letters, 2012
def estimate_LAD_MacArthurHorn(lidar_profile,n_ground_returns,layer_thickness,k,zero_nodata=True):
    n_layers = lidar_profile.size
    S = np.zeros(n_layers+1)
    S[1:]=np.cumsum(lidar_profile)
    S+=n_ground_returns
    if zero_nodata:
        S[S==0]=1 # This step is required to stop the base of the profile (final return)
                  # kicking out errors if there are no ground returns. This might be advisable
                  # if you were constructing a simple PAI map, but note that it will be
                  # affected by saturation in cases where the understory is not sampled

    S_in = S[1:]
    S_out= S[:-1]
    LAD_profile = np.log(S_in/S_out)/(k*layer_thickness)

    """
    # Shouldn't have any divide by zeros, but just in case...
    if np.sum(np.isfinite(LAD_profile)==False)>0:
        print np.sum(np.isfinite(LAD_profile)==False)
    LAD_profile[np.isfinite(LAD_profile)==False]==0
    """
    return LAD_profile

# Do some crunching to brute force the best fitting k for MacArther-Horn method.
def minimise_misfit_for_k(kmin,kmax,k_interval,subplot_LAIs,subplot_lidar_profiles,n_ground_returns,layer_thickness,minimum_height=0):

    n_layers = subplot_lidar_profiles.shape[1]
    heights = np.arange(1,n_layers+1)*layer_thickness
    mask = heights<=2
    # first of all, loop through the k values
    ks = np.arange(kmin,kmax+k_interval,k_interval)
    n_ks = ks.size
    misfit = np.zeros(n_ks)
    misfit_b = np.zeros(n_ks)
    n_subplots = subplot_LAIs.size
    for i in range(0,n_ks):
        # now loop through the subplots
        for j in range(0,n_subplots):
            LAD_profile = estimate_LAD_MacArthurHorn(subplot_lidar_profiles[j,:],n_ground_returns[j],layer_thickness,ks[i])
            LAD_profile[mask]=0
            misfit[i] += np.sqrt((np.sum(LAD_profile)-subplot_LAIs[j])**2)
    misfit/=float(n_subplots)
    best_k_LAD_profiles = np.zeros((subplot_lidar_profiles.shape))
    best_k = ks[misfit==np.min(misfit)]
    for j in range(0,n_subplots):
        best_k_LAD_profiles[j,:] = estimate_LAD_MacArtherHorn(subplot_lidar_profiles[j,:],n_ground_returns[j],layer_thickness,best_k)
    print("Field LAI: ", np.mean(subplot_LAIs), "+/-", np.std(subplot_LAIs),"; LiDAR LAI: ",np.mean(np.sum(best_k_LAD_profiles,axis=1)), "+/-", np.std(np.sum(best_k_LAD_profiles,axis=1)), "; best k: ", best_k, " m-1")

    return misfit, ks, best_k_LAD_profiles, best_k


def calculate_bestfit_LAD_profile(subplot_coordinate_file,LAI_file,las_file,Plot_name,minimum_height=0):
    subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
    field_LAI = aux.load_field_LAI(LAI_file)
    lidar_pts = lidar.load_lidar_data(las_file)

    n_subplots = subplot_polygons[Plot_name].shape[0]
    max_height = 80
    layer_thickness = 1
    n_layers = np.ceil(max_height/layer_thickness)
    subplot_lidar_profiles = np.zeros((n_subplots,n_layers))
    n_ground_returns = np.zeros(n_subplots)
    subplot_LAI = np.zeros(n_subplots)

    for i in range(0,n_subplots):
        print("Subplot: ", subplot_labels[Plot_name][i])
        sp_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,subplot_polygons[Plot_name][i,:,:])
        heights,subplot_lidar_profiles[i,:],n_ground_returns[i] = bin_returns(sp_pts, max_height, layer_thickness)
        subplot_LAI[i] = field_LAI['LAI'][np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)]

    kmin = 0.20
    kmax = 5.
    kinc = 0.005
    misfit, ks, best_k_LAD_profiles, best_k = minimise_misfit_for_k(kmin,kmax,kinc,subplot_LAI,subplot_lidar_profiles,n_ground_returns,layer_thickness,minimum_height)
    return heights, best_k_LAD_profiles


# A wrapper function that bins returns and then converts to MacArthur-Horn profile
def estimate_LAD_MacArthurHorn_full(sample_pts,max_height,layer_thickness,minimum_height=2):
    heights,first_return_profile,n_ground_returns = bin_returns(sample_pts, max_height, layer_thickness)
    LAD_MacArthurHorn = estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)
    mask = heights <= minimum_height
    LAD_MacArthurHorn[mask] = 0
    return heights, LAD_MacArthurHorn


# A function that returns the analytical solution to the maximum PAI that can be
# retrieved for a given point density for a specified sampling footprint using
# the MacArthur-Horn approach. As above, the constant k transforms the effective
# PAI into the actual PAI
# This should work with numbers or arrays of input point densities.
def calculate_analytical_limit(point_density,sample_footprint_area,k,theta=0):
    PAI_limit = (np.cos(theta)/k)*np.log(point_density*sample_footprint_area)
    return PAI_limit
