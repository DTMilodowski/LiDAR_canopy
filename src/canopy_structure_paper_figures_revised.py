###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
############################################################################################################### 
import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats

from matplotlib import rcParams
from datetime import datetime

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/')
import load_field_data as cen

import canopy_structure_plots as csp

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

colour = ['#46E900','#1A2BCE','#E0007F']
rgb = [[70,233,0],[26,43,206],[224,0,127]]
labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
stems_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_small_plot_tree_data.csv'

D1_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN04.csv'
D2_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN05.csv'
SAFE_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_plots/SAFE_SmallTreeCensus_year1only.csv'

# also define output directory (for saving figures)
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/FiguresRevised/'

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
#Plots = ['B North']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.
beta = 0.35 # morphological parameter in inventory-based canopy model
beta_min = 0.26 # 95% range reported for US tree species by Purves et al. (2007)
beta_max = 0.44 # 95% range reported for US tree species by Purves et al. (2007)
kappa = 0.76

heights = np.arange(0,max_height,layer_thickness)+layer_thickness

#------------------------------------------------------------------------------------
# DECLARATIONS
# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
MacArthurHorn_LAD_mean = {}
radiative_LAD = {}
radiative_DTM_LAD = {}
radiative_DTM_LAD_old = {}
radiative_LAD_mean = {}
radiative_DTM_LAD_mean = {}
inventory_LAD = {}
inventory_LAD_std = {}
lidar_profiles ={}
lidar_profiles_adjusted ={}
penetration_limit ={}

MacArthurHorn_LAI = {}
radiative_LAI = {}
radiative_DTM_LAI = {}
MacArthurHorn_LAI_mean = {}
radiative_LAI_mean = {}
radiative_DTM_LAI_mean = {}
inventory_LAI = {}

plot_point_cloud = {}
#------------------------------------------------------------------------------------
# LOADING DATA
# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = io.load_lidar_data(las_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
DBH_BAAD, H_BAAD, D_BAAD = field.load_BAAD_crown_allometry_data(allometry_file)

# load stem data
#stem_data = field.load_SAFE_small_plot_data(stems_file)
DC1_stem_data = field.load_Danum_stem_census(D1_stemcensus) # all subplots
DC2_stem_data = field.load_Danum_stem_census(D2_stemcensus) # all subplots
SAFE_stem_data = field.load_SAFE_small_stem_census(SAFE_stemcensus) # subset of subplots

#------------------------------------------------------------------------------------
# MAIN ANALYSIS
# LiDAR PROFILES LOOP
# loop through all plots to be analysed
plot_point_cloud= np.load('plot_point_clouds.npz')['arr_0'][()]
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    n_subplots = subplot_polygons[Plot_name].shape[0]

    #------------------------------------------------------------------------------------
    # CLIP DATA TO PLOT
    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)
    #plot_lidar_pts = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon,filter_by_first_return_location=True)    
    #plot_point_cloud[Plots[pp]]=plot_lidar_pts.copy()
    plot_lidar_pts=plot_point_cloud[Plots[pp]]
    print "canopy height = ", np.percentile(plot_lidar_pts[plot_lidar_pts[:,3]==1,2],99), "m"
    """
    # LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        print i, "Subplot: ", subplot_labels[Plot_name][i]
        subplot_index = subplot_labels[Plot_name][i]-1
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:],filter_by_first_return_location=True)
        pt_count += sp_pts.shape[0]
        #LAD2.calculate_LAD_rad_DTM_test_individual_scan_angles(sp_pts,80.,1.,2.,3.)
        LAD2.calculate_LAD_rad_DTM_test_ensemble(sp_pts,80.,1.,2.,3.)
    """

    #------------------------------------------------------------------------------------
    # SET UP ARRAYS TO HOST RESULTS
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]
    # set up some arrays to host the radiative transfer based profiles
    heights_rad = np.arange(0,max_height+layer_thickness,layer_thickness)
    LAD_rad = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM_old = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    #heights = np.arange(0,max_height)+1
    LAD_MH = np.zeros((n_subplots, heights.size))
    
    # set up array to host the lidar return profiles
    lidar_return_profiles = np.zeros((n_subplots, heights_rad.size, max_return))
    lidar_return_profiles_adj = np.zeros((n_subplots, heights_rad.size, max_return))

    penetration_lim = np.zeros((n_subplots,heights.size))

    #------------------------------------------------------------------------------------
    # LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        #print "Subplot: ", subplot_labels[Plot_name][i]
        subplot_index = subplot_labels[Plot_name][i]-1
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:],filter_by_first_return_location=True)
        pt_count += sp_pts.shape[0]

        # first loop through the return numbers to calculate the radiative LAD profiles
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad[subplot_index,:,rr]=u.copy()
        lidar_return_profiles[subplot_index,:,:n.shape[2]] = np.sum(n.copy(),axis=1)

        # now repeat but for adjusted profiles, accounting for imperfect penetration of LiDAR pulses into canopy
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM[subplot_index,:,rr]=u.copy()
            u,n,I,U = LAD2.calculate_LAD_DTM_old(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM_old[subplot_index,:,rr]=u.copy()
        lidar_return_profiles_adj[subplot_index,:,:n.shape[2]] = np.sum(n.copy(),axis=1)

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[subplot_index,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
        # get the LiDAR penetration limits for the subplots
        penetration_lim[subplot_index,np.cumsum(first_return_profile)==0]=1.

    print "average point density = ", pt_count/10.**4, " pts/m^2"

    #------------------------------------------------------------------------------------
    # AVERAGE 1 ha PROFILES
    LAD_MH_mean = np.nansum(LAD_MH,axis=0)/(np.sum(np.isfinite(LAD_MH),axis=0)).astype('float')
    LAD_rad_mean = np.nansum(LAD_rad,axis=0)/(np.sum(np.isfinite(LAD_rad),axis=0)).astype('float')
    LAD_rad_DTM_mean = np.nansum(LAD_rad_DTM,axis=0)/(np.sum(np.isfinite(LAD_rad_DTM),axis=0)).astype('float')
    
    # CLEANING AND STORING
    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    """
    # - set NaN values to zero
    LAD_rad[np.isnan(LAD_rad)]=0
    LAD_rad_DTM[np.isnan(LAD_rad_DTM)]=0
    LAD_MH[np.isnan(LAD_MH)]=0
    """
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=np.nan
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_rad[:,mask]=np.nan
    LAD_rad_DTM[:,mask]=np.nan
    LAD_rad_DTM_old[:,mask]=np.nan

    # now the profiles are ready to be stored into their relevant dictionaries
    lidar_profiles[Plot_name] = lidar_return_profiles.copy()
    lidar_profiles_adjusted[Plot_name] = lidar_return_profiles_adj.copy()
    penetration_limit[Plot_name] = penetration_lim.copy()
    
    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
    radiative_LAD[Plot_name] = LAD_rad.copy()
    radiative_DTM_LAD[Plot_name] = LAD_rad_DTM.copy()
    radiative_DTM_LAD_old[Plot_name] = LAD_rad_DTM_old.copy()
    MacArthurHorn_LAD_mean[Plot_name] = LAD_MH_mean.copy()
    radiative_LAD_mean[Plot_name] = LAD_rad_mean.copy()
    radiative_DTM_LAD_mean[Plot_name] = LAD_rad_DTM_mean.copy()

    # Also ready to store their respective LAI profiles
    MacArthurHorn_LAI[Plot_name] = np.nansum(LAD_MH,axis=1)*layer_thickness
    radiative_LAI[Plot_name] = np.nansum(LAD_rad,axis=1)*layer_thickness
    radiative_DTM_LAI[Plot_name] = np.nansum(LAD_rad_DTM,axis=1)*layer_thickness
    MacArthurHorn_LAI_mean[Plot_name] = np.nansum(LAD_MH,axis=1)*layer_thickness
    radiative_LAI_mean[Plot_name] = np.nansum(LAD_rad,axis=1)*layer_thickness
    radiative_DTM_LAI_mean[Plot_name] = np.nansum(LAD_rad_DTM,axis=1)*layer_thickness
#----------------------------------------------------------------------------
np.savez('plot_point_clouds.npz',plot_point_cloud)

#----------------------------------------------------------------------------
# FIELD INVENTORY ANALYSIS - MONTE-CARLO UNCERTAINTY PROPAGATION
# First get allometric relationships required
# 1) height-depth allometry
a, b, CF, r_sq, p, temp, PI_u, PI_l = field.log_log_linear_regression(H_BAAD,D_BAAD)
# 2) DBH-Crown area allometry
a_A, b_A, CF_A, r_sq_A, p_A, temp, PI_u_A, PI_l_A = field.log_log_linear_regression(field_data['DBH_field'],field_data['CrownArea'])
# 3) DBH-Height allometry
a_ht, b_ht, CF_ht, r_sq_ht, p_ht, temp, PI_u_ht, PI_l_ht = field.log_log_linear_regression(field_data['DBH_field'],field_data['Height_field'])

# Version one - no montecarlo routine
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    n_subplots = subplot_polygons[Plot_name].shape[0]

    # mask out dead and broken trees
    dead_mask = np.all((field_data['dead_flag1']==-1,field_data['dead_flag2']==-1,field_data['dead_flag3']==-1),axis=0)
    brokenlive_mask = field_data['brokenlive_flag']==-1
    mask = np.all((field_data['plot']==Plot_name,np.isfinite(field_data['DBH_field']),np.isfinite(field_data['Xfield']),dead_mask,brokenlive_mask),axis=0)

    Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
    #field_LAD_profiles[subplot_index,:], CanopyV = field.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)

    # now building the canopy model
    buff = 20
    xmin = np.min(field_data['Xfield'][mask])
    xmax = np.max(field_data['Xfield'][mask])
    ymin = np.min(field_data['Yfield'][mask])
    ymax = np.max(field_data['Yfield'][mask])
    x = np.arange(xmin-buff,xmax+buff,1.)+0.5
    y = np.arange(ymin-buff,ymax+buff,1.)+0.5
    z = np.arange(0,80.,1.)+0.5
    z=z[::-1]
    x0 = field_data['Xfield'][mask]
    y0 = field_data['Yfield'][mask]
    Z0 = 80 - Ht
    Zmax = Depth
    Rmax = np.sqrt(Area/np.pi)
    beta_list = np.ones(x0.size)*beta
    crown_model = field.generate_3D_canopy(x,y,z,x0,y0,Z0,Zmax,Rmax,beta_list)

    #condense into vertical profile
    field_LAD_profile = np.sum(np.sum(crown_model,axis=1),axis=0)/10.**4
    smallstem_LAD_profiles = np.zeros((n_subplots,heights.size))
    for i in range(0,n_subplots):
        #print "Subplot: ", subplot_labels[Plot_name][i]
        subplot_index = subplot_labels[Plot_name][i]-1
        # add small stem contributions
        if Plot_name == 'DC1':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC1_stem_data['dbh'],DC1_stem_data['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        elif Plot_name == 'DC2':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC2_stem_data['dbh'],DC2_stem_data['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        else:
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(SAFE_stem_data[Plot_name]['dbh'],SAFE_stem_data[Plot_name]['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

        smallstem_LAD_profiles[i,:] = field.calculate_LAD_profiles_from_stem_size_distributions(heights, Area, Depth, Ht, StemDensity, beta)

    field_LAD_profile+=np.mean(smallstem_LAD_profiles,axis=0)
    inventory_LAD[Plot_name] = field_LAD_profile.copy()
    #inventory_LAD_std[Plot_name] = np.std(field_LAD_profiles,axis=0)
    inventory_LAI[Plot_name] = np.sum(field_LAD_profile)


np.savez('canopy_profiles.npz',(inventory_LAD,MacArthurHorn_LAD,radiative_LAD,radiative_DTM_LAD))
"""
n_iter = 1000
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    n_subplots = subplot_polygons[Plot_name].shape[0]
    # set up array to host inventory profiles
    #field_LAD_profiles = np.zeros((n_iter,heights.size))

    # mask out dead and broken trees
    dead_mask = np.all((field_data['dead_flag1']==-1,field_data['dead_flag2']==-1,field_data['dead_flag3']==-1),axis=0)
    brokenlive_mask = field_data['brokenlive_flag']==-1
    mask = np.all((field_data['plot']==Plot_name,np.isfinite(field_data['DBH_field']),dead_mask,brokenlive_mask),axis=0)

    Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
    #field_LAD_profiles[subplot_index,:], CanopyV = field.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)

    # now building the canopy model
    buff = 20
    xmin = np.min(field_data['Xfield'][mask])
    xmax = np.max(field_data['Xfield'][mask])
    ymin = np.min(field_data['Yfield'][mask])
    ymax = np.max(field_data['Yfield'][mask])
    x = np.arange(xmin-buff,xmax+buff,1.)+0.5
    y = np.arange(ymin-buff,ymax+buff,1.)+0.5
    z = np.arange(0,80.,1.)+0.5
    z=z[::-1]
    x0 = field_data['Xfield'][mask]
    y0 = field_data['Yfield'][mask]
    Z0 = 80 - Ht
    Zmax = Depth
    Rmax = np.sqrt(Area)
    beta = np.ones(x0.size)*0.6
    crown_model = field.generate_3D_canopy(x,y,z,x0,y0,Z0,Zmax,Rmax,beta)
    
    # ITERATE MONTE-CARLO PROCEDURE
    #------------------------------------------------------------------------------------
    for ii in range(0,n_iter):
        # now get field inventory estimate - remove dead trees and broken live trees with no crown
        # Note that we only deal with the 1ha plot level estimates as errors relating stem based
        # vs. area based are problematic at subplot level
        Ht,Area,Depth = field.calculate_crown_dimensions_mc(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], DBH_BAAD, H_BAAD,D_BAAD, a_ht, b_ht, a_A, b_A, a, b)
        
        field_LAD_profiles[ii,:], CanopyV = field.calculate_LAD_profiles_generic_mc(heights, Area, Depth, Ht, beta_min,beta_max, subplot_area*25)

        # add small stem contributions
        if Plot_name == 'DC1':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions_mc(DC1_stem_data['dbh'].ravel(),DC1_stem_data['stem_density'].ravel(),a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        elif Plot_name == 'DC2':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions_mc(DC2_stem_data['dbh'].ravel(),DC2_stem_data['stem_density'].ravel(),a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        else:
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions_mc(SAFE_stem_data[Plot_name]['dbh'].ravel(),SAFE_stem_data[Plot_name]['stem_density'].ravel(),a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

        smallstem_LAD_profile = field.calculate_LAD_profiles_from_stem_size_distributions(heights, Area, Depth, Ht, StemDensity, beta)
        field_LAD_profiles[ii,:]+=smallstem_LAD_profile
    
        inventory_LAD[Plot_name] = np.mean(field_LAD_profiles,axis=0)
        inventory_LAD_std[Plot_name] = np.std(field_LAD_profiles,axis=0)
        inventory_LAI[Plot_name] = np.sum(field_LAD_profiles,axis=1)*layer_thickness
"""    
#----------------------------------------------------------------------------
# NOW MAKE PLOTS
# Figure 1 - Location map, with Hansen data and plot locations
figure_name = output_dir+'Fig1_Location_map.png'
figure_number = 1
csp.plot_location_map(figure_name,figure_number)

# Figure 2 - Transmission ratio figure; adapted to illustrate the variance in
# transmission ratios across the landscape

# Figure 3 - Sketch illustrating canopy model
# <Manual figure>

# Figure 4 - Allometric models; include confidence intervals, and add vertical band
# illustrating the 10 cm DBH cutoff
figure_name = output_dir + 'Fig4_allometric_relationships.png'
figure_number = 5
csp.plot_allometric_relationships(figure_name,figure_number,field_file,allometry_file)

# Figure 5 - Point clouds and profiles across degradation gradient
figure_name = output_dir + 'Fig5_pointclouds_and_profiles_test.png'
figure_number = 5
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
csp.plot_point_clouds_and_profiles(figure_name,figure_number, gps_pts_file,plot_point_cloud,heights,heights_rad, lidar_profiles,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,inventory_LAD)

# Figure 6 - PAI comparison between methods
# Add R-sq for both 20 m x 20 m and 1 ha estimates (to do)
figure_name = output_dir + 'Fig6_comparing_LiDAR_PAI_estimates.png'
figure_number = 6
csp.compare_LiDAR_PAI(figure_name,figure_number,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,layer_thickness=1)

# Figure 7 - Cross-plot canopy layers (new)
figure_name = output_dir + 'Fig7_crossplot_LiDAR_PAD_residual_profiles.png'
figure_number = 7
csp.plot_canopy_layer_residuals(figure_name,figure_number,heights,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean)


# Figure 8 - Summary plots from sensitivity analysis for old growth plot and logged plot
# Switch to split violin plots
# <see sensitivity analysis script>

# Figure 9 - Sensitivity analysis of vertical profiles to point density
# Comparison of OG vs Moderately Logged vs. Heavily Logged
# <see sensitivity analysis script>

# Figure 10 - Sensitivity analysis of vertical profiles to spatial resolution
# <see sensitivity analysis script>

# Figure 11 - Sensitivity analysis of unsampled voxels
# <see sensitivity analysis script>

# Figure 12 - Comparison against canopy volume estimates
figure_name = output_dir + 'Fig12_comparing_LiDAR_PAI_and_inventory_test.png'
figure_number = 12
csp.plot_LAI_vs_inventory(figure_name,figure_number,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,inventory_LAD,inventory_LAI)

# Figure 13 - Comparison against canopy volume distributions (analagous to Figure 7)
figure_name = output_dir + 'Fig7_crossplot_LiDAR_and_inventory_profiles_test.png'
figure_number = 13


gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
import plot_LAD_profiles as pprof
pprof.plot_subplot_LAD_profiles(radiative_DTM_LAD['Seraya'][:,:,-1],heights_rad[::-1],'g','test','test')
