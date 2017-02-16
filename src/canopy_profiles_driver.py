###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
############################################################################################################### 
import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/')
import statistics_tools as stats

# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_carbonplots_FieldMapcensus2016.csv'
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'

# also define output directory (for saving figures)
output_dir = './Figures/'

# define important parameters for canopy profile estimation
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.
beta = 0.6 # morphological parameter in inventory-based canopy model

# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
radiative_LAD = {}
radiative_DTM_LAD = {}
inventory_LAD = {}
lidar_profiles ={}
lidar_profiles_adjusted ={}

MacArthurHorn_LAI = {}
radiative_LAI = {}
radiative_DTM_LAI = {}
inventory_LAI = {}
Hemisfer_LAI = {}

# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = lidar.load_lidar_data(las_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
a, b, CF, r_sq, p, H, D = field.retrieve_crown_allometry_file(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)

# load LAI estimates from hemiphotos
field_LAI = aux.load_field_LAI(LAI_file)

# loop through all plots to be analysed
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)
    plot_lidar_pts = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon)
    
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]

    # set up some arrays to host the radiative transfer based profiles
    heights_rad = np.arange(0,max_height+1)
    LAD_rad = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    heights = np.arange(0,max_height)+1
    LAD_MH = np.zeros((n_subplots, heights.size))
    
    # set up array to host the lidar return profiles
    lidar_return_profiles = np.zeros((n_subplots, heights_rad.size, max_return))
    lidar_return_profiles_adj = np.zeros((n_subplots, heights_rad.size, max_return))

    # set up array to host inventory profiles
    field_LAD_profiles = np.zeros((n_subplots,heights.size))

    # set up array to hold the hemisfer LAI estimate
    LAI_hemisfer = np.zeros(n_subplots)

    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:])

        # first of all, loop through the return numbers to calculate the radiative LAD profiles
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad[i,:,rr]=u.copy()
        lidar_return_profiles[i,:,:] = np.sum(n.copy(),axis=1)

        # now repeat but for adjusted profiles, accounting for imperfect penetration of LiDAR pulses into canopy
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM[i,:,rr]=u.copy()
        lidar_return_profiles_adj[i,:,:] = np.sum(n.copy(),axis=1)

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[i,:] = LAD1.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)

        # now get field inventory estimate
        mask = np.all((field_data['plot']==Plot_name,field_data['subplot']==subplot_labels[Plot_name][i]),axis=0)
        Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        field_LAD_profiles[i,:], CanopyV = field.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)
        # now load in the LAI estimates from the hemispherical photographs
        Hemisfer_mask = np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)
        LAI_hemisfer[i] = field_LAI['LAI'][Hemisfer_mask]

    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - set NaN values to zero
    LAD_rad[np.isnan(LAD_rad)]=0
    LAD_rad_DTM[np.isnan(LAD_rad_DTM)]=0
    LAD_MH[np.isnan(LAD_MH)]=0
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=0
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_rad[:,mask]=0
    LAD_rad_DTM[:,mask]=0

    # now the profiles are ready to be stored into their relevant dictionaries
    lidar_profiles[Plot_name] = lidar_return_profiles.copy()
    lidar_profiles_adjusted[Plot_name] = lidar_return_profiles_adj.copy()

    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
    radiative_LAD[Plot_name] = LAD_rad.copy()
    radiative_DTM_LAD[Plot_name] = LAD_rad_DTM.copy()
    inventory_LAD[Plot_name] = field_LAD_profiles.copy()

    # Also ready to store their respective LAI profiles
    MacArthurHorn_LAI[Plot_name] = np.sum(LAD_MH,axis=1)
    radiative_LAI[Plot_name] = np.sum(LAD_rad,axis=1)
    radiative_DTM_LAI[Plot_name] = np.sum(LAD_rad_DTM,axis=1)
    inventory_LAI[Plot_name] = np.sum(field_LAD_profiles,axis=1)
    Hemisfer_LAI[Plot_name] = LAI_hemisfer.copy()

#########################################################################################
# Now plot some figures
#----------------------------------------------------------------------------------------
# Figure 1 - Transmittance ratios for the point cloud.
# This figure plots the number of vegetation returns that propagate through further into
# the canopy, defined as the transmittance ratio.
 
# A figure illustrating transmittance ratio between successive returns 
plt.figure(1, facecolor='White',figsize=[4,4])
ax21 = plt.subplot2grid((1,1),(0,0))
ax21.set_xlabel('return number')
ax21.set_ylabel('transmittance ratio')
for i in range(0,4):
    if i==0:
        ax21.plot(1,1,'o',color='blue')
    else:
        N_veg = float(np.all((all_lidar_pts[:,3]==i,all_lidar_pts[:,4]==1),axis=0).sum())
        N_i = float((all_lidar_pts[:,3]==i+1).sum())
        ax21.plot(i+1,N_i/N_veg,'o',color='blue')
ax21.set_ylim(0,1.1)
ax21.set_xlim(0,5)
plt.tight_layout()
plt.savefig('transmittance_ratios.png')
plt.show()
