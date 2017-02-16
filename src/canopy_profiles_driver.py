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
coordinate_file = 'BALI_plot_coordinates.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_carbonplots_FieldMapcensus2016.csv'
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'

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

# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
radiative_LAD_1st = {}
radiative_LAD_2nd = {}
radiative_LAD_3rd = {}
radiative_DTM_LAD_1st = {}
radiative_DTM_LAD_2nd = {}
radiative_DTM_LAD_3rd = {}
field_LAD = {}
lidar_profiles ={}

MacArthurHorn_LAI = {}
radiative_LAI_1st = {}
radiative_LAI_2nd = {}
radiative_LAI_3rd = {}
radiative_DTM_LAI_1st = {}
radiative_DTM_LAI_2nd = {}
radiative_DTM_LAI_3rd = {}
field_LAI = {}
Hemisfer_LAI = {}

# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_pts = lidar.load_lidar_data(las_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
a, b, CF, r_sq, p, H, D = field.retrieve_crown_allometry_file(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)

# load LAI estimates from hemiphotos
field_LAI = aux.load_field_LAI(LAI_file)


