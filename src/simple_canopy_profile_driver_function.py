import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import structural_metrics as structure
import plot_LAD_profiles as plot_LAD

# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'

# also define output directory (for saving figures)
output_dir = './Figures/'
pointcloud_dir = './subplot_pointclouds/'

# define important parameters for canopy profile estimation

# this set of parameters gets used by both MacArthur-Horn and the radiative transfer model
max_height = 80
layer_thickness = 1
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2. # ignore profiles <2m due to difficulties distinguishing ground return increasing error

# this second set of parameters is only used for the radiative transfer model
leaf_angle_dist = 'spherical'
max_return = 3
heights_rad = np.arange(0,max_height+1)

# load LiDAR point cloud and clip to polygon
lidar_pts = lidar.load_lidar_data(las_file)
coord_pairs = 
sample_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,polygon)

# MacArthur-Horn method (Stark et al., 2012)
heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
LAD_MacArthurHorn = LAD1.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)

#Radiative transfer approach (Milodowski building on Detto et al., 2015)
u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,leaf_angle_dist)
LAD_rad=u[::-1]

