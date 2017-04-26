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
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
#Plots = ['Belian']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1
layer_thickness_2m = 2
#layer_thickness_5m = 5
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.

# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
radiative_DTM_LAD = {}
MacArthurHorn_LAD_2m = {}
#MacArthurHorn_LAD_5m = {}
#radiative_LAD_5m = {}
radiative_LAD_2m = {}
#radiative_LAD_noS = {}
# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = lidar.load_lidar_data(las_file)

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
    heights_rad_2m = np.arange(0,max_height+1,2.)
    #heights_rad_5m = np.arange(0,max_height+1,5.)
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_2m = np.zeros((n_subplots,heights_rad_2m.size,max_return))
    #LAD_rad_5m = np.zeros((n_subplots,heights_rad_5m.size,max_return))
    #LAD_rad_2m_noS  = np.zeros((n_subplots,heights_rad_2m.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    heights = np.arange(0,max_height)+1
    heights_2m = np.arange(0,max_height,2.)+2
    #heights_5m = np.arange(0,max_height,5.)+5
    LAD_MH = np.zeros((n_subplots, heights.size))
    LAD_MH_2m = np.zeros((n_subplots, heights_2m.size))
    #LAD_MH_5m = np.zeros((n_subplots,heights_5m.size))
    
    """
    # write point cloud to file for testing
    out = open(pointcloud_dir + Plot_name+'_pts.csv','w')
    out.write('x,y,z,plot,subplot\n')
    """
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
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

