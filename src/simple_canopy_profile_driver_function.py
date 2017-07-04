import numpy as np
import LiDAR_tools as lidar
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
from matplotlib import pyplot as plt

# start by defining input files
las_file = '../data/Carbon_plot_point_cloud_buffer.las'

# define important parameters for canopy profile estimation

# define the xy coordinates of the location of interest
target_xy = [577601.951,526741.5143] # replace with x and y coordinates of site - this point is located in the middle of B_north
radius = 10. # this defines the neighbourhood radius for sampling the point cloud.  Suggest a value of 10 m to start with.

# this set of parameters gets used by both MacArthur-Horn and the radiative transfer model
max_height = 80
layer_thickness = 1
minimum_height = 2. # ignore profiles <2 m due to difficulties distinguishing ground return increasing error

# this second set of parameters is only used for the radiative transfer model
leaf_angle_dist = 'spherical' # other options include 'erectophile' and 'planophile'
max_return = 3

# load LiDAR point cloud and clip to neighbourhood around a specified point
lidar_pts = lidar.load_lidar_data(las_file)
sample_pts = lidar.filter_lidar_data_by_neighbourhood(lidar_pts,target_xy,radius)

 # sample_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,polygon) # there is also scope to clip point clouds using a polygon if preferred, but for camera trap, point centres are probably better options.

# MacArthur-Horn method (Stark et al., 2012)
heights, LAD_MacArthurHorn = LAD1.estimate_LAD_MacArthurHorn_full(sample_pts,max_height,layer_thickness,minimum_height)
LAI_MacArthurHorn = np.sum(LAD_MacArthurHorn)

#Radiative transfer approach (Milodowski building on Detto et al., 2015)
heights_rad,LAD_rad = LAD2.calculate_LAD_rad_DTM_full(sample_pts,max_height,layer_thickness,minimum_height,max_return,leaf_angle_dist)
LAI_rad = np.sum(LAD_rad)

# quick plot of the profiles for comparison
plt.plot(LAD_MacArthurHorn,heights,'-')
plt.plot(LAD_rad,heights_rad,'-')
plt.show()

## NOTES
# Of the two methods, the radiative transfer scheme gives systematically higher LAD values,
# although the profile characteristics are consistent.  Moreover the MacArthur Horn estimates
# are typically corrected empirically using a scalar to account for under-prediction of total
# leaf area.  For example Stark et al. (2012) found that LAI values predicted using the
# MacArthur-Horn method with airborne LiDAR predicted LAI values ~70% of true values for two
# Amazon sites.  I wrote a script to find the correction factor based on minimising misfit
# against subplot LAI derived from hemispherical photographs, but it turns out that hemiphotos
# seem to be terrible at predicting LAI in this landscape (paper draft forthcoming!)
#
# A second difference is that the radiative transfer scheme utilises multiple discrete returns.
# This may in part be driving the difference in the observed LAD profiles.  Happy to discuss
# further.
