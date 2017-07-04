import sys
import numpy as np
import LiDAR_tools as lidar
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
from matplotlib import pyplot as plt

# start by defining input files
#las_file = 'Carbon_plot_point_cloud_buffer.las'
las_file = sys.argv[3]
# define important parameters for canopy profile estimation
#print sys.argv
# define the xy coordinates of the location of interest
#target_xy = [577601.951,526741.5143] # replace with x and y coordinates of site - this point is located in the middle of LFE
target_xy=[float(sys.argv[1]),float(sys.argv[2])]

#radius = 10. # this defines the neighbourhood radius for sampling the point cloud.  Suggest a value of 10 m to start with.
radius = float(sys.argv[4])
# this set of parameters gets used by both MacArthur-Horn and the radiative transfer model
#max_height = 80
max_height = float(sys.argv[5])
#layer_thickness = 1
layer_thickness = float(sys.argv[7])
#minimum_height = 2. # ignore profiles <2 m due to difficulties distinguishing ground return increasing error
minimum_height = float(sys.argv[6])
# this second set of parameters is only used for the radiative transfer model
#leaf_angle_dist = 'spherical' # other options include 'erectophile' and 'planophile'
leaf_angle_dist = sys.argv[8]
#max_return = 3
max_return = int(sys.argv[9])

# choose method
method = sys.argv[10]

# output details
out_file = sys.argv[11]

# load LiDAR point cloud and clip to neighbourhood around a specified point
lidar_pts = lidar.load_lidar_data(las_file)
sample_pts = lidar_pts # we assume that all fitering has already been done.

#sample_pts = lidar.filter_lidar_data_by_neighbourhood(lidar_pts,target_xy,radius)
#sample_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,polygon) # there is also scope to clip point clouds using a polygon if preferred, but for camera trap, point centres are probably better options.



# MacArthur-Horn method (Stark et al., 2012)
if method == 'macarthur_horn':
    print "using MacArthur-Horn model"
    heights, LAD = LAD1.estimate_LAD_MacArthurHorn_full(sample_pts,max_height,layer_thickness,minimum_height)

#Radiative transfer approach (Milodowski building on Detto et al., 2015)
elif method == 'radiative_transfer':
    print "using radiative transfer model"
    heights,LAD = LAD2.calculate_LAD_rad_DTM_full(sample_pts,max_height,layer_thickness,minimum_height,max_return,leaf_angle_dist)
    
else:
    print "incorrect argument, so using MacArthur-Horn model"
    heights, LAD = LAD1.estimate_LAD_MacArthurHorn_full(sample_pts,max_height,layer_thickness,minimum_height)

out = open(out_file,'w')
N=heights.size
out.write('height, LAD\n')
for i in range(0,N):
    out.write(str(heights[i]) + ',' + str(LAD[i]) + '\n')
out.close()


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
