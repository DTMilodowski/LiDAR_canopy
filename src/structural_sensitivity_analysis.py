###############################################################################################################
# This driver function analyses the sensisitivity of the LiDAR-based metrics to spatial scale and point density
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


#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_file = 'Carbon_plot_point_cloud_buffer.las'

gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('S32','f16','f16','f16','f16')}
plot_coordinates = np.genfromtxt(gps_pts_file, skiprows = 1, delimiter = ',',dtype=datatype)

plot = 'Belian'
alt_plot = 'maliau_belian'

max_height = 80.
heights = np.arange(0.,max_height)+1
n_layers = heights.size
plot_width = 100.
sample_res = np.array([1.,2.,5.,10.,20.,25.,50.,100.])
keys = ['1','2','5','10','20','25','50','100']

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
mask = plot_coordinates['plot']==alt_plot
affine=lstsq.least_squares_affine_matrix(plot_coordinates['x'][mask],plot_coordinates['y'][mask],plot_coordinates['x_prime'][mask],plot_coordinates['y_prime'][mask])
plot_bbox = np.array(lstsq.apply_affine_transformation(plot_coordinates['x'][mask],plot_coordinates['y'][mask],affine)).transpose()

pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox)

PAD_profiles_MH = {}
PAD_profiles_rad1 = {}
PAD_profiles_rad2 = {}


for ss in range(0,sample_res.size):
    # Now create the subplot grids
    rows = int(plot_width/sample_res[ss])
    cols = int(plot_width/sample_res[ss])
    x=np.arange(0,plot_width+sample_res[ss],sample_res[ss])
    y=np.arange(0,plot_width+sample_res[ss],sample_res[ss])

    xv,yv=np.asarray(np.meshgrid(x,y))

    xv=xv.reshape(xv.size)
    yv=yv.reshape(yv.size)

    xy=np.asarray([xv,yv])

    n_points_in_grid = xv.size
    xi=np.zeros(n_points_in_grid)
    yi=np.zeros(n_points_in_grid)
    for ii in range(0,n_points_in_grid):
        Xi = np.ones((3,1))
        Xi[0]=xv[ii]
        Xi[1]=yv[ii]
        
        Xi_prime = np.dot(affine,Xi)
        xi[ii] = Xi_prime[0]
        yi[ii] = Xi_prime[1]

    x_prime = xi.reshape(rows+1,cols+1)
    y_prime = yi.reshape(rows+1,cols+1)

    count = 0
    subplot = []
    for i in range(0,rows):
        for j in range(0,cols):
            bbox = [ [x_prime[i,j], x_prime[i+1,j], x_prime[i+1,j+1], x_prime[i,j+1], x_prime[i,j]],
                     [y_prime[i,j], y_prime[i+1,j], y_prime[i+1,j+1], y_prime[i,j+1], y_prime[i,j]] ]
            subplot.append( np.asarray(bbox) )


    # for each of the subplots, clip point cloud and model PAD
    n_subplots = len(subplot)
    PAD_MH = np.zeros(n_layers,n_subplots+2)
    PAD_rad1 = np.zeros(n_layers,n_subplots+2)
    PAD_rad2 = np.zeros(n_layers,n_subplots+2)
    
    plot_lidar_pts = lidar.filter_lidar_data_by_polygon(pts,bbox_polygon)
    

