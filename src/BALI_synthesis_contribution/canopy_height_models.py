###############################################################################################################
# This driver function analyses the sensisitivity of the LiDAR-based metrics to spatial scale and point density
###############################################################################################################
import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams

import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import least_squares_fitting as lstsq
from scipy import stats
import time
import copy

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_file = 'Carbon_plot_point_cloud_buffer.las'

gps_pts_file = 'GPS_points_file_for_least_squares_fitting_sensitivity_analysis_version.csv'
datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('<U16','f16','f16','f16','f16')}
plot_coordinates = np.genfromtxt(gps_pts_file, skip_header = 0, delimiter = ',',dtype=datatype)

#plot = 'Belian'
#plot = 'E'
#plot = 'B North'
plot = sys.argv[1]
if plot == 'BNorth':
    plot = 'B North'
print(plot)
max_height = 80.
layer_thickness = 1.
heights = np.arange(0.,max_height)+1
heights_rad = np.arange(0,max_height+1)
n_layers = heights.size
plot_width = 100.
sample_res = 0.25

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
mask = plot_coordinates['plot']==plot
affine=lstsq.least_squares_affine_matrix(plot_coordinates['x'][mask],plot_coordinates['y'][mask],plot_coordinates['x_prime'][mask],plot_coordinates['y_prime'][mask])
plot_bbox = np.array(lstsq.apply_affine_transformation(np.array([0.,100.,100.,0.]),np.array([100.,100.,0.,0.]),affine)).transpose() # simple square bounding box applied for all sensitivity analyses

pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox,filter_by_first_return_location=True)
pts[pts[:,2]<0,2]=0
n_returns = pts.shape[0]


# Loop through all the spatial scales of interest
print("generating sample grid")
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
subplots = []
row_idx = []
col_idx = []
for i in range(0,rows):
    for j in range(0,cols):
        bbox = [ [x_prime[i,j], x_prime[i+1,j], x_prime[i+1,j+1], x_prime[i,j+1], x_prime[i,j]],
                 [y_prime[i,j], y_prime[i+1,j], y_prime[i+1,j+1], y_prime[i,j+1], y_prime[i,j]] ]
        subplots.append( np.asarray(bbox).transpose() )
        row_idx.append[i]
        col_idx.append[j]

n_subplots=len(subplot)
chm = np.zeros((rows,cols),dtype='float')*np.nan

#-----------------------------------------------------
# now get highest return in each grid cell to define the CHM
print("generating CHM")
#-------------------
starting_ids, trees = io.create_KDTree(pts)
# for each of the subplots, clip point cloud and model PAD and get the metrics
for pp in range(0,n_subplots):
    # query the tree to locate points of interest
    # note that we will only have one tree for the number of points in sensitivity analysis
    centre_x = np.mean(subplots[keys[ss]][pp][0:4,0])
    centre_y = np.mean(subplots[keys[ss]][pp][0:4,1])
    radius = np.sqrt(sample_res[ss]**2/2.)
    ids = trees[0].query_ball_point([centre_x,centre_y], radius)
    sp_pts = lidar.filter_lidar_data_by_polygon(pts[ids],subplots[pp],filter_by_first_return_location=False)
    #------
    if np.sum(sp_pts[:,2]>=0)>0: # check for returns within this column
        chm[row_idx[pp],col_idx[pp]] = np.max(sp_pts[:,2])
    else:
        # adaptive neighbourhood expansion to account for penetration
        # limitations, particularly for fine grids/low point densities
        nodata_test = np.isnan(chm[row_idx[pp],col_idx[pp]])

        while nodata_test:
            # expand neighbourhood for point cloud sample
            sp_pts_iter = pts[trees[0].query_ball_point([centre_x,centre_y], radius)]
            if np.sum(sp_pts_iter[:,2]>=0)>0: # check for returns within this column
                chm[row_idx[pp],col_idx[pp]] = np.max(sp_pts_iter[:,2])
            radius+=0.1
            nodata_test = np.isnan(chm[row_idx[pp],col_idx[pp]])

end_time = time.time()
print('\t loop time = ', end_time - start_time)

np.save("CHM_25cm_%s.npy" % plot, chm)
