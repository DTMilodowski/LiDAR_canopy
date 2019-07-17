###############################################################################################################
# This driver function analyses the sensisitivity of the LiDAR-based metrics to spatial scale and point density
###############################################################################################################
import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams

sys.path.append('../')
import LiDAR_io as io
import LiDAR_tools as lidar

#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
las_files = np.genfromtxt(las_list,delimiter=',',dtype='S256')
n_files = las_files.size

plot = 'OP'
print(plot)
plot_width = 100.
sample_res = 0.5

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
center = [550140.405,512920.137]
polygon = np.array([[center[0]-plot_width/2.,center[1]-plot_width/2.],
                      [center[0]+plot_width/2.,center[1]-plot_width/2.],
                      [center[0]+plot_width/2.,center[1]+plot_width/2.],
                      [center[0]-plot_width/2.,center[1]+plot_width/2.]]) # simple square bounding box applied for all sensitivity analyses

pts, starting_ids, trees = io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=False)
N_trees = len(trees)
pts[pts[:,2]<0,2]=0
n_returns = pts.shape[0]


# Loop through all the spatial scales of interest
print("generating sample grid")
# Now create the subplot grids
rows = int(plot_width/sample_res)
cols = int(plot_width/sample_res)
x=np.arange(0,plot_width+sample_res,sample_res)+[polygon[0,0]]
y=np.arange(0,plot_width+sample_res,sample_res)+[polygon[0,1]]

x_prime,y_prime=np.asarray(np.meshgrid(x,y))

count = 0
subplots = []
row_idx = []
col_idx = []
for i in range(0,rows):
    for j in range(0,cols):
        bbox = [ [x_prime[i,j], x_prime[i+1,j], x_prime[i+1,j+1], x_prime[i,j+1], x_prime[i,j]],
                 [y_prime[i,j], y_prime[i+1,j], y_prime[i+1,j+1], y_prime[i,j+1], y_prime[i,j]] ]
        subplots.append( np.asarray(bbox).transpose() )
        row_idx.append(i)
        col_idx.append(j)

n_subplots=len(subplots)
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
    centre_x = np.mean(subplots[pp][0:4,0])
    centre_y = np.mean(subplots[pp][0:4,1])
    radius = np.sqrt(sample_res**2/2.)
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
            radius+=0.5
            nodata_test = np.isnan(chm[row_idx[pp],col_idx[pp]])

plt.imshow(chm);plt.show()
np.save("CHM_50cm_%s.npy" % plot, chm)
