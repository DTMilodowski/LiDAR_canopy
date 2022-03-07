###############################################################################################################
# This driver function analyses the sensisitivity of the LiDAR-based metrics to spatial scale and point density
###############################################################################################################
import numpy as np
import sys

sys.path.append('../')
import LiDAR_io as io
import LiDAR_tools as lidar
import LiDAR_MacHorn_LAD_profiles as LAD1

#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
las_files_temp = np.genfromtxt(las_list,delimiter=',',dtype='S256')
las_files = []
for ii, filepath in enumerate(las_files_temp):
    las_files.append(filepath.decode('utf-8'))
n_files = len(las_files)

plot = 'OP'
print(plot)
plot_width = 100.
sample_res = 20.
max_height = 80
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
subplot_width=20.
kappa = 0.50

heights = np.arange(0,max_height,layer_thickness)+layer_thickness

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
center = [550140.405,512920.137]
polygon = np.array([[center[0]-plot_width/2.,center[1]-plot_width/2.],
                      [center[0]+plot_width/2.,center[1]-plot_width/2.],
                      [center[0]+plot_width/2.,center[1]+plot_width/2.],
                      [center[0]-plot_width/2.,center[1]+plot_width/2.]]) # simple square bounding box applied for all sensitivity analyses

plot_lidar_pts, starting_ids, trees = io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=False)
N_trees = len(trees)
plot_lidar_pts[plot_lidar_pts[:,2]<0,2]=0


# Loop through all the spatial scales of interest
print("generating sample grid")
# Now create the subplot grids
rows = int(plot_width/sample_res)
cols = int(plot_width/sample_res)
x=np.arange(0,plot_width+sample_res,sample_res)+[polygon[0,0]]
y=np.arange(0,plot_width+sample_res,sample_res)+[polygon[0,1]]

x_prime,y_prime=np.asarray(np.meshgrid(x,y))

count = 0
subplot_polygons = []
row_idx = []
col_idx = []
for i in range(0,rows):
    for j in range(0,cols):
        bbox = [ [x_prime[i,j], x_prime[i+1,j], x_prime[i+1,j+1], x_prime[i,j+1], x_prime[i,j]],
                 [y_prime[i,j], y_prime[i+1,j], y_prime[i+1,j+1], y_prime[i,j+1], y_prime[i,j]] ]
        subplot_polygons.append( np.asarray(bbox).transpose() )
        row_idx.append(i)
        col_idx.append(j)

subplots=np.asarray(subplot_polygons)
#-----------------------------------------------------
# now get highest return in each grid cell to define the CHM
print("generating canopy profiles")
n_subplots = subplots.shape[0]

#------------------------------------------------------------------------------------
# CLIP DATA TO PLOT
# clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
print("canopy height = ", np.percentile(plot_lidar_pts[plot_lidar_pts[:,3]==1,2],99), "m")

#------------------------------------------------------------------------------------
# SET UP ARRAYS TO HOST RESULTS
MacArthurHorn_LAD={}
LAD_MH = np.zeros((n_subplots, heights.size))

#------------------------------------------------------------------------------------
# LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
pt_count = 0.
# for each of the subplots, clip point cloud and model PAD and get the metrics
for pp in range(0,n_subplots):
    # query the tree to locate points of interest
    # note that we will only have one tree for the number of points in sensitivity analysis
    centre_x = np.mean(subplots[pp][0:4,0])
    centre_y = np.mean(subplots[pp][0:4,1])
    radius = np.sqrt(sample_res**2/2.)
    ids = trees[0].query_ball_point([centre_x,centre_y], radius)
    sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts[ids],subplots[pp],filter_by_first_return_location=False)
    pt_count += sp_pts.shape[0]
    # now get MacArthur-Horn profiles
    """
    heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
    LAD_MH[pp,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile,
                                                    n_ground_returns,
                                                    layer_thickness,
                                                    kappa,
                                                    zero_nodata=False)
    """
    heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts, max_height, layer_thickness)
    LAD_MH[pp,:] = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,
                                                    weighted_n_ground_returns,
                                                    layer_thickness,
                                                    kappa,
                                                    zero_nodata=False)
                                                                
    # Check for columns for which no pulses hit ground without interception.
    # Below the depth at which the canopy is not penetrated by first returns
    # some methods are infinite. Thus we need to expand the search radius
    # iteratively, so that we can build up a full vertical profile. Note that
    # this potentially gives rise to coarsening effective resolution down the
    # profile, but this seems preferable to more crude gap-filling schemes.
    nodata_test = np.any(~np.isfinite(LAD_MH[pp]))
    centre_x = np.mean(subplots[pp][0:4,0])
    centre_y = np.mean(subplots[pp][0:4,1])

    while nodata_test:
        # expand neighbourhood for point cloud sample
        ids = trees[0].query_ball_point([centre_x,centre_y], radius)
        sp_pts_iter = plot_lidar_pts[ids]

        # get MacArthur-Horn profiles
        nodata_gaps = ~np.isfinite(LAD_MH[pp])
        """
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts_iter, max_height, layer_thickness)
        LAD_MH[pp,nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile,
                                                                n_ground_returns,
                                                                layer_thickness,
                                                                kappa,
                                                                zero_nodata=False)[nodata_gaps]
        """
        heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts_iter, max_height, layer_thickness)
        LAD_MH[subplot_index,nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,
                                                                weighted_n_ground_returns,
                                                                layer_thickness,
                                                                kappa,
                                                                zero_nodata=False)[nodata_gaps]

        # update check
        radius+=1.
        nodata_test = np.any(~np.isfinite(LAD_MH[pp]))

print("average point density = ", pt_count/10.**4, " pts/m^2")

#------------------------------------------------------------------------------------
# CLEANING AND STORING
# now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
# to the LAD distributions
# - remove all profile values below minimum height prior to comparison
mask = heights <= 2.
LAD_MH[:,mask]=np.nan
MacArthurHorn_LAD[plot] = LAD_MH.copy()
#----------------------------------------------------------------------------
np.savez('lidar_canopy_profiles_adaptive_OP_for_synthesis.npz',(MacArthurHorn_LAD))
