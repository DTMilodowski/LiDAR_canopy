import numpy as np
from scipy import stats
import sys
from matplotlib import pyplot as plt
import fiona
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD

from matplotlib import rcParams


# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

# Some files
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_laz_files/laz_list.txt' ## CHANGE AS REQUIRED
pts_shp = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/CameraTraps/OM2015_16.shp' ## CHANGE AS REQUIRED
laz_files = True ## CHANGE AS REQUIRED

# Some parameters
radius = 10.
max_height = 80.
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7
layer_thickness = 1

# Open shapefile and get points & bounding box for ROI
shapefile= fiona.open(pts_shp)
N_pts = len(list(shapefile))
pts = np.zeros((N_pts,2))
for pp in range(0,N_pts):
    pts[pp,:]=np.asarray(list(shapefile)[pp]['geometry']['coordinates'])

W,S,E,N = shapefile.bounds
W-=radius
S-=radius
E+=radius
N+=radius

# Read in LiDAR points for region of interest
# - first find the tiles
polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(lastile_list,polygon, laz_files=True)
N_trees = len(trees)


sample_PAD = np.zeros((N_pts, heights.size))
sample_pts = []
# Now loop through the shapefile points and sample defined radius, before extracting canopy profiles
for pp in range(0,N_pts):
    for tt in range(0,N_trees):
        # retrieve point clouds samples
        ids = trees[tt].query_ball_point(pts[pp], radius)
        sample_pts.append(lidar_pts[ids+starting_ids_for_trees[tt]])
        
    # calculate PAD profile
    heights,first_return_profile,n_ground_returns = PAD.bin_returns(sample_pts[pp], max_height, layer_thickness)
    PAD[pp,:] = PAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
        
# Compute metrics...
