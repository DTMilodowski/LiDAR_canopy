import numpy as np
from scipy import stats
import sys
import os
from matplotlib import pyplot as plt
import fiona
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD
import structural_metrics as struct

from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

# Some files
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
pts_shp = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/CameraTraps/OM2015_16.shp' ## CHANGE AS REQUIRED
laz_files = False ## CHANGE AS REQUIRED

# Some parameters
min_PAD = 0.1
radius = 10.
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7

# Open shapefile and get points & bounding box for ROI
shapefile= fiona.open(pts_shp)

N_pts = len(shapefile)
features = list(shapefile)
pts = np.zeros((N_pts,2))
pt_camera = np.empty(N_pts,dtype='S16')
pt_scale = np.zeros(N_pts)
for pp in range(0,N_pts):
    pts[pp,:]=np.asarray(features[pp]['geometry']['coordinates'])
    pt_camera[pp]=np.asarray(features[pp]['properties']['id'])
    pt_scale[pp]=np.asarray(features[pp]['properties']['scale'])

cameras = np.unique(pt_camera)
N_cameras = cameras.size
structure_dict={}
for cc in range(0,N_cameras):
    camera_dict = {}
    mask = pt_camera==cameras[cc]
    pts_iter = pts[mask,:]
    N_pts_iter = pts_iter.shape[0]
    scale_iter = pt_scale[mask]

    #W,S,E,N = shapefile.bounds
    
    W=np.min(pts_iter[:,0])-radius
    S=np.min(pts_iter[:,1])-radius
    E=np.max(pts_iter[:,0])+radius
    N=np.max(pts_iter[:,1])+radius

    # Read in LiDAR points for region of interest
    # - first find the tiles
    polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
    lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(las_list,polygon, laz_files=laz_files)
    N_trees = len(trees)

    PAD = np.zeros((N_pts_iter, heights.size))
    PAI = np.zeros(N_pts_iter)
    pt_dens = np.zeros(N_pts_iter)
    shape_PAD = np.zeros(N_pts_iter)
    mean_ht = np.zeros(N_pts_iter)
    sd = np.zeros(N_pts_iter)
    skew = np.zeros(N_pts_iter)
    kurt = np.zeros(N_pts_iter)
    n_layers = np.zeros(N_pts_iter)
    frechet = np.zeros(N_pts_iter)
    height = np.zeros(N_pts_iter)
    Shannon = np.zeros(N_pts_iter)

    # Now loop through the shapefile points and sample defined radius, before extracting canopy profiles
    for pp in range(0,N_pts_iter):
        # retrieve point clouds samples        
        sample_pts = []
        for tt in range(0,N_trees):
            ids = trees[tt].query_ball_point(pts_iter[pp], radius)
            sample_pts.append(lidar_pts[ids+starting_ids_for_trees[tt]])
        sample_pts = np.asarray(sample_pts)
        # calculate PAD profile
        heights,first_return_profile,n_ground_returns = PAD.bin_returns(sample_pts, max_height, layer_thickness)
        PAD[pp,:] = PAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
        PAD[pp,heights<min_height]=0

        # Compute metrics...
        PAI[pp] = np.sum(PAD[pp,:])
        pt_dens[pp] = sample_pts.shape[0]/(np.pi*radius*2.)
        temp1, temp2, shape_PAD[pp] = struct.calculate_canopy_shape(heights,PAD[pp,:])
        mean_ht[pp], sd[pp], skew[pp], kurt[pp] = struct.calculate_moments_of_distribution(heights,PAD[pp,:])
        n_layers[pp] = struct.calculate_number_of_contiguous_layers(heights,PAD[pp,:],min_PAD)
        frechet[pp] = struct.calculate_mean_Frechet_distance()
        height[pp] = np.percentile(sample_pts[sample_pts[sample_pts[:,3]==1],2],99)
        Shannon[pp] = struct.calculate_Shannon_index(PAD[pp,:])

    # store into summary dictionary
    camera_dict['PAD'] = PAD.copy()
    camera_dict['PAI'] = PAI.copy()
    camera_dict['scale'] = scale_iter.copy()
    camera_dict['height'] = height.copy()
    camera_dict['pt_density'] = pt_dens.copy()
    camera_dict['shape'] = shape.copy()
    camera_dict['std_dev'] = sd.copy()
    camera_dict['skew'] = skew.copy()
    camera_dict['kurtosis'] = kurt.copy()
    camera_dict['n_layers'] = n_layers.copy()
    camera_dict['frechet'] = frechet.copy()
    camera_dict['ShannonIndex'] = Shannon.copy()
    structure_dict[cameras[cc]] = camera_dict
