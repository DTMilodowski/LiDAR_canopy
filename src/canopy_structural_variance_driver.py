import numpy as np
from scipy import stats
import sys
import os
from matplotlib import pyplot as plt
import fiona
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD1
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

# some output files
profile_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/Data/CameraTraps/PAD_profiles_for_'
stats_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/Data/CameraTraps/Structural_stats_for_SAFE_camera_traps.csv'

# Some parameters
min_PAD = 0.1
radius = 10.
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7

# Open shapefile and get points & bounding box for ROI
print 'Loading sample locations from shapefile'
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

print 'Starting loop through camera traps'
for cc in range(0,N_cameras):
    print '\tCamera trap: %s; number %.0f of %.0f' % (cameras[cc], cc+1, N_cameras)
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
    lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=laz_files)
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
        sample_pts = np.array([])
        for tt in range(0,N_trees):
            ids = trees[tt].query_ball_point(pts_iter[pp], radius)
            if len(ids)>0:
                if sample_pts.size==0:
                    sample_pts = lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]
                else:
                    sample_pts = np.concatenate((sample_pts,lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]),axis=0)

        if sample_pts.size == 0:
            PAD[pp,:] = -9999
            PAI[pp] = -9999
            pt_dens[pp] = 0.
            shape_PAD[pp] = -9999
            mean_ht[pp] = -9999
            sd[pp] = -9999
            skew[pp] = -9999
            kurt[pp] = -9999
            n_layers[pp] = -9999
            frechet[pp] = -9999
            height[pp] = -9999
            Shannon[pp] = -9999

        elif np.sum(sample_pts[:,3]==1) == 0:
            PAD[pp,:] = -9999
            PAI[pp] = -9999
            pt_dens[pp] = 0.
            shape_PAD[pp] = -9999
            mean_ht[pp] = -9999
            sd[pp] = -9999
            skew[pp] = -9999
            kurt[pp] = -9999
            n_layers[pp] = -9999
            frechet[pp] = -9999
            height[pp] = -9999
            Shannon[pp] = -9999            

        else:
            #print sample_pts.shape, np.sum(sample_pts[:,3]==1)
            # calculate PAD profile
            heights,first_return_profile,n_ground_returns = PAD1.bin_returns(sample_pts, max_height, layer_thickness)
            PAD[pp,:] = PAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
            PAD[pp,heights<min_height]=0

            # Compute metrics...
            PAI[pp] = np.sum(PAD[pp,:])
            pt_dens[pp] = sample_pts.shape[0]/(np.pi*radius**2.)

            if PAI[0] == -9999:
                frechet[pp] = -9999
            elif PAI[pp] == 0:
                frechet = -9999
            else:
                frechet[pp] = struct.calculate_mean_Frechet_distance(np.asarray([PAD[0,:],PAD[pp,:]]),heights)
            height[pp] = np.percentile(sample_pts[sample_pts[:,3]==1,2],99)

            if PAI[pp]==0:
                Shannon = -9999
                n_layers[pp] = -9999
                shape_PAD[pp] = -9999
                mean_ht[pp], sd[pp], skew[pp], kurt[pp] = -9999
            else:
                Shannon[pp] = struct.calculate_Shannon_index(PAD[pp,:])
                temp1, temp2, shape_PAD[pp] = struct.calculate_canopy_shape(heights,PAD[pp,:])
                mean_ht[pp], sd[pp], skew[pp], kurt[pp] = struct.calculate_moments_of_distribution(heights,PAD[pp,:])
                n_layers[pp] = struct.calculate_number_of_contiguous_layers(heights,PAD[pp,:],min_PAD)


        # Write PAD profiles to file
        f = open(profile_file+cameras[cc]+".csv","w") #opens file
        f.write("height, ")
        n_profiles = PAD.shape[0]
        for pp in range(0,n_profiles):
            f.write(cameras[cc]+'_'+str(pp+1).zfill(3)+', ')
        f.write("\n")

        for ll in range(0,heights.size):
            f.write('%.1f ' % heights[ll])
            for pp in range(0,n_profiles):
                f.write(', %.5f' % PAD[pp,ll])
            f.write("\n")
        f.close()

        sample_pts=None

    # store into summary dictionary
    #camera_dict['PAD'] = PAD.copy()
    camera_dict['PAI'] = PAI.copy()
    camera_dict['scale'] = scale_iter.copy()
    camera_dict['height'] = height.copy()
    camera_dict['pt_density'] = pt_dens.copy()
    camera_dict['shape'] = shape_PAD.copy()
    camera_dict['mean_ht'] = mean_ht.copy()
    camera_dict['std_dev'] = sd.copy()
    camera_dict['skew'] = skew.copy()
    camera_dict['kurtosis'] = kurt.copy()
    camera_dict['n_layers'] = n_layers.copy()
    camera_dict['Frechet'] = frechet.copy()
    camera_dict['Shannon'] = Shannon.copy()
    structure_dict[cameras[cc]] = camera_dict


# finally write two output files 1) the profiles; 2) the stats
"""
print "\t writing ", profile_file 
f = open(profile_file,"w") #opens file
f.write("height, ")
for cc in range(0,N_cameras):
    cam=cameras[cc]
    PAD = structure_dict[cam]['PAD']
    n_profiles = PAD.shape[0]
    for pp in range(0,n_profiles):
        f.write(cam+'_'+str(pp+1).zfill(3)+', ')
    f.write("\n")

for ll in range(0,heights.size):
    f.write('%.1f, ' % heights[ll])
    for cc in range(0,N_cameras):
        cam=cameras[cc]
        PAD = structure_dict[cam]['PAD']
        n_profiles = PAD.shape[0]
        for pp in range(0,n_profiles):
            f.write('%.5f' % PAD[pp,ll])
        f.write("\n")
f.close()
"""

print "\t writing ", stats_file 
f = open(stats_file,"w") #opens file
f.write("profile, camera, scale, point_density, canopy_height, PAI, n_layers, canopy_shape, mean_PAD, std_PAD, skew_PAD, kurtosis_PAD, Frechet_distance, Shannon_index\n")
for cc in range(0,N_cameras):
    cam=cameras[cc]
    PAI = structure_dict[cam]['PAI']
    scale = structure_dict[cam]['scale']
    height = structure_dict[cam]['height']
    density = structure_dict[cam]['pt_density']
    layers = structure_dict[cam]['n_layers']
    shape = structure_dict[cam]['shape']
    mean_ht = structure_dict[cam]['mean_ht']
    sd = structure_dict[cam]['std_dev']
    skew = structure_dict[cam]['skew']
    kurt = structure_dict[cam]['kurtosis']
    fr = structure_dict[cam]['Frechet']
    Shan = structure_dict[cam]['Shannon']
    n_profiles = PAI.size
    for pp in range(0,n_profiles):
        f.write('%s, %s, %.0f, %.5f, %.5f, %.5f, %.0f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f , %.5f\n' % (cam+'_'+str(pp+1).zfill(3), cam, scale[pp], density[pp], height[pp], PAI[pp], layers[pp], shape[pp], mean_ht[pp], sd[pp], skew[pp], kurt[pp], fr[pp], Shan[pp]))

f.close()
