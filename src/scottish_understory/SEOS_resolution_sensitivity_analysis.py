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
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats
import time
import copy
import SEOS_sensitivity_analysis_plots as splt

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_list = '../las_lists/carbomap_site1_lastiles.txt'

N=6206000.
S=6205900.
W=344520.
E=344620.

max_height = 40.
layer_thickness = 0.5
heights = np.arange(0.,max_height,layer_thickness)+0.5
heights_rad = np.arange(0,max_height+0.5,layer_thickness)
n_layers = heights.size
plot_width = 100.
sample_res = np.array([2.,5.,10.,20.,50.])
keys = ['2m','5m','10m','20m','50m']
kappa = 1.
max_k = 3
n_iter = 100

area = 100.*100.
target_pulse_density = 1000.
target_shots = area*target_pulse_density

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
plot_bbox = np.asarray([[W,N],[E,N],[E,S],[W,S]])
pts, starting_ids, trees = io.load_lidar_data_by_polygon(las_list,plot_bbox,laz_files=False,max_pts_per_tree = 5*10**5)
N_trees = len(trees)

pts[np.any((pts[:,4]==3,pts[:,4]==4,pts[:,4]==5),axis=0),4]=1
pts[pts[:,2]<0,2]=0
shots = np.unique(pts[:,6]) # gps times
n_returns = pts.shape[0]
n_shots=shots.size

subplots = {}
PAD_profiles_MH = {}
PAD_profiles_rad2 = {}
penetration_limit = {}

# Loop through all the spatial scales of interest
print("generating sample grid")
for ss in range(0,sample_res.size):
    print(sample_res[ss])
    # Now create the subplot grids
    rows = int(plot_width/sample_res[ss])
    cols = int(plot_width/sample_res[ss])
    x=np.arange(0,plot_width+sample_res[ss],sample_res[ss])+W
    y=np.arange(0,plot_width+sample_res[ss],sample_res[ss])+S

    xv,yv=np.asarray(np.meshgrid(x,y))

    count = 0
    subplot = []
    for i in range(0,rows):
        for j in range(0,cols):
            bbox = [ [xv[i,j], xv[i+1,j], xv[i+1,j+1], xv[i,j+1], xv[i,j]],
                     [yv[i,j], yv[i+1,j], yv[i+1,j+1], yv[i,j+1], yv[i,j]] ]
            subplot.append( np.asarray(bbox).transpose() )

    subplots[keys[ss]] = subplot

    n_subplots=len(subplot)
    PAD = np.zeros((n_iter,n_subplots,n_layers),dtype='float')

    # store all profiles in relevant dictionary
    PAD_profiles_MH[keys[ss]] = PAD.copy()
    PAD_profiles_rad2[keys[ss]] = PAD.copy()
    penetration_limit[keys[ss]] = PAD.copy()

PAD = None
temp_dic = None

#-----------------------------------------------------
# now do the sensitivity analysis
# for each of the subplots, clip point cloud and model PAD and get the metrics
# iterate through all iterations, so that we sample point cloud minimum number of times
for ii in range(0,n_iter):
    start_time = time.time()
    print('iteration ', ii+1,'/',n_iter)
    # subsample the point cloud - this is tricky because we are sampling with replacement, but for
    # every sample there could be up to kmax returns!  All these returns need to be pulled out so
    # that we are resampling by pulse, not point
    shots_iter = np.random.choice(shots,size=int(target_shots))
    mask = np.in1d(pts[:,6],shots_iter)
    pts_iter = pts[mask,:]
    #-------------------
    # this chunk of code makes sure that pulses sampled multiple times are properly accounted for
    temp_shots = shots_iter.copy()
    vals, idx_start= np.unique(temp_shots, return_index=True)
    temp_shots = np.delete(temp_shots,idx_start)
    count =1
    while temp_shots.size>0:
        count+=1
        mask = np.in1d(pts[:,6],temp_shots)
        pts_iter = np.concatenate((pts_iter,pts[mask,:]))
        vals, idx_start= np.unique(temp_shots, return_index=True)
        temp_shots = np.delete(temp_shots,idx_start)
    #-------------------
    starting_ids, trees = io.create_KDTree(pts_iter)
    N_trees = len(trees)
    # for each of the subplots, clip point cloud and model PAD and get the metrics
    for ss,res in enumerate(sample_res):
        n_subplots = len(subplots[keys[ss]])
        for pp in range(0,n_subplots):
            # query the tree to locate points of interest
            # note that we will only have one tree for the number of points in sensitivity analysis
            centre_x = np.mean(subplots[keys[ss]][pp][0:4,0])
            centre_y = np.mean(subplots[keys[ss]][pp][0:4,1])
            radius = np.sqrt(res**2/2.)
            # retrieve point clouds samples
            sp_pts = np.array([])
            for tt in range(0,N_trees):
                ids = trees[tt].query_ball_point([centre_x,centre_y], radius)
                if len(ids)>0:
                    if sp_pts.size==0:
                        sp_pts = lidar.filter_lidar_data_by_polygon(pts_iter[np.asarray(ids)+starting_ids[tt]],subplots[keys[ss]][pp])
                    else:
                        sp_iter = lidar.filter_lidar_data_by_polygon(pts_iter[np.asarray(ids)+starting_ids[tt]],subplots[keys[ss]][pp])
                        sp_pts = np.concatenate((sp_pts,sp_iter),axis=0)
                        sp_iter = None
            #------
            if sp_pts.size==0:
                PAD_profiles_MH[keys[ss]][ii,pp,:] = np.nan
                PAD_profiles_rad2[keys[ss]][ii,pp,:]=np.nan
                penetration_limit[keys[ss]][ii,pp,:] = 1.
            elif np.sum(sp_pts[:,3]==1)>0:
                heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
                PAD_profiles_MH[keys[ss]][ii,pp,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
                penetration_limit[keys[ss]][ii,pp,:] = np.cumsum(first_return_profile)==0
                #------
                u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
                PAD_profiles_rad2[keys[ss]][ii,pp,:]=u[::-1][1:].copy()
            else:
                PAD_profiles_MH[keys[ss]][ii,pp,:] = np.nan
                PAD_profiles_rad2[keys[ss]][ii,pp,:]=np.nan
                penetration_limit[keys[ss]][ii,pp,:] = 1.

    end_time = time.time()
    print('\t loop time = ', end_time - start_time)

np.save("SEOS_MH_sensitivity_resolution.npy", PAD_profiles_MH)
np.save("SEOS_rad2_sensitivity_resolution.npy", PAD_profiles_rad2)
np.save("SEOS_penetration_limit_resolution.npy", penetration_limit)

splt.plot_profile_sensitivity_resolution(1,'SEOS_sensitivity_analysis_resolution.png',heights,PAD_profiles_MH)
