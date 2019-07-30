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
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats
import time
import copy

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
#sample_res = np.array([2.,5.,10.,20.,25.,50.,100.])
#keys = ['2m','5m','10m','20m','25m','50m','100m']
sample_res = np.array([5.,10.,20.,50.,100.])
keys = ['5m','10m','20m','50m','100m']
kappa = 0.70
max_k = 2
n_iter = 250

area = 10.**4
#target_point_density = np.array([5., 10., 20., 30., 40.])
target_pulse_density = np.array([25])
keys_2 = ['25']
#target_points = area*target_point_density
target_shots = (area*target_pulse_density).astype('int')

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
mask = plot_coordinates['plot']==plot
affine=lstsq.least_squares_affine_matrix(plot_coordinates['x'][mask],plot_coordinates['y'][mask],plot_coordinates['x_prime'][mask],plot_coordinates['y_prime'][mask])
plot_bbox = np.array(lstsq.apply_affine_transformation(np.array([-20.,120.,120.,-20.]),np.array([120.,120.,-20.,-20.]),affine)).transpose() # simple square bounding box applied for all sensitivity analyses

pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox,
                                                filter_by_first_return_location=True)
pts[pts[:,2]<0,2]=0
n_returns = pts.shape[0]
shots = np.unique(pts[:,6])
#n_shots=shots.size
#target_shots = (np.ceil(target_points*n_shots/float(n_returns))).astype('int')

subplots = {}
PAD_profiles_MH = {}
PAD_profiles_MH_wt = {}
PAD_profiles_rad1 = {}
PAD_profiles_rad2 = {}
penetration_limit = {}

# Loop through all the spatial scales of interest
print("generating sample grid")
for ss in range(0,sample_res.size):
    print(sample_res[ss])
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
            subplot.append( np.asarray(bbox).transpose() )

    subplots[keys[ss]] = subplot

    n_subplots=len(subplot)
    PAD = np.zeros((n_iter,n_subplots,n_layers),dtype='float')
    # set up dictionaries to hold the profiles
    temp_dic={}
    #for dd in range(0,target_points.size):
    for dd in range(0,target_shots.size):
        temp_dic[keys_2[dd]] = PAD.copy()

    # store all profiles in relevant dictionary
    PAD_profiles_MH[keys[ss]] = copy.deepcopy(temp_dic)
    PAD_profiles_MH_wt[keys[ss]] = copy.deepcopy(temp_dic)
    PAD_profiles_rad1[keys[ss]] = copy.deepcopy(temp_dic)#temp_dic.copy()
    PAD_profiles_rad2[keys[ss]] = copy.deepcopy(temp_dic)#temp_dic.copy()
    penetration_limit[keys[ss]] = copy.deepcopy(temp_dic)#temp_dic.copy()

PAD = None
temp_dic = None

#-----------------------------------------------------
# now do the sensitivity analysis
# Loop through all the target point densities
print("sensitivity analysis")
#for dd in range(0,target_points.size):
for dd in range(0,target_shots.size):
    #print('target point density = ', target_point_density[dd])
    print('target pulse density = ', target_pulse_density[dd])
    # iterate through all iterations, so that we sample point cloud minimum number of times
    for ii in range(0,n_iter):
        start_time = time.time()
        #print('iteration ', ii+1,'/',n_iter, ' for point density ', dd+1,' of ', target_points.size)
        print('iteration ', ii+1,'/',n_iter, ' for point density ', dd+1,' of ', target_shots.size)
        # subsample the point cloud - this is tricky because we are sampling with replacement, but for
        # every sample there could be up to kmax returns!  All these returns need to be pulled out so
        # that we are resampling by pulse, not point

        shots_iter = np.random.choice(shots,size=target_shots[dd])
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

        #pts_iter=pts.copy()
        #-------------------
        starting_ids, trees = io.create_KDTree(pts_iter)
        # loop through each sampling resolution
        for ss in range(0,sample_res.size):
            print('\t - sample res = ', keys[ss])
            n_subplots = len(subplots[keys[ss]])
            # for each of the subplots, clip point cloud and model PAD and get the metrics
            for pp in range(0,n_subplots):
                # query the tree to locate points of interest
                # note that we will only have one tree for the number of points in sensitivity analysis
                centre_x = np.mean(subplots[keys[ss]][pp][0:4,0])
                centre_y = np.mean(subplots[keys[ss]][pp][0:4,1])
                radius = np.sqrt(sample_res[ss]**2/2.)
                ids = trees[0].query_ball_point([centre_x,centre_y], radius)
                sp_pts = lidar.filter_lidar_data_by_polygon(pts_iter[ids],subplots[keys[ss]][pp],
                                                                filter_by_first_return_location=True)
                #------
                if np.sum(sp_pts[:,3]==1)>0: # check for returns within this column
                    heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
                    mh_profile = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
                    pen_limit = np.cumsum(first_return_profile)==0
                    #------
                    heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts, max_height, layer_thickness)
                    mh_wt_profile = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,n_ground_returns,layer_thickness,kappa,zero_nodata=False)
                    #------
                    u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical',test_sensitivity=True)
                    rad1_profile=u[::-1][1:].copy()
                    #------
                    u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical',test_sensitivity=True)
                    rad2_profile=u[::-1][1:].copy()
                else:
                    mh_profile      = np.zeros(heights.size)*np.nan
                    mh_wt_profile   = np.zeros(heights.size)*np.nan
                    rad1_profile    = np.zeros(heights.size)*np.nan
                    rad2_profile    = np.zeros(heights.size)*np.nan
                    pen_limit       = np.ones(heights.size)

                # adaptive neighbourhood expansion to account for penetration
                # limitations, particularly for fine grids/low point densities
                nodata_test = np.any((np.any(~np.isfinite(mh_profile[2:])),
                                      np.any(~np.isfinite(rad1_profile[2:])),
                                      np.any(~np.isfinite(rad2_profile[2:]))))

                while np.all((nodata_test,radius<=50)):
                    # expand neighbourhood for point cloud sample
                    radius+=2.5
                    sp_pts_iter = pts_iter[trees[0].query_ball_point([centre_x,centre_y], radius)]

                    # recalculate profiles at coarser resolution and use these
                    # for gap-filling
                    nodata_gaps = np.isnan(rad1_profile)
                    u,n,I,U = LAD2.calculate_LAD(sp_pts_iter,heights_rad,max_k,'spherical',test_sensitivity=True)
                    rad1_profile[nodata_gaps]=u[::-1][1:][nodata_gaps]

                    # now repeat but for adjusted profiles, accounting for imperfect penetration of LiDAR pulses into canopy
                    nodata_gaps = np.isnan(rad2_profile)
                    if sample_res[ss]>=10:
                        u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts_iter,heights_rad,max_k,'spherical',test_sensitivity=True)
                    else:
                        u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts_iter,heights_rad,max_k,'spherical',test_sensitivity=True)

                    rad2_profile[nodata_gaps]=u[::-1][1:][nodata_gaps]

                    # now get MacArthur-Horn profiles
                    nodata_gaps = ~np.isfinite(mh_profile)
                    heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts_iter, max_height, layer_thickness)
                    mh_profile[nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile,
                                                                            n_ground_returns,
                                                                            layer_thickness,
                                                                            kappa,
                                                                            zero_nodata=False)[nodata_gaps]
                    # and do the same for weighted returns
                    nodata_gaps = ~np.isfinite(mh_wt_profile)
                    heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts_iter, max_height, layer_thickness)
                    mh_wt_profile[nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,
                                                                            n_ground_returns,
                                                                            layer_thickness,
                                                                            kappa,
                                                                            zero_nodata=False)[nodata_gaps]
                    # update check
                    nodata_test = np.any((np.any(~np.isfinite(mh_profile[2:])),
                                          np.any(~np.isfinite(rad1_profile[2:])),
                                          np.any(~np.isfinite(rad2_profile[2:]))))

                # Fill dictionaries
                PAD_profiles_MH[keys[ss]][keys_2[dd]][ii,pp,:] = mh_profile.copy()
                PAD_profiles_MH_wt[keys[ss]][keys_2[dd]][ii,pp,:] = mh_wt_profile.copy()
                PAD_profiles_rad1[keys[ss]][keys_2[dd]][ii,pp,:] = rad1_profile.copy()
                PAD_profiles_rad2[keys[ss]][keys_2[dd]][ii,pp,:] = rad2_profile.copy()
                penetration_limit[keys[ss]][keys_2[dd]][ii,pp,:] = pen_limit.copy()
                plt.plot(rad2_profile,'-',color='blue')
        end_time = time.time()
        print('\t loop time = ', end_time - start_time)

np.save("MH_resolution_sensitivity_adaptive_%s.npy" % plot, PAD_profiles_MH)
np.save("MH_resolution_wt_sensitivity_adaptive_%s.npy" % plot, PAD_profiles_MH_wt)
np.save("rad1_resolution_sensitivity_adaptive_%s.npy" % plot, PAD_profiles_rad1)
np.save("rad2_resolution_sensitivity_adaptive_%s.npy" % plot, PAD_profiles_rad2)
np.save("penetration_limit_resolution_adaptive_%s.npy" % plot, penetration_limit)
