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

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
las_file = 'Carbon_plot_point_cloud_buffer.las'

gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('S32','f16','f16','f16','f16')}
plot_coordinates = np.genfromtxt(gps_pts_file, skiprows = 1, delimiter = ',',dtype=datatype)

plot = 'Belian'
alt_plot = 'maliau_belian'

max_height = 80.
layer_thickness = 1.
heights = np.arange(0.,max_height)+1
heights_rad = np.arange(0,max_height+1)
n_layers = heights.size
plot_width = 100.
sample_res = np.array([5.,10.,20.,25.,50.,100.])
keys = ['5m','10m','20m','25m','50m','100m']
kappa = 0.72
max_k = 2
n_iter = 250

area = 10.**4
target_point_density = np.array([5., 10., 15., 20., 25., 30., 40.])
keys_2 = ['5','10','15','20','25','30','40']
target_points = area*target_point_density

#---------------------------------------------------------------------------------------------------------------
# Load the data etc.
mask = plot_coordinates['plot']==alt_plot
affine=lstsq.least_squares_affine_matrix(plot_coordinates['x'][mask],plot_coordinates['y'][mask],plot_coordinates['x_prime'][mask],plot_coordinates['y_prime'][mask])
plot_bbox = np.array(lstsq.apply_affine_transformation(plot_coordinates['x'][mask],plot_coordinates['y'][mask],affine)).transpose()

pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox)
n_returns = pts.shape[0]
shots = np.unique(pts[:,-1]) 
n_shots=shots.size
target_shots = (np.ceil(target_points*n_shots/float(n_returns))).astype('int')
shots = np.unique(pts[:,-1]) 

subplots = {}
PAD_profiles_MH = {}
PAD_profiles_rad1 = {}
PAD_profiles_rad2 = {}

# Loop through all the spatial scales of interest
for ss in range(0,sample_res.size):
    print sample_res[ss]
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
    PAD = np.zeros((n_iter,n_subplots,n_layers))
    # set up dictionaries to hold the profiles
    pt_dens_MH={}
    pt_dens_rad1={}
    pt_dens_rad2={}
    for dd in range(0,target_points.size):       
        pt_dens_MH[keys_2[dd]] = PAD.copy() 
        pt_dens_rad1[keys_2[dd]] = PAD.copy()
        pt_dens_rad2[keys_2[dd]] = PAD.copy()

    # store all profiles in relevant dictionary
    PAD_profiles_MH[keys[ss]] = pt_dens_MH.copy()
    PAD_profiles_rad1[keys[ss]] = pt_dens_rad1.copy()
    PAD_profiles_rad2[keys[ss]] = pt_dens_rad2.copy()
    

PAD = None

#-----------------------------------------------------
# now do the sensitivity analysis 
# Loop through all the target point densities
for dd in range(0,target_points.size):
    print 'target point density = ', target_point_density[dd]
    # iterate through all iterations, so that we sample point cloud minimum number of times
    for ii in range(0,n_iter):
        start_time = time.time()
        print 'iteration ', ii+1,'/',n_iter, ' for point density ', dd+1,' of ', target_points.size 
        # subsample the point cloud - this is tricky because we are sampling with replacement, but for
        # every sample there could be up to kmax returns!  All these returns need to be pulled out so
        # that we are resampling by pulse, not point
        shots_iter = np.random.choice(shots,size=target_shots[dd])
        mask = np.in1d(pts[:,-1],shots_iter)
        pts_iter = pts[mask,:]

        #-------------------
        # this chunk of code makes sure that pulses sampled multiple times are properly accounted for 
        temp_shots = shots_iter.copy()
        vals, idx_start= np.unique(temp_shots, return_index=True)
        temp_shots = np.delete(temp_shots,idx_start)
        count =1
        while temp_shots.size>0:
            count+=1
            mask = np.in1d(pts[:,-1],temp_shots)
            pts_iter = np.concatenate((pts_iter,pts[mask,:]))
            vals, idx_start= np.unique(temp_shots, return_index=True)
            temp_shots = np.delete(temp_shots,idx_start)
        #-------------------
        
        #print '\t\t', target_point_density[dd],'-' , pts_iter.shape[0]/100.**2 , count
        starting_ids, trees = io.create_KDTree(pts_iter) 
        # loop through each sampling resolution 
        for ss in range(0,sample_res.size):
            print '\t - sample res = ', keys[ss]
            n_subplots = len(subplots[keys[ss]])
            # for each of the subplots, clip point cloud and model PAD and get the metrics
            for pp in range(0,n_subplots):
                # query the tree to locate points of interest
                # note that we will only have one tree for number of points in sensitivity analysis  
                centre_x = np.mean(subplots[keys[ss]][pp][0:4,0])
                centre_y = np.mean(subplots[keys[ss]][pp][0:4,1])
                radius = np.sqrt(sample_res[ss]**2/2.)              
                ids = trees[0].query_ball_point([centre_x,centre_y], radius)
                sp_pts = lidar.filter_lidar_data_by_polygon(pts_iter[ids],subplots[keys[ss]][pp])
                #------
                heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
                PAD_profiles_MH[keys[ss]][keys_2[dd]][ii,pp,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
                #------
                u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical')
                PAD_profiles_rad1[keys[ss]][keys_2[dd]][ii,pp,:]=u[::-1][1:].copy()
                #------
                u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
                PAD_profiles_rad2[keys[ss]][keys_2[dd]][ii,pp,:]=u[::-1][1:].copy()
        end_time = time.time()
        print '\t loop time = ', end_time - start_time


np.save("MH_sensitivity.npy", PAD_profiles_MH)
np.save("rad1_sensitivity.npy", PAD_profiles_rad1)
np.save("rad2_sensitivity.npy", PAD_profiles_rad2)

# now plot the layers of interest           
# 1a) the mean profile from all iterations with 50% and 95% CI
PAD_profiles_MH = np.load("MH_sensitivity.npy")[()]
PAD_profiles_rad1 = np.load("rad1_sensitivity.npy")[()]
PAD_profiles_rad2 = np.load("rad2_sensitivity.npy")[()]


colour = ['#46E900','#1A2BCE','#E0007F']

plt.figure(1, facecolor='White',figsize=[9,9])
ax1a = plt.subplot2grid((2,5),(0,0))
ax1a.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1a.set_ylabel('height / m',fontsize=axis_size)
ax1a.annotate('a - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1b = plt.subplot2grid((2,5),(0,1),sharex=ax1a,sharey=ax1a)
ax1b.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1b.annotate('b - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1c = plt.subplot2grid((2,5),(0,2),sharex=ax1a,sharey=ax1a)
ax1c.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1c.annotate('c - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                       
ax1d = plt.subplot2grid((2,5),(0,3),sharex=ax1a,sharey=ax1a)
ax1d.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1d.annotate('d - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                       
ax1e = plt.subplot2grid((2,5),(0,4),sharex=ax1a,sharey=ax1a)
ax1e.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1e.annotate('e - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                            
# loop through the subplots
sp = [ax1a,ax1b,ax1c,ax1d,ax1e]
pkeys = ['5m','10m','20m','50m','100m']
for pp in range(0,len(sp)):
  rad2_sp = np.mean(PAD_profiles_rad2[pkeys[pp]]['40'],axis=1)  
  rad2 = np.median(rad2_sp,axis=0)
  rad2perc = np.percentile(rad2_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[3][2:],color=colour[2],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad2perc[1][2:],rad2perc[2][2:],color=colour[2],alpha=0.3)
 
  rad1_sp = np.mean(PAD_profiles_rad1[pkeys[pp]]['40'],axis=1)  
  rad1 = np.median(rad1_sp,axis=0)
  rad1perc = np.percentile(rad1_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[3][2:],color=colour[1],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad1perc[1][2:],rad1perc[2][2:],color=colour[1],alpha=0.3)
 
  MH_sp = np.mean(PAD_profiles_MH[pkeys[pp]]['40'],axis=1)  
  MH = np.median(MH_sp,axis=0)
  MHperc = np.percentile(MH_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[3][2:],color=colour[0],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],MHperc[1][2:],MHperc[2][2:],color=colour[0],alpha=0.3)

  if pp==2:
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)') 
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)
  else:                
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

# 1b) the widths of the 50% (solid) and 95% (dashed) confidence intervals of each vertical layer
# this gives an indication of the reliability of the profile at a given canopy depth
ax1f = plt.subplot2grid((2,5),(1,0),sharey=ax1a)
ax1f.set_xlabel('CI / %',fontsize=axis_size)
ax1f.set_ylabel('height / m',fontsize=axis_size)
ax1f.annotate('f - 5 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1g = plt.subplot2grid((2,5),(1,1),sharex=ax1f,sharey=ax1a)
ax1g.set_xlabel('CI / %',fontsize=axis_size)
ax1g.annotate('g - 10 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1h = plt.subplot2grid((2,5),(1,2),sharex=ax1f,sharey=ax1a)
ax1h.set_xlabel('CI / %',fontsize=axis_size)
ax1h.annotate('h - 20 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
                                                       
ax1i = plt.subplot2grid((2,5),(1,3),sharex=ax1f,sharey=ax1a)
ax1i.set_xlabel('CI / %',fontsize=axis_size)
ax1i.annotate('i - 50 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
                                                       
ax1j = plt.subplot2grid((2,5),(1,4),sharex=ax1f,sharey=ax1a)
ax1j.set_xlabel('CI / %',fontsize=axis_size)
ax1j.annotate('j - 100 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

sp = [ax1f,ax1g,ax1h,ax1i,ax1j]
for pp in range(0,len(sp)):
  MHl = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],2.5,axis=0)
  MH25 = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],25,axis=0)
  MH75 = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],75,axis=0)
  MHu = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],97.5,axis=0)
  MH_95CI_i = ((MHu-MHl)/np.mean(PAD_profiles_MH[pkeys[pp]]['40'],axis=0))*100.
  MH_50CI_i = ((MH75-MH25)/np.mean(PAD_profiles_MH[pkeys[pp]]['40'],axis=0))*100.
  MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
  MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
  sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
  sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
  if pp==2:
      sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='MacArthur-Horn 95% CI')
      sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='MacArthur-Horn 50% CI')
  else:
      sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
      sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
  
  rad1l = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],2.5,axis=0)
  rad125 = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],25,axis=0)
  rad175 = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],75,axis=0)
  rad1u = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],97.5,axis=0)
  rad1_95CI_i = ((rad1u-rad1l)/np.mean(PAD_profiles_rad1[pkeys[pp]]['40'],axis=0))*100.
  rad1_50CI_i = ((rad175-rad125)/np.mean(PAD_profiles_rad1[pkeys[pp]]['40'],axis=0))*100.
  rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
  rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
  if pp==2:
      sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
      sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
  else:
      sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
      sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

  rad2l = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],2.5,axis=0)
  rad225 = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],25,axis=0)
  rad275 = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],75,axis=0)
  rad2u = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],97.5,axis=0)
  rad2_95CI_i = ((rad2u-rad2l)/np.mean(PAD_profiles_rad2[pkeys[pp]]['40'],axis=0))*100.
  rad2_50CI_i = ((rad275-rad225)/np.mean(PAD_profiles_rad2[pkeys[pp]]['40'],axis=0))*100.
  rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
  rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
  if pp==2:
      sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
      sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
      lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
  else:
      sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
      sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)
                      
ax1a.set_ylim(0,80)
ax1a.set_xlim(0,0.7)
ax1f.set_xlim(0,520)
ax1a.locator_params(axis='x',nbins=5)
ax1b.locator_params(axis='x',nbins=5)
ax1c.locator_params(axis='x',nbins=5)
ax1d.locator_params(axis='x',nbins=5)
ax1e.locator_params(axis='x',nbins=5)
ax1f.locator_params(axis='x',nbins=5)
ax1g.locator_params(axis='x',nbins=5)
ax1h.locator_params(axis='x',nbins=5)
ax1i.locator_params(axis='x',nbins=5)
ax1j.locator_params(axis='x',nbins=5)

yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + ax1g.get_yticklabels() + ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.subplots_adjust(hspace=0.4, wspace = 0.1, bottom = 0.2)
plt.savefig('PAD_resolution_sensitivity.png')
plt.show()


#==============================================================================================            
# 2) Equivalent plots but this time using different point densities at target resolution          
pkeys=['5','10','15','20','25','30','40']
plt.figure(2, facecolor='White',figsize=[9,9])
ax2a = plt.subplot2grid((2,5),(0,0))
ax2a.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2a.set_ylabel('height / m',fontsize=axis_size)
ax2a.annotate('a - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2b = plt.subplot2grid((2,5),(0,1),sharex=ax2a,sharey=ax2a)
ax2b.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2b.annotate('b - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2c = plt.subplot2grid((2,5),(0,2),sharex=ax2a,sharey=ax2a)
ax2c.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2c.annotate('c - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                       
ax2d = plt.subplot2grid((2,5),(0,3),sharex=ax2a,sharey=ax2a)
ax2d.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2d.annotate('d - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                       
ax2e = plt.subplot2grid((2,5),(0,4),sharex=ax2a,sharey=ax2a)
ax2e.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2e.annotate('e - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
                                                            
# loop through the subplots
sp = [ax2a,ax2b,ax2c,ax2d,ax2e]
pkeys = ['5','10','20','30','40']
for pp in range(0,len(sp)):
  rad2_sp = np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=1)  
  rad2 = np.median(rad2_sp,axis=0)
  rad2perc = np.percentile(rad2_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[3][2:],color=colour[2],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad2perc[1][2:],rad2perc[2][2:],color=colour[2],alpha=0.3)
 
  rad1_sp = np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=1)  
  rad1 = np.median(rad1_sp,axis=0)
  rad1perc = np.percentile(rad1_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[3][2:],color=colour[1],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad1perc[1][2:],rad1perc[2][2:],color=colour[1],alpha=0.3)
 
  MH_sp = np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=1)  
  MH = np.median(MH_sp,axis=0)
  MHperc = np.percentile(MH_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[3][2:],color=colour[0],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],MHperc[1][2:],MHperc[2][2:],color=colour[0],alpha=0.3)

  if pp==2:
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)') 
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
  else:                
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

# 1b) the widths of the 50% (solid) and 95% (dashed) confidence intervals of each vertical layer
# this gives an indication of the reliability of the profile at a given canopy depth
ax2f = plt.subplot2grid((2,5),(1,0),sharey=ax2a)
ax2f.set_xlabel('CI / %',fontsize=axis_size)
ax2f.set_ylabel('height / m',fontsize=axis_size)
ax2f.annotate('f - 5 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2g = plt.subplot2grid((2,5),(1,1),sharex=ax2f,sharey=ax2a)
ax2g.set_xlabel('CI / %',fontsize=axis_size)
ax2g.annotate('g - 10 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2h = plt.subplot2grid((2,5),(1,2),sharex=ax2f,sharey=ax2a)
ax2h.set_xlabel('CI / %',fontsize=axis_size)
ax2h.annotate('h - 20 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
                                                       
ax2i = plt.subplot2grid((2,5),(1,3),sharex=ax2f,sharey=ax2a)
ax2i.set_xlabel('CI / %',fontsize=axis_size)
ax2i.annotate('i - 30 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
                                                       
ax2j = plt.subplot2grid((2,5),(1,4),sharex=ax2f,sharey=ax2a)
ax2j.set_xlabel('CI / %',fontsize=axis_size)
ax2j.annotate('j - 40 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

sp = [ax2f,ax2g,ax2h,ax2i,ax2j]
for pp in range(0,len(sp)):
    MHl = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],2.5,axis=0)
    MH25 = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],25,axis=0)
    MH75 = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],75,axis=0)
    MHu = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],97.5,axis=0)
    MH_95CI_i = ((MHu-MHl)/np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
    MH_50CI_i = ((MH75-MH25)/np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
    MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
    MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='MacArthur-Horn 95% CI')
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='MacArthur-Horn 50% CI')
    else:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
  
    rad1l = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],2.5,axis=0)
    rad125 = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],25,axis=0)
    rad175 = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],75,axis=0)
    rad1u = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],97.5,axis=0)
    rad1_95CI_i = ((rad1u-rad1l)/np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=0))*100.
    rad1_50CI_i = ((rad175-rad125)/np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=0))*100.
    rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
    rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
    else:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    rad2l = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],2.5,axis=0)
    rad225 = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],25,axis=0)
    rad275 = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],75,axis=0)
    rad2u = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],97.5,axis=0)
    rad2_95CI_i = ((rad2u-rad2l)/np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=0))*100.
    rad2_50CI_i = ((rad275-rad225)/np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=0))*100.
    rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
    rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
        sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
    else:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)
                      
ax2a.set_ylim(0,80)
ax2a.set_xlim(0,0.7)
ax2f.set_xlim(0,520)
ax2a.locator_params(axis='x',nbins=5)
ax2b.locator_params(axis='x',nbins=5)
ax2c.locator_params(axis='x',nbins=5)
ax2d.locator_params(axis='x',nbins=5)
ax2e.locator_params(axis='x',nbins=5)
ax2f.locator_params(axis='x',nbins=5)
ax2g.locator_params(axis='x',nbins=5)
ax2h.locator_params(axis='x',nbins=5)
ax2i.locator_params(axis='x',nbins=5)
ax2j.locator_params(axis='x',nbins=5)

yticklabels = ax2b.get_yticklabels() + ax2c.get_yticklabels() + ax2d.get_yticklabels() + ax2e.get_yticklabels() + ax2g.get_yticklabels() + ax2h.get_yticklabels() + ax2i.get_yticklabels() + ax2j.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.subplots_adjust(hspace=0.4, wspace = 0.1, bottom = 0.2)
plt.savefig('PAD_point_density_sensitivity.png')
plt.show()

                                                                                                  
# 3) PAI vs resolution at different shot spacings
N_res = len(keys)
plt.figure(3, facecolor='White',figsize=[9,9])
ax3a = plt.subplot2grid((3,3),(0,0))
ax3a.set_ylabel('PAI',fontsize=axis_size)
ax3a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3b = plt.subplot2grid((3,3),(0,1),sharex=ax3a)
ax3b.set_title('Point density = 5 pts m$^{-2}$', fontsize=10)
ax3b.annotate('b - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3c = plt.subplot2grid((3,3),(0,2),sharex=ax3a)
ax3c.annotate('c - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

boxplot_a=[]
boxplot_b=[]
boxplot_c=[]
for kk in range(0,N_res):
    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['5'],axis=1)[:,2:],axis=1)  
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['5'],axis=1)[:,2:],axis=1)  
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['5'],axis=1)[:,2:],axis=1)  
    boxplot_a.append(PAI_MH.copy())
    boxplot_b.append(PAI_rad1.copy())
    boxplot_c.append(PAI_rad2.copy())

bpa = ax3a.boxplot(boxplot_a, patch_artist=True)
bpb = ax3b.boxplot(boxplot_b, patch_artist=True)
bpc = ax3c.boxplot(boxplot_c, patch_artist=True)

bps=[bpa,bpb,bpc]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)



ax3d = plt.subplot2grid((3,3),(1,0),sharex=ax3a,sharey=ax3a)
ax3d.set_ylabel('PAI',fontsize=axis_size)
ax3d.annotate('d - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3e = plt.subplot2grid((3,3),(1,1),sharex=ax3a,sharey=ax3b)
ax3e.set_title('Point density = 20 pts m$^{-2}$', fontsize=10)
ax3e.annotate('e - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3f = plt.subplot2grid((3,3),(1,2),sharex=ax3a,sharey=ax3c)
ax3f.annotate('f - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

boxplot_d=[]
boxplot_e=[]
boxplot_f=[]
for kk in range(0,N_res):

    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['20'],axis=1)[:,2:],axis=1)  
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['20'],axis=1)[:,2:],axis=1)  
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['20'],axis=1)[:,2:],axis=1)  
    boxplot_d.append(PAI_MH.copy())
    boxplot_e.append(PAI_rad1.copy())
    boxplot_f.append(PAI_rad2.copy())

bpd = ax3d.boxplot(boxplot_d, patch_artist=True)
bpe = ax3e.boxplot(boxplot_e, patch_artist=True)
bpf = ax3f.boxplot(boxplot_f, patch_artist=True)

bps=[bpd,bpe,bpf]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)

ax3g = plt.subplot2grid((3,3),(2,0),sharex=ax3a,sharey=ax3a)
ax3g.set_ylabel('PAI',fontsize=axis_size)
ax3g.annotate('g - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3g.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

ax3h = plt.subplot2grid((3,3),(2,1),sharex=ax3a,sharey=ax3b)
ax3h.set_title('Point density = 40 pts m$^{-2}$', fontsize=10)
ax3h.annotate('h - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3h.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

ax3i = plt.subplot2grid((3,3),(2,2),sharex=ax3a,sharey=ax3c)
ax3i.annotate('i - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3i.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

boxplot_g=[]
boxplot_h=[]
boxplot_i=[]
for kk in range(0,N_res):

    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['40'],axis=1)[:,2:],axis=1)  
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['40'],axis=1)[:,2:],axis=1)  
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['40'],axis=1)[:,2:],axis=1)  
    boxplot_g.append(PAI_MH.copy())
    boxplot_h.append(PAI_rad1.copy())
    boxplot_i.append(PAI_rad2.copy())

bpg = ax3g.boxplot(boxplot_g, patch_artist=True)
bph = ax3h.boxplot(boxplot_h, patch_artist=True)
bpi = ax3i.boxplot(boxplot_i, patch_artist=True)

bps=[bpg,bph,bpi]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)

x_locs = [1,2,3,4,5,6]
ax3c.set_xticks(x_locs)
xticks=ax3c.get_xticks().tolist()
xticks=keys
ax3c.set_xticklabels(xticks,fontsize=axis_size)

xticklabels = ax3a.get_xticklabels() + ax3b.get_xticklabels() + ax3c.get_xticklabels() + ax3d.get_xticklabels() + ax3e.get_xticklabels() + ax3f.get_xticklabels()

yticklabels = ax3h.get_yticklabels() + ax3b.get_yticklabels() + ax3c.get_yticklabels() + ax3i.get_yticklabels() + ax3e.get_yticklabels() + ax3f.get_yticklabels()
ax3a.set_ylim(5,11)
ax3b.set_ylim(2,10)
ax3c.set_ylim(6,20)

plt.setp(xticklabels,visible=False)
plt.subplots_adjust(wspace = 0.4,hspace=0.2)
plt.savefig('PAI_sensitivity.png')

plt.show()
