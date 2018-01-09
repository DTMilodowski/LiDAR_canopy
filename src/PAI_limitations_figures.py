#===============================================================================
# PAI_limitations_figures.py
# D.T.Milodowski, November 2017
#-------------------------------------------------------------------------------
# This function contains the scripts used to produce the figures in the paper:
# "Point density imposes limitations on PAI estimation from discrete return
# LiDAR"
#-------------------------------------------------------------------------------
import numpy as np
import fiona

# import plotting libraries
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.ticker as plticker

import sys
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.set_cmap(cmaps.viridis)

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 16
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

# code to get trio of nice colourblind friendly colours 
cmap = cm.get_cmap('plasma')
scale = np.arange(0.,3.)
scale /=2.5
colour = cmap(scale)

#colour = ['#46E900','#1A2BCE','#E0007F']

# import required LiDAR libaries
import LiDAR_MacHorn_LAD_profiles as MH
import LiDAR_io as lidar_io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import raster_io as io

# Directory listings
SAVEDIR = '~/'

# Load files
dens_file = 'arrays_for_lombok/SAFE_pointcloud_metrics_10m_point_density_data.tif'
PAI_file = 'arrays_for_lombok/SAFE_pointcloud_metrics_10m_pai_data.tif'

dens, geo, coord = io.load_GeoTIFF_band_and_georeferencing(dens_file)
PAI, geo, coord = io.load_GeoTIFF_band_and_georeferencing(PAI_file)
dens[np.isnan(PAI)]=np.nan

coords = np.array([[565042,522612],[572092,520634],[549913,514144]])
labels = ['A','B','C']

rows,cols = PAI.shape
N = geo[3]
S = geo[3] + (rows+1)*geo[5]
W = geo[0]
E = geo[0] + (cols+1)*geo[1]
Y, X = np.mgrid[slice(S,N,-geo[5]),slice(W,E,geo[1])]


#-------------------------------------------------------------------------------
# Get analytical solution
k = 0.7
theta_deg = np.asarray([0, 5, 10])
theta_rad = theta_deg*np.pi/180.
radius = np.asarray([10.,20.,05.])
area = np.pi*radius**2
dens_a = np.arange(0.01,np.nanmax(dens),0.01)

PAImax = MH.calculate_analytical_limit(dens,area[0],k,theta_rad[0])

PAImax_10m_00deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[0])
PAImax_20m_00deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[0])
PAImax_05m_00deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[0])
"""
PAImax_10m_05deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[1])
PAImax_20m_05deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[1])
PAImax_30m_05deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[1])

PAImax_10m_10deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[2])
PAImax_20m_10deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[2])
PAImax_30m_10deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[2])
"""
#-------------------------------------------------------------------------------
# Figure 2 - this figure illustrates the analytical solution presented in this
# paper describing the threshold PAI that can be detected using discrete return
# LiDAR
fig = plt.figure(2, facecolor='White',figsize=[12,6])
ax1a= plt.subplot2grid((1,7),(0,0),colspan=3)
ax1a.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2)
ax1a.set_xlabel('point density / pts m$^{-2}$',fontsize = axis_size)
ax1a.set_ylabel('PAI$_{max}$',fontsize = axis_size)

ax1a.plot(dens_a,PAImax_20m_00deg,dashes=[8, 5],c=colour[2],label = '%.3f ha' % (np.pi*radius[1]**2/10.**4))
ax1a.plot(dens_a,PAImax_10m_00deg,dashes=[16, 5],c=colour[1],label = '%.3f ha' % (np.pi*radius[0]**2/10.**4))
ax1a.plot(dens_a,PAImax_05m_00deg,'-',c=colour[0],label = '%.3f ha' % (np.pi*radius[2]**2/10.**4))
ax1a.legend(loc='lower right')

ax1b= plt.subplot2grid((1,7),(0,3),colspan=3,sharex=ax1a,sharey=ax1a)
ax1b.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2, color='white')
ax1b.set_xlabel('point density / pts m$^{-2}$',fontsize = axis_size)
ax1b.set_ylabel('PAI',fontsize = axis_size)

hb = ax1b.hexbin(dens.reshape(dens.size), PAI.reshape(PAI.size), gridsize=(1000,200), bins='log', cmap='plasma')

ax1c= plt.subplot2grid((1,7),(0,6))
cb = fig.colorbar(hb, cax=ax1c)
cb.set_label('log$_{10}$(Number of grid cells)',fontsize = axis_size)

ax1b.plot(dens_a,PAImax_10m_00deg,'-',c='white',linewidth=2)
ax1b.set_xlim(0,30)
ax1b.set_ylim(0,np.nanmax(PAI))
plt.tight_layout()

plt.savefig('Fig2_SAFE_point_density_vs_PAI.png')
plt.savefig('Fig2_SAFE_point_density_vs_PAI.pdf')

#-------------------------------------------------------------------------------
# Figure 1 - this figure presents maps of PAI and point density across the SAFE
# landscape

# Load shapefiles
shapefile_dir = '../Data/Fig1_Shapefiles/'
land_cover_file = 'HCS_Stratification_DOI.shp'
vjr_file = 'VJR_prj.shp'
ea_file = 'EA_prj.shp'

land_cover = fiona.open(shapefile_dir+land_cover_file)
vjr = fiona.open(shapefile_dir+vjr_file)             
ea = fiona.open(shapefile_dir+ea_file)

patches = []
colours = ['#30A000','#85EE58','#E756A8','0.75']

for poly in land_cover:
  color_iter = colours[poly['properties']['HCS_CLass']-1]
  if poly['geometry']['type']=='MultiPolygon':
    Npoly = len(poly['geometry']['coordinates'][1])
    for nn in range(0,Npoly):
      polygon = Polygon(poly['geometry']['coordinates'][1][nn], True,ec='None',fc=color_iter)
      patches.append(polygon)
  else:  
    polygon = Polygon(poly['geometry']['coordinates'][0], True,ec='None',fc=color_iter)
    patches.append(polygon)
    

VJR_poly = vjr[0]
polygon = Polygon(VJR_poly['geometry']['coordinates'][0], True,ec='#1A2CCE',fc='None')
patches.append(polygon)


ea_poly = ea[0]
polygon = Polygon(ea_poly['geometry']['coordinates'][0], True,ec='#1A2CCE',fc='None')
patches.append(polygon)

coords_trans = coords-np.array([W,S])/float(geo[1])+np.array([W,S])

fig = plt.figure(2, facecolor='White',figsize=[12,12]) 
loc_x = plticker.MultipleLocator(base=10**4) 
loc_y = plticker.MultipleLocator(base=10**4) 
loc_cb = plticker.MultipleLocator(base=10) 
loc_cc = plticker.MultipleLocator(base=4) 
loc_cd = plticker.MultipleLocator(base=0.05) 

ax2a= plt.subplot2grid((2,2),(0,0))
ax2a.yaxis.set_major_locator(loc_y)
ax2a.xaxis.set_major_locator(loc_x)
p2a = PatchCollection(patches, match_original=True)
ax2a.annotate('a - Land cover', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
ax2a.add_collection(p2a) 
ax2a.set_aspect('equal', adjustable='box-forced')     
ax2a.set_xlim(xmin=W,xmax=E)
ax2a.set_ylim(ymin=S,ymax=N)
divider = make_axes_locatable(ax2a)
ax_cb = divider.new_horizontal(size="5%", pad=0.05, pack_start=False)
ax2a.annotate('VJR', xy=(560437,516426), xycoords='data',backgroundcolor='none',horizontalalignment='center', verticalalignment='center')
ax2a.annotate('SAFE', xy=(564805,520225), xycoords='data',backgroundcolor='none',horizontalalignment='center', verticalalignment='center')
for tick in ax2a.get_yticklabels():
    tick.set_rotation(90)
#ax2a.set_xticklabels([])       

for pp in range(0,3):
  ax2a.plot(coords[pp,0],coords[pp,1],'o',c='white')
  ax2a.annotate(labels[pp], xy=coords[pp]+250, xycoords='data',color='black')

    
ax2b= plt.subplot2grid((2,2),(0,1))
ax2b.yaxis.set_major_locator(loc_y)
ax2b.xaxis.set_major_locator(loc_x)
ax2b.annotate('b - Point density', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
ax2b.set_xlim(xmin=W,xmax=E)
ax2b.set_ylim(ymin=S,ymax=N)
im2b=ax2b.imshow(dens,vmin=0,vmax=30,cmap='plasma',origin='lower',extent=[W,E,S,N])   
#densm = np.ma.masked_where(np.isnan(dens),dens)
#im2b = ax2b.pcolormesh(X,Y,densm,vmin=0,vmax=30,cmap='plasma')                     
ax2b.axis('image')    
#ax2b.set_xticklabels([])
#ax2b.set_yticklabels([])           
ax2b.yaxis.set_major_locator(loc_y)
ax2b.xaxis.set_major_locator(loc_x)                
for tick in ax2b.get_yticklabels():
    tick.set_rotation(90)

divider2b = make_axes_locatable(ax2b)
cax2b = divider2b.append_axes("right", size="5%", pad=0.05)
cbar2b=plt.colorbar(im2b, cax=cax2b)
cbar2b.ax.set_ylabel('point density / pts m$^{-2}$',fontsize = axis_size)
cbar2b.solids.set_edgecolor("face")             
cbar2b.locator = loc_cb
cbar2b.update_ticks()
"""
for pp in range(0,3):
  ax2b.plot(coords[pp,0],coords[pp,1],'o',c='white')
  ax2b.annotate(labels[pp], xy=coords_trans[pp], xycoords='data', xytext=(2, 2), textcoords='offset points',color='white')
"""
ax2c= plt.subplot2grid((2,2),(1,0))
ax2c.yaxis.set_major_locator(loc_y)
ax2c.xaxis.set_major_locator(loc_x)
ax2c.annotate('c - PAI', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
im2c = ax2c.imshow(PAI,cmap='viridis',origin='lower',extent=[W,E,S,N])
#PAIm = np.ma.masked_where(np.isnan(PAI),PAI)
#im2c = ax2c.pcolormesh(X,Y,PAIm,cmap='viridis')
ax2c.axis('image')                                 
for tick in ax2c.get_yticklabels():
    tick.set_rotation(90)

divider2c = make_axes_locatable(ax2c)
cax2c = divider2c.append_axes("right", size="5%", pad=0.05)
cbar2c=plt.colorbar(im2c, cax=cax2c)
cbar2c.ax.set_ylabel('PAI',fontsize = axis_size)
cbar2c.solids.set_edgecolor("face")             
cbar2c.locator = loc_cc
cbar2c.update_ticks()
"""
for pp in range(0,3):
  ax2c.plot(coords[pp,0],coords[pp,1],'o',c='white')
  ax2c.annotate(labels[pp], xy=coords_trans[pp], xycoords='data', xytext=(2, 2), textcoords='offset points',color='white')
"""

ax2d= plt.subplot2grid((2,2),(1,1))
ax2d.yaxis.set_major_locator(loc_y)
ax2d.xaxis.set_major_locator(loc_x)
ax2d.annotate('d - PAI/PAI$_{max}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
im2d = ax2d.imshow(PAI/PAImax,vmin = 0.85, vmax=1, cmap='plasma',origin='lower',extent=[W,E,S,N])
#proximity = np.ma.masked_where(np.isnan(PAI),PAI/PAImax)
#im2d = ax2d.pcolormesh(X,Y,proximity,vmin=0.85,vmax=1,cmap='plasma')
ax2d.axis('image')
#ax2d.set_yticklabels([])  

divider2d = make_axes_locatable(ax2d)
cax2d = divider2d.append_axes("right", size="5%", pad=0.05)
cbar2d=plt.colorbar(im2d, cax=cax2d)
cbar2d.ax.set_ylabel('PAI/PAI$_{max}$',fontsize = axis_size)
cbar2d.solids.set_edgecolor("face")           
cbar2d.locator = loc_cd
cbar2d.update_ticks()                
for tick in ax2d.get_yticklabels():
    tick.set_rotation(90)
    
plt.tight_layout()
plt.savefig('Fig1_SAFE_point_density_PAI_maps.png')
plt.savefig('Fig1_SAFE_point_density_PAI_maps.pdf')
plt.show()

# Figure 3 - investigating role of point density on PAI estimates made at specific points
# (controlling for other sources of variation)

las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt'
laz_files = False
# Some parameters
min_PAD = 0.1
radius = 10.
area = np.pi*radius**2
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7
n_iterations = 100

target_dens = np.arange(40,0,-0.2)
target_points = (np.ceil(area*target_dens)).astype('int')
n_dens = target_dens.size
PAI_iter = np.zeros((3,n_dens,n_iterations))

sample_pts_collated = []

for pp in range(0,3):
    print "point %i, x = %.0f, y = %.0f" % (pp+1, coords[pp,0], coords[pp,1])
    # define a bounding box around target points to load in point cloud around area of interest
    E = coords[pp,0]+500
    N = coords[pp,1]+500
    W = coords[pp,0]-500
    S = coords[pp,1]-500
    # Read in LiDAR points for region of interest
    polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
    lidar_pts, starting_ids_for_trees, trees = lidar_io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=laz_files)
    N_trees = len(trees)
    # retrieve point clouds samples        
    sample_pts = np.array([])
    for tt in range(0,N_trees):
      ids = trees[tt].query_ball_point(coords[pp], radius)
      if len(ids)>0:
        if sample_pts.size==0:
          sample_pts = lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]
        else:
          sample_pts = np.concatenate((sample_pts,lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]),axis=0)

    sample_pts_collated.append(sample_pts.copy())
    sample_pts = None
          
# Now loop through the points again pulling out the metrics
for pp in range(0,3):
    # If we have the returns, then calculate metric of interest - in
    # this case the PAI
    sample_pts = sample_pts_collated[pp].copy()
    if sample_pts.size > 0:
      # keep only first returns
      sample_pts=sample_pts[sample_pts[:,3]==1,:]
      if sample_pts.size > 0:
        sample_ids = np.arange(sample_pts.shape[0])
        for dd in range(0,n_dens):
          for ii in range(0,n_iterations):
            sample_pts_iter = sample_pts[np.random.choice(sample_ids,size=target_points[dd]),:]
            # calculate PAD profile
            heights,first_return_profile,n_ground_returns = MH.bin_returns(sample_pts_iter, max_height, layer_thickness)
            PADprof = MH.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
            
            # remove lowermost portion of profile
            PAD_iter = PADprof.copy()
            PAD_iter[heights<min_height]=0
            
            PAI_iter[pp,dd,ii] = np.sum(PAD_iter)
              
      else:
        print "no first returns in neighbourhood"
    else:
      print "no returns in neighbourhood"

# Now plot up the results
PAI_A = np.mean(PAI_iter[0,:,:],axis=1)
ulim_A = np.mean(PAI_iter[0,:,:],axis=1)+np.std(PAI_iter[0,:,:],axis=1)
llim_A =  np.mean(PAI_iter[0,:,:],axis=1)-np.std(PAI_iter[0,:,:],axis=1)

PAI_B = np.mean(PAI_iter[1,:,:],axis=1)
ulim_B = np.mean(PAI_iter[1,:,:],axis=1)+np.std(PAI_iter[1,:,:],axis=1)
llim_B =  np.mean(PAI_iter[1,:,:],axis=1)-np.std(PAI_iter[1,:,:],axis=1)

PAI_C = np.mean(PAI_iter[2,:,:],axis=1)
ulim_C = np.mean(PAI_iter[2,:,:],axis=1)+np.std(PAI_iter[2,:,:],axis=1)
llim_C =  np.mean(PAI_iter[2,:,:],axis=1)-np.std(PAI_iter[2,:,:],axis=1)

fig = plt.figure(3, facecolor='White',figsize=[6,6])
ax3= plt.subplot2grid((1,1),(0,0))
ax3.plot(target_dens,PAI_A,'-',c=colour[2],label = 'A')
ax3.fill_between(target_dens,llim_A,ulim_A,color=colour[2],alpha=0.25)
ax3.plot(target_dens,PAI_B,'-',c=colour[1],label = 'B')
ax3.fill_between(target_dens,llim_B,ulim_B,color=colour[1],alpha=0.25)
ax3.plot(target_dens,PAI_C,'-',c=colour[0],label = 'C')
ax3.fill_between(target_dens,llim_C,ulim_C,color=colour[0],alpha=0.25)
ax3.plot(dens_a,PAImax_10m_00deg,dashes=[8,5],c='k',linewidth=2,label = 'limit')
ax3.legend(loc='lower right')
ax3.set_xlim(xmin=0,xmax=40)
ax3.set_xlabel('point density / pts m$^{-2}$', fontsize = axis_size)
ax3.set_ylabel('PAI', fontsize = axis_size)
plt.tight_layout()
plt.savefig('Figure3_point_density_vs_PAI_pointwise.png')
plt.savefig('Figure3_point_density_vs_PAI_pointwise.pdf')
plt.show()

# Figure 4 - investigating role of resolution on PAI estimates made at specific points
# (controlling for other sources of variation)

# first construct sampling grid
ha = 100
div = np.array([10.,8.,6.,5.,4.,3.,2.,1.])
keys = ['10m','12.5m','16.7m','20m','25m','33m','50m','100m']
res = ha/div

coords_00 = coords-ha/2
n_grids = res.size
subplots = []
subplots.append({})
subplots.append({})
subplots.append({})
PAI_res = {}

PAI_mean = np.zeros((3,div.size))
PAI_sd = np.zeros((3,div.size))

for pp in range(0,3):
  for ss in range(0,n_grids):
    x = np.arange(coords_00[pp,0],coords_00[pp,0]+ha+1,res[ss])
    y = np.arange(coords_00[pp,1],coords_00[pp,1]+ha+1,res[ss])
    
    xv,yv=np.asarray(np.meshgrid(x,y))
    
    rr,cc = xv.shape
    rr-=1
    cc-=1
    subplot = []
    for i in range(0,rr):
        for j in range(0,cc):
            bbox = [ [xv[i,j], xv[i+1,j], xv[i+1,j+1], xv[i,j+1], xv[i,j]],
                     [yv[i,j], yv[i+1,j], yv[i+1,j+1], yv[i,j+1], yv[i,j]] ]
            subplot.append( np.asarray(bbox).transpose() )

    subplots[pp][keys[ss]] = subplot
    
    n_subplots=len(subplot)
    PAI_res[keys[ss]] = np.zeros((3,n_subplots))

# now get point clouds
sample_pts_collated = []
starting_ids_collated = []
trees_collated = []
for pp in range(0,3):
    print "point %i, x = %.0f, y = %.0f" % (pp+1, coords[pp,0], coords[pp,1])
    # define a bounding box around target points to load in point cloud around area of interest
    E = coords[pp,0]+500
    N = coords[pp,1]+500
    W = coords[pp,0]-500
    S = coords[pp,1]-500
    # Read in LiDAR points for region of interest
    polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
    lidar_pts, starting_ids_for_trees, trees = lidar_io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=laz_files)

    sample_pts_collated.append(lidar_pts.copy())
    starting_ids_collated.append(starting_ids_for_trees.copy())
    trees_collated.append(np.asarray(trees))

# Now loop through the subplots and sample the point cloud
for pp in range(0,3):
    print "point %i, x = %.0f, y = %.0f" % (pp+1, coords[pp,0], coords[pp,1])
    # loop through each sampling resolution
    lidar_pts = sample_pts_collated[pp].copy()
    starting_ids_for_trees = starting_ids_collated[pp].copy()
    N_trees = trees_collated[pp].size
    
    for ss in range(0,res.size):
      print '\t - sample res = ', keys[ss]
      n_subplots = len(subplots[pp][keys[ss]])
      rad_ss = np.sqrt(res[ss]**2/2.)              
      # for each of the subplots, clip point cloud and model PAD and get the metrics
      for sp in range(0,n_subplots):
        # query the tree to locate points of interest
        # note that we will only have one tree for number of points in sensitivity analysis  
        centre_x = np.mean(subplots[pp][keys[ss]][sp][0:4,0])
        centre_y = np.mean(subplots[pp][keys[ss]][sp][0:4,1])
        radius = np.sqrt(res[ss]**2/2.)
        # retrieve point clouds samples        
        sample_pts = np.array([])
        for tt in range(0,N_trees):
          ids = trees_collated[pp][tt].query_ball_point(np.array([centre_x,centre_y]), radius)
          if len(ids)>0:
            if sample_pts.size==0:
              sample_pts = lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]
            else:
              sample_pts = np.concatenate((sample_pts,lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]),axis=0)
              
        # keep only first returns
        sample_pts=sample_pts[sample_pts[:,3]==1,:]
        sp_pts = lidar.filter_lidar_data_by_polygon(sample_pts,subplots[pp][keys[ss]][sp])
        #------
        heights,first_return_profile,n_ground_returns = MH.bin_returns(sp_pts, max_height, layer_thickness)
        PADprof = MH.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, k)
        # remove lowermost portion of profile
        PADprof[heights<min_height]=0
        PAI_res[keys[ss]][pp,sp] = PADprof.sum()
      PAI_mean[pp,ss] = np.mean(PAI_res[keys[ss]][pp,:])
      PAI_sd[pp,ss] = np.std(PAI_res[keys[ss]][pp,:]) 

PAI_serr = np.zeros(PAI_sd.shape)
for pp in range(0,3):
  for ss in range(0,res.size):
    PAI_serr[pp,ss]=PAI_sd[pp,ss] / np.sqrt(PAI_res[keys[ss]][pp,:].size)


# Now want to get spatial scaling of canopy variance
# First create array for 10 m resolution case
PAI_array_10m = np.zeros((3,10,10))
for pp in range(0,3):
  sp = 0
  for rr in range(0,10):
    for cc in range(0,10):
      PAI_array_10m[pp,rr,cc]=PAI_res['10m'][pp,sp]
      sp+=1

test_res = np.arange(1,11)
e = np.zeros((3,len(test_res)))
bias = np.zeros((3,len(test_res)))

for tt in range(0,len(test_res)):
  for pp in range(0,3):
    temp_host = np.zeros((10-test_res[tt]+1,10-test_res[tt]+1))
    temp_host2 = np.zeros((10-test_res[tt]+1,10-test_res[tt]+1))
    for rr in range(0,10-test_res[tt]+1):
      for cc in range(0,10-test_res[tt]+1):
        sample_PAI = PAI_array_10m[pp,rr:rr+test_res[tt],cc:cc+test_res[tt]]
        sample_E = sample_PAI-np.mean(sample_PAI)
        temp_host2[rr,cc]= np.mean(sample_E)
        temp_host[rr,cc] = -(1/k)*np.log(np.mean(np.exp(-k*sample_E)))
        #temp_host[rr,cc] = np.std(PAI_array_10m[pp,rr:rr+test_res[tt]+1,cc:cc+test_res[tt]+1])
    #print temp_host.shape
    #PAI_std_scaling[pp,tt] = np.mean(temp_host)
    bias[pp,tt] = np.mean(temp_host)
    e[pp,tt]=np.mean(temp_host2)
    print pp,tt,bias[pp,tt], e[pp,tt]
                                      
# now use linear interpolation to estimate bias at each of the resolutions used
# in this analysis  
bias_interpolated = np.zeros((3,res.size))    
for pp in range(0,3):
  for rr in range(0,res.size):
    res1 = test_res[test_res*10<=res[rr]][-1]
    res2 = test_res[test_res*10>=res[rr]][0]
    bias1 = bias[pp,test_res==res1][0]
    bias2 = bias[pp,test_res==res2][0]
    print res[rr],res1,res2
    if res1!=res2:
      bias_interpolated[pp,rr] = bias1+(bias2-bias1)*(res[rr]/10-res1)/(res2-res1)
    else:
      bias_interpolated[pp,rr] = bias1
PAI_corrected = PAI_mean-bias_interpolated

# Now plot up the results
fig = plt.figure(4, facecolor='White',figsize=[7,6])
ax4a= plt.subplot2grid((1,5),(0,0),colspan=4)
ax4a.errorbar(res,PAI_mean[0,:],yerr=2*PAI_serr[0,:],marker='o',c=colour[2],label = 'A',linestyle='none')
ax4a.plot(res,PAI_corrected[0,:],marker='^',c=colour[2],linestyle='none')
ax4a.axhline(PAI_mean[0,0],c=colour[2],linestyle=':')
ax4a.errorbar(res,PAI_mean[1,:],yerr=2*PAI_serr[1,:],marker='o',c=colour[1],label = 'B',linestyle='none')
ax4a.plot(res,PAI_corrected[1,:],marker='^',c=colour[1],linestyle='none')
ax4a.axhline(PAI_mean[1,0],c=colour[1],linestyle=':')
ax4a.errorbar(res,PAI_mean[2,:],yerr=2*PAI_serr[2,:],marker='o',c=colour[0],label = 'C',linestyle='none')
ax4a.plot(res,PAI_corrected[2,:],marker='^',c=colour[0],linestyle='none')
ax4a.axhline(PAI_mean[2,0],c=colour[0],linestyle=':')
ax4a.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax4a.set_xlabel('spatial resolution / m', fontsize = axis_size)
ax4a.set_ylabel('PAI', fontsize = axis_size)

plt.tight_layout()
plt.savefig('Figure4_resolution_vs_PAI_pointwise.png')
plt.savefig('Figure4_resolution_vs_PAI_pointwise.pdf')
plt.show()
