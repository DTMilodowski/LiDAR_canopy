#===============================================================================
# PAI_limitations_figures.py
# D.T.Milodowski, November 2017
#-------------------------------------------------------------------------------
# This function contains the scripts used to produce the figures in the paper:
# "Point density imposes limitations on PAI estimation from discrete return
# LiDAR"
#-------------------------------------------------------------------------------
import numpy as np

# import plotting libraries
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
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
axis_size = rcParams['font.size']+4
colour = ['#46E900','#1A2BCE','#E0007F']

# import required LiDAR libaries
import LiDAR_MacHorn_LAD_profiles as MH
import raster_io as io

# Directory listings
SAVEDIR = '~/'

# Load files
dens_file = 'arrays_for_lombok/SAFE_pointcloud_metrics_10m_point_density_data.tif'
PAI_file = 'arrays_for_lombok/SAFE_pointcloud_metrics_10m_pai_data.tif'

dens, geo, coord = io.load_GeoTIFF_band_and_georeferencing(dens_file)
PAI, geo, coord = io.load_GeoTIFF_band_and_georeferencing(PAI_file)
dens[np.isnan(PAI)]=np.nan

rows,cols = PAI.shape
N = geo[3]
S = geo[3] + (rows+1)*geo[5]
W = geo[0]
E = geo[0] + (cols+1)*geo[1]
#-------------------------------------------------------------------------------
# Get analytical solution
k = 0.7
theta_deg = np.asarray([0, 5, 10])
theta_rad = theta_deg*np.pi/180.
radius = np.asarray([10.,20.,30.])
area = np.pi*radius**2
dens_a = np.arange(0.01,np.nanmax(dens),0.01)

PAImax = MH.calculate_analytical_limit(dens,area[0],k,theta_rad[0])

PAImax_10m_00deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[0])
PAImax_20m_00deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[0])
PAImax_30m_00deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[0])
"""
PAImax_10m_05deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[1])
PAImax_20m_05deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[1])
PAImax_30m_05deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[1])

PAImax_10m_10deg = MH.calculate_analytical_limit(dens_a,area[0],k,theta_rad[2])
PAImax_20m_10deg = MH.calculate_analytical_limit(dens_a,area[1],k,theta_rad[2])
PAImax_30m_10deg = MH.calculate_analytical_limit(dens_a,area[2],k,theta_rad[2])
"""
#-------------------------------------------------------------------------------
# Figure 1 - this figure illustrates the analytical solution presented in this
# paper describing the threshold PAI that can be detected using discrete return
# LiDAR
fig = plt.figure(1, facecolor='White',figsize=[12,6])
ax1a= plt.subplot2grid((1,7),(0,0),colspan=3)
ax1a.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2)
ax1a.set_xlabel('point density / pts m$^{-2}$',fontsize = axis_size)
ax1a.set_ylabel('PAI$_{max}$',fontsize = axis_size)

ax1a.plot(dens_a,PAImax_10m_00deg,'-',c=colour[0],label = '%.3f ha' % (np.pi*radius[0]**2/10.**4))
ax1a.plot(dens_a,PAImax_20m_00deg,'-',c=colour[1],label = '%.3f ha' % (np.pi*radius[1]**2/10.**4))
ax1a.plot(dens_a,PAImax_30m_00deg,'-',c=colour[2],label = '%.3f ha' % (np.pi*radius[2]**2/10.**4))
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

plt.savefig('Fig1_SAFE_point_density_vs_PAI.png')
plt.savefig('Fig1_SAFE_point_density_vs_PAI.pdf')

#-------------------------------------------------------------------------------
# Figure 2 - this figure presents maps of PAI and point density across the SAFE
# landscape

# Load shapefiles
shapefile_dir = '../Data/Fig1_Shapefiles/'
land_cover_file = 'HCS_Stratification_DOI.shp'
vjr_file = 'VJR_prj.shp'
ea_file = 'EA_prj.shp'

land_cover = shp.Reader(shapefile_dir+land_cover_file)
ea = shp.Reader(shapefile_dir+ea_file)
vjr = shp.Reader(shapefile_dir+vjr_file)
"""
for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y)
"""

fig = plt.figure(2, facecolor='White',figsize=[12,12])

ax2b= plt.subplot2grid((2,2),(0,1))
ax2b.annotate('b - Point density', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2)
im2b=ax2b.imshow(dens,vmin=0,vmax=30,cmap='plasma',origin='lower',extent=[W,E,N,S])
ax2b.axis('image')

divider2b = make_axes_locatable(ax2b)
cax2b = divider2b.append_axes("right", size="5%", pad=0.05)
cbar2b=plt.colorbar(im2b, cax=cax2b)
cbar2b.ax.set_ylabel('point density / pts m$^{-2}$',fontsize = axis_size)
cbar2b.solids.set_edgecolor("face")

ax2c= plt.subplot2grid((2,2),(1,0))
ax2c.annotate('c - PAI', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2)
im2c = ax2c.imshow(PAI,cmap='viridis',origin='lower',extent=[W,E,N,S])
ax2c.axis('image')

divider2c = make_axes_locatable(ax2c)
cax2c = divider2c.append_axes("right", size="5%", pad=0.05)
cbar2c=plt.colorbar(im2c, cax=cax2c)
cbar2c.ax.set_ylabel('PAI',fontsize = axis_size)
cbar2c.solids.set_edgecolor("face")

ax2d= plt.subplot2grid((3,1),(2,0))
ax2d.annotate('d - PAI/PAI$_{max}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size+2)
im2d = ax2d.imshow(PAI/PAImax,vmin = 0.85, vmax=1, cmap='plasma',origin='lower',extent=[W,E,N,S])
ax2d.axis('image')

divider2d = make_axes_locatable(ax2d)
cax2d = divider2d.append_axes("right", size="5%", pad=0.05)
cbar2d=plt.colorbar(im2d, cax=cax2d)
cbar2d.ax.set_ylabel('PAI/PAI$_{max}$',fontsize = axis_size)
cbar2d.solids.set_edgecolor("face")

plt.tight_layout()
plt.savefig('Fig2_SAFE_point_density_PAI_maps.png')
plt.savefig('Fig2_SAFE_point_density_PAI_maps.pdf')
plt.show()
