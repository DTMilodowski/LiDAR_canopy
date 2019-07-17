###############################################################################################################
# This driver function analyses the sensisitivity of the LiDAR-based metrics to spatial scale and point density
###############################################################################################################
import numpy as np
import sys
sys.path.append('../')
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats
import least_squares_fitting as lstsq
import time
import seaborn as sns
sns.set()
#---------------------------------------------------------------------------------------------------------------
# Some filenames & params
plots_b = [b'Belian',b'E',b'B North',b'OP']
plots = ['Belian','E','B North','OP']
chm_files = ['CHM_50cm_Belian.npy','CHM_50cm_E.npy','CHM_50cm_BNorth.npy','CHM_50cm_OP.npy']
chm = {}
x = np.arange(0,100,0.5)+.25
y = np.arange(0,100,0.5)+.25
xx,yy=np.meshgrid(x,y)
for pp,plot in enumerate(plots):
    chm[plot]=np.load(chm_files[pp])

profile_file = 'lidar_canopy_profiles_adaptive.npz'
canopy_profiles = np.load(profile_file)['arr_0'][()]
op_file = 'lidar_canopy_profiles_adaptive_OP.npz'
canopy_profiles[b'OP']=np.load(op_file)['arr_0'][()]['OP']

gps_pts_file = '../GPS_points_file_for_least_squares_fitting_sensitivity_analysis_version.csv'
datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('<U16','f16','f16','f16','f16')}
plot_coordinates = np.genfromtxt(gps_pts_file, skip_header = 0, delimiter = ',',dtype=datatype)

heights = np.arange(0.,80)+1
plot_width = 100.

#---------------------------------------------------------------------------------------------------------------
# Generate North arrows
n_arrow={}
n_arrow_text={}
for plot in (plots):
    if plot in plot_coordinates['plot']:
        mask = plot_coordinates['plot']==plot
        affine=lstsq.least_squares_affine_matrix(plot_coordinates['x_prime'][mask],plot_coordinates['y_prime'][mask],
                                                plot_coordinates['x'][mask],plot_coordinates['y'][mask])
        n_arrow[plot] = np.array(lstsq.apply_affine_transformation(np.array([0.,0.]),
                                            np.array([0.,15.]),affine)).transpose() # simple square bounding box applied for all sensitivity analyses
        n_arrow_text[plot] = np.array(lstsq.apply_affine_transformation(np.array([0.,0.]),
                                            np.array([0.,18.]),affine)).transpose() # simple square bounding box applied for all sensitivity analyses
    else:
        n_arrow[plot] = np.array(([0.,0.],[0.,15.])) # simple square bounding box applied for all sensitivity analyses
        n_arrow_text[plot] = np.array(([0.,0.],[0.,18]))

    shift_x = np.max(n_arrow[plot][:,0])-90
    shift_y = np.max(n_arrow[plot][:,1])-90
    n_arrow[plot][:,0]-=shift_x
    n_arrow[plot][:,1]-=shift_y
    n_arrow_text[plot][:,0]-=shift_x
    n_arrow_text[plot][:,1]-=shift_y

# Plot up the figure
labels = ['Old growth (MLA-01)','Moderately logged (SAF-03)','Heavily logged (SAF-02)','Oil Palm (OP)']
fig,axes = plt.subplots(nrows=2,ncols=2,figsize=[6,6],subplot_kw={'aspect':'equal',
                        'xlim':(0,100),'ylim':(0,100)})
for pp,ax in enumerate(axes.ravel()):
    mean_prof = np.mean(canopy_profiles[plots_b[pp]],axis=0)
    sem_prof = stats.sem(canopy_profiles[plots_b[pp]],axis=0)
    # chm
    im=ax.pcolormesh(xx,yy,chm[plots[pp]],cmap='viridis',vmin=0,vmax=80)
    # profile
    ax.plot(mean_prof*100,heights,'-',color='white')
    ax.fill_betweenx(heights,(mean_prof-sem_prof)*100,(mean_prof+sem_prof)*100,color='white',alpha=0.5,linewidth=0)
    # N arrow
    N = n_arrow[plots[pp]]
    N_txt = n_arrow_text[plots[pp]]
    ax.annotate('',xytext=N[0],xy=N[1],arrowprops=dict(arrowstyle='-|>',color='white'))
    ax.annotate('N',xy=N_txt[1],color='white',va='center',ha='center',fontsize=10)
    # scale bar
    ax.annotate('20m',xy=(92.5,15),color='white',va='center',ha='right',fontsize=10)
    ax.plot((95,95),(5,25),'-',color='white',lw=1)
    ax.plot((92.5,97.5),(25,25),'-',color='white',lw=1)
    ax.plot((92.5,97.5),(5,5),'-',color='white',lw=1)
    # axis label
    ax.set_title(labels[pp])
    ax.set_xticklabels(ax.get_xticklabels(),visible=False)
    ax.set_yticklabels(ax.get_yticklabels(),visible=False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.84,0.245,0.04,0.51])
fig.colorbar(im,cax=cbar_ax, label='canopy height / m')
fig.savefig('canopy_structure_across_gradient.png')
fig.show()
