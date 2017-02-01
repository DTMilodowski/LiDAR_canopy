import numpy as np
import laspy as las
import LiDAR_tools as LiDAR

# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This function
# returns True or False.  The algorithm is called
# the "Ray Casting Method".
# the point_in_poly algorithm was found here:
# http://geospatialpython.com/2011/01/point-in-polygon.html
def point_in_poly(x,y,poly):

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

# This one is my own version of the ray-trace algorithm which utilises the numpy arrays so that a list of x and y coordinates can be processed in one call and only points inside polygon are returned alongside the indices in case required for future referencing.  This saves a fair bit of looping.
def points_in_poly(x,y,poly):
    n = len(poly)
    inside=np.zeros(x.size,dtype=bool)
    xints=np.zeros(x.size)

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y=poly[i % n]
        if p1y!=p2y:
            xints[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)] = (y[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
        if p1x==p2x:
            inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)] = np.invert(inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)])
        else:
            inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x),x<=xints],axis=0)] = np.invert(inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x),x<=xints],axis=0)])
        p1x,p1y = p2x,p2y

    return x[inside],y[inside], inside
        


# This function loads the subplot coordinates from a csv file.  File columns should be as follows:
# Plot Subplot 'X0', 'Y0', 'X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3'
def load_boundaries(coordinate_list):
    
    datatype = {'names': ('Plot', 'Subplot', 'X0', 'Y0', 'X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3'), 'formats': ('S32','i8','f16','f16','f16','f16','f16','f16','f16','f16')}
    coords = np.genfromtxt(coordinate_list, skiprows = 1, delimiter = ',',dtype=datatype)
    plot_name=np.unique(coords['Plot'])
    coordinate_dict = {}
    subplot_dict = {}
    for pp in range(0,plot_name.size):
        n_subplots = np.sum(coords['Plot']==plot_name[pp])
        subplot_polygons = np.zeros((n_subplots,5,2))
        subplot_list = np.zeros(n_subplots)
        for i in range(0,n_subplots):
            subplot_polygons[i,0,0]=coords['X0'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,0,1]=coords['Y0'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,1,0]=coords['X1'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,1,1]=coords['Y1'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,2,0]=coords['X2'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,2,1]=coords['Y2'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,3,0]=coords['X3'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,3,1]=coords['Y3'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,4,0]=coords['X0'][coords['Plot']==plot_name[pp]][i]
            subplot_polygons[i,4,1]=coords['Y0'][coords['Plot']==plot_name[pp]][i]
            subplot_list[i]=coords['Subplot'][coords['Plot']==plot_name[pp]][i]
        coordinate_dict[plot_name[pp]]=subplot_polygons.copy()
        subplot_dict[plot_name[pp]]=subplot_list.copy()
    return coordinate_dict, subplot_dict

# get bounding box for list of coordinates
def get_bounding_box(coordinate_list):
    N = coordinate_list.shape[0]
    bbox = np.zeros((4,2))
    top = coordinate_list[:,1].max()
    bottom = coordinate_list[:,1].min()
    left = coordinate_list[:,0].min()
    right = coordinate_list[:,0].max()
    bbox[0,0]=left
    bbox[0,1]=top
    bbox[1,0]=right
    bbox[1,1]=top
    bbox[2,0]=right
    bbox[2,1]=bottom
    bbox[3,0]=left
    bbox[3,1]=bottom
    return bbox

# Load lidar data => x,y,z,return,class
def load_lidar_data(las_file):#,subplot_coords,max_height,bin_width):
    lasFile = las.file.File(las_file,mode='r')
    #pts = np.vstack((lasFile.x*lasFile.header.scale[0]+lasFile.header.offset[0], lasFile.y*lasFile.header.scale[1]+lasFile.header.offset[1], lasFile.z*lasFile.header.scale[2]+lasFile.header.offset[2], lasFile.return_num, lasFile.classification)).transpose()
    pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num, lasFile.classification, lasFile.scan_angle_rank)).transpose()
    pts = pts[pts[:,2]>=0,:]
    print "loaded ", pts[:,0].size, " points"
    return pts

# filter lidar wth polygon
def filter_lidar_data_by_polygon(in_pts,polygon):
    x,y,inside = points_in_poly(in_pts[:,0],in_pts[:,1],polygon)
    pts = in_pts[inside,:]
    return pts

# bin lidar returns 
def bin_returns(pts, max_height, layer_thickness):
    # filter the points to leave only first returns - not sure why this is necessary, but Cam group suggest - point for discussion
    print pts.shape
    pts=pts[pts[:,3]==1,:]
    
    # calculate n ground points
    n_ground_returns = np.sum(pts[:,4]==2)
    # filter to consider only veg returns
    can_pts = pts[pts[:,4]==1,:]
    
    # now set up bins
    lower_lims=np.arange(0,max_height,layer_thickness)
    upper_lims = lower_lims + layer_thickness
    n_bins = lower_lims.size
    profile = np.zeros(n_bins)
    heights = lower_lims+layer_thickness
    n_returns = can_pts.shape[0]

    # bin data
    bin = can_pts[:,2]//layer_thickness

    for i in range(0,n_returns):
        if bin[i]<n_bins and bin[i]>=0: # in case there are some spuriously high returns outside the permitted range
            profile[bin[i]]+=1

    return heights,profile,n_ground_returns

# get bounding box from las file
def get_lasfile_bbox(las_file):
    lasFile = las.file.File(las_file,mode='r')
    max_xyz = lasFile.header.max
    min_xyz = lasFile.header.min
    UR = np.asarray[max_xyz[0],max_xyz[1]]
    LR = np.asarray[max_xyz[0],min_xyz[1]]
    UL = np.asarray[min_xyz[0],max_xyz[1]]
    LL = np.asarray[min_xyz[0],min_xyz[1]]
    
    return UR, LR, UL, LL

# find all las files from a list that are located within a specified polygon
def find_las_files_by_polygon(file_list,polygon):
    las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
    keep = []
    n_files = las_files.size
    for i in range(0,n_files):
        UR, LR, UL, LL = get_lasfile_bbox(las_files[i])
        in_pts = np.asarray(UR,LR,UL,LL)
        x,y,inside = points_in_poly(in_pts[:,0],in_pts[:,1],polygon)
        if inside.sum()>0:
            keep.append(las_files[i])
    print keep
    return keep

# load all lidar points from multiple las files witin specified polygon.  The file list needs to have either the full or relative path to the files included.
def load_lidar_data_by_polygon(file_list,polygon):
    keep_files = find_las_files_by_polygon(file_list,polygon)
    n_files = len(keep_files)
    if n_files == 0:
        print 'WARNING: No files within specified polygon - try again'
    else:
        tile_pts = load_lidar_data(keep_files[0])
        pts = filter_lidar_data_by_polygon(tile_pts,polygon)
        for i in range(1,n_files):
            tile_pts = load_lidar_data(keep_files[i])
            tile_pts_filt = filter_lidar_data_by_polygon(tile_pts,polygon)
            pts = np.concatenate(pts,tile_pts_filt,axis=0)

    print "loaded ", pts[:,0].size, " points"
    return pts


# Use MacArther-Horn method to estimate LAD profile from the lidar return profile.  See methods described by Stark et al., Ecology Letters, 2012
def estimate_LAD_MacArtherHorn(lidar_profile,n_ground_returns,layer_thickness,k):
    n_layers = lidar_profile.size    
    S = np.zeros(n_layers+1)
    S[1:]=np.cumsum(lidar_profile)
    S+=n_ground_returns
    #S+=(2*n_ground_returns) # Harding et al., 2001 correction to account for the fact that ground reflectance is typically lower than canopy
    S[S==0]=1 # This step is required to stop the base of the profile (final return) kicking out errors if there are no ground returns
    S_in = S[1:]
    S_out= S[:-1]
    LAD_profile = np.log(S_in/S_out)/(k*layer_thickness)
    
    # Shouldn't have any divide by zeros, but just in case...
    if np.sum(np.isfinite(LAD_profile)==False)>0:
        print np.sum(np.isfinite(LAD_profile)==False)
    LAD_profile[np.isfinite(LAD_profile)==False]==0

    return LAD_profile

# Load spreadsheet of LAI derived from hemispherical photographs (courtesy of Terhi Riutta at Oxford).  LAI estimated using Hemisfer.
def load_field_LAI(LAI_file):
    datatype = {'names': ('ForestType','Plot', 'Subplot', 'LAI'), 'formats': ('S32','S32','i8','f16')}
    hemisfer_LAI = np.genfromtxt(LAI_file, skiprows = 1, delimiter = ',',dtype=datatype)

    return hemisfer_LAI

# Do some crunching to brute force the best fitting k for MacArther-Horn method.
def minimise_misfit_for_k(kmin,kmax,k_interval,subplot_LAIs,subplot_lidar_profiles,n_ground_returns,layer_thickness,minimum_height=0):

    n_layers = subplot_lidar_profiles.shape[1]
    heights = np.arange(1,n_layers+1)*layer_thickness
    mask = heights<=2
    # first of all, loop through the k values
    ks = np.arange(kmin,kmax+k_interval,k_interval)
    n_ks = ks.size
    misfit = np.zeros(n_ks)
    misfit_b = np.zeros(n_ks)
    n_subplots = subplot_LAIs.size
    for i in range(0,n_ks):
        # now loop through the subplots
        for j in range(0,n_subplots):
            LAD_profile = estimate_LAD_MacArtherHorn(subplot_lidar_profiles[j,:],n_ground_returns[j],layer_thickness,ks[i])
            LAD_profile[mask]=0
            misfit[i] += np.sqrt((np.sum(LAD_profile)-subplot_LAIs[j])**2)
    misfit/=float(n_subplots)
    best_k_LAD_profiles = np.zeros((subplot_lidar_profiles.shape))
    best_k = ks[misfit==np.min(misfit)]
    for j in range(0,n_subplots):
        best_k_LAD_profiles[j,:] = estimate_LAD_MacArtherHorn(subplot_lidar_profiles[j,:],n_ground_returns[j],layer_thickness,best_k)
    print "Field LAI: ", np.mean(subplot_LAIs), "+/-", np.std(subplot_LAIs),"; LiDAR LAI: ",np.mean(np.sum(best_k_LAD_profiles,axis=1)), "+/-", np.std(np.sum(best_k_LAD_profiles,axis=1)), "; best k: ", best_k, " m-1"

    return misfit, ks, best_k_LAD_profiles, best_k


def calculate_bestfit_LAD_profile(subplot_coordinate_file,LAI_file,las_file,Plot_name,minimum_height=0):
    subplot_polygons, subplot_labels = load_boundaries(subplot_coordinate_file)
    field_LAI = load_field_LAI(LAI_file)
    lidar_pts = load_lidar_data(las_file)
    
    n_subplots = subplot_polygons[Plot_name].shape[0]
    max_height = 80
    layer_thickness = 1
    n_layers = np.ceil(max_height/layer_thickness)
    subplot_lidar_profiles = np.zeros((n_subplots,n_layers))
    n_ground_returns = np.zeros(n_subplots)
    subplot_LAI = np.zeros(n_subplots)

    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        sp_pts = filter_lidar_data_by_polygon(lidar_pts,subplot_polygons[Plot_name][i,:,:])
        heights,subplot_lidar_profiles[i,:],n_ground_returns[i] = bin_returns(sp_pts, max_height, layer_thickness)
        subplot_LAI[i] = field_LAI['LAI'][np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)]
        
    kmin = 0.20
    kmax = 5.
    kinc = 0.005
    misfit, ks, best_k_LAD_profiles, best_k = minimise_misfit_for_k(kmin,kmax,kinc,subplot_LAI,subplot_lidar_profiles,n_ground_returns,layer_thickness,minimum_height)
    return heights, best_k_LAD_profiles


if __name__ == "__main__":
    LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'
    las_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/C_plots.las'
    las_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/Danum_plots.las'
    las_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/Plot_E.las'
    subplot_coordinate_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/src/Cplots_subplot_coordinates.csv'
    subplot_coordinate_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/src/Danum_subplot_coordinates.csv'
    subplot_coordinate_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/LAI_profiles/LAD_estimation/src/E_subplot_coordinates.csv'
    Plot_name = 'Belian'

    subplot_polygons, subplot_labels = LiDAR.load_boundaries(subplot_coordinate_file)
    field_LAI = LiDAR.load_field_LAI(LAI_file)
    lidar_pts = LiDAR.load_lidar_data(las_file)
    
    n_subplots = subplot_polygons[Plot_name].shape[0]
    max_height = 80
    layer_thickness = 1
    minimum_height = 2
    n_layers = np.ceil(max_height/layer_thickness)
    subplot_lidar_profiles = np.zeros((n_subplots,n_layers))
    n_ground_returns = np.zeros(n_subplots)
    subplot_LAI = np.zeros(n_subplots)

    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        sp_pts = LiDAR.filter_lidar_data_by_polygon(lidar_pts,subplot_polygons[Plot_name][i,:,:])
        heights,subplot_lidar_profiles[i,:],n_ground_returns[i] = LiDAR.bin_returns(sp_pts, max_height, layer_thickness)
        subplot_LAI[i] = field_LAI['LAI'][np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)]
        
    kmin = 0.20
    kmax = 5.
    kinc = 0.005
    misfit, ks, best_k_LAD_profiles, best_k = LiDAR.minimise_misfit_for_k(kmin,kmax,kinc,subplot_LAI,subplot_lidar_profiles,n_ground_returns,layer_thickness,minimum_height)
    # option here would be to monte-carlo the position estimates, then create misfit vs. ks "envelopes", but runtime for one iteration is ~1 minute, so this is likely to be prohibitive for more than ~500.

    plt.figure(1, facecolor='White',figsize=[7,7])
    ax1 = plt.subplot2grid((2,2),(0,0))
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12) 
    ax1.plot(ks,misfit,'-k')
    ax1.set_ylim(0,6)
    ax1.set_xlim(0,5)
    ax1.set_ylabel('Average subplot misfit in LAI')
    ax1.set_xlabel('k / m$^{-1}$')

    ax2 = plt.subplot2grid((2,2),(0,1),rowspan=2)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax2.annotate(Plot_name, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=12) 
    for i in range(0,n_subplots):
        ax2.plot(best_k_LAD_profiles[i,:],heights,'-',c='g',alpha=0.2,linewidth=1)
    ax2.plot(np.mean(best_k_LAD_profiles,axis=0),heights,'-',c='k',linewidth=1)
    ax2.set_ylim(0,80)
    ax2.set_ylabel('Height / m')
    ax2.set_xlabel('LAD / m$^2$m$^{-1}$')

    ax3 = plt.subplot2grid((2,2),(1,0))
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12) 
    ax3.annotate('k='+str(best_k[0])+ 'm$^{-1}$', xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=12)
    ax3.plot([0,10],[0,10],'-',color='0.5')
    for i in range(0,n_subplots):
        ax3.plot(subplot_LAI[i],np.sum(best_k_LAD_profiles[i,:]),'.',c='k')
    ax3.set_ylim(0,10)
    ax3.set_xlim(0,10)
    ax3.set_xlabel('LAI from Hemisfer')
    ax3.set_ylabel('LAI from LiDAR')


    plt.tight_layout()
    plt.savefig(Plot_name+'_LAD_profile.png')
    plt.show()
