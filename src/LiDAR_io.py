## LiDAR_io.py
## This library contains the code to read LiDAR point clouds from .las files.
## It then creates kd-trees associated with each tile. This will subsequently
## allow v. efficient spatial searches to be done.  Moreover, the when dealing
## with >10^6 points (build time on our linux server ~0.5 s), we split the data
## into multiple trees due to the non-linear increase in build time with the
## number of points. 
import numpy as np
import laspy as las
from scipy import spatial
import LiDAR_tools as lidar

# get bounding box from las file
def get_lasfile_bbox(las_file):
    lasFile = las.file.File(las_file,mode='r')
    max_xyz = lasFile.header.max
    min_xyz = lasFile.header.min
    UR = np.asarray([max_xyz[0],max_xyz[1]])
    LR = np.asarray([max_xyz[0],min_xyz[1]])
    UL = np.asarray([min_xyz[0],max_xyz[1]])
    LL = np.asarray([min_xyz[0],min_xyz[1]])
    
    return UR, LR, UL, LL

# Load lidar data => x,y,z,return,class, scan angle
# Returns: - a numpy array containing the points
def load_lidar_data(las_file):
    lasFile = las.file.File(las_file,mode='r')
    pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num, lasFile.classification, lasFile.scan_angle_rank)).transpose().astype('float32')
    pts = pts[pts[:,2]>=0,:]
    print "loaded ", pts[:,0].size, " points"
    return pts

# a similar script, but now only loading points within bbox into memory
def load_lidar_data_by_bbox(las_file,N,S,E,W):
    lasFile = las.file.File(las_file,mode='r')

    # conditions for points to be included
    X_valid = np.logical_and((lasFile.x <= E), (lasFile.x >= W))
    Y_valid = np.logical_and((lasFile.y <= N), (lasFile.y >= S))
    Z_valid = lasFile.z >= 0
    ii = np.where(np.logical_and(X_valid, Y_valid, Z_valid))

    pts = np.vstack((lasFile.x[ii], lasFile.y[ii], lasFile.z[ii], lasFile.return_num[ii], lasFile.classification[ii], lasFile.scan_angle_rank[ii])).transpose().astype('float32')
    print "loaded ", pts[:,0].size, " points"
    return pts

#---------------------------------------------------------------------------
# KD-Trees :-)
# Creates kd-trees to host data.
# RETURNS  - a second  array containing the starting indices of the points
#            associated with a given tree for cross checking against the
#            point cloud 
#          - a list of trees
def create_KDTree(pts,max_pts_per_tree = 10**6):
    npts = pts.shape[0]
    ntrees = int(np.ceil(npts/float(max_pts_per_tree)))
    print npts,ntrees, int(ntrees)
    trees = []
    starting_ids = []
    
    for tt in range(0,ntrees):
        i0=tt*max_pts_per_tree
        i1 = (tt+1)*max_pts_per_tree
        if i1 < pts.shape[0]:
            trees.append(spatial.cKDTree(pts[i0:i1,0:2],leafsize=32,balanced_tree=True))
        else:
            trees.append(spatial.cKDTree(pts[i0:,0:2],leafsize=32,balanced_tree=True))
        starting_ids.append(i0)

    return np.asarray(starting_ids,dtype='int'), trees


#----------------------------------------------------------------------------
## Now we have some more involved reading scripts that only load in the data
## that satisfy certain neighbourhood criteria specified either as polygon or
## a neighbourhood

#----------------------------------------------------------------------------
# First of all here are scripts that use a polygon to load subset
#----------------------------------------------------------------------------

# find all las files from a list that are located within a specified polygon
def find_las_files_by_polygon(file_list,polygon):
    las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
    keep = []
    n_files = las_files.size
    for i in range(0,n_files):
        UR, LR, UL, LL = get_lasfile_bbox(las_files[i])
        las_box = np.asarray([UR,LR,LL,UL])
        x,y,inside = lidar.points_in_poly(polygon[:,0],polygon[:,1],las_box)
        if inside.sum()>0:
            keep.append(las_files[i])
    print 'las tiles to load in:', len(keep)
    for ll in range(0,len(keep)):
        print keep[ll]
    return keep

# load all lidar points from multiple las files witin specified polygon.  The file list needs to have either the full or relative path to the files included.
# polygon is a 2D array with N_pts*rows and two cols (x,y)
def load_lidar_data_by_polygon(file_list,polygon,max_pts_per_tree = 10**6):
    keep_files = find_las_files_by_polygon(file_list,polygon)
    n_files = len(keep_files)
    trees = []
    starting_ids = np.asarray([])

    # first case scenario that no points in ROI
    if n_files == 0:
        print 'WARNING: No files within specified polygon - try again'
        pts = np.array([])
                              
    # otherwise, we have work to do!
    else:
        W = polygon[:,0].min()
        E = polygon[:,0].max()
        S = polygon[:,1].min()
        N = polygon[:,1].max()
        tile_pts = load_lidar_data_by_bbox(keep_files[0],N,S,E,W)
        pts = lidar.filter_lidar_data_by_polygon(tile_pts,polygon)
                
        # now repeat for subsequent tiles
        for i in range(1,n_files):
            tile_pts = load_lidar_data(keep_files[i])
            pts_ = lidar.filter_lidar_data_by_polygon(tile_pts,polygon)
            pts = np.concatenate((pts,pts_),axis=0)
    
    # now create KDTrees
    starting_ids, trees = create_KDTree(pts)

    print "loaded ", pts.shape[0], " points"
    return pts, starting_ids, trees

# equivalent file but for a single las file
def load_lidar_file_by_polygon(lasfile,polygon,max_pts_per_tree = 10**6):
    W = polygon[:,0].min()
    E = polygon[:,0].max()
    S = polygon[:,1].min()
    N = polygon[:,1].max()
    tile_pts = load_lidar_data_by_bbox(lasfile,N,S,E,W)
    pts = lidar.filter_lidar_data_by_polygon(tile_pts,polygon)
    # now create KDTrees
    starting_ids, trees = create_KDTree(pts)

    print "loaded ", pts.shape[0], " points"
    return pts, starting_ids, trees


#----------------------------------------------------------------------------
# Next here are scripts that use a point and circular neighbourhood instead
#----------------------------------------------------------------------------
# Similar to above but using a focal point (specified by xy) and neighbourhood (specified by radius) to find .las tiles rather than using an input polygon
def find_las_files_by_neighbourhood(file_list,xy,radius):
    polygon = np.asarray([[xy[0]+radius,xy[1]+radius], [xy[0]+radius,xy[1]-radius], [xy[0]-radius,xy[1]-radius], [xy[0]-radius,xy[1]+radius]])
    keep = find_las_files_by_polygon(file_list,polygon)
    return keep

def load_lidar_data_by_neighbourhood(file_list,xy,radius,max_pts_per_tree = 10**6):
    keep_files = find_las_files_by_neighbourhood(file_list,xy,radius)
    n_files = len(keep_files)
    trees = []
    starting_ids = np.asarray([])

    if n_files == 0:
        print 'WARNING: No files within specified neighbourhood - try again'
        pts = np.array([])

    else:
        # first tile
        tile_pts = load_lidar_data(keep_files[0])
        pts = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
        # and loop through the remaining tiles
        for i in range(1,n_files):
            tile_pts = load_lidar_data(keep_files[i])
            pts_ = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
            pts = np.concatenate((pts,pts_),axis=0)

    # now create KDTrees
    starting_ids, trees = create_KDTree(pts)

    print "loaded ", pts[:,0].size, " points"
    return pts,starting_ids, trees



##------------------------------------------------------------------------------
## Functions to write point cloud files
##------------------------------------------------------------------------------

# This function writes a set of lidar returns into a csv file, so that the same 
# point cloud samples can be loaded into different software packages
def points_to_csv(pts,outfile):
    n_pts,temp = pts.shape
    f = open(outfile,"w") #opens file
    f.write("X, Y, Z, k, Class, A\n")
    for i in range(0,n_pts):
        f.write(str(pts[i,0])+','+str(pts[i,1])+','+str(pts[i,2])+','+str(pts[i,3])+','+str(pts[i,4])+','+str(pts[i,5])+'\n')
    f.close()
