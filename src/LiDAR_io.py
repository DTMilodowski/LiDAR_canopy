## LiDAR_io.py
## This library contains the code to read LiDAR point clouds from .las files.
## It then creates kd-trees associated with each tile. This will subsequently
## allow v. efficient spatial searches to be done.  Moreover, the when dealing
## with >10^6 points (build time on our linux server ~0.5 s), we split the data
## into multiple trees due to the non-linear increase in build time with the
## number of points.
import numpy as np
import laspy as las
import os
from scipy import spatial
import LiDAR_tools as lidar
import time

# get bounding box from las file
def get_lasfile_bbox(las_file):
    lasFile = las.file.File(las_file,mode='r-')
    max_xyz = lasFile.header.max
    min_xyz = lasFile.header.min
    UR = np.asarray([max_xyz[0],max_xyz[1]])
    LR = np.asarray([max_xyz[0],min_xyz[1]])
    UL = np.asarray([min_xyz[0],max_xyz[1]])
    LL = np.asarray([min_xyz[0],min_xyz[1]])
    lasFile.close()
    return UR, LR, UL, LL

# get bounding box from many las files or laz files
def get_bbox_of_multiple_tiles(file_list,laz_files=False,return_zlim=False):
    las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
    n_files = las_files.size

    if laz_files:
        temp_file = 'temp_%i.las' % np.round(np.random.random()*10**9).astype(int)
        os.system("las2las %s %s" % (las_files[0],temp_file))
        lasFile = las.file.File('%s' % temp_file,mode='r-')
        max_xyz = lasFile.header.max
        min_xyz = lasFile.header.min
        xmin = min_xyz[0]
        ymin = min_xyz[1]
        xmax = max_xyz[0]
        ymax = max_xyz[1]
        zmin = min_xyz[2]
        zmax = max_xyz[2]
        lasFile.close()
        os.system("rm %s" % temp_file)

        for i in range(1,n_files):
            temp_file = 'temp_%i.las' % np.round(np.random.random()*10**9).astype(int)
            os.system("las2las %s %s" % (las_files[i],temp_file))
            lasFile = las.file.File('%s' % temp_file,mode='r-')
            max_xyz = lasFile.header.max
            min_xyz = lasFile.header.min
            xmin = min(xmin,min_xyz[0])
            ymin = min(ymin,min_xyz[1])
            xmax = max(xmax,max_xyz[0])
            ymax = max(ymax,max_xyz[1])
            zmin = min(zmin,max_xyz[2])
            zmax = max(zmax,max_xyz[2])
            lasFile.close()
            os.system("rm %s" % temp_file)

    else:

        lasFile = las.file.File(las_files[0],mode='r-')
        max_xyz = lasFile.header.max
        min_xyz = lasFile.header.min
        xmin = min_xyz[0]
        ymin = min_xyz[1]
        xmax = max_xyz[0]
        ymax = max_xyz[1]
        zmin = min_xyz[2]
        zmax = max_xyz[2]
        lasFile.close()

        for i in range(1,n_files):
            lasFile = las.file.File(las_files[i],mode='r-')
            max_xyz = lasFile.header.max
            min_xyz = lasFile.header.min
            xmin = min(xmin,min_xyz[0])
            ymin = min(ymin,min_xyz[1])
            xmax = max(xmax,max_xyz[0])
            ymax = max(ymax,max_xyz[1])
            zmin = min(zmin,min_xyz[2])
            zmax = max(zmax,max_xyz[2])
            lasFile.close()

    UR = np.asarray([xmax,ymax])
    LR = np.asarray([xmax,ymin])
    UL = np.asarray([xmin,ymax])
    LL = np.asarray([xmin,ymin])
    Zlim = np.asarray([zmin,zmax])
    if return_zlim:
        return UR, LR, UL, LL, Zlim
    else:
        return UR, LR, UL, LL


# Load lidar data => x,y,z,return,class, scan angle, gps_time
# Returns: - a numpy array containing the points
# Added checks to see what scan angle information is available
def load_lidar_data(las_file,print_npts=True):
    lasFile = las.file.File(las_file,mode='r')

    try:
        pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num,
                    lasFile.classification, lasFile.scan_angle,
                    lasFile.gps_time, lasFile.num_returns)).transpose()
    except:
        try:
            pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num,
                    lasFile.classification, lasFile.scan_angle_rank,
                    lasFile.gps_time, lasFile.num_returns)).transpose()
        except:
            print('warning - no scan angle information, assuming zero for all points')
            pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num,
                    lasFile.classification, np.zeros(lasFile.x.size),
                    lasFile.gps_time, lasFile.num_returns)).transpose()
        else:
            print('warning - scan angle not present, returning scan angle rank')

    pts = pts[pts[:,2]>=0,:]
    if print_npts:
        print("loaded ", pts[:,0].size, " points")
    lasFile.close()
    return pts

# a similar script, but now only loading points within bbox into memory
def load_lidar_data_by_bbox(las_file,N,S,E,W,print_npts=True):
    lasFile = las.file.File(las_file,mode='r')

    # conditions for points to be included
    X_valid = np.logical_and((lasFile.x <= E), (lasFile.x >= W))
    Y_valid = np.logical_and((lasFile.y <= N), (lasFile.y >= S))
    Z_valid = lasFile.z >= 0
    ii = np.where(np.logical_and(X_valid, Y_valid, Z_valid))

    try:
        pts = np.vstack((lasFile.x[ii], lasFile.y[ii], lasFile.z[ii], lasFile.return_num[ii],
                    lasFile.classification[ii], lasFile.scan_angle[ii],
                    lasFile.gps_time[ii], lasFile.num_returns[ii])).transpose()
    except:
        try:
            pts = np.vstack((lasFile.x[ii], lasFile.y[ii], lasFile.z[ii], lasFile.return_num[ii],
                    lasFile.classification[ii], lasFile.scan_angle_rank[ii],
                    lasFile.gps_time[ii], lasFile.num_returns[ii])).transpose()
        except:
            print('warning - no scan angle information, assuming zero for all points')
            pts = np.vstack((lasFile.x[ii], lasFile.y[ii], lasFile.z[ii], lasFile.return_num[ii],
                    lasFile.classification[ii], np.zeros(lasFile.x[ii].size),
                    lasFile.gps_time[ii], lasFile.num_returns[ii])).transpose()
        else:
            print('warning - scan angle not present, returning scan angle rank')

    #pts = np.vstack((lasFile.x[ii], lasFile.y[ii], lasFile.z[ii],
    #                lasFile.return_num[ii], lasFile.classification[ii],
    #                lasFile.scan_angle_rank[ii], lasFile.gps_time[ii],
    #                lasFile.num_returns[ii])).transpose()

    if print_npts:
        print("loaded ", pts[:,0].size, " points")
    lasFile.close()
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
def find_las_files_by_polygon(file_list,polygon,print_keep=False):
    las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
    keep = []
    n_files = las_files.size
    for i in range(0,n_files):
        UR, LR, UL, LL = get_lasfile_bbox(las_files[i])
        las_box = np.asarray([UR,LR,LL,UL])
        x,y,inside = lidar.points_in_poly(las_box[:,0],las_box[:,1],polygon) # fix this test
        if inside.sum()>0:
            keep.append(las_files[i])
        else:
            x,y,inside = lidar.points_in_poly(polygon[:,0],polygon[:,1],las_box)
            if inside.sum()>0:
                keep.append(las_files[i])
    if print_keep:
        print('las tiles to load in:', len(keep))
        for ll in range(0,len(keep)):
            print(keep[ll])
    return keep

# load all lidar points from multiple las files witin specified polygon.  The file list needs to have either the full or relative path to the files included.
# polygon is a 2D array with N_pts*rows and two cols (x,y).
# I've added a fudge to deal with laz files, which uses las2las from lastools to create a temporary .las file before reading with laspy. Undoubtedly not the most elegant solution, but it works :-)
def load_lidar_data_by_polygon(file_list,polygon,max_pts_per_tree = 10**6, laz_files=False,print_keep=False):
    W = polygon[:,0].min()
    E = polygon[:,0].max()
    S = polygon[:,1].min()
    N = polygon[:,1].max()
    if laz_files:
        keep_files = find_laz_files_by_polygon(file_list,polygon,print_keep)
    else:
        keep_files = find_las_files_by_polygon(file_list,polygon,print_keep)

    n_files = len(keep_files)
    trees = []
    starting_ids = np.asarray([])

    # first case scenario that no points in ROI
    print('\t\tloading %.0f tiles...' % n_files)
    start=time.time()
    if n_files == 0:
        print('WARNING: No files within specified polygon - try again')
        pts = np.array([])

    # otherwise, we have work to do!
    else:
        if laz_files:
            temp_file = 'temp_%i.las' % np.round(np.random.random()*10**9).astype(int)
            os.system("las2las %s %s" % (keep_files[0],temp_file))
            #os.system("las2las %s temp.las" % keep_files[0])
            tile_pts = load_lidar_data_by_bbox(temp_file,N,S,E,W,print_npts=False)
            os.system("rm %s" % temp_file)
        else:
            tile_pts = load_lidar_data_by_bbox(keep_files[0],N,S,E,W,print_npts=False)

        pts = lidar.filter_lidar_data_by_polygon(tile_pts,polygon)

        # now repeat for subsequent tiles
        for i in range(1,n_files):
            if laz_files:
                temp_file = 'temp_%i.las' % np.round(np.random.random()*10**9).astype(int)
                os.system("las2las %s %s" % (keep_files[i],temp_file))
                tile_pts = load_lidar_data_by_bbox(temp_file,N,S,E,W,print_npts=False)
                os.system("rm %s" % temp_file)
            else:
                tile_pts = load_lidar_data_by_bbox(keep_files[i],N,S,E,W,print_npts=False)

            pts_ = lidar.filter_lidar_data_by_polygon(tile_pts,polygon)
            pts = np.concatenate((pts,pts_),axis=0)

    end=time.time()
    print('\t\t\t...%.3f s' % (end-start))
    # now create KDTrees
    print('\t\tbuilding KD-trees...')
    start=time.time()
    starting_ids, trees = create_KDTree(pts)
    end=time.time()
    print('\t\t\t...%.3f s' % (end-start))

    print("loaded ", pts.shape[0], " points into ", len(trees), " KDTrees")
    return pts, starting_ids, trees

# equivalent file but for a single las file
def load_lidar_file_by_polygon(lasfile,polygon,max_pts_per_tree = 10**6,print_keep=False,filter_by_first_return_location=False):
    W = polygon[:,0].min()
    E = polygon[:,0].max()
    S = polygon[:,1].min()
    N = polygon[:,1].max()
    tile_pts = load_lidar_data_by_bbox(lasfile,N,S,E,W,print_npts=False)
    pts = lidar.filter_lidar_data_by_polygon(tile_pts,polygon,filter_by_first_return_location)
    # now create KDTrees
    starting_ids, trees = create_KDTree(pts)

    print("loaded ", pts.shape[0], " points")
    return pts, starting_ids, trees

# equivalent scripts for laz files - calls las2las to transform to .las files
# before reading in with laspy
def find_laz_files_by_polygon(file_list,polygon,print_keep=False):
    laz_files = np.genfromtxt(file_list,delimiter=',',dtype='unicode')
    keep = []
    n_files = laz_files.size
    for i in range(0,n_files):
        temp_file = 'temp_%i.las' % np.round(np.random.random()*10**9).astype(int)
        os.system("las2las %s %s" % (laz_files[i],temp_file))
        UR, LR, UL, LL = get_lasfile_bbox(temp_file)
        las_box = np.asarray([UR,LR,LL,UL])
        x,y,inside = lidar.points_in_poly(las_box[:,0],las_box[:,1],polygon)
        if inside.sum()>0:
            keep.append(laz_files[i])
        os.system("rm %s" % temp_file)
    if print_keep:
        print('las tiles to load in:', len(keep))
        for ll in range(0,len(keep)):
            print(keep[ll])
    return keep

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
        print('WARNING: No files within specified neighbourhood - try again')
        pts = np.array([])

    else:
        # first tile
        tile_pts = load_lidar_data(keep_files[0],print_npts==False)
        pts = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
        # and loop through the remaining tiles
        for i in range(1,n_files):
            tile_pts = load_lidar_data(keep_files[i])
            pts_ = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
            pts = np.concatenate((pts,pts_),axis=0)

    # now create KDTrees
    starting_ids, trees = create_KDTree(pts)

    print("loaded ", pts[:,0].size, " points")
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
