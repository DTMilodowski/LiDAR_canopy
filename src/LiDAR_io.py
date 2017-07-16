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


# Load lidar data => x,y,z,return,class, scan angle
# Also creates kd-trees to host data.
# Returns: - a numpy array containing the points
#          - a second  array containing the starting indices of the points
#            associated with a given tree for cross checking against the
#            point cloud 
#          - a list of trees
def load_lidar_data(las_file,max_pts_per_tree = 10**6):
    lasFile = las.file.File(las_file,mode='r')
    pts = np.vstack((lasFile.x, lasFile.y, lasFile.z, lasFile.return_num, lasFile.classification, lasFile.scan_angle_rank),dtype='float16').transpose()
    pts = pts[pts[:,2]>=0,:]

    # now create kdtrees :-)
    npts = pts.size
    ntrees = np.ceil(npts/float(max_pts_per_tree));
    trees = []
    starting_ids = np.zeros(n_trees,dtype='int')
    
    for tt in range(0,ntrees):
        i0=tt*max_pts_per_tree
        i1 = (tt+1)*max_pts_per_tree
        trees.append(spatial.cKDTree(pts[i0:i1,0:2],leafsize=32,balanced_tree=True))
        starting_ids[tt] = i0

    print "loaded ", pts[:,0].size, " points"
    return pts, starting_ids, 


