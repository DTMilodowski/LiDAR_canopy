import numpy as np
import laspy as las

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
        
# This retrieves all points within circular neighbourhood,  Terget point is the location around which the neighbourhood search is conducted, for a specified search radius.  x and y are vectors with the x and y coordinates of the test points
def points_in_radius(x,y,target_x, target_y,radius):
    inside=np.zeros(x.size,dtype=bool)
    d2=(x-target_x)**2+(y-target_y)**2
    inside = d2<=radius**2
    return x[inside],y[inside], inside
        



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
    pts = np.zeros((0,in_pts.shape[1]))
    if in_pts.shape[0]>0:
        x,y,inside = points_in_poly(in_pts[:,0],in_pts[:,1],polygon)
        pts = in_pts[inside,:]
    else:
        print "\t\t\t no points in polygon"
    return pts

# filter lidar by circular neighbourhood
def filter_lidar_data_by_neighbourhood(in_pts,target_xy,radius):
    pts = np.zeros((0,in_pts.shape[1]))
    if in_pts.shape[0]>0:
        x,y,inside =  points_in_radius(in_pts[:,0],in_pts[:,1],target_xy[0],target_xy[1],radius)
        pts = in_pts[inside,:]
    else:
        print "\t\t\t no points in neighbourhood"
    return pts


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

# find all las files from a list that are located within a specified polygon
def find_las_files_by_polygon(file_list,polygon):
    las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
    keep = []
    n_files = las_files.size
    for i in range(0,n_files):
        UR, LR, UL, LL = get_lasfile_bbox(las_files[i])
        las_box = np.asarray([UR,LR,LL,UL])
        x,y,inside = points_in_poly(polygon[:,0],polygon[:,1],las_box)
        if inside.sum()>0:
            keep.append(las_files[i])
    print 'las tiles to load in:', len(keep)
    for ll in range(0,len(keep)):
        print keep[ll]
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
            pts = np.concatenate((pts,tile_pts_filt),axis=0)

    print "loaded ", pts[:,0].size, " points"
    return pts

# Similar to above but using a focal point (specified by xy) and neighbourhood (specified by radius) to find .las tiles rather than using an input polygon
def find_las_files_by_neighbourhood(file_list,xy,radius):
    polygon = np.asarray([[xy[0]+radius,xy[1]+radius], [xy[0]+radius,xy[1]-radius], [xy[0]-radius,xy[1]-radius], [xy[0]-radius,xy[1]+radius]])
    keep = find_las_files_by_polygon(file_list,polygon)
    return keep

def load_lidar_data_by_neighbourhood(file_list,xy,radius):
    polygon = np.asarray([[xy[0]+radius,xy[1]+radius], [xy[0]+radius,xy[1]-radius], [xy[0]-radius,xy[1]-radius], [xy[0]-radius,xy[1]+radius]])

    keep_files = find_las_files_by_polygon(file_list,polygon)
    n_files = len(keep_files)
    if n_files == 0:
        print 'WARNING: No files within specified neighbourhood - try again'
    else:
        tile_pts = load_lidar_data(keep_files[0])
        pts = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
        for i in range(1,n_files):
            tile_pts = load_lidar_data(keep_files[i])
            tile_pts_filt = filter_lidar_data_by_neighbourhood(tile_pts,xy,radius)
            pts = np.concatenate((pts,tile_pts_filt),axis=0)

    print "loaded ", pts[:,0].size, " points"
    return pts
