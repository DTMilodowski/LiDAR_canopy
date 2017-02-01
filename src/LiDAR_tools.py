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

