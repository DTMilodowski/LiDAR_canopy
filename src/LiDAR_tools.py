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

# filter lidar wth polygon
# This function has been updated to include an option to filter by first return location.
# The reason for this is so full collections of returns associated with each LiDAR pulse
# can be retrieved, which can be an issue at edges in multi-return analyses
def filter_lidar_data_by_polygon(in_pts,polygon,filter_by_first_return_location = False):
    pts = np.zeros((0,in_pts.shape[1]))
    if in_pts.shape[0]>0:
        if filter_by_first_return_location:
            # find first returns
            mask = in_pts[:,3]==1
            x_temp, y_temp, inside_temp = points_in_poly(in_pts[mask,0],in_pts[mask,1],polygon)
            shots = np.unique(in_pts[mask,6][inside_temp]) # index 6 refers to GPS time
            inside = np.in1d(in_pts[:,6],shots) # this function retrieves all points corresponding to this GPS time
            x = in_pts[inside,0]
            y = in_pts[inside,1]
            x_temp=None
            y_temp=None
            inside_temp=None
        else:
            x,y,inside = points_in_poly(in_pts[:,0],in_pts[:,1],polygon)
        pts = in_pts[inside,:]
    else:
        print("\t\t\t no points in polygon")
    return pts

# filter lidar by circular neighbourhood
def filter_lidar_data_by_neighbourhood(in_pts,target_xy,radius):
    pts = np.zeros((0,in_pts.shape[1]))
    if in_pts.shape[0]>0:
        x,y,inside =  points_in_radius(in_pts[:,0],in_pts[:,1],target_xy[0],target_xy[1],radius)
        pts = in_pts[inside,:]
    else:
        print( "\t\t\t no points in neighbourhood")
    return pts
