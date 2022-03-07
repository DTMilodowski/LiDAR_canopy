import numpy as np
from scipy import stats
###################################################################################
# This set of functions provides some auxilliary tools needed for some of the
# LAD profiles analysis
# @author DTMilodowski
###################################################################################


# Load spreadsheet of LAI derived from hemispherical photographs (courtesy of Terhi Riutta at Oxford).  LAI estimated using Hemisfer.
def load_field_LAI(LAI_file):
    datatype = {'names': ('ForestType','Plot', 'Subplot', 'LAI'), 'formats': ('S32','S32','i8','f16')}
    hemisfer_LAI = np.genfromtxt(LAI_file, skip_header = 1, delimiter = ',',dtype=datatype)

    return hemisfer_LAI

# This function loads the subplot coordinates from a csv file.  File columns should be as follows:
# Plot Subplot 'X0', 'Y0', 'X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3'
def load_boundaries(coordinate_list):

    datatype = {'names': ('Plot', 'Subplot', 'X0', 'Y0', 'X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3'), 'formats': ('U32','i8','f16','f16','f16','f16','f16','f16','f16','f16')}
    coords = np.genfromtxt(coordinate_list, skip_header = 1, delimiter = ',',dtype=datatype)
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

# This function is similar to above but now will read an arbitrary number of vertices so that generic polygons can be used
def load_generic_boundaries(coordinate_list):

    f = open(coordinate_list,'r')
    lines=f.readlines()[0:]
    N = len(lines)

    coordinate_dict = {}
    for i in range(1,N):
        line = lines[i].strip().split(",")
        Ncols = len(line)
        plot = line[0]
        vertices = []#np.zeros((Ncols-1,2))
        vertices.append([])
        vertices.append([])
        for j in range(0,(Ncols-1)/2):
            if len(line[2*j+1])>0:
                vertices[0].append(float(line[2*j+1]))
                vertices[1].append(float(line[2*j+2]))
        coordinate_dict[plot]=np.asarray(vertices).transpose()
    f.close()
    return coordinate_dict

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

# get bounding box for list of coordinates
def get_bounding_box_with_buffer(coordinate_list,buffer_width):
    N = coordinate_list.shape[0]
    bbox = np.zeros((4,2))
    top = coordinate_list[:,1].max()+buffer_width
    bottom = coordinate_list[:,1].min()-buffer_width
    left = coordinate_list[:,0].min()-buffer_width
    right = coordinate_list[:,0].max()+buffer_width
    bbox[0,0]=left
    bbox[0,1]=top
    bbox[1,0]=right
    bbox[1,1]=top
    bbox[2,0]=right
    bbox[2,1]=bottom
    bbox[3,0]=left
    bbox[3,1]=bottom
    return bbox


# write R-squared and p-value
def get_rsquared_annotation(x,y):

    m, c, r, p, err = stats.linregress(x,y)
    p_str = ''
    if  p <= 0.001:
        p_str='***'
    elif  p <= 0.01:
        p_str='**'
    elif  p <= 0.05:
        p_str='*'
    elif p <= 0.1:
        p_str='$^.$'

    annotation = '$R^2$ = %.3f' % r**2
    annotation = annotation + p_str
    return annotation
