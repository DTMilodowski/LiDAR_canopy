# Some functions to read/write geoTiffs
import numpy as np
from osgeo import gdal
import os
import osr
import sys


# Function to write geotiff from an array
def write_raster_to_GeoTiff(array,geoTrans, OUTFILE_prefix, EPSG_CODE='4326', north_up=True):
    #-----------------------------------
    # some dummy variables to be updated
    NBands = 1
    NRows = 0
    NCols = 0
    #-----------------------------------

    #-----------------------------------
    # Check orientation of array before writing to raster to ensure 
    # compatibility with GIS platforms.
    if north_up:
        # for north_up array, need the n-s resolution (element 5) to be negative
        if geoTrans[5]>0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2: 
            print 'array has less than two dimensions! Unable to write to raster'
            sys.exit(1)  
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print 'array has too many dimensions! Unable to write to raster'
            sys.exit(1)  

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2: 
            print 'array has less than two dimensions! Unable to write to raster'
            sys.exit(1)  
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print 'array has too many dimensions! Unable to write to raster'
            sys.exit(1)  
    
    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2: 
        print 'array has less than two dimensions! Unable to write to raster'
        sys.exit(1)  
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print 'array has too many dimensions! Unable to write to raster'
        sys.exit(1)  
    #-----------------------------------

    #-----------------------------------
    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'_data.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    dataset.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset.GetRasterBand(1).WriteArray( array )
    dataset = None
    return 0
    #-----------------------------------
