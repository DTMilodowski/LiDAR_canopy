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
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2:
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2:
        print('array has less than two dimensions! Unable to write to raster')
        sys.exit(1)
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print('array has too many dimensions! Unable to write to raster')
        sys.exit(1)
    #-----------------------------------

    #-----------------------------------
    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    if NBands==1:
        dataset.GetRasterBand(1).SetNoDataValue( -9999 )
        dataset.GetRasterBand(1).WriteArray( array )
    else:
        for bb in range(0,NBands):
            dataset.GetRasterBand(bb+1).SetNoDataValue( -9999 )
            dataset.GetRasterBand(bb+1).WriteArray( array[:,:,bb] )
    dataset = None
    return 0
    #-----------------------------------


# Function to write geotiff from an array
def write_raster_to_GeoTiff_UTM(array,geoTrans, OUTFILE_prefix, utm_zone, northern_hemisphere = True, datum = "WGS84", north_up=True):

    if northern_hemisphere == True:
        hem = "northern"
        bnorth = 1
    else:
        hem = "southern"
        bnorth = 0

    prjCS = "UTM %i (%s) in %s hemisphere." % (utm_zone, datum, hem)

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
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            #array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            #for i in range(0,NBands):
            #    array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2:
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2:
        print('array has less than two dimensions! Unable to write to raster')
        sys.exit(1)
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print('array has too many dimensions! Unable to write to raster')
        sys.exit(1)
    #-----------------------------------

    #-----------------------------------
    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetProjCS( prjCS )
    srs.SetWellKnownGeogCS( datum )
    srs.SetUTM( utm_zone, bnorth )
    #srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    if NBands==1:
        dataset.GetRasterBand(1).SetNoDataValue( -9999 )
        dataset.GetRasterBand(1).WriteArray( array )
    else:
        for bb in range(0,NBands):
            dataset.GetRasterBand(bb+1).SetNoDataValue( -9999 )
            dataset.GetRasterBand(bb+1).WriteArray( array[:,:,bb] )
    dataset = None
    return 0
    #-----------------------------------

# Simple script to load GeoTIFF
def load_GeoTIFF_band_and_georeferencing(File,band_number=1):

    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        ds = gdal.Open(File)
    except(RuntimeError, e):
        print('unable to open ' + File)
        print(e)
        sys.exit(1)

    source_band = ds.GetRasterBand(band_number)
    if source_band is None:
        print("BAND MISSING")
        sys.exit(1)

    array = np.array(ds.GetRasterBand(band_number).ReadAsArray(),dtype=np.float64)
    geoTrans = ds.GetGeoTransform()
    coord_sys = ds.GetProjectionRef()

    if geoTrans[-1] < 0:
        array=np.flipud(array)

    return array, geoTrans, coord_sys



# Function to write an array to a geoTIFF
def write_array_to_GeoTiff_with_coordinate_system(array,geoTrans,coord_sys,OUTFILE,north_up=True):
    NBands = 1
    NRows = 0
    NCols = 0

    if north_up:
        # for north_up array, need the n-s resolution (element 5) to be negative
        if geoTrans[5]>0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2:
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2:
            print('array has less than two dimensions! Unable to write to raster')
            sys.exit(1)
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print('array has too many dimensions! Unable to write to raster')
            sys.exit(1)

    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2:
        print('array has less than two dimensions! Unable to write to raster')
        sys.exit(1)
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print('array has too many dimensions! Unable to write to raster')
        sys.exit(1)

    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE, NCols, NRows, NBands, gdal.GDT_Float32 )
    print("Generating raster %s with %i bands" % ((OUTFILE),NBands))
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference(wkt=coord_sys)
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    if len(array.shape) == 2:
        dataset.GetRasterBand(1).SetNoDataValue( -9999 )
        dataset.GetRasterBand(1).WriteArray( array )
    elif len(array.shape) == 3:
        for bb in range(0,NBands):
            dataset.GetRasterBand(bb+1).SetNoDataValue( -9999 )
            dataset.GetRasterBand(bb+1).WriteArray( array[:,:,bb] )
    dataset = None
    return 0
