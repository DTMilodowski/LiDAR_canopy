rm(list= ls())

library(lidR)
library(raster)
library(sp)
library(rgdal)
library(rgeos)

# Load wrapper functions:
source("R/Functions/CanProf_funs.R")
source("R/Functions/CanProf_poly2.R")

# Set file paths: ----
las_path<-"C:/Data/Malaysia/Processed/SAFE/SAFE_CHM_laz"
LAStools_path<-"C:/Users/Forecol/Documents/LAStools/bin"
LiDAR_canopy_path<-"LiDAR_canopy/src"
py_path<-"/Python27/python"
out_path<- "CanProf_output"
#las_path<-tmp_path

# Set global variables:----
id<-"Cam_Loc__I"
radius <- 20

# Load points and covert to polygons:----
aoi_points<-readOGR("Data/shapefiles/CT_Riparian_prj.shp")
aoi_points<-aoi_points[-111,] # remove one strange point
loc_id<-aoi_points@data[,id]

plot(aoi_points)
aoi_points<-gBuffer(aoi_points, width = radius, byid=TRUE, capStyle = "ROUND")
aoi_data<-data.frame(location = loc_id)
rownames(aoi_data)<-row.names(aoi_points)
aoi_poly<-SpatialPolygonsDataFrame(aoi_points, 
                                   data = aoi_data)

# Run canopy profile extraction:----
CanProf_poly2(aoi_poly,
             radius =  radius,
             las_path=las_path,
             LAStools_path = LAStools_path,
             py_path = py_path,
             LiDAR_canopy_path = LiDAR_canopy_path,
             remove_buffer = TRUE)

