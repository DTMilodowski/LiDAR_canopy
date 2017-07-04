# ==============================================================
#
# CanProf_poly is a function that enables a canopy profile
# to be calculated for a SpatialPolygonsDataFrame
#
# It is essentially a wrapper for several functions in 
# LAStools as well as David Milodowski's python based 
# canopy profiler. As such the functions requires file paths
# to the various applications, including lastools and python2.7.
#
# The remaining arguments are mostly for setting up the canopy
# profiler but remove_buffer is present to deal with the case
# of las tiles that contain buffers. These are dropped prior 
# to clipping to the area of interest. The cores argument is
# takes an integer referring to the number of cores used for 
# buffer removal using lastile.
#
# show.output.on.console enables the user to view the output from
# the console, which is very useful for identifying errors.
#
# Tom Swinfield
# 17-07-04
#
# ==============================================================

# x<-aoi_poly
# h.min <- 2
# h.max<-80
# layer.thickness <-1
# max.return<-3
#cores<-19

# aoi<-x[x$location=="R60-A-T29 (RD)",]
# las_path<-"C:/Temp/CanProf_tmp2"

CanProf_poly2 <-
  function(x,
           radius = 10, 
           h.min = 2,
           h.max = 80,
           layer.thickness = 1,
           leaf.angle.dist = "spherical",
           max.return=3,
           method= "macarthur_horn",
           las_path,
           LAStools_path,
           py_path,
           LiDAR_canopy_path,
           tmp_path ="C:/Temp/CanProf_tmp",
           out_path = "CanProf_output",
           out_overwrite = FALSE,
           remove_buffer = TRUE,
           cores = parallel::detectCores()-1,
           show.output.on.console = FALSE){
    radius<-as.character(radius)
    h.min<-as.character(h.min)
    h.max<-as.character(h.max)
    layer.thickness<-as.character(layer.thickness)
    max.return<-as.character(max.return)
    
    cat("================================================\n")
    cat("Extraction of canopy profile metrics\n")
    cat("Initiated: ", date(), "\n")
    
    if(dir.exists(tmp_path))
      file.remove(list.files(tmp_path, full.name = TRUE))
    else
      dir.create(tmp_path)
    
    if(!dir.exists(out_path))
      dir.create(out_path)
    else{
      if(length(list.files(out_path) != 0 & out_overwrite))
        stop(out_path, "folder already exists, please choose alternative directory")
    }
    
    # Load the catalogue of las files: -----
    ctg <- catalog(las_path)
    #rownames(ctg)<-paste("id_", rownames(ctg), sep ="")
    ctg_bbox<-ctg_sp(ctg) # convert to spatialpolygons
    proj4string(ctg_bbox)<-CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                               +towgs84=0,0,0")
    plot(ctg_bbox)
    ctg_df<-as.data.frame(ctg)
    
    # Find the tiles that match with the polygons:
    tile_match<-lapply(1:nrow(x), function(i) {
      gIntersects(x[i,], ctg_bbox, byid=TRUE)
    })
    # Test if any tiles match:
    intersect_test<-sapply(tile_match, any)
    # List the matching tiles for each polygon:
    poly_match<-lapply(tile_match, function(X){
      rownames(X)[X]
    })
    # Create the complete list of matching tiles without duplicates:
    tile_match<-unlist(poly_match)
    tile_match<-tile_match[!duplicated(tile_match)]
    
    if(remove_buffer){
      # Drop buffers from tiles: ----
      cat("------ Removing buffers -----\n")
      args<-paste(file.path(LAStools_path, "lastile"),
                  "-i",
                  paste(tile_match, collapse=" "),
                  "-remove_buffer",
                  "-odir",
                  tmp_path,
                  "-cores",
                  cores)
      system(args, show.output.on.console = show.output.on.console)
      cat("------ Buffers removed -----\n")
    }
    
    for(i in 1:nrow(x)){  
      aoi<-x[i,]
      id<-aoi$location
      if(show.output.on.console)
        cat("id: ", as.character(id), "\n")
      id<-gsub("([[:punct:]])", "_", id) # remove punctuation
      id<-gsub(" ", "", id) # remove spaces
      
      # save the aoi shapefile:
      writeOGR(aoi, dsn = tmp_path, "aoi",
               driver = "ESRI Shapefile", 
               overwrite_layer = TRUE)
      
      # Run the data extraction and processing only for polygons
      # with lidar data:
      if(intersect_test[i]){
        # point the file path to the right place for the tiles to use:
        if(remove_buffer){
          tile_path<-paste(file.path(tmp_path, basename(poly_match[[i]])),
                           collapse = " ")
          tile_path<-gsub("laz", "las", tile_path)
        }
        else
          tile_path<-paste(poly_match[[i]], collapse = " ")
      
        # Clip out area of interest: ----
        args<-paste(file.path(LAStools_path, "lasclip"),
                    "-i",
                    tile_path,
                    "-merged",
                    "-o",
                    file.path(tmp_path, "aoi.las"),
                    "-poly",
                    file.path(tmp_path, "aoi.shp"),
                    "-v")
        system(args, show.output.on.console = show.output.on.console)
        
        if(file.exists(file.path(tmp_path, "aoi.las"))){
          # Run DM lidar canopy profile driver: -----
          args<-paste(py_path,
                      file.path(LiDAR_canopy_path,
                                "simple_canopy_profile_driver_command_line_function.py"),
                      "0", "0",
                      file.path(tmp_path, "aoi.las"), 
                      radius, h.max, h.min, layer.thickness,
                      leaf.angle.dist,
                      max.return,
                      method,
                      paste(out_path, "/", id, ".csv", sep="")
                      )
          system(args, show.output.on.console = show.output.on.console)
        }
      }
      if(show.output.on.console)
        cat("\n")
    }
    system(paste("rm -rf", tmp_path)) # remove temporary folder
    
    cat("Complete: ", date(), "\n")
    cat("================================================\n")
    return(0)
  }
