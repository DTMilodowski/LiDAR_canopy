# ==============================================================
#
#
#
#
#
#
# ==============================================================

CanProf_poly <-
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
           LiDAR_canopy_path) {
  radius<-as.character(radius)
  h.min<-as.character(h.min)
  h.max<-as.character(h.max)
  layer.thickness<-as.character(layer.thickness)
  max.return<-as.character(max.return)
  
  cat("================================================\n")
  cat("Extraction of canopy profile metrics\n")
  cat("Initiated: ", date(), "\n")
  
  if(dir.exists("C:/Temp/CanProf_tmp"))
    file.remove(list.files("C:/Temp/CanProf_tmp", full.name = TRUE))
  else
    dir.create("C:/Temp/CanProf_tmp")
  
  if(!dir.exists("CanProf_output"))
    dir.create("CanProf_output")
  
  # Load the catalogue of las files: -----
  ctg <- catalog(las_path)
  #rownames(ctg)<-paste("id_", rownames(ctg), sep ="")
  ctg_bbox<-ctg_sp(ctg) # convert to spatialpolygons
  proj4string(ctg_bbox)<-CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                             +towgs84=0,0,0")
  plot(ctg_bbox)
  ctg_df<-as.data.frame(ctg)
  
  for(i in 1:nrow(x)){  
    aoi<-x[i,]
    id<-aoi$location
    cat("id: ", as.character(id), "\n")
    #id<-gsub("[[:punct:]]", " ", id) # remove punctuation
    
    writeOGR(aoi, dsn = "C:/Temp/CanProf_tmp", "aoi",
             driver = "ESRI Shapefile", 
             overwrite_layer = TRUE)
    
    # Check for intersection: -----
    tile_match<-gIntersects(aoi, ctg_bbox, byid=TRUE)
    tile_match<-rownames(tile_match)[tile_match]
    
    if(length(tile_match)>1){
    # Drop buffers from tiles: ----
      args<-paste(file.path(LAStools_path, "lastile"),
                  "-i",
                  paste(tile_match, collapse=" "),
                  "-remove_buffer",
                  "-odir",
                  "C:/Temp/CanProf_tmp",
                  "-cores",
                  "1")
      system(args)
    
      # Clip out area of interest: ----
      args<-paste(file.path(LAStools_path, "lasclip"),
                  "-i",
                  paste(list.files("C:/Temp/CanProf_tmp",
                                   pattern = "*.(las|laz)",
                                   full.names = TRUE),
                        collapse = " "),
                  "-merged",
                  "-o",
                  "C:/Temp/CanProf_tmp/aoi.las",
                  "-poly",
                  "C:/Temp/CanProf_tmp/aoi.shp",
                  "-v")
      
    }
    else{
      args<-paste(file.path(LAStools_path, "lasclip"),
                  "-i",
                  tile_match,
                  "-o",
                  "C:/Temp/CanProf_tmp/aoi.las",
                  "-poly",
                  "C:/Temp/CanProf_tmp/aoi.shp",
                  "-v")
    }
    # Clip out area of interest: ----
    system(args)
      
    # Run DM lidar canopy profile driver: -----
    
    args<-paste("/Python27/python",
                file.path(LiDAR_canopy_path,
                          "simple_canopy_profile_driver_command_line_function.py"),
                "0", "0",
                "C:/Temp/CanProf_tmp/aoi.las", 
                radius, h.max, h.min, layer.thickness,
                leaf.angle.dist,
                max.return,
                method,
                paste("CanProf_output/", id, ".csv", sep="")
    )
    system(args)
    file.remove(list.files("C:/Temp/CanProf_tmp", full.name = TRUE))
    cat("\n")
  }
  system("rm -rf C:/Temp/CanProf_tmp") # remove temporary folder
  
  
  cat("Complete: ", date(), "\n")
  cat("================================================\n")
  return(0)
}
