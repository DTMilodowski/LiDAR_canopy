bbox_make<-function(min.x, max.x, min.y, max.y, ID){
  coords<-matrix(0, nrow =5, ncol=2, dimnames = list(NULL, c("x", "y")))
  coords[,"x"]<-c(min.x, min.x, max.x, max.x, min.x)
  coords[,"y"]<-c(min.y, max.y, max.y, min.y, min.y)
  coords<-Polygons(list(Polygon(coords)), ID = ID)
  return(coords)
}

ctg_sp<-function(x){
  x_ext<-as.data.frame(x[,c("Min.X", "Max.X", "Min.Y", "Max.Y")])
  ID<-x$filename
  #  ID<-rownames(x)
  x_bbox<-lapply(1:nrow(x_ext), function(id) {
    e<-as.numeric(x_ext[id,])
    bbox<-bbox_make(e[1], e[2], e[3], e[4], 
                    ID=ID[id])
    return(bbox)
  })
  x_bbox<-sp::SpatialPolygons(x_bbox)
  return(x_bbox)
}