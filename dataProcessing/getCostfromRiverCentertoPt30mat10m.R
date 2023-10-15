#Calculate cost of getting to all points outside the river 
#from the closest point in the middle of the river
#Use middle of river and outside points defined by 10m DEM
#calculate slope using 30m DEM

rm(list=ls())

library(sf)
library(terra)

################################################################################
#30m 

#load Schuylkill River shapefile
setwd("C:/Users/svr5482/Reification/Philly/data")
schuylkill_shp <- read_sf(dsn = ".", layer = "schuylkill10")

dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
crs(dem10)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m<- xyFromCell(dem10,1:ncell(dem10))

dem30<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_30m_new.asc")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
elev.30mat10m<- extract(dem30,coords.10m)

dem30at10<- dem10
values(dem30at10)<- elev.30mat10m
writeRaster(dem30at10,file="C:/Users/svr5482/Downloads/NewDEMS/norristown_30mat10m_new.tif")

mask10 <- rasterize(schuylkill_shp, dem10)
terra::plot(mask10)
coords.mask.10<- xyFromCell(mask10,1:ncell(mask10))
maskvals.10m<- extract(mask10,coords.mask.10)

riverbds.inds<- which(maskvals.10m$layer==1)

library(raster)
library(sp)

dem30<- raster::raster("C:/Users/svr5482/Downloads/NewDEMS/norristown_30mat10m_new.tif")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/middle_river_coords.10m.RData")

costs_RiverCentertoOutsidePts<- list()

for(i in 1:length(unique(coords.10m[,1]))){
  xval<- unique(coords.10m[,1])[i]
  x.inds<- which(coords.10m[,1]==xval)
  riverbds_atx.inds<- intersect(x.inds,riverbds.inds)
  yval<- coords.10m[riverbds_atx.inds,2]
  aboveriver_atx<- intersect(x.inds,which(coords.10m[,2]>max(yval)))
  belowriver_atx<- intersect(x.inds,which(coords.10m[,2]<min(yval)))
  outsideriver_atx<- c(aboveriver_atx,belowriver_atx) 
  #now we have the indices at this x value that are above or below the river
  
  A <- sp::SpatialPoints(cbind(middle_river_coords.10m[i,1], middle_river_coords.10m[i,2]))
  
  #costs<- rep(NA,length(outsideriver_atx))
  
  #for(j in 1:length(outsideriver_atx)){
  j=1
  B <- sp::SpatialPoints(cbind(coords.10m[outsideriver_atx[j],1],
                               coords.10m[outsideriver_atx[j],2]))
  cost30<- movecost(dtm= dem30, origin = A, destin = B, graph.out = FALSE)
  #costs[j]<- cost30$accumulated.cost.raster@data@values[outsideriver_atx[j]]
  costs<- cost30$accumulated.cost.raster@data@values[outsideriver_atx]
  #}
  costs_RiverCentertoOutsidePts[[i]]<- costs
  #print(i)
}

save(costs_RiverCentertoOutsidePts,
     file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_30mat10m.RData")


coords.30m<- xyFromCell(dem30,1:ncell(dem30))

min(coords.30m[,1]);max(coords.30m[,1])
min(coords.10m[,1]);max(coords.10m[,1])

min(coords.30m[,2])
min(coords.10m[,2])

max(coords.30m[,2])
max(coords.10m[,2])
