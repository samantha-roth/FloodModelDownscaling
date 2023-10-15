#Calculate cost of getting to all points outside the river 
#from the closest point in the middle of the river

rm(list=ls())

library(sf)
library(terra)

################################################################################
#10m 

#load Schuylkill River shapefile
setwd("C:/Users/svr5482/Reification/Philly/data")
schuylkill_shp <- read_sf(dsn = ".", layer = "schuylkill10")

dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
crs(dem10)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m<- xyFromCell(dem10,1:ncell(dem10))
elev.10m<- extract(dem10,coords.10m)
ncell.10m<- length(elev.10m$norristown_10m)

mask10 <- rasterize(schuylkill_shp, dem10)
terra::plot(mask10)
coords.mask.10<- xyFromCell(mask10,1:ncell(mask10))
maskvals.10m<- extract(mask10,coords.mask.10)

riverbds.inds<- which(maskvals.10m$layer==1)

library(raster)
library(sp)
library(movecost)

dem10<- raster::raster("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
crs(dem10)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/middle_river_coords.10m.RData")

costs_RiverCentertoOutsidePts<- list()
coords_OutsideRiverbyX<- list()

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
    cost10<- movecost(dtm= dem10, origin = A, destin = B, graph.out = FALSE)
    #costs[j]<- cost10$accumulated.cost.raster@data@values[outsideriver_atx[j]]
    costs<- cost10$accumulated.cost.raster@data@values[outsideriver_atx]
  #}
  costs_RiverCentertoOutsidePts[[i]]<- costs
  #print(i)
  
  coords_OutsideRiverbyX[[i]]<- cbind(outsideriver_atx,coords.10m[outsideriver_atx,])
}

save(costs_RiverCentertoOutsidePts,
     file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_10m.RData")

save(coords_OutsideRiverbyX,
     file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/coords_OutsideRiverbyX_10m.RData")


