rm(list=ls())

nHWMs=7
nRuns=200

library(sf)
library(terra)

#load Schuylkill River shapefile
setwd("C:/Users/svr5482/Reification/Philly/data")
schuylkill_shp <- read_sf(dsn = ".", layer = "schuylkill3")


Run3m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs3m/Norristown/nCh/test45/Extent/Run_45.asc")
crs(Run3m)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.3m<- xyFromCell(Run3m,1:ncell(Run3m))
elev.3m<- extract(Run3m,coords.3m)
ncell.3m<- length(elev.3m$norristown_3m)


mask5 <- rasterize(schuylkill_shp, Run3m)
#plot(mask5)
coords.mask.5<- xyFromCell(mask5,1:ncell(mask5))
maskvals.3m<- extract(mask5,coords.mask.5)

writeRaster(mask5, filename= "C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/river3m.tif", overwrite=TRUE)

#m5<- rast("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/river3m.tif")

#find the grid cells where every cell above it in its column is NA and it is NA

yvals.3m<- unique(coords.3m[,2])
xvals.3m<- unique(coords.3m[,1])

above.river.inds<- c()

for(x in 1:length(xvals.3m)){
  xcoords<- which(coords.3m[,1]==xvals.3m[x])
  maskvals.x<- maskvals.3m$layer[xcoords]
  coords.x<- coords.3m[xcoords,]
  one.inds<- which(maskvals.x==1)
  max.y<- max(coords.x[one.inds,2])
  x.above.river.inds<- intersect(xcoords,which(coords.3m[,2]>max.y))
  above.river.inds<- c(above.river.inds,x.above.river.inds)
}

below.river.inds<- c()

for(x in 1:length(xvals.3m)){
  xcoords<- which(coords.3m[,1]==xvals.3m[x])
  maskvals.x<- maskvals.3m$layer[xcoords]
  coords.x<- coords.3m[xcoords,]
  one.inds<- which(maskvals.x==1)
  min.y<- min(coords.x[one.inds,2])
  x.below.river.inds<- intersect(xcoords,which(coords.3m[,2]<min.y))
  below.river.inds<- c(below.river.inds,x.below.river.inds)
}


not.river.inds3m<- c(above.river.inds,below.river.inds)

plot(coords.3m[not.river.inds3m,1],coords.3m[not.river.inds3m,2])

save(not.river.inds3m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/not.river.inds3m.RData")
