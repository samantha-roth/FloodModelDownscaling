
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

#min.x.inds<- which(coords.10m[,1]==min(coords.10m[,1]))
#furthestleftriver.inds<- intersect(min.x.inds,riverbds.inds)

#y.inds<- coords.10m[furthestleftriver.inds,2]

#middle_river_coords.10m<- coords.10m[which.min(abs(coords.10m[,2]-mean(y.inds))),]

middle_river_coords.10m<- matrix(NA, nrow=length(unique(coords.10m[,1])), ncol=2)

for(i in 1:length(unique(coords.10m[,1]))){
  xval<- unique(coords.10m[,1])[i]
  x.inds<- which(coords.10m[,1]==xval)
  riverbds_atx.inds<- intersect(x.inds,riverbds.inds)
  yval<- coords.10m[riverbds_atx.inds,2]
  if(length(yval)==2){
    middle_river_coords.10m[i,2]<- coords.10m[which.min(abs(coords.10m[,2]-mean(yval))),2]
  }
  if(length(yval>2)){
    pick2<- order(yval)[1:2]
    middle_river_coords.10m[i,2]<- coords.10m[which.min(abs(coords.10m[,2]-mean(yval[pick2]))),2]
  }
  middle_river_coords.10m[i,1]<- xval
  
  #print(x.inds)
}


save(middle_river_coords.10m, file= "C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/middle_river_coords.10m.RData")


plot(mask10)
points(x=middle_river_coords.10m[,1],y=middle_river_coords.10m[,2])
################################################################################
#30m 

#load Schuylkill River shapefile
setwd("C:/Users/svr5482/Reification/Philly/data")
schuylkill_shp <- read_sf(dsn = ".", layer = "schuylkill30")

dem30<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_30m_new.asc")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.30m<- xyFromCell(dem30,1:ncell(dem30))
elev.30m<- extract(dem30,coords.30m)
ncell.30m<- length(elev.30m$norristown_30m)


mask30 <- rasterize(schuylkill_shp, dem30)
terra::plot(mask30)
coords.mask.30<- xyFromCell(mask30,1:ncell(mask30))
maskvals.30m<- extract(mask30,coords.mask.30)

riverbds.inds<- which(maskvals.30m$layer==1)

#min.x.inds<- which(coords.30m[,1]==min(coords.30m[,1]))
#furthestleftriver.inds<- intersect(min.x.inds,riverbds.inds)

#y.inds<- coords.30m[furthestleftriver.inds,2]

#middle_river_coords.30m<- coords.30m[which.min(abs(coords.30m[,2]-mean(y.inds))),]

middle_river_coords.30m<- matrix(NA, nrow=length(unique(coords.30m[,1])), ncol=2)

for(i in 1:length(unique(coords.30m[,1]))){
  xval<- unique(coords.30m[,1])[i]
  x.inds<- which(coords.30m[,1]==xval)
  riverbds_atx.inds<- intersect(x.inds,riverbds.inds)
  yval<- coords.30m[riverbds_atx.inds,2]
  if(length(yval)==2){
    middle_river_coords.30m[i,2]<- coords.30m[which.min(abs(coords.30m[,2]-mean(yval))),2]
  }
  if(length(yval>2)){
    pick2<- order(yval)[1:2]
    middle_river_coords.30m[i,2]<- coords.30m[which.min(abs(coords.30m[,2]-mean(yval[pick2]))),2]
  }
  middle_river_coords.30m[i,1]<- xval
  
  #print(x.inds)
}


save(middle_river_coords.30m, file= "C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/middle_river_coords.30m.RData")


plot(mask30)
points(x=middle_river_coords.30m[,1],y=middle_river_coords.30m[,2])
################################################################################
#50m 

#load Schuylkill River shapefile
setwd("C:/Users/svr5482/Reification/Philly/data")
schuylkill_shp <- read_sf(dsn = ".", layer = "schuylkill50")

dem50<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_50m_new.asc")
crs(dem50)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.50m<- xyFromCell(dem50,1:ncell(dem50))
elev.50m<- extract(dem50,coords.50m)
ncell.50m<- length(elev.50m$norristown_50m)


mask50 <- rasterize(schuylkill_shp, dem50)
terra::plot(mask50)
coords.mask.50<- xyFromCell(mask50,1:ncell(mask50))
maskvals.50m<- extract(mask50,coords.mask.50)

riverbds.inds<- which(maskvals.50m$layer==1)

#min.x.inds<- which(coords.50m[,1]==min(coords.50m[,1]))
#furthestleftriver.inds<- intersect(min.x.inds,riverbds.inds)

#y.inds<- coords.50m[furthestleftriver.inds,2]

#middle_river_coords.50m<- coords.50m[which.min(abs(coords.50m[,2]-mean(y.inds))),]

middle_river_coords.50m<- matrix(NA, nrow=length(unique(coords.50m[,1])), ncol=2)

for(i in 1:length(unique(coords.50m[,1]))){
  xval<- unique(coords.50m[,1])[i]
  x.inds<- which(coords.50m[,1]==xval)
  riverbds_atx.inds<- intersect(x.inds,riverbds.inds)
  yval<- coords.50m[riverbds_atx.inds,2]
  if(length(yval)==2){
    middle_river_coords.50m[i,2]<- coords.50m[which.min(abs(coords.50m[,2]-mean(yval))),2]
  }
  if(length(yval>2)){
    pick2<- order(yval)[1:2]
    middle_river_coords.50m[i,2]<- coords.50m[which.min(abs(coords.50m[,2]-mean(yval[pick2]))),2]
  }
  middle_river_coords.50m[i,1]<- xval
  
  #print(x.inds)
}


save(middle_river_coords.50m, file= "C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/middle_river_coords.50m.RData")


plot(mask50)
points(x=middle_river_coords.50m[,1],y=middle_river_coords.50m[,2])
