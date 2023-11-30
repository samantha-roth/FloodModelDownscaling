#trim the shp file containing the rivers in PA
rm(list=ls())

library(sf)
library(terra)

#load DEMs at 5m, 10m, 30m, and 50m resolution
################################################################################

dem3<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_3m_resample.asc")
crs(dem3)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.3m<- xyFromCell(dem3,1:ncell(dem3))
elev.3m<- extract(dem3,coords.3m)
ncell.3m<- length(elev.3m$norristown_3m)

################################################################################

dem5<- rast("C:/Users/svr5482/Downloads/norristown/norristown/norristown_5m.asc")
crs(dem5)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.5m<- xyFromCell(dem5,1:ncell(dem5))
elev.5m<- extract(dem5,coords.5m)
ncell.5m<- length(elev.5m$norristown_5m)

################################################################################

dem10<- rast("C:/Users/svr5482/Downloads/norristown_10m_agg.asc")
crs(dem10)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m<- xyFromCell(dem10,1:ncell(dem10))
elev.10m<- extract(dem10,coords.10m)
ncell.10m<- length(elev.10m$norristown_10m)

################################################################################

dem30<- rast("C:/Users/svr5482/Downloads/norristown_30m_agg.asc")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.30m<- xyFromCell(dem30,1:ncell(dem30))
elev.30m<- extract(dem30,coords.30m)
ncell.30m<- length(elev.30m$norristown_30m)

################################################################################

dem50<- rast("C:/Users/svr5482/Downloads/norristown_50m_agg.asc")
crs(dem50)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.50m<- xyFromCell(dem50,1:ncell(dem50))
elev.50m<- extract(dem50,coords.50m)
ncell.50m<- length(elev.50m$norristown_50m)

################################################################################

#load rivers shapefile
setwd("C:/Users/svr5482/Reification/Philly/data/majrivrs")
allPArivers <- read_sf(dsn = ".", layer = "majrivrs")

#plot(allPArivers)

# Transform the CRS of the shapefile to match that of the dem10 SpatRaster object
#allPArivers <- st_transform(allPArivers, crs(dem10))

################################################################################
#find the maxes of the mins and the mins of the maxes

c(min(coords.50m[,1]),min(coords.30m[,1]),min(coords.10m[,1]))
c(min(coords.50m[,2]),min(coords.30m[,2]),min(coords.10m[,2]))

keepX10m<- intersect(which(coords.10m[,1]>=min(coords.50m[,1])),which(coords.10m[,1]<=max(coords.50m[,1])))
keepY10m<- intersect(which(coords.10m[,2]>=min(coords.50m[,2])),which(coords.10m[,2]<=max(coords.50m[,2])))
keep10m<- intersect(keepX10m,keepY10m)

keepX30m<- intersect(which(coords.30m[,1]>=min(coords.50m[,1])),which(coords.30m[,1]<=max(coords.50m[,1])))
keepY30m<- intersect(which(coords.30m[,2]>=min(coords.50m[,2])),which(coords.30m[,2]<=max(coords.50m[,2])))
keep30m<- intersect(keepX30m,keepY30m)

coordsIwant.10m<- coords.10m[keep10m,]
coordsIwant.30m<- coords.30m[keep30m,]

################################################################################
allPArivers3 <- st_transform(allPArivers, crs(dem3))

# Define the spatial domain to crop to (a rectangle in this example)

bbox3 <- st_bbox(dem3)

schuylkill_shp3<- st_crop(allPArivers3, bbox3)

st_write(schuylkill_shp3, "C:/Users/svr5482/Reification/Philly/data/schuylkill3.shp",append=FALSE)

################################################################################

# Define the spatial domain to crop to (a rectangle in this example)

bbox5 <- st_bbox(dem5)

schuylkill_shp5<- st_crop(allPArivers, bbox5)
#plot(schuylkill_shp5)

st_write(schuylkill_shp5, "C:/Users/svr5482/Reification/Philly/data/schuylkill5.shp",append=FALSE)

################################################################################

# Define the spatial domain to crop to (a rectangle in this example)

bbox10 <- st_bbox(dem10)

schuylkill_shp10<- st_crop(allPArivers, bbox10)
#plot(schuylkill_shp10)

st_write(schuylkill_shp10, "C:/Users/svr5482/Reification/Philly/data/schuylkill10.shp",append=FALSE)
################################################################################

# Define the spatial domain to crop to (a rectangle in this example)

bbox30 <- st_bbox(dem30)

schuylkill_shp30<- st_crop(allPArivers, bbox30)
#plot(schuylkill_shp30)

st_write(schuylkill_shp30, "C:/Users/svr5482/Reification/Philly/data/schuylkill30.shp",append=FALSE)
################################################################################

# Define the spatial domain to crop to (a rectangle in this example)

bbox50 <- st_bbox(dem50)

schuylkill_shp50<- st_crop(allPArivers, bbox50)
#plot(schuylkill_shp50)

st_write(schuylkill_shp50, "C:/Users/svr5482/Reification/Philly/data/schuylkill50.shp",append=FALSE)
