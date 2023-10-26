
rm(list=ls())

library(terra)


setwd("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD")


################################################################################

dem10<- rast("norristown_10m_new.asc")
crs(dem10)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m<- xyFromCell(dem10,1:ncell(dem10))
elev.10m<- extract(dem10,coords.10m)

################################################################################

dem30<- rast("norristown_30m_new.asc")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
elev.30m<- extract(dem30,coords.10m)

################################################################################

dem50<- rast("norristown_50m_new.asc")
crs(dem50)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
elev.50m<- extract(dem50,coords.10m)

###############################################################################################
#10m
FlowAcc_n10m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/FlowAcc_n10m.tif")
FlowAcc.10m<- extract(FlowAcc_n10m,coords.10m)

TRI_n10m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/TRI_n10m.tif")
TRI.10m<- extract(TRI_n10m,coords.10m)

Slope_n10m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/Slope_n10m.tif")
Slope.10m<- extract(Slope_n10m,coords.10m)

###############################################################################################
#30m
FlowAcc_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/FlowAcc_n30m.tif")
FlowAcc.30m<- extract(FlowAcc_n30m,coords.10m)

TRI_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/TRI_n30m.tif")
TRI.30m<- extract(TRI_n30m,coords.10m)

Slope_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFlooding/Slope_n30m.tif")
Slope.30m<- extract(Slope_n30m,coords.10m)

###############################################################################################
#50m

FlowAcc_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/FlowAcc_n50m.tif")
FlowAcc.50m<- extract(FlowAcc_n50m,coords.10m)

TRI_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/TRI_n50m.tif")
TRI.50m<- extract(TRI_n50m,coords.10m)

Slope_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFlooding/Slope_n50m.tif")
Slope.50m<- extract(Slope_n50m,coords.10m)

###############################################################################################

DEMinfo_10m<- cbind(coords.10m,elev.10m$norristown_10m,Slope.10m$Slope_n10m,TRI.10m$Band_1,FlowAcc.10m$FlowAcc_n10m)
colnames(DEMinfo_10m)<- c("x","y","elev","slope","TRI","flowacc")
save(DEMinfo_10m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_30m<- cbind(coords.10m,elev.30m$norristown_30m,Slope.30m$Slope_n30m,TRI.30m$Band_1,FlowAcc.30m$FlowAcc_n30m)
colnames(DEMinfo_30m)<- c("x","y","elev","slope","TRI","flowacc")
save(DEMinfo_30m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30mat10m.RData")

DEMinfo_50m<- cbind(coords.10m,elev.50m$norristown_50m,Slope.50m$Slope_n50m,TRI.50m$Band_1,FlowAcc.50m$FlowAcc_n50m)
colnames(DEMinfo_50m)<- c("x","y","elev","slope","TRI","flowacc")
save(DEMinfo_50m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat10m.RData")
