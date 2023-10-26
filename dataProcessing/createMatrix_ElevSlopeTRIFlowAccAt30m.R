#save elevation, slope, TRI, and flow accumulation for each DEM in matrices

rm(list=ls())

library(terra)

setwd("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD")

#load DEMs
################################################################################

dem30<- rast("norristown_30m_new.asc")
crs(dem30)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.30m<- xyFromCell(dem30,ncell(dem30))
elev.30m<- extract(dem30,coords.30m)

################################################################################

dem50<- rast("norristown_50m_new.asc")
crs(dem50)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
elev.50m<- extract(dem50,coords.30m)

###############################################################################################
#30m
FlowAcc_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/FlowAcc_n30m.tif")
coords.FlowAcc_n30m<- xyFromCell(FlowAcc_n30m,1:ncell(FlowAcc_n30m))
FlowAcc.30m<- extract(FlowAcc_n30m,coords.30m)
ncell.FlowAcc_n30m<- length(FlowAcc.30m$FlowAcc_n30m)

TRI_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/TRI_n30m.tif")
coords.TRI_n30m<- xyFromCell(TRI_n30m,1:ncell(TRI_n30m))
TRI.30m<- extract(TRI_n30m,coords.30m)
ncell.TRI_n30m<- length(TRI.30m$Band_1)

Slope_n30m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFlooding/Slope_n30m.tif")
coords.Slope_n30m<- xyFromCell(Slope_n30m,1:ncell(Slope_n30m))
Slope.30m<- extract(Slope_n30m,coords.30m)
ncell.Slope_n30m<- length(Slope.30m$Slope_n30m)

###############################################################################################
#50m
FlowAcc_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/FlowAcc_n50m.tif")
coords.FlowAcc_n50m<- xyFromCell(FlowAcc_n50m,1:ncell(FlowAcc_n50m))
FlowAcc.50m<- extract(FlowAcc_n50m,coords.30m)
ncell.FlowAcc_n50m<- length(FlowAcc.50m$FlowAcc_n50m)

TRI_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/TRI_n50m.tif")
coords.TRI_n50m<- xyFromCell(TRI_n50m,1:ncell(TRI_n50m))
TRI.50m<- extract(TRI_n50m,coords.30m)
ncell.TRI_n50m<- length(TRI.50m$Band_1)

Slope_n50m<- rast("C:/Users/svr5482/Documents/ArcGIS/Projects/PhillyFloodNewDEMs/Slope_n50m.tif")
coords.Slope_n50m<- xyFromCell(Slope_n50m,1:ncell(Slope_n50m))
Slope.50m<- extract(Slope_n50m,coords.30m)
ncell.Slope_n50m<- length(Slope.50m$Slope_n50m)

###############################################################################################

DEMinfo_30m<- cbind(coords.30m,elev.30m$norristown_30m,Slope.30m$Slope_n30m,TRI.30m$Band_1,FlowAcc.30m$FlowAcc_n30m)
colnames(DEMinfo_30m)<- c("x","y","elev","slope","TRI","flowacc")
save(DEMinfo_30m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30m.RData")

DEMinfo_50m<- cbind(coords.30m,elev.50m$norristown_50m,Slope.50m$Slope_n50m,TRI.50m$Band_1,FlowAcc.50m$FlowAcc_n50m)
colnames(DEMinfo_50m)<- c("x","y","elev","slope","TRI","flowacc")
save(DEMinfo_50m,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat30m.RData")
