rm(list=ls())

################################################################################
#10m

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

#load calibrated runs
Run10m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/Norristown/nCh/simplecal/Extent/Run_1.asc")
coords.10m<- xyFromCell(Run10m,1:ncell(Run10m))
flood.10m<- extract(Run10m,coords.10m)

binFlood.10m<- ifelse(flood.10m$Run_1>0,1,0)

plot(DEMinfo_10m[,"elev"],binFlood.10m)
floodInds<- which(binFlood.10m==1)

maxElevFlood10m<- max(DEMinfo_10m[floodInds,"elev"])
minElevNoFlood10m<- min(DEMinfo_10m[-floodInds,"elev"])

save(floodInds,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/floodInds10m.RData")

################################################################################
#30m

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30m.RData")

#load calibrated runs
Run30m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs30m/Norristown/nCh/simplecal/Extent/Run_1.asc")
coords.30m<- xyFromCell(Run30m,1:ncell(Run30m))
flood.30m<- extract(Run30m,coords.30m)

binFlood.30m<- ifelse(flood.30m$Run_1>0,1,0)

plot(DEMinfo_30m[,"elev"],binFlood.30m)
floodInds<- which(binFlood.30m==1)

maxElevFlood30m<- max(DEMinfo_30m[floodInds,"elev"])
minElevNoFlood30m<- min(DEMinfo_30m[-floodInds,"elev"])

save(floodInds,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/floodInds30m.RData")

################################################################################
#50m

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50m.RData")

#load calibrated runs
Run50m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/simplecal/Extent/Run_1.asc")
coords.50m<- xyFromCell(Run50m,1:ncell(Run50m))
flood.50m<- extract(Run50m,coords.50m)

binFlood.50m<- ifelse(flood.50m$Run_1>0,1,0)
plot(DEMinfo_50m[,"elev"],binFlood.50m)
floodInds<- which(binFlood.50m==1)

maxElevFlood50m<- max(DEMinfo_50m[floodInds,"elev"])
minElevNoFlood50m<- min(DEMinfo_50m[-floodInds,"elev"])

save(floodInds,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/floodInds50m.RData")

