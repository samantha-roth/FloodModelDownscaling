rm(list=ls())

library(terra)

nRuns=200
################################################################################

dem5<- rast("C:/Users/svr5482/Downloads/norristown/norristown/norristown_5m.asc")
crs(dem5)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.5m<- xyFromCell(dem5,1:ncell(dem5))
elev.5m<- extract(dem5,coords.5m)
ncell.5m<- length(elev.5m$norristown_5m)

################################################################################
#load the indices that are not in the river at 5m resolution
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/not.river.inds5m.RData")

coordsIwant.5m<- coords.5m[not.river.inds5m,]

SS<-nrow(coordsIwant.5m)

################################################################################

closest10mto5m_noRiver<- matrix(NA, ncol= nRuns, nrow= SS)

for(i in 1:nRuns){
  run<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/Norristown/nCh/Extent/Run_",i,".asc"))
  vals<- extract(run,coordsIwant.5m)
  closest10mto5m_noRiver[,i]<- vals$max
}

save(closest10mto5m_noRiver,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest10mto5m_noRiver.RData")

################################################################################

closest30mto5m_noRiver<- matrix(NA, ncol= nRuns, nrow= SS)

for(i in 1:nRuns){
  run<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs30m/Norristown/nCh/Extent/Run_",i,".asc"))
  vals<- extract(run,coordsIwant.5m)
  closest30mto5m_noRiver[,i]<- vals$max
}

save(closest30mto5m_noRiver,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest30mto5m_noRiver.RData")

################################################################################

closest50mto5m_noRiver<- matrix(NA, ncol= nRuns, nrow= SS)

for(i in 1:nRuns){
  run<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/Extent/Run_",i,".asc"))
  vals<- extract(run,coordsIwant.5m)
  closest50mto5m_noRiver[,i]<- vals$max
}

save(closest50mto5m_noRiver,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest50mto5m_noRiver.RData")
