rm(list=ls())

library(terra)

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/not.river.inds5m.RData")

runST<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs5m/Norristown/nCh/ST/Extent/Run_1.asc")

crs(runST)<- "+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.5m<- xyFromCell(runST,1:ncell(runST))
floodST<- extract(runST,coords.5m)
ncell.5m<- length(floodST$max)

floodInds<- which(floodST>0)
noFloodInds<- which(floodST==0)

indsErr<-intersect(floodInds,not.river.inds5m)
obsErr<- rnorm(length(indsErr),0,.03)

obsWE_5m<- floodST$max
obsWE_5m[indsErr]<- floodST$max[indsErr]+obsErr
obsWE_OR_5m<- obsWE_5m[not.river.inds5m]

save(obsWE_OR_5m,indsErr,obsErr,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/obsWE_OR_5m.RData")
