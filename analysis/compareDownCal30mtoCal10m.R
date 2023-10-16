#see how well calibrated 30m flood heights predict calibrated 10m flood heights

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

library(terra)

dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
dem30<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_30m_new.asc")

Run10m<- rast("C:/Users/svr5482/Downloads/Run_10m.asc")
Run30m<- rast("C:/Users/svr5482/Downloads/Run_30m.asc")

coords.10m<- xyFromCell(Run10m,1:ncell(Run10m))
flood.10m<- extract(Run10m,coords.10m)

flood.30mat10m<- extract(Run30m,coords.10m)
goodinds<- which(!is.na(flood.30mat10m))

elev.10m<- extract(dem10,coords.10m)
elev.30mat10m<- extract(dem30,coords.10m)

WSE30mat10m<- flood.30mat10m$Run_30m[goodinds] + elev.30mat10m$norristown_30m_new[goodinds]
WSE10m<- flood.10m$Run_10m[goodinds]+elev.10m$norristown_10m_new[goodinds]

elevdiff<- elev.30mat10m$norristown_30m_new[goodinds]- elev.10m$norristown_10m_new[goodinds]
wsediff<- WSE30mat10m - WSE10m

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSE30mVS10m.jpeg",width = 500, height = 500)
plot(WSE30mat10m,WSE10m,main="Calibrated WSE at 30m VS 10m",
     xlab="Calibrated WSE at 30m", ylab= "Calibrated WSE at 10m")
dev.off()

plot(elevdiff,wsediff,main="30m - 10m: Difference in elev VS difference in calibrated WSE")


fit<- lm(WSE10m~ poly(WSE30mat10m,2))
summary(fit)

plot(WSE10m,fit$residuals)

plot(elevdiff,fit$residuals)

fitres<- lm(fit$residuals~elevdiff)
summary(fitres)

residRun<- Run10m
values(residRun)<- fit$residuals

vals<- extract(residRun,coords.10m)

#all equal

plot(residRun)
summary(fit$residuals)
#75% of calibrated 30m flood heights are 
#within .46 m of the calibrate 10m flood height


######################################################################

coords.30m<- xyFromCell(Run30m,1:ncell(Run30m))

min.x30<- min(coords.30m[,1])
min.y30<- min(coords.30m[,2])
max.x30<- max(coords.30m[,1])
max.y30<- max(coords.30m[,2])

xIndsIwant1<- which(coords.10m[,1]>=min.x30)
xIndsIwant2<- which(coords.10m[,1]<=max.x30)
yIndsIwant1<- which(coords.10m[,2]>=min.y30)
yIndsIwant2<- which(coords.10m[,2]<=max.y30)

yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)

coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
coordsIwant.10m<- coords.10m[coordsIwantInds,]

elevIwant.10m<- extract(dem10,coordsIwant.10m)
floodIwant.10m<- extract(Run10m,coordsIwant.10m)
WSEIwant.10m = elevIwant.10m+ floodIwant.10m

y10m<- unique(coordsIwant.10m[,2])
x10m<- unique(coordsIwant.10m[,1])
y30m<- unique(coords.30m[,2])
x30m<- unique(coords.30m[,1])

nx10m<- length(x10m)
ny10m<- length(y10m)
nx30m<- length(x30m)
ny30m<- length(y30m)

library(akima)

vals.30m<- extract(Run30m,coords.30m)
z1= matrix(vals.30m$Run_30m, nrow= ny30m, ncol= nx30m, byrow= TRUE)

z= matrix(NA,nrow= ny30m, ncol= nx30m)
for(j in 1:nx30m){ z[,j]<- rev(z1[,j]) }

test<- bilinear(x= rev(y30m), y= x30m, z= z, 
                x0= rev(coordsIwant.10m[,2]), y0= coordsIwant.10m[,1])

z2= matrix(test$z, nrow= ny10m, ncol= nx10m, byrow= TRUE)

#translate back to original coordinates
downscaled.z<- matrix(NA,nrow= ny10m, ncol= nx10m)
for(k in 1:nx10m){downscaled.z[,k]<- rev(z2[,k])}
downscaled.z.vec<- c(t(downscaled.z))

plot(downscaled.z.vec,WSEIwant.10m$norristown_10m_new)

plot(downscaled.z.vec,WSE10m,main="Calibrated, downscaled 30m WSE VS calibrated 10m WSE",
     xlab="Calibrated, downscaled 30m WSE", ylab= "Calibrated 10m WSE")
