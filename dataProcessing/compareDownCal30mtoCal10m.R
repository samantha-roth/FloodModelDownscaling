#explore whether elevation, slope, TRI, flow accumulation, and cost at 10m
#and differences between these variables at 30m and 10m impact residuals
#from either model WSE10m ~ poly(WSE30mat10m,2)
#OR if they impact downscaled.z.vec- WSE10m or WSE30mat10m-WSE10m 

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

n10<- 104

library(terra)

#load DEMs
dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
dem30<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_30m_new.asc")

#load calibrated runs
Run10m<- rast("C:/Users/svr5482/Downloads/Run_10m.asc")
Run30m<- rast("C:/Users/svr5482/Downloads/Run_30m.asc")
ncell30m<- ncell(Run30m)

coords.10m<- xyFromCell(Run10m,1:ncell(Run10m))
flood.10m<- extract(Run10m,coords.10m)

flood.30mat10m<- extract(Run30m,coords.10m)
goodinds<- which(!is.na(flood.30mat10m))

elev.10m<- extract(dem10,coords.10m)
elev.30mat10m<- extract(dem30,coords.10m)

#load cost of getting to each point outside the river
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_10m.RData")
costs10<- costs_RiverCentertoOutsidePts
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/coords_OutsideRiverbyX_10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_30mat10m.RData")
costs30at10<- costs_RiverCentertoOutsidePts

#compute differences in cost
costdiffs<- rep(NA,nrow(coords.10m))
allInds<- NA
for(i in 1:n10){
  c10<- costs10[[i]]
  c30<- costs30at10[[i]]
  coords<- coords_OutsideRiverbyX[[i]]
  inds<- coords[,1]
  allInds<- c(allInds,inds)
  cdiff<- c30-c10
  costdiffs[inds]<- cdiff
}
allInds<- allInds[-1]

length(which(!is.na(costdiffs)))==length(allInds) 

keepInds<- intersect(which(!is.na(costdiffs)),goodinds)

#compute WSE at each resolution on 10m grid

WSE30mat10m<- flood.30mat10m$Run_30m[keepInds] + elev.30mat10m$norristown_30m_new[keepInds]
WSE10m<- flood.10m$Run_10m[keepInds]+elev.10m$norristown_10m_new[keepInds]

#compute difference in calibrated WSEs at 30m on 10m grid vs 10m
wsediff<- WSE30mat10m - WSE10m

WSH30mat10m<- flood.30mat10m$Run_30m[keepInds] + 
  elev.30mat10m$norristown_30m_new[keepInds] - elev.10m$norristown_10m_new[keepInds]
WSH10m<- flood.10m$Run_10m[keepInds] 

#compute difference in calibrated WSEs at 30m on 10m grid vs 10m
wshdiff<- WSH30mat10m - WSH10m

######################################################################
#Next we use bilinear interpolation to get the WSH at 30m on 10m grid

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
WSHIwant.10m = floodIwant.10m

y10m<- unique(coordsIwant.10m[,2])
x10m<- unique(coordsIwant.10m[,1])
y30m<- unique(coords.30m[,2])
x30m<- unique(coords.30m[,1])

nx10m<- length(x10m)
ny10m<- length(y10m)
nx30m<- length(x30m)
ny30m<- length(y30m)

library(akima)

elev.30m<- extract(dem30,coords.30m)
wsh.30m<- extract(Run30m,coords.30m)
vals.30m<- wsh.30m + elev.30m 
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
downscaled.z.vec<- downscaled.z.vec- elevIwant.10m$norristown_10m_new

coordstoCompare<- intersect(keepInds,coordsIwantInds)

downscale30m<- rep(NA,length(flood.10m$Run_10m))
downscale30m[coordsIwantInds]<- downscaled.z.vec

downscale30m<- downscale30m[coordstoCompare]

WSHtocompare<- flood.10m$Run_10m[coordstoCompare] 

plot(downscale30m,WSHtocompare,main="Calibrated, downscaled 30m WSH VS calibrated 10m WSH",
     xlab="Calibrated, downscaled 30m WSH", ylab= "Calibrated 10m WSH")

#how well can we predict calibrated WSH at 10m using calibrated WSH at 30m on 10m grid?
fit<- lm(WSHtocompare~ downscale30m)
summary(fit)
resids<- fit$residuals
#incredibly damn well, better than when not bilinearly interpolated!

plot(WSHtocompare,resids)
plot(fit$fitted.values,WSHtocompare)

summary(resids)
#50% of calibrated 30m flood heights are 
#within (-.4667, .4113) m of the calibrate 10m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(coords.10m[coordstoCompare,],resids))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals in space- WSHtocompare~ downscale30m") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")
dev.off()

#highish residuals near bottom of river
#lowish residuals near middle and for below river, near upper
#most extreme values above river near top of region

################################################################################

#load DEM info at 30m and 10m resolution on 10m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30mat10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_10m<- DEMinfo_10m[coordstoCompare,]
DEMinfo_30m<- DEMinfo_30m[coordstoCompare,]

DEMinfo_diff<- DEMinfo_30m-DEMinfo_10m

plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

DEMinfo_diff[,6]<- DEMinfo_30m[,6]/ncell30m- DEMinfo_10m[,6]/nrow(coords.10m)

plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_30m[,6]/ncell30m/DEMinfo_10m[,6]/nrow(coords.10m))

plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

#IF ANYTHING: more negative elevation difference --> positive residual,
# positive elevation difference --> negative residual

bin_res<- ifelse(resids>0,1,0)
bin_elev.diff<- ifelse(DEMinfo_diff[,3]>0,1,0)

fit.sign.res<- glm(bin_res~bin_elev.diff[coordstoCompare],family="binomial")
summary(fit.sign.res)
#not even close to significant

################################################################################

#now assume we can only predict WSH at 10m with downscale 30m

resids<- WSHtocompare- downscale30m

summary(resids)
#50% of calibrated 30m flood heights are 
#within (-.4667, .4113) m of the calibrate 10m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(coords.10m[coordstoCompare,],resids))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals in space- WSHtocompare~ downscale30m") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")
dev.off()

#highish residuals near bottom of river
#lowish residuals near middle and for below river, near upper
#most extreme values above river near top of region

################################################################################

#load DEM info at 30m and 10m resolution on 10m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30mat10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_10m<- DEMinfo_10m[coordstoCompare,]
DEMinfo_30m<- DEMinfo_30m[coordstoCompare,]

DEMinfo_diff<- DEMinfo_30m-DEMinfo_10m

plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

DEMinfo_diff[,6]<- DEMinfo_30m[,6]/ncell30m- DEMinfo_10m[,6]/nrow(coords.10m)

plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_30m[,6]/ncell30m/DEMinfo_10m[,6]/nrow(coords.10m))

plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

#IF ANYTHING: more negative elevation difference --> positive residual,
# positive elevation difference --> negative residual

bin_res<- ifelse(resids>0,1,0)
bin_elev.diff<- ifelse(DEMinfo_diff[,3]>0,1,0)

fit.sign.res<- glm(bin_res~bin_elev.diff[coordstoCompare],family="binomial")
summary(fit.sign.res)
#not even close to significant
