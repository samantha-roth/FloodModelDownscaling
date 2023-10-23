#explore whether elevation, slope, TRI, flow accumulation, and cost at 10m
#and differences between these variables at 50m and 50m impact residuals
#from either model WSE10m ~ poly(WSE50mat10m,2)
#OR if they impact downscaled.z.vec- WSE10m or WSE50mat10m-WSE10m 

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

n10<- 104

library(terra)

#load DEMs
dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
dem50<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_50m_new.asc")

#load calibrated runs
Run10m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/Norristown/nCh/simplecalQs/Extent/Run_1.asc")
Run50m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/simplecalQs/Extent/Run_1.asc")
ncell50m<- ncell(Run50m)

coords.10m<- xyFromCell(Run10m,1:ncell(Run10m))
flood.10m<- extract(Run10m,coords.10m)

flood.50mat10m<- extract(Run50m,coords.10m)
goodinds<- which(!is.na(flood.50mat10m))

elev.10m<- extract(dem10,coords.10m)
elev.50mat10m<- extract(dem50,coords.10m)

#load cost of getting to each point outside the river
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_10m.RData")
costs10<- costs_RiverCentertoOutsidePts
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/coords_OutsideRiverbyX_10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_50mat10m.RData")
costs50at10<- costs_RiverCentertoOutsidePts

#compute differences in cost
costdiffs<- rep(NA,nrow(coords.10m))
allInds<- NA
for(i in 1:n10){
  c10<- costs10[[i]]
  c50<- costs50at10[[i]]
  coords<- coords_OutsideRiverbyX[[i]]
  inds<- coords[,1]
  allInds<- c(allInds,inds)
  cdiff<- c50-c10
  costdiffs[inds]<- cdiff
}
allInds<- allInds[-1]

length(which(!is.na(costdiffs)))==length(allInds) 

keepInds<- intersect(which(!is.na(costdiffs)),goodinds)

#compute WSE at each resolution on 10m grid

WSE50mat10m<- flood.50mat10m$Run50m_simplecal[keepInds] + elev.50mat10m$norristown_50m_new[keepInds]
WSE10m<- flood.10m$Run10m_simplecal[keepInds]+elev.10m$norristown_10m_new[keepInds]

#compute difference in calibrated WSEs at 50m on 10m grid vs 10m
wsediff<- WSE50mat10m - WSE10m

WSH50mat10m<- flood.50mat10m$Run50m_simplecal[keepInds] + 
  elev.50mat10m$norristown_50m_new[keepInds] - elev.10m$norristown_10m_new[keepInds]
WSH10m<- flood.10m$Run10m_simplecal[keepInds] 

#compute difference in calibrated WSEs at 50m on 10m grid vs 10m
wshdiff<- WSH50mat10m - WSH10m

######################################################################
#Next we use bilinear interpolation to get the WSH at 50m on 10m grid

coords.50m<- xyFromCell(Run50m,1:ncell(Run50m))

min.x50<- min(coords.50m[,1])
min.y50<- min(coords.50m[,2])
max.x50<- max(coords.50m[,1])
max.y50<- max(coords.50m[,2])

xIndsIwant1<- which(coords.10m[,1]>=min.x50)
xIndsIwant2<- which(coords.10m[,1]<=max.x50)
yIndsIwant1<- which(coords.10m[,2]>=min.y50)
yIndsIwant2<- which(coords.10m[,2]<=max.y50)

yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)

coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
coordsIwant.10m<- coords.10m[coordsIwantInds,]

elevIwant.10m<- extract(dem10,coordsIwant.10m)
floodIwant.10m<- extract(Run10m,coordsIwant.10m)
WSHIwant.10m = floodIwant.10m

y10m<- unique(coordsIwant.10m[,2])
x10m<- unique(coordsIwant.10m[,1])
y50m<- unique(coords.50m[,2])
x50m<- unique(coords.50m[,1])

nx10m<- length(x10m)
ny10m<- length(y10m)
nx50m<- length(x50m)
ny50m<- length(y50m)

library(akima)

elev.50m<- extract(dem50,coords.50m)
wsh.50m<- extract(Run50m,coords.50m)
vals.50m<- wsh.50m + elev.50m 
z1= matrix(vals.50m$Run50m_simplecal, nrow= ny50m, ncol= nx50m, byrow= TRUE)

z= matrix(NA,nrow= ny50m, ncol= nx50m)
for(j in 1:nx50m){ z[,j]<- rev(z1[,j]) }

test<- bilinear(x= rev(y50m), y= x50m, z= z, 
                x0= rev(coordsIwant.10m[,2]), y0= coordsIwant.10m[,1])

z2= matrix(test$z, nrow= ny10m, ncol= nx10m, byrow= TRUE)

#translate back to original coordinates
downscaled.z<- matrix(NA,nrow= ny10m, ncol= nx10m)
for(k in 1:nx10m){downscaled.z[,k]<- rev(z2[,k])}
downscaled.z.vec<- c(t(downscaled.z))
downscaled.z.vec<- downscaled.z.vec- elevIwant.10m$norristown_10m_new

indstoCompare<- intersect(keepInds,coordsIwantInds)

downscale50m<- rep(NA,length(flood.10m$Run10m_simplecal))
downscale50m[coordsIwantInds]<- downscaled.z.vec

downscale50m<- downscale50m[indstoCompare]

WSHtocompare<- flood.10m$Run10m_simplecal[indstoCompare] 

plot(downscale50m,WSHtocompare,main="Calibrated, downscaled 50m WSH VS calibrated 10m WSH",
     xlab="Calibrated, downscaled 50m WSH", ylab= "Calibrated 10m WSH")

#how well can we predict calibrated WSH at 10m using calibrated WSH at 50m on 10m grid?
fit<- lm(WSHtocompare~ downscale50m)
summary(fit)
resids<- fit$residuals
#incredibly damn well, better than when not bilinearly interpolated!

plot(WSHtocompare,resids)
plot(fit$fitted.values,WSHtocompare)

summary(resids)
#50% of calibrated 50m flood heights are 
#within (-.18356, .15672) m of the calibrate 10m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(coords.10m[indstoCompare,],resids))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals in space- WSHtocompare~ downscale50m") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")
dev.off()

#highish residuals near bottom of river
#lowish residuals near middle and for below river, near upper
#most extreme values above river near top of region

################################################################################

#load DEM info at 50m and 10m resolution on 10m grid
#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat10m.RData")
#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

#DEMinfo_10m<- DEMinfo_10m[indstoCompare,]
#DEMinfo_50m<- DEMinfo_50m[indstoCompare,]

#DEMinfo_diff<- DEMinfo_50m-DEMinfo_10m

#plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

#plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

#plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

#DEMinfo_diff[,6]<- DEMinfo_50m[,6]/ncell50m- DEMinfo_10m[,6]/nrow(coords.10m)

#plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

#DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_50m[,6]/ncell50m/DEMinfo_10m[,6]/nrow(coords.10m))

#plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

##IF ANYTHING: more negative elevation difference --> positive residual,
## positive elevation difference --> negative residual

#bin_res<- ifelse(resids>0,1,0)
#bin_elev.diff<- ifelse(DEMinfo_diff[,3]>0,1,0)

#fit.sign.res<- glm(bin_res~bin_elev.diff[indstoCompare],family="binomial")
#summary(fit.sign.res)
##not even close to significant

################################################################################

#now assume we can only predict WSH at 10m with downscale 50m


length(which(downscale50m<0))

downscale50m[which(downscale50m<0)]<- 0

resids<- WSHtocompare- downscale50m

summary(resids)

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(indstoCompare,coords.10m[indstoCompare,],resids))
colnames(residsbyloc.df)<-c("inds","x","y","value")

save(residsbyloc.df,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH10-downscale50.RData")

library(ggplot2)

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnscl,cal50-cal10resids_spatial.jpeg", width = 800, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=1.6) + ggtitle("Residuals in space- WSH10m- downscale50m") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()

#highish residuals near bottom of river
#lowish residuals near middle and for below river, near upper
#most extreme values above river near top of region

#make a dataframe to plot residuals in space
elevbyloc.df<- as.data.frame(cbind(coords.10m[indstoCompare,],
                                   elev.10m$norristown_10m_new[indstoCompare]))
colnames(elevbyloc.df)<-c("x","y","value")

ggplot(elevbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("elev in space") +theme_bw()+
  scale_color_gradient(low="blue", high="red")



################################################################################

#load DEM info at 50m and 10m resolution on 10m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_10m<- DEMinfo_10m[indstoCompare,]
DEMinfo_50m<- DEMinfo_50m[indstoCompare,]

DEMinfo_diff<- DEMinfo_50m-DEMinfo_10m

plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

DEMinfo_diff[,6]<- DEMinfo_50m[,6]/ncell50m- DEMinfo_10m[,6]/nrow(coords.10m)

plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_50m[,6]/ncell50m/DEMinfo_10m[,6]/nrow(coords.10m))

plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

#IF ANYTHING: more negative elevation difference --> positive residual,
# positive elevation difference --> negative residual

bin_res<- ifelse(resids>0,1,0)
bin_elev.diff<- ifelse(DEMinfo_diff[,3]>0,1,0)

fit.sign.res<- glm(bin_res~bin_elev.diff[indstoCompare],family="binomial")
summary(fit.sign.res)
#not even close to significant

################################################################################

#use cost from center of river to point

plot(costdiffs[indstoCompare],resids,xlab="Cost Difference",ylab="resid")

################################################################################
#now just compare to DEM info at 10m

plot(DEMinfo_10m[,3],resids,xlab="elev",ylab="resid")

plot(DEMinfo_10m[,4],resids,xlab="slope",ylab="resid")

plot(DEMinfo_10m[,5],resids,xlab="TRI",ylab="resid")

plot(DEMinfo_10m[,6],resids,xlab="Flow Acc",ylab="resid")
