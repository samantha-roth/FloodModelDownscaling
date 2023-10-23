#explore whether elevation, slope, TRI, flow accumulation, and cost at 30m
#and differences between these variables at 50m and 50m impact residuals
#from either model WSE30m ~ poly(WSE50mat30m,2)
#OR if they impact downscaled.z.vec- WSE30m or WSE50mat30m-WSE30m 

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

n30<- 35

library(terra)

#load DEMs
dem30<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_30m_new.asc")
dem50<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_50m_new.asc")

#load calibrated runs
Run30m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs30m/Norristown/nCh/simplecalQs/Extent/Run_1.asc")
Run50m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/simplecalQs/Extent/Run_1.asc")
ncell50m<- ncell(Run50m)

coords.30m<- xyFromCell(Run30m,1:ncell(Run30m))
flood.30m<- extract(Run30m,coords.30m)

flood.50mat30m<- extract(Run50m,coords.30m)
goodinds<- which(!is.na(flood.50mat30m))

elev.30m<- extract(dem30,coords.30m)
elev.50mat30m<- extract(dem50,coords.30m)

#load cost of getting to each point outside the river
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_30m.RData")
costs30<- costs_RiverCentertoOutsidePts
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/coords_OutsideRiverbyX_30m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/costs_RiverCentertoOutsidePts_50mat30m.RData")
costs50at30<- costs_RiverCentertoOutsidePts

#compute differences in cost
costdiffs<- rep(NA,nrow(coords.30m))
allInds<- NA
for(i in 1:n30){
  c30<- costs30[[i]]
  c50<- costs50at30[[i]]
  coords<- coords_OutsideRiverbyX[[i]]
  inds<- coords[,1]
  allInds<- c(allInds,inds)
  cdiff<- c50-c30
  costdiffs[inds]<- cdiff
}
allInds<- allInds[-1]

length(which(!is.na(costdiffs)))==length(allInds) 

keepInds<- intersect(which(!is.na(costdiffs)),goodinds)

#compute WSE at each resolution on 30m grid

WSE50mat30m<- flood.50mat30m$Run50m_simplecal[keepInds] + elev.50mat30m$norristown_50m_new[keepInds]
WSE30m<- flood.30m$Run30m_simplecal[keepInds]+elev.30m$norristown_30m_new[keepInds]

#compute difference in calibrated WSEs at 50m on 30m grid vs 30m
wsediff<- WSE50mat30m - WSE30m

WSH50mat30m<- flood.50mat30m$Run50m_simplecal[keepInds] + 
  elev.50mat30m$norristown_50m_new[keepInds] - elev.30m$norristown_30m_new[keepInds]
WSH30m<- flood.30m$Run30m_simplecal[keepInds] 

#compute difference in calibrated WSEs at 50m on 30m grid vs 30m
wshdiff<- WSH50mat30m - WSH30m

######################################################################
#Next we use bilinear interpolation to get the WSH at 50m on 30m grid

coords.50m<- xyFromCell(Run50m,1:ncell(Run50m))

min.x50<- min(coords.50m[,1])
min.y50<- min(coords.50m[,2])
max.x50<- max(coords.50m[,1])
max.y50<- max(coords.50m[,2])

xIndsIwant1<- which(coords.30m[,1]>=min.x50)
xIndsIwant2<- which(coords.30m[,1]<=max.x50)
yIndsIwant1<- which(coords.30m[,2]>=min.y50)
yIndsIwant2<- which(coords.30m[,2]<=max.y50)

yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)

coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
coordsIwant.30m<- coords.30m[coordsIwantInds,]

elevIwant.30m<- extract(dem30,coordsIwant.30m)
floodIwant.30m<- extract(Run30m,coordsIwant.30m)
WSHIwant.30m = floodIwant.30m

y30m<- unique(coordsIwant.30m[,2])
x30m<- unique(coordsIwant.30m[,1])
y50m<- unique(coords.50m[,2])
x50m<- unique(coords.50m[,1])

nx30m<- length(x30m)
ny30m<- length(y30m)
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
                x0= rev(coordsIwant.30m[,2]), y0= coordsIwant.30m[,1])

z2= matrix(test$z, nrow= ny30m, ncol= nx30m, byrow= TRUE)

#translate back to original coordinates
downscaled.z<- matrix(NA,nrow= ny30m, ncol= nx30m)
for(k in 1:nx30m){downscaled.z[,k]<- rev(z2[,k])}
downscaled.z.vec<- c(t(downscaled.z))
downscaled.z.vec<- downscaled.z.vec- elevIwant.30m$norristown_30m_new

indstoCompare<- intersect(keepInds,coordsIwantInds)

downscale50m<- rep(NA,length(flood.30m$Run30m_simplecal))
downscale50m[coordsIwantInds]<- downscaled.z.vec

downscale50m<- downscale50m[indstoCompare]

WSHtocompare<- flood.30m$Run30m_simplecal[indstoCompare] 

plot(downscale50m,WSHtocompare,main="Calibrated, downscaled 50m WSH VS calibrated 30m WSH",
     xlab="Calibrated, downscaled 50m WSH", ylab= "Calibrated 30m WSH")

#how well can we predict calibrated WSH at 30m using calibrated WSH at 50m on 30m grid?
fit<- lm(WSHtocompare~ downscale50m)
summary(fit)
resids<- fit$residuals
#incredibly damn well, better than when not bilinearly interpolated!

plot(WSHtocompare,resids)
plot(fit$fitted.values,WSHtocompare)

summary(resids)
#50% of calibrated 50m flood heights are 
#within (-.18356, .15672) m of the calibrate 30m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(coords.30m[indstoCompare,],resids))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals in space- WSHtocompare~ downscale50m") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red") +
  theme(plot.title = element_text(size=24),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 18))
dev.off()

#highish residuals near bottom of river
#lowish residuals near middle and for below river, near upper
#most extreme values above river near top of region

################################################################################

##load DEM info at 50m and 30m resolution on 30m grid
#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat30m.RData")
#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30m.RData")

#DEMinfo_30m<- DEMinfo_30m[indstoCompare,]
#DEMinfo_50m<- DEMinfo_50m[indstoCompare,]

#DEMinfo_diff<- DEMinfo_50m-DEMinfo_30m

#plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

#plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

#plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

#DEMinfo_diff[,6]<- DEMinfo_50m[,6]/ncell50m- DEMinfo_30m[,6]/nrow(coords.30m)

#plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

#DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_50m[,6]/ncell50m/DEMinfo_30m[,6]/nrow(coords.30m))

#plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

##IF ANYTHING: more negative elevation difference --> positive residual,
## positive elevation difference --> negative residual

#bin_res<- ifelse(resids>0,1,0)
#bin_elev.diff<- ifelse(DEMinfo_diff[,3]>0,1,0)

#fit.sign.res<- glm(bin_res~bin_elev.diff[indstoCompare],family="binomial")
#summary(fit.sign.res)
##not even close to significant

################################################################################

#now assume we can only predict WSH at 30m with downscale 50m

resids<- WSHtocompare- downscale50m

summary(resids)
#50% of calibrated 50m flood heights are 
#within (-.4667, .4113) m of the calibrate 30m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(indstoCompare,coords.30m[indstoCompare,],resids))
colnames(residsbyloc.df)<-c("inds","x","y","value")

save(residsbyloc.df,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH30-downscale50.RData")

library(ggplot2)

#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnsclno50to10.RData")

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/dnscl,cal50-cal30resids_spatial.jpeg", width = 800, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=5) + ggtitle("Residuals in space- WSH30m- downscale50m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
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
elevbyloc.df<- as.data.frame(cbind(coords.30m[indstoCompare,],
                                   elev.30m$norristown_30m_new[indstoCompare]))
colnames(elevbyloc.df)<-c("x","y","value")

ggplot(elevbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("elev in space") +theme_bw()+
  scale_color_gradient(low="blue", high="red")



################################################################################

#load DEM info at 50m and 30m resolution on 30m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info50mat30m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30m.RData")

DEMinfo_30m<- DEMinfo_30m[indstoCompare,]
DEMinfo_50m<- DEMinfo_50m[indstoCompare,]

DEMinfo_diff<- DEMinfo_50m-DEMinfo_30m

plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

DEMinfo_diff[,6]<- DEMinfo_50m[,6]/ncell50m- DEMinfo_30m[,6]/nrow(coords.30m)

plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_50m[,6]/ncell50m/DEMinfo_30m[,6]/nrow(coords.30m))

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
#now just compare to DEM info at 30m

plot(DEMinfo_30m[,3],resids,xlab="elev",ylab="resid")

plot(DEMinfo_30m[,4],resids,xlab="slope",ylab="resid")

plot(DEMinfo_30m[,5],resids,xlab="TRI",ylab="resid")

plot(DEMinfo_30m[,6],resids,xlab="Flow Acc",ylab="resid")
