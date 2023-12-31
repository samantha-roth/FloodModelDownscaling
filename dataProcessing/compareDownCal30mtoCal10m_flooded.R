
rm(list=ls())

#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/floodInds30m.RData")
#floodInds30m<- floodInds
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/floodInds30mat10m.RData")
floodInds30mat10m<- floodInds; rm(floodInds)

#explore whether elevation, slope, TRI, flow accumulation, and cost at 10m
#and differences between these variables at 30m and 10m impact residuals
#from either model WSE10m ~ poly(WSE30mat10m,2)
#OR if they impact downscaled.z.vec- WSE10m or WSE30mat10m-WSE10m 

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots/30to10")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots/30to10")
}

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds")
}

n10<- 104

library(terra)

#load DEMs
dem10<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_10m_new.asc")
dem30<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_30m_new.asc")

#load calibrated runs
Run10m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/Norristown/nCh/simplecal/Extent/Run_1.asc")
Run30m<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs30m/Norristown/nCh/simplecal/Extent/Run_1.asc")
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

oldkeepInds<- intersect(which(!is.na(costdiffs)),goodinds)
keepInds<- intersect(oldkeepInds,floodInds30mat10m)

#compute WSE at each resolution on 10m grid

WSE30mat10m<- flood.30mat10m$Run_1[keepInds] + elev.30mat10m$norristown_30m_new[keepInds]
WSE10m<- flood.10m$Run_1[keepInds]+elev.10m$norristown_10m_new[keepInds]

#compute difference in calibrated WSEs at 30m on 10m grid vs 10m
wsediff<- WSE30mat10m - WSE10m

WSH30mat10m<- flood.30mat10m$Run_1[keepInds] + 
  elev.30mat10m$norristown_30m_new[keepInds] - elev.10m$norristown_10m_new[keepInds]
WSH10m<- flood.10m$Run_1[keepInds] 

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
z1= matrix(vals.30m$Run_1, nrow= ny30m, ncol= nx30m, byrow= TRUE)

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

indstoCompare<- intersect(keepInds,coordsIwantInds)

downscale30m<- rep(NA,length(flood.10m$Run_1))
downscale30m[coordsIwantInds]<- downscaled.z.vec

downscale30m<- downscale30m[indstoCompare]

WSHtocompare<- flood.10m$Run_1[indstoCompare] 

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30VS10_floodInds.jpeg",
     width = 800, height = 700)
plot(downscale30m,WSHtocompare,main="Calibrated, downscaled 30m WSH VS calibrated 10m WSH",
     xlab="Calibrated, downscaled 30m WSH", ylab= "Calibrated 10m WSH")
dev.off()

#predict WSH at 10m with downscale 30m

#length(which(downscale30m<0))

#downscale30m[which(downscale30m<0)]<- 0

resids<- WSHtocompare- downscale30m

summary(resids)
#30% of calibrated 30m flood heights are 
#within (-.4667, .4113) m of the calibrate 10m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(indstoCompare,coords.10m[indstoCompare,],resids))
colnames(residsbyloc.df)<-c("inds","x","y","value")

save(residsbyloc.df,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH10-downscale30_floodInds.RData")

library(ggplot2)

#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")
#load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnsclno30to10.RData")

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl,cal30-cal10resids_spatialfloodInds.jpeg", width = 800, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=1.8) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()


#make a dataframe to plot residuals in space
elevbyloc.df<- as.data.frame(cbind(coords.10m[indstoCompare,],
                                   elev.10m$norristown_10m_new[indstoCompare]))
colnames(elevbyloc.df)<-c("x","y","value")

ggplot(elevbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("elev in space") +theme_bw()+
  scale_color_gradient(low="blue", high="red")



################################################################################

#load DEM info at 30m and 10m resolution on 10m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30mat10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_10m<- DEMinfo_10m[indstoCompare,]
DEMinfo_30m<- DEMinfo_30m[indstoCompare,]

DEMinfo_diff<- DEMinfo_30m-DEMinfo_10m

plot(DEMinfo_diff[,3],resids,xlab="elev diff",ylab="resid")

plot(DEMinfo_diff[,4],resids,xlab="slope diff",ylab="resid")

plot(DEMinfo_diff[,5],resids,xlab="TRI diff",ylab="resid")

DEMinfo_diff[,6]<- DEMinfo_30m[,6]/ncell30m- DEMinfo_10m[,6]/nrow(coords.10m)

plot(DEMinfo_diff[,6],resids,xlab="Flow Acc Pct diff",ylab="resid")

DEMinfo_diff<- cbind(DEMinfo_diff,DEMinfo_30m[,6]/ncell30m/DEMinfo_10m[,6]/nrow(coords.10m))

plot(DEMinfo_diff[,7],resids,xlab="Flow Acc ratio",ylab="resid")

################################################################################

#use cost from center of river to point

plot(costdiffs[indstoCompare],resids,xlab="Cost Difference",ylab="resid")

################################################################################
#now just compare to DEM info at 10m

plot(DEMinfo_10m[,3],resids,xlab="elev",ylab="resid")

#jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/slope10mVSdnscl,cal30-cal10resids.jpeg", width = 800, height = 700)
plot(DEMinfo_10m[,4],resids,xlab="slope",ylab="resid",main="downscale from 30m to 10m resid vs 10m slope")
#dev.off()

plot(DEMinfo_10m[,5],resids,xlab="TRI",ylab="resid")

plot(DEMinfo_10m[,6],resids,xlab="Flow Acc",ylab="resid")
###############################################################
#two way interactions

residsbylocandDEM.df<- as.data.frame(cbind(indstoCompare,resids,DEMinfo_10m))
colnames(residsbylocandDEM.df)[1:2]<-c("inds","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30-10res_elevVSslope_floodInds.jpeg", width = 800, height = 700)
ggplot(residsbylocandDEM.df, aes(elev, slope, color=value))+
  geom_point(size=1.5) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30-10res_elevVSTRI_floodInds.jpeg", width = 800, height = 700)
ggplot(residsbylocandDEM.df, aes(elev, TRI, color=value))+
  geom_point(size=1.5) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30-10res_elevVSflowacc_floodInds.jpeg", width = 800, height = 700)
ggplot(residsbylocandDEM.df, aes(elev, flowacc, color=value))+
  geom_point(size=1.5) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30-10res_slopeVSTRI_floodInds.jpeg", width = 800, height = 700)
ggplot(residsbylocandDEM.df, aes(slope, TRI, color=value))+
  geom_point(size=1.5) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/30to10/floodInds/dnscl30-10res_slopeVSflowacc_floodInds.jpeg", width = 800, height = 700)
ggplot(residsbylocandDEM.df, aes(slope, flowacc, color=value))+
  geom_point(size=1.5) + ggtitle("Residuals in space- WSH10m- downscale30m") +theme_bw()+
  #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
  theme(plot.title = element_text(size=24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18))
dev.off()
