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

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSE30mVS10m.jpeg",width = 500, height = 500)
plot(WSE30mat10m,WSE10m,main="Calibrated WSE at 30m VS 10m",
     xlab="Calibrated WSE at 30m", ylab= "Calibrated WSE at 10m")
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSH30mVS10m.jpeg",width = 500, height = 500)
plot(WSH30mat10m,WSH10m,main="Calibrated WSH at 30m VS 10m",
     xlab="Calibrated WSH at 30m", ylab= "Calibrated WSH at 10m")
dev.off()

#how well can we predict calibrated WSH at 10m using calibrated WSH at 30m on 10m grid?
fit<- lm(WSH10m~ WSH30mat10m)
summary(fit)
resids<- fit$residuals
#pretty damn well!

plot(WSH10m,resids)
plot(fit$fitted.values,WSH10m)

summary(resids)
#50% of calibrated 30m flood heights are 
#within about .5 m of the calibrate 10m flood height

#make a dataframe to plot residuals in space
residsbyloc.df<- as.data.frame(cbind(coords.10m[keepInds,],resids))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSHmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals by location") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")
dev.off()
#the residuals appear to be larger in absolute value above than below the river
#below the river the residuals appear to be more negative close to the upper portion of the river

flood.30mat10m_tokeep<- flood.30mat10m$Run_30m[keepInds]
zerothirtyInds<- which(flood.30mat10m_tokeep==0)

residRun<- Run10m
values(residRun)[keepInds]<- fit$residuals
values(residRun)[-keepInds]<- NA
vals<- extract(residRun,coords.10m)
plot(residRun)

################################################################################

#load DEM info at 30m and 10m resolution on 10m grid
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info30mat10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/DEM_info10m.RData")

DEMinfo_10m<- DEMinfo_10m[keepInds,]
DEMinfo_30m<- DEMinfo_30m[keepInds,]

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
bin_elev.diff<- ifelse(elev.diff>0,1,0)

fit.sign.res<- glm(bin_res~bin_elev.diff[keepInds],family="binomial")
summary(fit.sign.res)

################################################################################

#elev.10m<- extract(dem10,coords.10m)
#elev.30mat10m<- extract(dem30,coords.10m)

#compute difference in elevation at 30m vs 10m
#elevdiff<- elev.30mat10m$norristown_30m_new[keepInds]- elev.10m$norristown_10m_new[keepInds]

#plot(elevdiff,wsediff,main="30m - 10m: Difference in elev VS difference in calibrated WSE")

#plot(elevdiff,fit$residuals)

#fitres<- lm(fit$residuals~elevdiff)
#summary(fitres)
