#compare costs of water transport to differences between WSEs

#see how well calibrated 30m flood heights predict calibrated 10m flood heights

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

n10<- 104

library(terra)

dem10<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_10m_new.asc")
dem30<- rast("C:/Users/svr5482/Downloads/NewDEMS/norristown_30m_new.asc")

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

WSE30mat10m<- flood.30mat10m$Run_30m[keepInds] + elev.30mat10m$norristown_30m_new[keepInds]
WSE10m<- flood.10m$Run_10m[keepInds]+elev.10m$norristown_10m_new[keepInds]

wsediff<- WSE30mat10m - WSE10m

fit<- lm(WSE10m~ poly(WSE30mat10m,2))
summary(fit)

resids<- fit$residuals

residsbyloc.df<- as.data.frame(cbind(coords.10m[keepInds,],
                                         fit$residuals))
colnames(residsbyloc.df)<-c("x","y","value")

library(ggplot2)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSEmodresids_spatial.jpeg", width = 700, height = 700)
ggplot(residsbyloc.df, aes(x, y, color=value))+
  geom_point(shape="square",size=2) + ggtitle("Residuals by location") +theme_bw()+
  scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")
dev.off()


flood.30mat10m_tokeep<- flood.30mat10m$Run_30m[keepInds]
zerothirtyInds<- which(flood.30mat10m_tokeep==0)


elev.diff<- elev.30mat10m$norristown_30m_new-elev.10m$norristown_10m_new


jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/costdiffVSWSEmodresids.jpeg",
     width = 700, height = 700)
plot(costdiffs[keepInds],resids,main="Cost differences VS residuals")
dev.off()


jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/elevdiffVSWSEmodresids.jpeg",
     width = 700, height = 700)
plot(elev.diff[keepInds],fit$residuals,main="Elevation differences VS residuals")
dev.off()




#good, all points outside river are represented

wsediffVScostdiff.df<- data.frame("WSE_diff"= wsediff[keepInds],
                                  "cost_diff"= costdiffs[keepInds],
                                  "elev_diff"= elev.diff[keepInds])



#fit<- lm(WSE_diff~ cost_diff, wsediffVScostdiff.df)

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/costdiffVSWSEdiff.jpeg",
     width = 500, height = 500)
plot(wsediffVScostdiff.df$cost_diff,
     wsediffVScostdiff.df$WSE_diff,
     xlab= "Cost difference", ylab= "WSE difference")
dev.off()

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/elevdiffVSWSEdiff.jpeg",
     width = 500, height = 500)
plot(wsediffVScostdiff.df$elev_diff,
     wsediffVScostdiff.df$WSE_diff,
     xlab= "Elevation difference", ylab= "WSE difference")
dev.off()

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
