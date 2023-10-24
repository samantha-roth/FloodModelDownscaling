#explore whether elevation, slope, TRI, flow accumulation, and cost at 30m
#and differences between these variables at 50m and 30m impact residuals
#from either model WSE30m ~ poly(WSE50mat30m,2)
#OR if they impact downscaled.z.vec- WSE30m or WSE50mat30m-WSE30m 

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

library(terra)

#load DEMs
dem30<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_30m_new.asc")
dem50<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_50m_new.asc")

quantiles<- c(.01,.05,.1,.25,.5,.75,.9,.95,.99)
coords.30m<- xyFromCell(dem30,1:ncell(dem30))

n30<- length(unique(coords.30m))

for(i in 1:9){
  #load calibrated runs
  Run30m<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs30m/Norristown/nCh/simplecalQs/Extent/Run_",i,".asc"))
  Run50m<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/simplecalQs/Extent/Run_",i,".asc"))
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
  for(j in 1:n30){
    c30<- costs30[[j]]
    c50<- costs50at30[[j]]
    coords<- coords_OutsideRiverbyX[[j]]
    inds<- coords[,1]
    allInds<- c(allInds,inds)
    cdiff<- c50-c30
    costdiffs[inds]<- cdiff
  }
  allInds<- allInds[-1]
  
  length(which(!is.na(costdiffs)))==length(allInds) 
  
  keepInds<- intersect(which(!is.na(costdiffs)),goodinds)
  
  #compute WSE at each resolution on 30m grid
  
  WSE50mat30m<- flood.50mat30m$max[keepInds] + elev.50mat30m$norristown_50m_new[keepInds]
  WSE30m<- flood.30m$max[keepInds]+elev.30m$norristown_30m_new[keepInds]
  
  #compute difference in calibrated WSEs at 50m on 30m grid vs 30m
  wsediff<- WSE50mat30m - WSE30m
  
  WSH50mat30m<- flood.50mat30m$max[keepInds] + 
    elev.50mat30m$norristown_50m_new[keepInds] - elev.30m$norristown_30m_new[keepInds]
  WSH30m<- flood.30m$max[keepInds] 
  
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
  z1= matrix(vals.50m$max, nrow= ny50m, ncol= nx50m, byrow= TRUE)
  
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
  
  downscale50m<- rep(NA,length(flood.30m$max))
  downscale50m[coordsIwantInds]<- downscaled.z.vec
  
  downscale50m<- downscale50m[indstoCompare]
  
  WSHtocompare<- flood.30m$max[indstoCompare] 
  
  #plot(downscale50m,WSHtocompare,main="Calibrated, downscaled 50m WSH VS calibrated 30m WSH",
  #     xlab="Calibrated, downscaled 50m WSH", ylab= "Calibrated 30m WSH")
  
  #how well can we predict calibrated WSH at 30m using calibrated WSH at 50m on 30m grid?
  fit<- lm(WSHtocompare~ downscale50m)
  #summary(fit)
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
  
  jpeg(filename=paste0("C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial_q",quantiles[i],".jpeg"), width = 800, height = 700)
  print(ggplot(residsbyloc.df, aes(x, y, color=value))+
          geom_point(shape="square",size=1.8) + ggtitle(paste0("Q",quantiles[i],"- Residuals in space- WSHtocompare~ downscale50m")) +theme_bw()+
          scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red"))
  dev.off()
  
  #now assume we can only predict WSH at 30m with downscale 50m
  
  length(which(downscale50m<0))
  
  downscale50m[which(downscale50m<0)]<- 0
  
  resids<- WSHtocompare- downscale50m
  
  summary(resids)
  #50% of calibrated 50m flood heights are 
  #within (-.4667, .4113) m of the calibrate 30m flood height
  
  #make a dataframe to plot residuals in space
  residsbyloc.df<- as.data.frame(cbind(indstoCompare,coords.30m[indstoCompare,],resids))
  colnames(residsbyloc.df)<-c("inds","x","y","value")
  
  save(residsbyloc.df,file=paste0("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH30-downscale50_q",quantiles[i],".RData"))
  
  #load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")
  #load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnsclno50to30.RData")
  
  jpeg(filename=paste0("C:/Users/svr5482/Reification/Philly/plots/dnscl,cal50-cal30resids_spatial_q",quantiles[i],".jpeg"), width = 800, height = 700)
  print(ggplot(residsbyloc.df, aes(x, y, color=value))+
          geom_point(shape="square",size=1.8) + ggtitle(paste0("Q",quantiles[i],"- Residuals in space- WSH30m- downscale50m")) +theme_bw()+
          #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
          scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
          theme(plot.title = element_text(size=24),
                axis.title.x = element_text(size = 24),
                axis.text.x = element_text(size = 18),
                axis.title.y = element_text(size = 24),
                axis.text.y = element_text(size = 18)))
  dev.off()
  
}
