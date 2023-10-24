#explore whether elevation, slope, TRI, flow accumulation, and cost at 10m
#and differences between these variables at 50m and 10m impact residuals
#from either model WSE10m ~ poly(WSE50mat10m,2)
#OR if they impact downscaled.z.vec- WSE10m or WSE50mat10m-WSE10m 

rm(list=ls())

if(dir.exists("C:/Users/svr5482/Reification/Philly/plots")==F){
  dir.create("C:/Users/svr5482/Reification/Philly/plots")
}

n10<- 104

library(terra)



#load DEMs
dem10<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_10m_new.asc")
dem50<- rast("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD/norristown_50m_new.asc")

quantiles<- c(.01,.05,.1,.25,.5,.75,.9,.95,.99)

for(i in 1:9){
  #load calibrated runs
  Run10m<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/Norristown/nCh/simplecalQs/Extent/Run_",i,".asc"))
  Run50m<- rast(paste0("C:/Users/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/Norristown/nCh/simplecalQs/Extent/Run_",i,".asc"))
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
  for(j in 1:n10){
    c10<- costs10[[j]]
    c50<- costs50at10[[j]]
    coords<- coords_OutsideRiverbyX[[j]]
    inds<- coords[,1]
    allInds<- c(allInds,inds)
    cdiff<- c50-c10
    costdiffs[inds]<- cdiff
  }
  allInds<- allInds[-1]
  
  length(which(!is.na(costdiffs)))==length(allInds) 
  
  keepInds<- intersect(which(!is.na(costdiffs)),goodinds)
  
  #compute WSE at each resolution on 10m grid
  
  WSE50mat10m<- flood.50mat10m$max[keepInds] + elev.50mat10m$norristown_50m_new[keepInds]
  WSE10m<- flood.10m$max[keepInds]+elev.10m$norristown_10m_new[keepInds]
  
  #compute difference in calibrated WSEs at 50m on 10m grid vs 10m
  wsediff<- WSE50mat10m - WSE10m
  
  WSH50mat10m<- flood.50mat10m$max[keepInds] + 
    elev.50mat10m$norristown_50m_new[keepInds] - elev.10m$norristown_10m_new[keepInds]
  WSH10m<- flood.10m$max[keepInds] 
  
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
  z1= matrix(vals.50m$max, nrow= ny50m, ncol= nx50m, byrow= TRUE)
  
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
  
  downscale50m<- rep(NA,length(flood.10m$max))
  downscale50m[coordsIwantInds]<- downscaled.z.vec
  
  downscale50m<- downscale50m[indstoCompare]
  
  WSHtocompare<- flood.10m$max[indstoCompare] 
  
  #plot(downscale50m,WSHtocompare,main="Calibrated, downscaled 50m WSH VS calibrated 10m WSH",
  #     xlab="Calibrated, downscaled 50m WSH", ylab= "Calibrated 10m WSH")
  
  #how well can we predict calibrated WSH at 10m using calibrated WSH at 50m on 10m grid?
  fit<- lm(WSHtocompare~ downscale50m)
  #summary(fit)
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
  
  jpeg(filename=paste0("C:/Users/svr5482/Reification/Philly/plots/dnsclWSHmodresids_spatial_q",quantiles[i],".jpeg"), width = 800, height = 700)
  print(ggplot(residsbyloc.df, aes(x, y, color=value))+
          geom_point(shape="square",size=1.8) + ggtitle(paste0("Q",quantiles[i],"- Residuals in space- WSHtocompare~ downscale50m")) +theme_bw()+
          scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red"))
  dev.off()
  
  #now assume we can only predict WSH at 10m with downscale 50m
  
  length(which(downscale50m<0))
  
  downscale50m[which(downscale50m<0)]<- 0
  
  resids<- WSHtocompare- downscale50m
  
  summary(resids)
  #50% of calibrated 50m flood heights are 
  #within (-.4667, .4113) m of the calibrate 10m flood height
  
  #make a dataframe to plot residuals in space
  residsbyloc.df<- as.data.frame(cbind(indstoCompare,coords.10m[indstoCompare,],resids))
  colnames(residsbyloc.df)<-c("inds","x","y","value")
  
  save(residsbyloc.df,file=paste0("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH10-downscale50_q",quantiles[i],".RData"))
  
  #load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")
  #load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnsclno50to10.RData")
  
  jpeg(filename=paste0("C:/Users/svr5482/Reification/Philly/plots/dnscl,cal50-cal10resids_spatial_q",quantiles[i],".jpeg"), width = 800, height = 700)
  print(ggplot(residsbyloc.df, aes(x, y, color=value))+
          geom_point(shape="square",size=1.8) + ggtitle(paste0("Q",quantiles[i],"- Residuals in space- WSH10m- downscale50m")) +theme_bw()+
          #scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red",limits=c(minres,maxres))+
          scale_color_gradient2(midpoint=0, low="blue",mid="gray", high="red")+
          theme(plot.title = element_text(size=24),
                axis.title.x = element_text(size = 24),
                axis.text.x = element_text(size = 18),
                axis.title.y = element_text(size = 24),
                axis.text.y = element_text(size = 18)))
  dev.off()
  
}
