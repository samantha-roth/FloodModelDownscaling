
rm(list=ls())

library(ggplot2)

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH10-downscale30.RData")
residsbyloc.df10_30<- residsbyloc.df
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH10-downscale50.RData")
residsbyloc.df10_50<- residsbyloc.df
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/WSH30-downscale50.RData")
residsbyloc.df30_50<- residsbyloc.df

overlap.inds<- which(residsbyloc.df10_30$inds%in%residsbyloc.df10_50$inds)
residsbyloc.df10_30<- residsbyloc.df10_30[overlap.inds,]

jpeg(filename="C:/Users/svr5482/Reification/Philly/plots/WSH10m-dwnsclWSH50mVSWSH10m-dwnsclWSH30m.jpeg", 
     width = 700, height = 700)
plot(residsbyloc.df10_50$value,residsbyloc.df10_30$value,
     xlab= "WSH10m - downscaled WSH50m", ylab= "WSH30m - downscaled WSH50m")
dev.off()

summary(residsbyloc.df10_30$value)
summary(residsbyloc.df10_50$value)

minres<- min(c(residsbyloc.df10_30$value,residsbyloc.df10_50$value,residsbyloc.df30_50$value))

maxres<- max(c(residsbyloc.df10_30$value,residsbyloc.df10_50$value,residsbyloc.df30_50$value))

save(minres,maxres,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnscl.RData")

minres<- min(c(residsbyloc.df10_30$value,residsbyloc.df30_50$value))

maxres<- max(c(residsbyloc.df10_30$value,residsbyloc.df30_50$value))

save(minres,maxres,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/simplecal/MinandMaxWSH-Dnsclno50to10.RData")
