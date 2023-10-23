
rm(list=ls())

setwd("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m/ST30") 
load("chain_MH_ch_2ds2.RData")
res10m<- res[2:nrow(res),]; rm(res)

setwd("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/30m/ST30") 
load("chain_MH_ch_2ds2.RData")
res30m<- res[2:nrow(res),]; rm(res)

setwd("/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/50m/ST30") 
load("chain_MH_ch_2ds2.RData")
res50m<- res[2:nrow(res),]; rm(res)

q10m<- quantile(res10m[,1], probs = c(.01,.05,.1,.25,.5,.75,.9,.95,.99))
q30m<- quantile(res30m[,1], probs = c(.01,.05,.1,.25,.5,.75,.9,.95,.99))
q50m<- quantile(res50m[,1], probs = c(.01,.05,.1,.25,.5,.75,.9,.95,.99))


q30m-q10m
#there is more of a difference between the quantiles of these two posterior samples for larger quantiles
#but the differences are very small. How much of a difference does this make in terms of flood projections?
mean(res30m[,1])-mean(res10m[,1])

q50m-q30m
#there is no pattern in the differences between quantiles of these two posterior samples
mean(res50m[,1])-mean(res30m[,1])

save(q10m,file="/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/10m/ST30/q10m.RData")
save(q30m,file="/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/30m/ST30/q30m.RData")
save(q50m,file="/storage/work/svr5482/Reification/Philly/data/Norristown/nCh/calibration/50m/ST30/q50m.RData")
