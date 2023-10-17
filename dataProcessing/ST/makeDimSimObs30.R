#get 30 randomly sampled locations from the simulated observation

rm(list=ls())

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/obsWE_OR_5m.RData")

floodInds<- which(obsWE_OR_5m>0)

SS=30

set.seed(888)
RS<- sample(floodInds,SS)

obsWE_RS<- obsWE_OR_5m[RS]

save(obsWE_RS,RS,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/obsWE_RS.RData")

load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/5m_closest_10m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/5m_closest_30m.RData")
load("C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/5m_closest_50m.RData")

closest10m_RS<- hwm_closest.10m[RS,]
closest30m_RS<- hwm_closest.30m[RS,]
closest50m_RS<- hwm_closest.50m[RS,]

save(closest10m_RS,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest10m_RS.RData")
save(closest30m_RS,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest30m_RS.RData")
save(closest50m_RS,file="C:/Users/svr5482/Reification/Philly/data/Norristown/nCh/closest50m_RS.RData")

