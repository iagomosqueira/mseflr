#library(debug)
#mtrace(setPopDiffuse)
d.<-diffussion
d.[,  ,]<-0
d.[,1,1]<-1
d.[,2,2]<-1
d.[,3,3]<-1
d.[,4,4]<-1

pop.1<-setPopDiffuse(simPop,run.[1,c("F1","F2","F3","F4")],run.[1,c("SR1","SR2","SR3","SR4")],sel,avail,srs,srDeviates,projPeriod,diffussion,unlist(run.[1,c("Shift","yrShift","fShift","srShift.a","srShift.b")]))
pop.2<-setPopDiffuse(simPop,run.[1,c("F1","F2","F3","F4")],run.[1,c("SR1","SR2","SR3","SR4")],sel,avail,srs,srDeviates,projPeriod,d.,        unlist(run.[1,c("Shift","yrShift","fShift","srShift.a","srShift.b")]))

plot(window(FLStocks(pop.1[[1]],pop.2[[1]]),start=2004))

plot(window(pop.1,start=2004))
windows()
plot(window(pop.2,start=2004))

plot(window(FLStocks(pop.1[[1]],pop.2[[1]]),start=2004))
plot(window(FLStocks(pop.1[[2]],pop.2[[2]]),start=2004))
plot(window(FLStocks(pop.1[[3]],pop.2[[3]]),start=2004))
plot(window(FLStocks(pop.1[[4]],pop.2[[4]]),start=2004))
