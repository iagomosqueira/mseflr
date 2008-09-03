library(FLCore)
library(FLAssess)
library(FLash)
library(FLBRP)
library(FLXSA)
library(DBI)
library(RSQLite)
library(debug)

source("C:/FLR/packages/FLBRP/R/FLBRP-class.R")
source("C:/FLR/packages/FLBRP/R/FLBRP-methods.R")

my.dir<-"c:/papers/HerringLinking"
#my.dir<-"z:/pc2900/Papers/HerringLinking/"

## convert steepness & virgin biomass parameters to alpha and beta
BHSVPar<-function(sr) return(FLPar(array(sv2ab(params(sr)[,"s"],params(sr)[,"v"],params(sr)[,"spr0"],"bevholt"),dim=c(dim(params(sr))[1],2))))

## Read in historic populations
w    <-read.csv(paste(my.dir,"/inputs/mass.csv",sep=""))
f    <-read.csv(paste(my.dir,"/inputs/f.csv",   sep=""))
n    <-read.csv(paste(my.dir,"/inputs/n.csv",   sep=""))

p1<-FLStock(stock.n     =FLQuant(unlist(n[n[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            harvest     =FLQuant(unlist(f[f[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004),units="f"),
            m           =FLQuant(c(1.0,0.3,0.2,rep(0.1,6)),             dimnames=list(age=1:9,year=1970:2004)),
            stock.wt    =FLQuant(unlist(w[w[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            catch.wt    =FLQuant(unlist(w[w[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            landings.wt =FLQuant(unlist(w[w[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            discards.wt =FLQuant(unlist(w[w[,"popln"]==1,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            mat         =FLQuant(c(0.1,0.5,0.8,rep(1.0,6)),             dimnames=list(age=1:9,year=1970:2004)),
            m.spwn      =FLQuant(0.75,                                  dimnames=list(age=1:9,year=1970:2004)),
            harvest.spwn=FLQuant(0.75,                                  dimnames=list(age=1:9,year=1970:2004)))

stock(     p1)  <-computeStock(p1)
catch.n(   p1)  <-stock.n(p1)*harvest(p1)/(harvest(p1)+m(p1))*(1-exp(-harvest(p1)-m(p1)))
landings.n(p1)  <-catch.n(p1)
discards.n(p1)[]<-0
catch(     p1)  <-computeCatch(   p1)
landings(  p1)  <-computeLandings(p1)
discards(  p1)  <-computeDiscards(p1)
range(     p1,c("minfbar","maxfbar"))<-c(4,6) 

p2<-FLStock(stock.n     =FLQuant(unlist(n[n[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            harvest     =FLQuant(unlist(f[f[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004),units="f"),
            m           =FLQuant(c(1.0,0.3,0.2,rep(0.1,6)),             dimnames=list(age=1:9,year=1970:2004)),
            stock.wt    =FLQuant(unlist(w[w[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            catch.wt    =FLQuant(unlist(w[w[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            landings.wt =FLQuant(unlist(w[w[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            discards.wt =FLQuant(unlist(w[w[,"popln"]==2,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            mat         =FLQuant(c(0.1,0.5,0.8,rep(1.0,6)),             dimnames=list(age=1:9,year=1970:2004)),
            m.spwn      =FLQuant(0.25,                                  dimnames=list(age=1:9,year=1970:2004)),
            harvest.spwn=FLQuant(0.25,                                  dimnames=list(age=1:9,year=1970:2004)))

stock(     p2)  <-computeStock(p2)
catch.n(   p2)  <-stock.n(p2)*harvest(p2)/(harvest(p2)+m(p2))*(1-exp(-harvest(p2)-m(p2)))
landings.n(p2)  <-catch.n(p2)
discards.n(p2)[]<-0
catch(     p2)  <-computeCatch(   p2)
landings(  p2)  <-computeLandings(p2)
discards(  p2)  <-computeDiscards(p2)
range(     p2,c("minfbar","maxfbar"))<-c(4,6) 

p3<-FLStock(stock.n     =FLQuant(unlist(n[n[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            harvest     =FLQuant(unlist(f[f[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004),units="f"),
            m           =FLQuant(c(1.0,0.3,0.2,rep(0.1,6)),             dimnames=list(age=1:9,year=1970:2004)),
            stock.wt    =FLQuant(unlist(w[w[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            catch.wt    =FLQuant(unlist(w[w[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            landings.wt =FLQuant(unlist(w[w[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            discards.wt =FLQuant(unlist(w[w[,"popln"]==3,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            mat         =FLQuant(c(0.2,0.7,0.8,rep(1.0,6)),             dimnames=list(age=1:9,year=1970:2004)),
            m.spwn      =FLQuant(0.75,                                  dimnames=list(age=1:9,year=1970:2004)),
            harvest.spwn=FLQuant(0.75,                                  dimnames=list(age=1:9,year=1970:2004)))

stock(     p3)  <-computeStock(p3)
catch.n(   p3)  <-stock.n(p3)*harvest(p3)/(harvest(p3)+m(p3))*(1-exp(-harvest(p3)-m(p3)))
landings.n(p3)  <-catch.n(p3)
discards.n(p3)[]<-0
catch(     p3)  <-computeCatch(   p3)
landings(  p3)  <-computeLandings(p3)
discards(  p3)  <-computeDiscards(p3)
range(     p3,c("minfbar","maxfbar"))<-c(4,6) 

p4<-FLStock(stock.n     =FLQuant(unlist(n[n[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            harvest     =FLQuant(unlist(f[f[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004),units="f"),
            m           =FLQuant(c(1.0,0.3,0.2,rep(0.1,6)),             dimnames=list(age=1:9,year=1970:2004)),
            stock.wt    =FLQuant(unlist(w[w[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            catch.wt    =FLQuant(unlist(w[w[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            landings.wt =FLQuant(unlist(w[w[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            discards.wt =FLQuant(unlist(w[w[,"popln"]==4,][,c(-1,-2)]), dimnames=list(age=1:9,year=1970:2004)),
            mat         =FLQuant(c(0.2,0.8,rep(1.0,7)),                 dimnames=list(age=1:9,year=1970:2004)),
            m.spwn      =FLQuant(0.25,                                  dimnames=list(age=1:9,year=1970:2004)),
            harvest.spwn=FLQuant(0.25,                                  dimnames=list(age=1:9,year=1970:2004)))

stock(     p4)  <-computeStock(p4)
catch.n(   p4)  <-stock.n(p4)*harvest(p4)/(harvest(p4)+m(p4))*(1-exp(-harvest(p4)-m(p4)))
landings.n(p4)  <-catch.n(p4)
discards.n(p4)[]<-0
catch(     p4)  <-computeCatch(   p4)
landings(  p4)  <-computeLandings(p4)
discards(  p4)  <-computeDiscards(p4)
range(     p4,c("minfbar","maxfbar"))<-c(4,6) 

poplns<-FLStocks("1"=p1,"2"=p2,"3"=p3,"4"=p4)

rm(p1,p2,p3,p4,f,n,w)

plot(poplns,layout=c(4,1))
savePlot(paste(my.dir,"/figs/Poplns.jpeg",sep=""),type="jpeg")

xyplot(data~year|qname,data=lapply(poplns,ssb),  type="l",layout=c(4,1),scale="free", main="SSB")
savePlot(paste(my.dir,"/figs/SSB.jpeg",sep=""),  type="jpeg")
xyplot(data~year|qname,data=lapply(poplns,fbar), type="l",layout=c(4,1),scale="free", main="F")
savePlot(paste(my.dir,"/figs/fbar.jpeg",sep=""), type="jpeg")
xyplot(data~year|qname,data=lapply(poplns,catch),type="l",layout=c(4,1),scale="free", main="Yield")
savePlot(paste(my.dir,"/figs/catch.jpeg",sep=""),type="jpeg")
xyplot(data~year|qname,data=lapply(poplns,rec),  type="l",layout=c(4,1),scale="free", main="Recruits")
savePlot(paste(my.dir,"/figs/rec.jpeg",sep=""),  type="jpeg")

## Selection pattern
xyplot(data~age,groups=qname,data=lapply(poplns,function(x) apply(sweep(harvest(x)[,ac(1981:2004)],2,fbar(x)[,ac(1981:2004)],"/"),1,mean)),type="l",ylab="Selection-at-age",xlab="Age")
savePlot(paste(my.dir,"/figs/Selection.jpeg",sep=""),type="jpeg")

####### Stock Recruitment relationships ########################################
#### Bevholt steepness=0.75 ####################################################
## SR
sr.bh.75.1<-as.FLSR(poplns[[1]])
model(sr.bh.75.1)<-bevholt.sv()
sr.bh.75.1<-mle(sr.bh.75.1,fixed=list(s=0.75,spr0=spr0(FLBRP(poplns[[1]],mnYrs=35)))) #,nyrs=35))))

sr.bh.75.2<-as.FLSR(poplns[[2]])
model(sr.bh.75.2)<-bevholt.sv()
sr.bh.75.2<-mle(sr.bh.75.2,fixed=list(s=0.75,spr0=spr0(FLBRP(poplns[[2]],mnYrs=35)))) #,nyrs=35))))

sr.bh.75.3<-as.FLSR(poplns[[3]])
model(sr.bh.75.3)<-bevholt.sv()
sr.bh.75.3<-mle(sr.bh.75.3,fixed=list(s=0.75,spr0=spr0(FLBRP(poplns[[3]],mnYrs=35)))) #,nyrs=35))))

sr.bh.75.4<-as.FLSR(poplns[[4]])
model(sr.bh.75.4)<-bevholt.sv()
sr.bh.75.4<-mle(sr.bh.75.4,fixed=list(s=0.75,spr0=spr0(FLBRP(poplns[[4]],mnYrs=35)))) #,nyrs=35))))

#### Bevholt steepness=0.90 ####################################################
## SR
sr.bh.90.1<-as.FLSR(poplns[[1]])
model(sr.bh.90.1)<-bevholt.sv()
sr.bh.90.1<-mle(sr.bh.90.1,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[1]],mnYrs=35)))) #,nyrs=35))))

sr.bh.90.2<-as.FLSR(poplns[[2]])
model(sr.bh.90.2)<-bevholt.sv()
sr.bh.90.2<-mle(sr.bh.90.2,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[2]],mnYrs=35)))) #,nyrs=35))))

sr.bh.90.3<-as.FLSR(poplns[[3]])
model(sr.bh.90.3)<-bevholt.sv()
sr.bh.90.3<-mle(sr.bh.90.3,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[3]],mnYrs=35)))) #,nyrs=35))))

sr.bh.90.4<-as.FLSR(poplns[[4]])
model(sr.bh.90.4)<-bevholt.sv()
sr.bh.90.4<-mle(sr.bh.90.4,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[4]],mnYrs=35)))) #,nyrs=35))))
  
srs<-cbind(popln=rep(1:4,2),steepness=rep(c(.75,.90),each=4),
     rbind(BHSVPar(sr.bh.75.1)@.Data[1,,drop=T],
           BHSVPar(sr.bh.75.2)@.Data[1,,drop=T],
           BHSVPar(sr.bh.75.3)@.Data[1,,drop=T],
           BHSVPar(sr.bh.75.4)@.Data[1,,drop=T],
           BHSVPar(sr.bh.90.1)@.Data[1,,drop=T],
           BHSVPar(sr.bh.90.2)@.Data[1,,drop=T],
           BHSVPar(sr.bh.90.3)@.Data[1,,drop=T],
           BHSVPar(sr.bh.90.4)@.Data[1,,drop=T]))
           
####### Biological Reference Points ############################################ 
## steepness =0.75 #############################################################
brp.bh.75.1<-FLBRP(poplns[[1]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.75.1<-FLBRP(poplns[[1]],sr=sr.bh.75.1,fbar=seq(0,0.5,length.out=100))
brp.bh.75.1<-brp(brp.bh.75.1)

brp.bh.75.2<-FLBRP(poplns[[2]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.75.2<-FLBRP(poplns[[1]],sr=sr.bh.75.2,fbar=seq(0,0.5,length.out=100))
brp.bh.75.2<-brp(brp.bh.75.2)

brp.bh.75.3<-FLBRP(poplns[[3]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.75.3<-FLBRP(poplns[[1]],sr=sr.bh.75.3,fbar=seq(0,0.5,length.out=100))
brp.bh.75.3<-brp(brp.bh.75.3)

brp.bh.75.4<-FLBRP(poplns[[4]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.75.4<-FLBRP(poplns[[1]],sr=sr.bh.75.4,fbar=seq(0,0.5,length.out=100))
brp.bh.75.4<-brp(brp.bh.75.4)

## steepness =0.90 #############################################################
brp.bh.90.1<-FLBRP(poplns[[1]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.90.1<-FLBRP(poplns[[1]],sr=sr.bh.90.1,fbar=seq(0,0.5,length.out=100))
brp.bh.90.1<-brp(brp.bh.90.1)

brp.bh.90.2<-FLBRP(poplns[[2]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.90.2<-FLBRP(poplns[[1]],sr=sr.bh.90.2,fbar=seq(0,0.5,length.out=100))
brp.bh.90.2<-brp(brp.bh.90.2)

brp.bh.90.3<-FLBRP(poplns[[3]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.90.3<-FLBRP(poplns[[1]],sr=sr.bh.90.3,fbar=seq(0,0.5,length.out=100))
brp.bh.90.3<-brp(brp.bh.90.3)

brp.bh.90.4<-FLBRP(poplns[[4]],fbar=seq(0,0.5,length.out=100),mnYrs=35) #,nyrs=35)
brp.bh.90.4<-FLBRP(poplns[[1]],sr=sr.bh.90.4,fbar=seq(0,0.5,length.out=100))
brp.bh.90.4<-brp(brp.bh.90.4)

dimnames(refpts(brp.bh.75.1))$refpt[5]<-"FCrash"
refpts(brp.bh.75.1)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.75.1)<-computeRefpts(brp.bh.75.1)

dimnames(refpts(brp.bh.75.2))$refpt[5]<-"FCrash"
refpts(brp.bh.75.2)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.75.2)<-computeRefpts(brp.bh.75.2)

dimnames(refpts(brp.bh.75.3))$refpt[5]<-"FCrash"
refpts(brp.bh.75.3)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.75.3)<-computeRefpts(brp.bh.75.3)

dimnames(refpts(brp.bh.75.4))$refpt[5]<-"FCrash"
refpts(brp.bh.75.4)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.75.4)<-computeRefpts(brp.bh.75.4)

dimnames(refpts(brp.bh.90.1))$refpt[5]<-"FCrash"
refpts(brp.bh.90.1)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.90.1)<-computeRefpts(brp.bh.90.1)

dimnames(refpts(brp.bh.90.2))$refpt[5]<-"FCrash"
refpts(brp.bh.90.2)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.90.2)<-computeRefpts(brp.bh.90.2)

dimnames(refpts(brp.bh.90.3))$refpt[5]<-"FCrash"
refpts(brp.bh.90.3)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.90.3)<-computeRefpts(brp.bh.90.3)

dimnames(refpts(brp.bh.90.4))$refpt[5]<-"FCrash"
refpts(brp.bh.90.4)[1,     ,"FCrash"]<-c(NA,NA,NA,0,NA)
refpts(brp.bh.90.4)<-computeRefpts(brp.bh.90.4)

brps<-rbind(cbind(rbind(refpts(brp.bh.75.1)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.75.2)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.75.3)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.75.4)@.Data[1,c(1,2,3,4),"msy"]),
                  rbind(refpts(brp.bh.75.1)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.75.2)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.75.3)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.75.4)@.Data[1,1,         "FCrash"])),
            cbind(rbind(refpts(brp.bh.90.1)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.90.2)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.90.3)@.Data[1,c(1,2,3,4),"msy"],
                        refpts(brp.bh.90.4)@.Data[1,c(1,2,3,4),"msy"]),
                  rbind(refpts(brp.bh.90.1)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.90.2)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.90.3)@.Data[1,1,         "FCrash"],
                        refpts(brp.bh.90.4)@.Data[1,1,         "FCrash"])))

dimnames(brps)[[2]]<-c("Fmsy","MSY","Rmsy","Bmsy","Fcrash")
brps<-cbind(popln=rep(1:4,2),steepness=rep(c(.75,.90),each=4),brps)

brps12 <-brps[brps[,"popln"] %in% 1:2,]
brps123<-brps[brps[,"popln"] %in% 1:3,]

brps<-rbind(brps,
            cbind(12, c(0.75,0.9),NA,tapply(brps12[, "MSY"], brps12[, 2],sum),
                                     tapply(brps12[, "Rmsy"],brps12[, 2],sum),
                                     tapply(brps12[, "Bmsy"],brps12[, 2],sum),NA),
            cbind(123,c(0.75,0.9),NA,tapply(brps123[,"MSY"], brps123[,2],sum),
                                     tapply(brps123[,"Rmsy"],brps123[,2],sum),
                                     tapply(brps123[,"Bmsy"],brps123[,2],sum),NA))
dimnames(brps)[[1]]<-1:12
rm(brps12,brps123)

brps.rel<-brps
brps.rel[,"Bmsy"]<-brps[,"Bmsy"]/max(brps[,"Bmsy"])
brps.rel[,"MSY"] <-brps[, "MSY"]/max(brps[,"MSY"])
              
rm(brp.bh.75.1,brp.bh.75.2,brp.bh.75.3,brp.bh.75.4,
   brp.bh.90.1,brp.bh.90.2,brp.bh.90.3,brp.bh.90.4)

rm(sr.bh.75.1, sr.bh.75.2, sr.bh.75.3, sr.bh.75.4,
   sr.bh.90.1, sr.bh.90.2, sr.bh.90.3, sr.bh.90.4)

rm(BHSVPar)

## Set up simulations ##########################################################
histMinYr   <-1995
histMaxYr   <-2024
futureMaxYr <-2044
histPeriod  <-ac(histMinYr:histMaxYr)
projPeriod  <-ac(histMinYr:futureMaxYr)
nits        <-1000

simPop<-lapply(poplns,stf,nyrs=41)
simPop<-lapply(simPop,propagate,iter=nits)

## MC Stuff
srDeviates  <-FLQuant(exp(rnorm(nits*length(projPeriod)*4,0,.5))/exp((0.5^2)/2),dimnames=list(age=1,year=projPeriod,unit=1:4,iter=1:nits))
cpueDeviates<-FLQuant(exp(rnorm(9*length(projPeriod)*4*nits,0,0.3))/exp((0.3^2)/2),dimnames=list(age=1:9,year=projPeriod,unit=1:4,iter=1:nits))
sel         <-lapply(poplns,function(x) { FLQuant(c(sweep(harvest(x),2,fbar(x),"/")[,ac(sample(1981:2004,nits*20,replace=T))]),,dimnames=list(age=1:9,year=projPeriod,iter=1:nits))})

## Availability
avail      <-array(read.csv(paste(my.dir,"/inputs/mixing.csv",sep=""),header=F)[,1],c(9,4,4),dimnames=list(age=1:9,popln=1:4,fishery=1:4))

## Diffusion vectors
diffussion <-array(read.csv(paste(my.dir,"/inputs/mixing.csv",sep=""),header=F)[,2],c(9,4,4),dimnames=list(age=1:9,popln=1:4,popln=1:4))

save(poplns,file=paste(my.dir,"/data/poplns.RData",sep=""))

save(futureMaxYr,histMaxYr,histMinYr,histPeriod,projPeriod,
     avail,diffussion,
     simPop,
     srs,brps,brps.rel,
     cpueDeviates,srDeviates,sel,file=paste(my.dir,"/data/conditioning.RData",sep=""))
