library(lattice)
lattice.options(default.theme = canonical.theme(color = FALSE))

## SRRs
par(mfrow=c(2,2))
par(mar=c(4,4,3,1))
SSB=FLQuant(seq(0,max(ssb(sr.bh.75.1)),length.out=100))
plot(rec(sr.bh.75.1)~ssb(sr.bh.75.1),pch=19,col="black",xlim=c(0,max(SSB)),xlab="SSB",ylab="Recruits",main="Population 1")
lines(predict(sr.bh.75.1,ssb=SSB)~SSB,col="black",lwd=2)
lines(predict(sr.bh.90.1,ssb=SSB)~SSB,col="black",lwd=1)

SSB=FLQuant(seq(0,max(ssb(sr.bh.75.2)),length.out=100))
plot(rec(sr.bh.75.2)~ssb(sr.bh.75.2),pch=19,col="black",xlim=c(0,max(SSB)),xlab="SSB",ylab="Recruits",main="Population 2")
lines(predict(sr.bh.75.2,ssb=SSB)~SSB,col="black",lwd=2)
lines(predict(sr.bh.90.2,ssb=SSB)~SSB,col="black",lwd=1)

SSB=FLQuant(seq(0,max(ssb(sr.bh.75.3)),length.out=100))
plot(rec(sr.bh.75.3)~ssb(sr.bh.75.3),pch=19,col="black",xlim=c(0,max(SSB)),xlab="SSB",ylab="Recruits",main="Population 3")
lines(predict(sr.bh.75.3,ssb=SSB)~SSB,col="black",lwd=2)
lines(predict(sr.bh.90.3,ssb=SSB)~SSB,col="black",lwd=1)

SSB=FLQuant(seq(0,max(ssb(sr.bh.75.4)),length.out=100))
plot(rec(sr.bh.75.4)~ssb(sr.bh.75.4),pch=19,col="black",xlim=c(0,max(SSB)),xlab="SSB",ylab="Recruits",main="Population 4")
lines(predict(sr.bh.75.4,ssb=SSB)~SSB,col="black",lwd=2)
lines(predict(sr.bh.90.4,ssb=SSB)~SSB,col="black",lwd=1)

savePlot("C:/Papers/HerringLinking/figs/SR.jpeg",type="jpeg")


library(FLCore)
library(RSQLite)

Con<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/Final2.dbf")
xTS<-dbGetQuery(Con, "SELECT * FROM TimeSeries") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
xTS[,"popln"]       <-factor(xTS[,"popln"],levels=c(1,2,3,4,12,123),labels=c("population/stock 1","population/stock 2","population/stock 3","population/stock 4","population/stock 1&2","population/stock 1,2&3"))
xTS[,"scenario"]    <-factor(xTS[,"scenario"],levels=c(1:10),labels=c("scenario 1","scenario 2","scenario 3","scenario 4","scenario 5","scenario 6","scenario 7","scenario 8","scenario 9","scenario 10"))
xTS[xTS[,"year__1"]==1995,c("fbr_MP","ssb_MP","yld_MP","rec_MP","fbr_OM","ssb_OM","yld_OM","rec_OM")]<-NA

t1                      <-xTS[,2:19]
t2                      <-xTS[,c(2:15,20:23)]
dimnames(t1)[[2]][13:18]<-c("year","iter","Fbar","SSB","Recruits","Yield")
dimnames(t2)[[2]][13:18]<-c("year","iter","Fbar","SSB","Recruits","Yield")
xTS                     <-rbind(cbind(type="MP",t1),cbind(type="OM",t2))

#sp <- list(superpose.line=list(col=c("blue","red","black","red"), lwd=rep(c(2,1),2), lty=1))
sp <- list(superpose.line=list(col=c("blue","red","black","red"), lwd=c(2,2,2,2), lty=1))

#TS <-xTS[as.integer(xTS[,"popln"]) %in% c(1:2,5) & xTS[,"iter"]==3 & xTS[,"diffuse"]==0,]
TS <-xTS[as.integer(xTS[,"popln"]) %in% c(1:2,5) & xTS[,"iter"]==3,]
#TS <-xTS[xTS[,"iter"]==3,]

xyplot(SSB~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 2 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))

savePlot("C:/Papers/HerringLinking/figs/Scenario2SSB.jpeg",type="jpeg")

xyplot(SSB~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))

savePlot("C:/Papers/HerringLinking/figs/Scenario3SSB.jpeg",type="jpeg")

xyplot(SSB~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3 & TS[,"diffuse"]==1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))

savePlot("C:/Papers/HerringLinking/figs/Scenario3SSBDiffuse.jpeg",type="jpeg")

xyplot(SSB~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 4 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))

savePlot("C:/Papers/HerringLinking/figs/Scenario4SSB.jpeg",type="jpeg")

xyplot(Fbar~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 2 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab="Fishing Mortality",xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))
savePlot("C:/Papers/HerringLinking/figs/Scenario2Fbar.jpeg",type="jpeg")

xyplot(Fbar~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab="Fishing Mortality",xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))
savePlot("C:/Papers/HerringLinking/figs/Scenario3Fbar.jpeg",type="jpeg")

xyplot(Fbar~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3 & TS[,"diffuse"]==1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))

savePlot("C:/Papers/HerringLinking/figs/Scenario3FbarDiffuse.jpeg",type="jpeg")

xyplot(Fbar~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 4 & TS[,"diffuse"]==0,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab="Fishing Mortality",xlab="Year",

 key = list(lines = list(alpha=1,col=c("red","blue","black"),lty=rep(1,3),lwd=c(2,2,2)),
                  space="top",
                  text = list(lab = c("Actual","VPA Estimate with Fishery Index","VPA Estimate with Spawning Index")),
                  columns = 3))
savePlot("C:/Papers/HerringLinking/figs/Scenario4Fbar.jpeg",type="jpeg")

ylim.<-c(0,1.0)
xyplot(SSB~year | popln+scenario, data=TS[as.integer(TS[,"scenario"]) %in% c(2,3),],
       groups = as.factor(paste(cpue,type,diffuse)),
#       layout = c(5,3),
#       index.cond=list(1:5, 3:1),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
      auto.key = list(columns=4, points=F, lines=TRUE),
      main="FBar",scale="free" ,ylim=list(ylim.,ylim.,ylim.,ylim.)) #[1:6,1:3] #,scale="free")

xyplot(SSB~year | popln+scenario, data=TS[as.integer(TS[,"scenario"])>1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       layout = c(5,3),
       index.cond=list(1:5, 1:3),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
       auto.key = list(columns=4, points=F, lines=TRUE),
       #main="SSB") #,scale="free") #,ylim=list(c(0,.6),c(0,0.6),c(0,.6),c(0,.6),c(0,.6))) #[1:6,1:3] #,scale="free")



xyplot(Fbar~year | scenario+popln, data=xTS[as.integer(xTS[,"scenario"])>=5 & as.integer(xTS[,"popln"]) %in% 3:4 & xTS[,"iter"]==3 & xTS[,"cpue"]!=1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       layout = c(5,2),
       index.cond=list(1:5, 1:2),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
       auto.key = list(columns=4, points=F, lines=TRUE),main="FBar") #[1:6,1:3] #,scale="free")
windows()
xyplot(Recruits~year | scenario+popln, data=xTS[as.integer(xTS[,"scenario"])>=5 & as.integer(xTS[,"popln"]) %in% 3:4 & xTS[,"iter"]==3 & xTS[,"cpue"]!=1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       layout = c(5,2),
       index.cond=list(1:5, 1:2),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
       auto.key = list(columns=4, points=F, lines=TRUE),main="Recruits") #[1:6,1:3] #,scale="free")
windows()
xyplot(SSB~year | scenario+popln, data=xTS[as.integer(xTS[,"scenario"])>=5 & as.integer(xTS[,"popln"]) %in% 3:4 & xTS[,"iter"]==3 & xTS[,"cpue"]!=1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       layout = c(5,2),
       index.cond=list(1:5, 1:2),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
       auto.key = list(columns=4, points=F, lines=TRUE),main="SSB") #[1:6,1:3] #,scale="free")
windows()
xyplot(Yield~year | scenario+popln, data=xTS[as.integer(xTS[,"scenario"])>=5 & as.integer(xTS[,"popln"]) %in% 3:4 & xTS[,"iter"]==3 & xTS[,"cpue"]!=1,],
       groups = as.factor(paste(cpue,type,diffuse)),
       layout = c(5,2),
       index.cond=list(1:5, 1:2),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
#         panel.superpose(x, y, type="p",...)
       },
       par.settings = sp,
       auto.key = list(columns=4, points=F, lines=TRUE),main="Yield") #[1:6,1:3] #,scale="free")


sp <- list(superpose.line=list(col=1:4, lwd=c(1,2), lty=rep(1,2,4)))

xyplot(SSB~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       auto.key = list(columns=4, points=F, lines=TRUE),
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year")
       
xyplot(Fbar~year | popln, data=TS[as.integer(TS[,"scenario"]) %in% 3,],
       groups = as.factor(paste(cpue,type,diffuse)),
       panel = function(x, y, type, ...) {
         panel.superpose(x, y, type="l",...)
       },
       par.settings = sp,ylab=expression(SSB/B[MSY]),xlab="Year")