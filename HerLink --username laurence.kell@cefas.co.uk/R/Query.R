library(FLCore)
library(RSQLite)

Con<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/histNoDif.dbf")
x  <- dbGetQuery(Con, "SELECT * FROM fishery") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table


Con1<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist1.dbf")
Con2<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist131.dbf")
Con3<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist261.dbf")
Con4<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist391.dbf")
Con5<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist521.dbf")
Con6<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist651.dbf")
Con7<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist781.dbf")
Con8<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist911.dbf")
Con9<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist1041.dbf")
Con10<-dbConnect(dbDriver("SQLite"), dbname="C:/Papers/HerringLinking/results/hist1171.dbf")

x1   <- dbGetQuery(Con1, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x2   <- dbGetQuery(Con2, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x3   <- dbGetQuery(Con3, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x4   <- dbGetQuery(Con4, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x5   <- dbGetQuery(Con5, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x6   <- dbGetQuery(Con6, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x7   <- dbGetQuery(Con7, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x8   <- dbGetQuery(Con8, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x9   <- dbGetQuery(Con9, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table
x10  <- dbGetQuery(Con10, "SELECT * FROM historic") # WHERE steep = 0.75 AND stock = "her"")[c(-1)] #query table

x<-rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
rm(      x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
dbDisconnect(Con1)
dbDisconnect(Con2)
dbDisconnect(Con3)
dbDisconnect(Con4)
dbDisconnect(Con5)
dbDisconnect(Con6)
dbDisconnect(Con7)
dbDisconnect(Con8)
dbDisconnect(Con9)
dbDisconnect(Con10)

#densityplot( ~ data | refpt, data = refpts(MPBrp)[, "harvest",c(1:4,6) ], layout = c(1, 5), xlab = "Fishing mortality", scale="free", xlim=c(0.1,.4))

load("C:\\Papers\\HerringLinking\\data\\conditioning.RData")
x<-merge(x, brps, by.x="fishery", by.y="popln")
attach(x)
biasFMSY<-(FSPR30-Fmsy)/Fmsy
biasBMSY<-(BSPR30-Bmsy)/Bmsy
biasMSY <-(YSPR30-MSY)/MSY
#biasFbr <-(fbr_MP-fbr_OM)/Fmsy
#biasSSB <-(ssb_MP-ssb_OM)/Bmsy
biasFbr <-(x[,17]-fbr_OM)/Fmsy
biasSSB <-(x[,18]-ssb_OM)/Bmsy
biasYld <-(yld_MP-yld_OM)/MSY

par(mfrow=c(2,3))
hist(biasFMSY,xlim=c(-2.5,2.5), breaks=seq(-2,10,  length.out=120),  main=expression(F[MSY]))
hist(biasBMSY,xlim=c(-2.5,2.5), breaks=seq(-2.5,20,length.out=240),  main=expression(B[MSY]))
hist(biasMSY, xlim=c(-2.5,2.5), breaks=seq(-2,10,  length.out=120),  main=expression(MSY))
hist(biasFbr, xlim=c(-2.5,2.5), breaks=seq(-6,2.5, length.out= 80),  main="Fbar")
hist(biasSSB, xlim=c(-2.5,2.5), breaks=seq(-5,30.0,length.out=360),  main="SSB")
hist(biasYld, xlim=c(-2.5,2.5),                                      main="Yield")

t.<-apply( biasYld, cbind(x[,"F1"],x[,"F2"],x[,"F3"],x[,"F4"],x[,"cpue"]),median,na.rm=T) 
bwplot(biasYld~(x[,"F1"]+x[,"F2"]+x[,"F3"]+x[,"F4"])|x[,"cpue"]) 

t.<-tapply( biasMSY, x[,c("F1","F2","F3","F4","cpue")],median,na.rm=T) 
t.<-cbind(data=c(t.),expand.grid(dimnames(t.)))
bwplot(t.[,"data"]~t.[,"F3"]|t.[,"cpue"]) 

t.<-tapply( biasFbr, x[,c("F1","F2","F3","F4","cpue")],quantile, 0.75,na.rm=T) 
t.<-cbind(data=c(t.),expand.grid(dimnames(t.)))
bwplot(t.[,"data"]~t.[,"F1"]|t.[,"F2"]*t.[,"F3"]*t.[,"cpue"])