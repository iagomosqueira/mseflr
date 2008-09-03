#### Ricker ####################################################################
sr.rk.1<-as.FLSR(poplns[[1]])
model(sr.rk.1)<-ricker()
sr.rk.1<-mle(sr.rk.1)
plot(sr.rk.1)
savePlot("C:/Papers/HerringLinking/figs/SRRk1.jpeg",type="jpeg")

sr.rk.2<-as.FLSR(poplns[[2]])
model(sr.rk.2)<-ricker()
sr.rk.2<-mle(sr.rk.2)
plot(sr.rk.2)
savePlot("C:/Papers/HerringLinking/figs/SRRk2.jpeg",type="jpeg")

sr.rk.3<-as.FLSR(poplns[[3]])
model(sr.rk.3)<-ricker()
sr.rk.3<-mle(sr.rk.3)
plot(sr.rk.3)
savePlot("C:/Papers/HerringLinking/figs/SRRk43.jpeg",type="jpeg")

sr.rk.4<-as.FLSR(poplns[[4]])
model(sr.rk.4)<-ricker()
sr.rk.4<-mle(sr.rk.4)
plot(sr.rk.4)
savePlot("C:/Papers/HerringLinking/figs/SRRk4.jpeg",type="jpeg")

##### BRPs #####################################################################
brp.rk.1<-brp(FLBRP(poplns[[1]],sr=sr.rk.1,fbar=seq(0,0.5,length.out=100)))
plot(brp.rk.1,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRP1.jpeg",type="jpeg")

brp.rk.2<-brp(FLBRP(poplns[[2]],sr=sr.rk.2,fbar=seq(0,0.5,length.out=100)))
plot(brp.rk.2,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRP2.jpeg",type="jpeg")

brp.rk.3<-brp(FLBRP(poplns[[3]],sr=sr.rk.3,fbar=seq(0,0.5,length.out=100)))
plot(brp.rk.3,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRP3.jpeg",type="jpeg")

brp.rk.4<-brp(FLBRP(poplns[[4]],sr=sr.rk.4,fbar=seq(0,0.5,length.out=100)))
plot(brp.rk.4,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRP4.jpeg",type="jpeg")

#### Bevholt steepness=0.60 ####################################################
sr.bh.60.1<-as.FLSR(poplns[[1]])
model(sr.bh.60.1)<-bevholt.sv()
sr.bh.60.1<-mle(sr.bh.60.1,fixed=list(s=0.60,spr0=spr0(FLBRP(poplns[[1]]))))
plot(sr.bh.60.1)
savePlot("C:/Papers/HerringLinking/figs/SRBH751.jpeg",type="jpeg")

sr.bh.60.2<-as.FLSR(poplns[[2]])
model(sr.bh.60.2)<-bevholt.sv()
model(sr.2)<-bevholt.sv()
sr.bh.60.2<-mle(sr.bh.60.2,fixed=list(s=0.60,spr0=spr0(FLBRP(poplns[[2]]))))
plot(sr.bh.60.2)
savePlot("C:/Papers/HerringLinking/figs/SRBH752.jpeg",type="jpeg")

sr.bh.60.3<-as.FLSR(poplns[[3]])
model(sr.bh.60.3)<-bevholt.sv()
model(sr.3)<-bevholt.sv()
sr.bh.60.3<-mle(sr.bh.60.3,fixed=list(s=0.60,spr0=spr0(FLBRP(poplns[[3]]))))
plot(sr.bh.60.3)
savePlot("C:/Papers/HerringLinking/figs/SRBH753.jpeg",type="jpeg")

sr.bh.60.4<-as.FLSR(poplns[[4]])
model(sr.bh.60.4)<-bevholt.sv()
model(sr.4)<-bevholt.sv()
sr.bh.60.4<-mle(sr.bh.60.4,fixed=list(s=0.60,spr0=spr0(FLBRP(poplns[[4]]))))
plot(sr.bh.60.4)
savePlot("C:/Papers/HerringLinking/figs/SRBH754.jpeg",type="jpeg")

##### BRPs #####################################################################
brp.bh.60.1<-FLBRP(poplns[[1]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.60.1)<-BHSVPar(sr.bh.60.1)
sr.model( brp.bh.60.1)<-bevholt()$model
#brp.1<-FLBRP(poplns[[1]],sr=sr.bh.60.1,fbar=seq(0,0.5,length.out=100))
brp.bh.60.1<-brp(brp.bh.60.1)
plot(brp.bh.60.1,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH751.jpeg",type="jpeg")

brp.bh.60.2<-FLBRP(poplns[[2]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.60.2)<-BHSVPar(sr.bh.60.2)
sr.model( brp.bh.60.2)<-bevholt()$model
brp.bh.60.2<-brp(brp.bh.60.2)
plot(brp.bh.60.2,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH752.jpeg",type="jpeg")

brp.bh.60.3<-FLBRP(poplns[[3]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.3)<-BHSVPar(sr.bh.60.3)
sr.model( brp.bh.60.3)<-bevholt()$model
brp.bh.60.3<-brp(brp.bh.60.3)
plot(brp.bh.60.3,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH753.jpeg",type="jpeg")

brp.bh.60.4<-FLBRP(poplns[[4]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.60.4)<-BHSVPar(sr.bh.60.4)
sr.model( brp.bh.60.4)<-bevholt()$model
brp.bh.60.4<-brp(brp.bh.60.4)
plot(brp.bh.60.4,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH754.jpeg",type="jpeg")

#### Bevholt steepness=0.9 #####################################################
sr.bh.90.1<-as.FLSR(poplns[[1]])
model(sr.bh.90.1)<-bevholt.sv()
sr.bh.90.1<-mle(sr.bh.90.1,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[1]]))))
plot(sr.bh.90.1)
savePlot("C:/Papers/HerringLinking/figs/SRBH751.jpeg",type="jpeg")

sr.bh.90.2<-as.FLSR(poplns[[2]])
model(sr.bh.90.2)<-bevholt.sv()
model(sr.2)<-bevholt.sv()
sr.bh.90.2<-mle(sr.bh.90.2,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[2]]))))
plot(sr.bh.90.2)
savePlot("C:/Papers/HerringLinking/figs/SRBH752.jpeg",type="jpeg")

sr.bh.90.3<-as.FLSR(poplns[[3]])
model(sr.bh.90.3)<-bevholt.sv()
model(sr.3)<-bevholt.sv()
sr.bh.90.3<-mle(sr.bh.90.3,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[3]]))))
plot(sr.bh.90.3)
savePlot("C:/Papers/HerringLinking/figs/SRBH753.jpeg",type="jpeg")

sr.bh.90.4<-as.FLSR(poplns[[4]])
model(sr.bh.90.4)<-bevholt.sv()
model(sr.4)<-bevholt.sv()
sr.bh.90.4<-mle(sr.bh.90.4,fixed=list(s=0.90,spr0=spr0(FLBRP(poplns[[4]]))))
plot(sr.bh.90.4)
savePlot("C:/Papers/HerringLinking/figs/SRBH754.jpeg",type="jpeg")


##### BRPs #####################################################################
brp.bh.90.1<-FLBRP(poplns[[1]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.90.1)<-BHSVPar(sr.bh.90.1)
sr.model( brp.bh.90.1)<-bevholt()$model
#brp.1<-FLBRP(poplns[[1]],sr=sr.bh.90.1,fbar=seq(0,0.5,length.out=100))
brp.bh.90.1<-brp(brp.bh.90.1)
plot(brp.bh.90.1,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH751.jpeg",type="jpeg")

brp.bh.90.2<-FLBRP(poplns[[2]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.90.2)<-BHSVPar(sr.bh.90.2)
sr.model( brp.bh.90.2)<-bevholt()$model
brp.bh.90.2<-brp(brp.bh.90.2)
plot(brp.bh.90.2,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH752.jpeg",type="jpeg")

brp.bh.90.3<-FLBRP(poplns[[3]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.3)<-BHSVPar(sr.bh.90.3)
sr.model( brp.bh.90.3)<-bevholt()$model
brp.bh.90.3<-brp(brp.bh.90.3)
plot(brp.bh.90.3,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH753.jpeg",type="jpeg")

brp.bh.90.4<-FLBRP(poplns[[4]],fbar=seq(0,0.5,length.out=100))
sr.params(brp.bh.90.4)<-BHSVPar(sr.bh.90.4)
sr.model( brp.bh.90.4)<-bevholt()$model
brp.bh.90.4<-brp(brp.bh.90.4)
plot(brp.bh.90.4,obs=T)
savePlot("C:/Papers/HerringLinking/figs/BRPBH754.jpeg",type="jpeg")

#### Ricker depensation ########################################################
sr.1<-as.FLSR(poplns[[1]])
model(sr.1)<-ricker.d()
sr.rk.d.1<-mle(sr.1)
plot(sr.rk.d.1)
savePlot("C:/Papers/HerringLinking/figs/SRRkD1.jpeg",type="jpeg")

sr.2<-as.FLSR(poplns[[2]])
model(sr.2)<-ricker.d()
sr.rk.d.2<-mle(sr.2)
plot(sr.rk.d.2)
savePlot("C:/Papers/HerringLinking/figs/SRRkD2.jpeg",type="jpeg")

sr.3<-as.FLSR(poplns[[3]])
model(sr.3)<-ricker.d()
sr.rk.d.3<-mle(sr.3)
plot(sr.rk.d.3)
savePlot("C:/Papers/HerringLinking/figs/SRRkD3.jpeg",type="jpeg")

sr.4<-as.FLSR(poplns[[4]])
model(sr.4)<-ricker.d()
sr.rk.d.4<-mle(sr.4)
plot(sr.rk.d.4)
savePlot("C:/Papers/HerringLinking/figs/SRRkD4.jpeg",type="jpeg")


#### Bevholt depensation ########################################################
sr.1<-as.FLSR(poplns[[1]])
model(sr.1)<-bevholt.d()
sr.bh.d.1<-mle(sr.1)
plot(sr.bh.d.1)
savePlot("C:/Papers/HerringLinking/figs/SRBHD1.jpeg",type="jpeg")

sr.2<-as.FLSR(poplns[[2]])
model(sr.2)<-bevholt.d()
sr.bh.d.2<-mle(sr.2)
plot(sr.bh.d.2)
savePlot("C:/Papers/HerringLinking/figs/SRBHD2.jpeg",type="jpeg")

sr.3<-as.FLSR(poplns[[3]])
model(sr.3)<-bevholt.d()
sr.bh.d.3<-mle(sr.3)
plot(sr.bh.d.3)
savePlot("C:/Papers/HerringLinking/figs/SRBHD3.jpeg",type="jpeg")

sr.4<-as.FLSR(poplns[[4]])
model(sr.4)<-bevholt.d()
sr.bh.d.4<-mle(sr.4)
plot(sr.bh.d.4)
savePlot("C:/Papers/HerringLinking/figs/SRBHD4.jpeg",type="jpeg")

par(mfrow=c(2,2))
SSB=FLQuant(seq(0,max(ssb(sr.rk.1)),length.out=100))
plot(rec(sr.rk.1)~ssb(sr.rk.1),xlim=c(0,max(c(SSB))),ylim=c(0,max(c(rec(sr.rk.1)))),pch=19,xlab="SSB",ylab="Recruits",main="Popln 1")
lines(predict(sr.rk.1,   ssb=SSB)~SSB,col="blue")
lines(predict(sr.bh.90.1,ssb=SSB)~SSB,col="cyan")
lines(predict(sr.bh.60.1,ssb=SSB)~SSB,col="brown")
lines(predict(sr.bh.d.1, ssb=SSB)~SSB,col="red")
lines(predict(sr.rk.d.1, ssb=SSB)~SSB,col="green")
legend(0.0,9,c("Ricker", "Bev. Holt s=0.9", "Bev. Holt s=0.6", "Ricker + depen.", "Bev Holt. + depen."),
           col =c("cyan","blue","brown","red","green"),
           lty = 1, pch = -1,
       merge = TRUE, bg = 'gray90')

SSB=FLQuant(seq(0,max(ssb(sr.rk.2)),length.out=100))
plot(rec(sr.rk.2)~ssb(sr.rk.2),,xlim=c(0,max(c(SSB))),ylim=c(0,max(c(rec(sr.rk.2)))),pch=19,xlab="SSB",ylab="Recruits",main="Popln 3")
lines(predict(sr.rk.2,   ssb=SSB)~SSB,col="blue")
lines(predict(sr.bh.90.2,ssb=SSB)~SSB,col="cyan")
lines(predict(sr.bh.60.2,ssb=SSB)~SSB,col="brown")
lines(predict(sr.bh.d.2, ssb=SSB)~SSB,col="red")
lines(predict(sr.rk.d.2, ssb=SSB)~SSB,col="green")

SSB=FLQuant(seq(0,max(ssb(sr.rk.3)),length.out=100))
plot(rec(sr.rk.3)~ssb(sr.rk.3),,xlim=c(0,max(c(SSB))),ylim=c(0,max(c(rec(sr.rk.3)))),pch=19,xlab="SSB",ylab="Recruits",main="Popln 3")
lines(predict(sr.rk.3,   ssb=SSB)~SSB,col="blue")
lines(predict(sr.bh.90.3,ssb=SSB)~SSB,col="cyan")
lines(predict(sr.bh.60.3,ssb=SSB)~SSB,col="brown")
lines(predict(sr.bh.d.3, ssb=SSB)~SSB,col="red")
lines(predict(sr.rk.d.3, ssb=SSB)~SSB,col="green")

SSB=FLQuant(seq(0,max(ssb(sr.rk.4)),length.out=100))
plot(rec(sr.rk.4)~ssb(sr.rk.4),,xlim=c(0,max(c(SSB))),ylim=c(0,max(c(rec(sr.rk.4)))),pch=19,xlab="SSB",ylab="Recruits",main="Popln 4")
lines(predict(sr.rk.4,   ssb=SSB)~SSB,col="blue")
lines(predict(sr.bh.90.4,ssb=SSB)~SSB,col="cyan")
lines(predict(sr.bh.60.4,ssb=SSB)~SSB,col="brown")
lines(predict(sr.bh.d.4, ssb=SSB)~SSB,col="red")
lines(predict(sr.rk.d.4, ssb=SSB)~SSB,col="green")
savePlot("C:/Papers/HerringLinking/figs/SRs.jpeg",type="jpeg")

sr.sh.1<-as.FLSR(poplns[[1]])
model(sr.sh.1)<-shepherd()
sr.sh.1<-mle(sr.sh.1,fixed=list(c=2.5))
plot(sr.sh.1)
savePlot("C:/Papers/HerringLinking/figs/SRSh1.jpeg",type="jpeg")

sr.sh.2<-as.FLSR(poplns[[2]])
model(sr.sh.2)<-shepherd()
sr.sh.2<-mle(sr.sh.2,fixed=list(c=2.5))
plot(sr.sh.2)
savePlot("C:/Papers/HerringLinking/figs/SRSh2.jpeg",type="jpeg")

sr.sh.3<-as.FLSR(poplns[[3]])
model(sr.sh.3)<-shepherd()
sr.sh.3<-mle(sr.sh.3,fixed=list(c=2.5))
plot(sr.sh.3)
savePlot("C:/Papers/HerringLinking/figs/SRSh3.jpeg",type="jpeg")

sr.sh.4<-as.FLSR(poplns[[4]])
model(sr.sh.4)<-shepherd()
sr.sh.4<-mle(sr.sh.4,fixed=list(c=2.5))
plot(sr.sh.4)
savePlot("C:/Papers/HerringLinking/figs/SRSh4.jpeg",type="jpeg")
