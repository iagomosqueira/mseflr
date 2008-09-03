library(FLCore)
library(FLAssess)
library(FLash)
library(FLBRP)
library(FLXSA)

regimeShift<-function(ref){

    fr1 <- function(x,rfp) {  
         x1 <- x[1]
         x2 <- x[2]
      
         refMSY              <-refpts(rfp)@.Data[1,"yield","msy"]
         refBMSY             <-refpts(rfp)@.Data[1,  "ssb","msy"]
         sr.params(rfp)[,"a"]<-sr.params(rfp)[,"a"]*x1
         sr.params(rfp)[,"b"]<-sr.params(rfp)[,"b"]*x2
         res                 <-(refMSY *.6-computeRefpts(rfp)@.Data[1,"yield","msy"])^2+
                               (refBMSY   -computeRefpts(rfp)@.Data[1,  "ssb","msy"])^2
         
         return(res)
         }
    
    fr2 <- function(x,rfp) {  
         x1 <- x[1]
         x2 <- x[2]
      
         refMSY              <-refpts(rfp)@.Data[1,   "yield",   "msy"]
         refFMSY             <-refpts(rfp)@.Data[1, "harvest",   "msy"]
         sr.params(rfp)[,"a"]<-sr.params(rfp)[,"a"]*x1
         sr.params(rfp)[,"b"]<-sr.params(rfp)[,"b"]*x2
         res                 <-(refMSY *.6-computeRefpts(rfp)@.Data[1,   "yield","msy"])^2+
                               (refFMSY   -computeRefpts(rfp)@.Data[1, "harvest","msy"])^2
         
         return(res)
         }
    
    res1  <-optim(c(0.6,0.6), fr1, rfp=ref)
    bmsy.5<-ref
    sr.params(bmsy.5)[,"a"]<-sr.params(ref)[,"a"]*res1$par[1]
    sr.params(bmsy.5)[,"b"]<-sr.params(ref)[,"b"]*res1$par[2]
    refpts(bmsy.5)[1,,"FCrash"]<-c(NA,NA,NA,0,NA)
    bmsy.5<-brp(bmsy.5)
    
    res2  <-optim(c(.8,1.8), fr2, rfp=ref)
    fmsy.5<-ref
    sr.params(fmsy.5)[,"a"]<-sr.params(ref)[,"a"]*res2$par[1]
    sr.params(fmsy.5)[,"b"]<-sr.params(ref)[,"b"]*res2$par[2]
    refpts(fmsy.5)[1,,"FCrash"]<-c(NA,NA,NA,0,NA)
    fmsy.5<-brp(fmsy.5)
    
    par(mfrow=c(1,2))
    plot( yield( ref)  ~ssb(  ref), col="black",type="l",xlab="SSB",ylab="Yield")
    lines(yield(bmsy.5)~ssb(bmsy.5),col="red")
    lines(yield(fmsy.5)~ssb(fmsy.5),col="blue")
    
    plot( yield( ref)  ~fbar(  ref), col="black",type="l",xlab="Fishing Mortality",ylab="Yield")
    lines(yield(bmsy.5)~fbar(bmsy.5),col="red")
    lines(yield(fmsy.5)~fbar(fmsy.5),col="blue")
    
    return(c(res1$par,res2$par,refpts(bmsy.5)@.Data[1,c("ssb","harvest","yield"),"msy"],refpts(fmsy.5)@.Data[1,c("ssb","harvest","yield"),"msy"]))
    }
    
srsRegimeShift<-cbind(steepness =rep(c(0.75,0.9),each=4),
                      population=rep(1:4,             2),
                      rbind(regimeShift(brp.bh.75.1),
                            regimeShift(brp.bh.75.2),
                            regimeShift(brp.bh.75.3),
                            regimeShift(brp.bh.75.4),
                            regimeShift(brp.bh.90.1),
                            regimeShift(brp.bh.90.2),
                            regimeShift(brp.bh.90.3),
                            regimeShift(brp.bh.90.4)))
              
dimnames(srsRegimeShift)[[2]]<-c("steepness","population","a.1","b.1","a.2","b.2","bmsy.1","fmsy.1","msy.1","bmsy.2","fmsy.2","msy.2")
save(srsRegimeShift,file="C:/Papers/HerringLinking/data/srsRegimeShift.RData")
