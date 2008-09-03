library(FLCore)
library(FLAssess)
library(FLash)
library(FLBRP)
library(FLXSA)
library(DBI)
library(RSQLite)
#library(debug)

diffuse<-function(iPop,iYr,Pops,diffussion)
   {
   res<- Pops[[1]]@stock.n[,ac(iYr)]*diffussion[,1,iPop]+
         Pops[[2]]@stock.n[,ac(iYr)]*diffussion[,2,iPop]+
         Pops[[3]]@stock.n[,ac(iYr)]*diffussion[,3,iPop]+
         Pops[[4]]@stock.n[,ac(iYr)]*diffussion[,4,iPop]

   return(res)
   }

setPopDiffuse<-function(Pops,Fs,SR,sel,avail,srs,srDeviates,yrs,diffussion,shift="missing")
    {
    F1<-c(Fs[,"F1"])
    F2<-c(Fs[,"F2"])
    F3<-c(Fs[,"F3"])
    F4<-c(Fs[,"F4"])

    harvest(Pops[[1]])[,yrs ]<-setSel(1,F1,F2,F3,F4,sel,avail,yrs)
    harvest(Pops[[2]])[,yrs ]<-setSel(2,F1,F2,F3,F4,sel,avail,yrs)
    harvest(Pops[[3]])[,yrs ]<-setSel(3,F1,F2,F3,F4,sel,avail,yrs)
    harvest(Pops[[4]])[,yrs ]<-setSel(4,F1,F2,F3,F4,sel,avail,yrs)

    Ftrgt<-FLQuant(NA,dimnames=list(age="all",year=yrs,unit=1:4,iter=1))
    Ftrgt[,,1]<-apply(fbar(Pops[[1]][,ac(yrs)]),2,mean)
    Ftrgt[,,2]<-apply(fbar(Pops[[2]][,ac(yrs)]),2,mean)
    Ftrgt[,,3]<-apply(fbar(Pops[[3]][,ac(yrs)]),2,mean)
    Ftrgt[,,4]<-apply(fbar(Pops[[4]][,ac(yrs)]),2,mean)

    srPars<-list()
    srPars[[1]]<-FLPar(array(c(srs[srs[,"popln"]==1 & srs[,"steepness"]==c(SR["SR1"]),3:4]),c(1,2,length(yrs)),dimnames=list(iter=1,params=c("a","b"),year=yrs)))
    srPars[[2]]<-FLPar(array(c(srs[srs[,"popln"]==2 & srs[,"steepness"]==c(SR["SR2"]),3:4]),c(1,2,length(yrs)),dimnames=list(iter=1,params=c("a","b"),year=yrs)))
    srPars[[3]]<-FLPar(array(c(srs[srs[,"popln"]==3 & srs[,"steepness"]==c(SR["SR3"]),3:4]),c(1,2,length(yrs)),dimnames=list(iter=1,params=c("a","b"),year=yrs)))
    srPars[[4]]<-FLPar(array(c(srs[srs[,"popln"]==4 & srs[,"steepness"]==c(SR["SR4"]),3:4]),c(1,2,length(yrs)),dimnames=list(iter=1,params=c("a","b"),year=yrs)))

#    srPar<-array(NA,c(1,2,length(yrs),4),dimnames=list(iter=1,params=c("a","b"),year=yrs,unit=1:4))

#    srPar[,,,1]<-array(c(srs[srs[,"popln"]==1 & srs[,"steepness"]==c(SR["SR1"]),3:4]),c(1,2,length(yrs)))
#    srPar[,,,2]<-array(c(srs[srs[,"popln"]==2 & srs[,"steepness"]==c(SR["SR2"]),3:4]),c(1,2,length(yrs)))
#    srPar[,,,3]<-array(c(srs[srs[,"popln"]==3 & srs[,"steepness"]==c(SR["SR3"]),3:4]),c(1,2,length(yrs)))
#    srPar[,,,4]<-array(c(srs[srs[,"popln"]==4 & srs[,"steepness"]==c(SR["SR4"]),3:4]),c(1,2,length(yrs)))

    if (!(missing(shift) | is.null(shift)))
       {
       Ftrgt[ ,    as.numeric(yrs)>=shift["yrShift"],shift["Shift"]]  <-Ftrgt[,     as.numeric(yrs)>=shift["yrShift"],shift["Shift"]]*shift["fShift"]
       srPars[[shift["Shift"]]][,"a",as.numeric(yrs)>=shift["yrShift"]]<-srPars[[shift["Shift"]]][,"a",as.numeric(yrs)>=shift["yrShift"]]*shift["srShift.a"]
       srPars[[shift["Shift"]]][,"b",as.numeric(yrs)>=shift["yrShift"]]<-srPars[[shift["Shift"]]][,"b",as.numeric(yrs)>=shift["yrShift"]]*shift["srShift.a"]
       }

    dmns     <-dimnames(stock.n(Pops[[1]]))
    dmns$unit<-1:4
    for (iYr in as.numeric(yrs))
       {
       ## project
       for (iPop in 1:4)
          {
          trgt                      <-fwdTarget(year=iYr,quantity="f",value=c(Ftrgt[1,ac(iYr),iPop]))
          Pops[[iPop]][,ac(iYr+0:1)]<-fwd(Pops[[iPop]],target=trgt,sr.model="bevholt",sr.params=srPars[[iPop]][,,ac(iYr)],sr.residuals=srDeviates[,,iPop])[,ac(iYr+0:1)]
          }
          
       ## Diffuse at end of year
       if (iYr<max(as.numeric(yrs)))
          {
          dmns$year<-ac(iYr+1)
          stkn     <-FLQuant(NA,dimnames=dmns)
          for (iPop in 1:4)
             stkn[,,iPop]<-diffuse(iPop,iYr+1,Pops,diffussion)
          for (iPop in 1:4)
             stock.n(Pops[[iPop]])[,ac(iYr+1)]<-stkn[,ac(iYr+1),iPop]
          }
       }
         
    return(Pops)
    }

setPop<-function(iPop,F1,F2,F3,F4,sel,avail,Pop,srs,srDeviates,yrs,shift="missing")
    {
    ## First period
    srPar <-FLPar(array(srs,c(1,2),dimnames=list(iter=1,params=c("a","b"))))
    
    harvest(Pop)[,yrs ]<-setSel(iPop,F1,F2,F3,F4,sel,avail,yrs)
    if (!(missing(shift) | is.null(shift)))
       {
       projYrs<-yrs
       yrs    <-yrs[as.numeric(yrs)<shift["yrShift"]]
       }
       
    trgt<-fwdTarget(year=as.numeric(yrs),quantity="f",value=c(apply(fbar(Pop)[,yrs],2,mean)))

    Pop <-fwd(Pop,target=trgt,sr.model="bevholt",sr.params=srPar,sr.residuals=srDeviates[,,1])

    ## Second Period
    if (!(missing(shift) | is.null(shift)))
       {
       srPar[,"a"] <-srPar[,"a"]*shift["srShift.a"]
       srPar[,"b"] <-srPar[,"b"]*shift["srShift.b"]
       
       yrs <-projYrs[as.numeric(projYrs)>=shift["yrShift"]]
       
#       harvest(Pop)[,yrs] <-setSel(iPop,F1,F2,F3,F4,sel,avail,yrs)
       trgt               <-fwdTarget(year=as.numeric(yrs),quantity="f",value=c(apply(fbar(Pop)[,yrs],2,mean))*shift["fShift"])
       Pop                <-fwd(Pop,target=trgt,sr.model="bevholt",sr.params=srPar,sr.residuals=srDeviates[,,1])  
       }
       
    return(Pop)
    }

setSel<-function(iPop,F1,F2,F3,F4,sel,avail,yrs)
    {
    res <-F1*sweep(sel[[1]][,yrs],1,avail[,iPop,1],"*")+
          F2*sweep(sel[[2]][,yrs],1,avail[,iPop,2],"*")+
          F3*sweep(sel[[3]][,yrs],1,avail[,iPop,3],"*")+
          F4*sweep(sel[[4]][,yrs],1,avail[,iPop,4],"*")

    units(res)<-"f"
    return(res)
    }

getCatch<-function(iFish,f,sel,avail,Pop,yrs,shift="missing")
    {
    f1<-array(f,dim=length(yrs),dimnames=yrs)
    f2<-f1
    f3<-f2
    f4<-f3
    
    if (!missing(shift))
       {
       if (shift["Shift"]==1) f1[as.numeric(yrs)>=shift["yrShift"]]<-f1[as.numeric(yrs)>=shift["yrShift"]]*shift["fShift"] else
       if (shift["Shift"]==1) f1[as.numeric(yrs)>=shift["yrShift"]]<-f1[as.numeric(yrs)>=shift["yrShift"]]*shift["fShift"] else
       if (shift["Shift"]==1) f1[as.numeric(yrs)>=shift["yrShift"]]<-f1[as.numeric(yrs)>=shift["yrShift"]]*shift["fShift"] else
       if (shift["Shift"]==1) f1[as.numeric(yrs)>=shift["yrShift"]]<-f1[as.numeric(yrs)>=shift["yrShift"]]*shift["fShift"]
       }
            
    res<- sweep(catch.n(Pop[[1]])[,yrs],2,f1,"*")*sweep(sel[[iFish]][,yrs],1,avail[,1,iFish],"*")/harvest(Pop[[1]][,yrs])+
          sweep(catch.n(Pop[[2]])[,yrs],2,f2,"*")*sweep(sel[[iFish]][,yrs],1,avail[,2,iFish],"*")/harvest(Pop[[2]][,yrs])+
          sweep(catch.n(Pop[[3]])[,yrs],2,f3,"*")*sweep(sel[[iFish]][,yrs],1,avail[,3,iFish],"*")/harvest(Pop[[3]][,yrs])+
          sweep(catch.n(Pop[[4]])[,yrs],2,f4,"*")*sweep(sel[[iFish]][,yrs],1,avail[,4,iFish],"*")/harvest(Pop[[4]][,yrs])

    return(res)
    }

getStock<-function(iFish,avail,Pop)
    {
    res<- avail[,1,iFish]*stock.n(Pop[[1]])+
          avail[,2,iFish]*stock.n(Pop[[2]])+
          avail[,3,iFish]*stock.n(Pop[[3]])+
          avail[,4,iFish]*stock.n(Pop[[4]])

    return(res)
    }

getCatchWt<-function(iFish,f,sel,avail,Pop,yrs)
    {
    wt1 <-catch.wt(Pop[[1]])[,yrs]*catch.n(Pop[[1]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,1,iFish],"*")/harvest(Pop[[1]][,yrs])
    wt2 <-catch.wt(Pop[[2]])[,yrs]*catch.n(Pop[[2]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,2,iFish],"*")/harvest(Pop[[2]][,yrs])
    wt3 <-catch.wt(Pop[[3]])[,yrs]*catch.n(Pop[[3]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,3,iFish],"*")/harvest(Pop[[3]][,yrs])
    wt4 <-catch.wt(Pop[[4]])[,yrs]*catch.n(Pop[[4]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,4,iFish],"*")/harvest(Pop[[4]][,yrs])

    wt1[is.na(wt1)]<-0
    wt2[is.na(wt2)]<-0
    wt3[is.na(wt3)]<-0
    wt4[is.na(wt4)]<-0

    num1<-catch.n(Pop[[1]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,1,iFish],"*")/harvest(Pop[[1]][,yrs])
    num2<-catch.n(Pop[[2]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,2,iFish],"*")/harvest(Pop[[2]][,yrs])
    num3<-catch.n(Pop[[3]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,3,iFish],"*")/harvest(Pop[[3]][,yrs])
    num4<-catch.n(Pop[[4]])[,yrs]*f*sweep(sel[[iFish]][,yrs],1,avail[,4,iFish],"*")/harvest(Pop[[4]][,yrs])

    num1[is.na(num1)]<-0
    num2[is.na(num2)]<-0
    num3[is.na(num3)]<-0
    num4[is.na(num4)]<-0

    res <-(wt1+wt2+wt3+wt4)/(num1+num2+num3+num4)

    res[is.na(res)]<-0

    return(res)
    }

getStockWt<-function(iFish,avail,Pop,yrs)
    {
    wt <- stock.n(Pop[[1]])[,yrs]*avail[,1,iFish]*stock.n(Pop[[1]])[,yrs]+
          stock.n(Pop[[2]])[,yrs]*avail[,2,iFish]*stock.n(Pop[[2]])[,yrs]+
          stock.n(Pop[[3]])[,yrs]*avail[,3,iFish]*stock.n(Pop[[3]])[,yrs]+
          stock.n(Pop[[4]])[,yrs]*avail[,4,iFish]*stock.n(Pop[[4]])[,yrs]

    num<- avail[,1,iFish]*stock.n(Pop[[1]])[,yrs]+
          avail[,2,iFish]*stock.n(Pop[[2]])[,yrs]+
          avail[,3,iFish]*stock.n(Pop[[3]])[,yrs]+
          avail[,4,iFish]*stock.n(Pop[[4]])[,yrs]

    res <-wt/num

    return(res)
    }

MP<-function(iFish,run,Pop,Stock,cpue,xsa.ctrl,key,diffuse=0)
   {
   Stock<-Stock+FLXSA(Stock,cpue,xsa.ctrl,diag.flag=FALSE)

   print(plot(lapply(FLStocks(Pop,Stock),window,start=2000,end=2044),layout=c(4,1), prob=c(0.25,0.5,0.75), xlab="Year"))
   savePlot(paste(my.dir,"figs/",scen,key,".jpeg",sep=""),type="jpeg")

   if (FALSE)
      {
      brp.           <-FLBRP(Stock)
      sr.params(brp.)<-FLPar(apply(rec(Stock),6,function(x) exp(mean(log(x))))[,1,1,1,1,,drop=T])
      refpts(brp.)   <-refpts(brp.)[,,"spr.30"]
      refpts(brp.)   <-computeRefpts(brp.)

      # summary stats
      SmryStats<-cbind(SmryStats=iFish,F1 =run[1,"F1"],  F2=run[1, "F2"], F3=run[1, "F3"], F4=run[1, "F4"],
                                   SR1=run[1,"SR1"],SR2=run[1,"SR2"],SR3=run[1,"SR3"],SR4=run[1,"SR4"],
                                   diffuse=diffuse,as.data.frame(refpts(brp.)@.Data[,c("harvest","ssb","yield"),"spr.30"]))
      SmryStats<-cbind(SmryStats,fbr.MP=c(fbar(Stock)[,ac(histMaxYr)]),ssb.MP=c(ssb(Stock)[,ac(histMaxYr)]),yld.MP=c(computeCatch(Stock)[,ac(histMaxYr)]))
      SmryStats<-cbind(SmryStats,fbr.OM=c(fbar(Pop)[,ac(histMaxYr)]),ssb.OM=c(ssb(Pop)[,ac(histMaxYr)]),yld.OM=c(computeCatch(Pop,na.rm=T)[,ac(histMaxYr)]))

      names(SmryStats)[11:13]<-c("FSPR30","BSPR30","YSPR30")

      return(invisible(SmryStats))
      }
   else
      {
      # time series
      rps<-brps[brps[,"steepness"]==run[1,"SR1"] & brps[,"popln"]==iFish,c("MSY","Bmsy","Rmsy")]

      MP<-model.frame(FLQuants(fbr.MP=quantile(fbar(        Stock),na.rm=T),
                               ssb.MP=quantile(ssb(         Stock),na.rm=T)/rps["Bmsy"],
                               rec.MP=quantile(rec(         Stock),na.rm=T)/rps["Rmsy"],
                               yld.MP=quantile(computeCatch(Stock),na.rm=T)/rps["MSY"]))[,c(2,6:10)]

      yrs<-dimnames(m(Stock))$year
      OM <-model.frame(FLQuants(fbr.OM=quantile(fbar(          Pop),na.rm=T)[,yrs],
                                ssb.OM=quantile(ssb(           Pop),na.rm=T)[,yrs]/rps["Bmsy"],
                                rec.OM=quantile(rec(           Pop),na.rm=T)[,yrs]/rps["Rmsy"],
                                yld.OM=quantile(computeCatch(  Pop),na.rm=T)[,yrs]/rps["MSY"]))[,c(7:10)]
    
      DBTimeSeries<-cbind(scenario=scenario,popln=iFish,diffuse=diffuse,
                          F1 =run[1,"F1"],  F2=run[1, "F2"], F3=run[1, "F3"], F4=run[1, "F4"],
                          SR1=run[1,"SR1"],SR2=run[1,"SR2"],SR3=run[1,"SR3"],SR4=run[1,"SR4"],MP,OM)

      return(invisible(DBTimeSeries))
      }
   }

RunFunc<-function(run,yrs,Pop,srs,srDeviates,sel,cpueDeviates,avail,diffussion="missing")
  {
  ## set Fs
  iF1=run[1,"F1"]
  iF2=run[1,"F2"]
  iF3=run[1,"F3"]
  iF4=run[1,"F4"]

  ## set SRs
  iSR1=srs[srs[,"popln"]==1 & run[1,"SR1"]==srs[,"steepness"],3:4]
  iSR2=srs[srs[,"popln"]==2 & run[1,"SR2"]==srs[,"steepness"],3:4]
  iSR3=srs[srs[,"popln"]==3 & run[1,"SR3"]==srs[,"steepness"],3:4]
  iSR4=srs[srs[,"popln"]==4 & run[1,"SR4"]==srs[,"steepness"],3:4]

  cat(iF1,"\t",iF2,"\t",iF3,"\t",iF4,"\t",run[1,"SR1"],"\t",run[1,"SR2"],"\t",run[1,"SR3"],"\t",run[1,"SR4"],"\n")

  ## Make stock projections on a stock basis
  shift =unlist(run[1,c("Shift","yrShift","fShift","srShift.a","srShift.b")])

  if (missing(diffussion))
     {
     Pop[[1]]<-setPop(1,iF1,iF2,iF3,iF4,sel,avail,Pop[[1]],iSR1,srDeviates,yrs,if (shift["Shift"]==1) shift)
     Pop[[2]]<-setPop(2,iF1,iF2,iF3,iF4,sel,avail,Pop[[2]],iSR2,srDeviates,yrs,if (shift["Shift"]==2) shift)
     Pop[[3]]<-setPop(3,iF1,iF2,iF3,iF4,sel,avail,Pop[[3]],iSR3,srDeviates,yrs,if (shift["Shift"]==3) shift)
     Pop[[4]]<-setPop(4,iF1,iF2,iF3,iF4,sel,avail,Pop[[4]],iSR4,srDeviates,yrs,if (shift["Shift"]==4) shift)
  }else
     Pop<-setPopDiffuse(Pop,run[1,c("F1","F2","F3","F4")],run[1,c("SR1","SR2","SR3","SR4")],sel,avail,srs,srDeviates,yrs,diffussion,shift)

  ## Stocks
  Stock<-lapply(Pop,trim,year=as.numeric(yrs))
  cpue <-lapply(Stock,function(x) as(x,"FLIndex"))

  xsa.ctrl<-FLXSA.control(maxit=200,tsrange=100,tspower=0)

  #### Splitters
  key=0
  for (i in 1:4){
    iF<-switch(i, "1"=iF1, "2"=iF2, "3"=iF3, "4"=iF4)

    landings.n( Stock[[i]])[,yrs] <-getCatch(  i,iF,sel,avail,Pop,yrs,unlist(run[1,c("Shift","yrShift","fShift","srShift.a","srShift.b")]))
    catch.wt(   Stock[[i]])[,yrs] <-getCatchWt(i,iF,sel,avail,Pop,yrs)
    landings.wt(Stock[[i]])[,yrs] <-getCatchWt(i,iF,sel,avail,Pop,yrs)
    stock.wt(   Stock[[i]])[,yrs] <-getStockWt(i,       avail,Pop,yrs)
    catch.n(    Stock[[i]])[,yrs] <-landings.n(Stock[[i]])[,yrs]
    catch(      Stock[[i]])       <-computeCatch(     Stock[[i]],'all')
    landings(   Stock[[i]])       <-computeLandings(  Stock[[i]],'all')

    ## Index ~ population at spwaning time
    key<-key+1
    range(cpue[[i]],c("startf","endf"))<-c(0,0.01)
    cpue[[i]]@index[,yrs]              <-stock.n(Pop[[i]])[,yrs]*cpueDeviates[,,i][,yrs]
    DBTimeSeries                       <-MP(i,run,Pop[[i]],Stock[[i]],cpue[[i]],xsa.ctrl,key,ifelse(missing(diffussion),0,1))

    range(cpue[[i]],c("startf","endf"))<-c(0.0,0.01)
    cpue[[i]]@index[,yrs]              <-getStock(i,avail,Pop)[,yrs]*cpueDeviates[,,i][,yrs]

    key<-key+1
    DBTimeSeries<-rbind(cbind(cpue=1,DBTimeSeries),
                        cbind(cpue=2,MP(i,run,Pop[[i]],Stock[[i]],cpue[[i]],xsa.ctrl,key,ifelse(missing(diffussion),0,1))))
    dbWriteTable(TSCon, "TimeSeries", DBTimeSeries, append=T)
    }

  #### Lumpers
  ## 1&2
  stk12  <-Stock[[1]]+Stock[[2]]
  pop12  <-Pop[[  1]]+Pop[[  2]]
  harvest(pop12)<-calcF(m(pop12),catch.n(pop12),stock.n(pop12))
  cpue <-as(pop12,"FLIndex")
  range(cpue,c("startf","endf"))<-c(0,0.01)
  key<-key+1
  DBTimeSeries<-MP(12,run,pop12,stk12,cpue,xsa.ctrl,key,ifelse(missing(diffussion),0,1))

  cpue <-as(stk12,"FLIndex")
  range(cpue,c("startf","endf"))<-c(0,0.01)
  key<-key+1
  DBTimeSeries<-rbind(cbind(cpue=1,DBTimeSeries),
                  cbind(cpue=2,MP(12,run,pop12,stk12,cpue,xsa.ctrl,key,ifelse(missing(diffussion),0,1))))
  dbWriteTable(TSCon, "TimeSeries", DBTimeSeries, append=T)
                        
  ## 1,2&3
  stk123<-Stock[[1]]+Stock[[2]]+Stock[[3]]
  pop123<-Pop[[  1]]+Pop[[  2]]+Pop[[  3]]
  harvest(pop123)<-calcF(m(pop123),catch.n(pop123),stock.n(pop123))
  cpue <-as(pop123,"FLIndex")
  range(cpue,c("startf","endf"))<-c(0,0.01)
  key<-key+1
  DBTimeSeries<-MP(123,run,pop123,stk123,cpue,xsa.ctrl,key,ifelse(missing(diffussion),0,1))
  cpue <-as(stk123,"FLIndex")
  range(cpue,c("startf","endf"))<-c(0,0.01)
 
  key<-key+1
  DBTimeSeries<-rbind(cbind(cpue=1,DBTimeSeries),
                  cbind(cpue=2,MP(123,run,pop123,stk123,cpue,xsa.ctrl,key,ifelse(missing(diffussion),0,1))))
  dbWriteTable(TSCon, "TimeSeries", DBTimeSeries, append=T)
  }