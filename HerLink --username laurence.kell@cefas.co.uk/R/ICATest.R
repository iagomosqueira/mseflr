
xsa.ctrl<-FLXSA.control(maxit=60,shk.ages=3,tspower=0,tsrange=100)
index(cpue[[1]])[,ac(1995:2024)] <- getStock(1,        avail,histPop)[,ac(1995:2024)] #*cpueDeviates[,,1]
range(cpue[[1]],c("startf","endf"))<-c(0.0,0.1)
hSn<-histStock[[1]]+FLXSA(histStock[[1]],cpue[[1]],xsa.ctrl,diag.flag=FALSE)
catch(hSn)<-computeCatch(hSn)
landings(hSn)<-computeLandings(hSn)

index(cpue[[1]])[,ac(2005:2024)] <- getCatch(1,iF1,sel,avail,histPop)[,ac(2005:2024)] #*cpueDeviates[,,1]
range(cpue[[1]],c("startf","endf"))<-c(0.0,1.0)
hSc<-histStock[[1]]+FLXSA(histStock[[1]],cpue[[1]],xsa.ctrl,diag.flag=FALSE)
catch(hSc)<-computeCatch(hSc)
landings(hSc)<-computeLandings(hSc)

catch(   histPop[[1]])<-computeCatch(   histPop[[1]])
landings(histPop[[1]])<-computeLandings(histPop[[1]])

catch(   histStock[[1]])<-computeCatch(   histStock[[1]])
landings(histStock[[1]])<-computeLandings(histStock[[1]])

plot(lapply(FLStocks(histPop[[1]],histStock[[1]],hSn,hSc),window,start=1995))


cpue     <-lapply(histPop,function(x) as(x,"FLIndex"))
xsa.ctrl<-FLXSA.control(maxit=30,shk.ages=3,tspower=0,tsrange=100,vpa=T)
index(cpue[[1]])[,ac(1970:2004)]   <-stock.n(histPop[[1]])[,ac(1970:2004)] #*cpueDeviates[,,1]
range(cpue[[1]],c("startf","endf"))<-c(0.0,0.1)
hPn<-histPop[[1]]+FLXSA(histPop[[1]],cpue[[1]],xsa.ctrl,diag.flag=FALSE)
plot(FLStocks(hPn,histPop[[1]]))

harvest(hPn)[1,]
harvest(histPop[[1]])[1,]

apply(sweep(harvest(hPn,c(2,6),fbar(hPn),"/"),c(1,6),mean)
apply(sweep(harvest(histPop[[1]]),c(2,6),fbar(histPop[[1]]),"/"),c(1,6),mean)

apply(sweep(harvest(hPn)[1,].fbar(hPn),"/"),c(1,6),mean)

stock.n(histPop[[1]])*harvest(histPop[[1]])/(harvest(histPop[[1]])+m(histPop[[1]]))*(1-exp(-harvest(-histPop[[1]])-m(histPop[[1]])))

window(stock.n(histPop[[1]])*harvest(histPop[[1]])/(harvest(histPop[[1]])+m(histPop[[1]]))*(1-exp(-harvest(histPop[[1]])-m(histPop[[1]]))),start=1995,end=2004)/catch.n(histPop[[1]])[,ac(1995:2004)]

sep.2(      control)<-as.integer(NA)
sep.gradual(control)<-TRUE
sr(         control)<-FALSE
sr.age(     control)<-as.integer(1)
lambda.age( control)<-rep(1,9)
lambda.yr(  control)<-rep(1,5)
lambda.sr(  control)<- 0.1
index.model(control)<-"l"
index.cor(  control)<-0
sep.nyr(    control)<-as.integer(5)
sep.age(    control)<-as.integer(4)
sep.sel(    control)<-as.double(1)

index.var(cpue[[1]])[]<-1
ica                   <- FLICA(iter(histPop[[1]],1),FLIndices(iter(cpue[[1]],1)), control)

hP       <-window(iter(histPop[[1]],1),end=2004)
hC       <-as(hP,"FLIndex")
index(hC)<-stock.n(hP)
range(hC,c("startf","endf"))<-c(0.0,0.1)

xsa.ctrl<-FLXSA.control(maxit=100)
hX<-hP+FLXSA(hP,hC,xsa.ctrl,diag.flag=FALSE)

plot(FLStocks(hP,hX))

cC<-function(x) stock.n(x)*harvest(x)/(harvest(x)+m(x))*(1-exp(-harvest(x)-m(x)))
cC(hX)/cC(hP)
