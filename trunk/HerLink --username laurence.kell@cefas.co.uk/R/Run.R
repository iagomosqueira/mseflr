memory.limit(size=4000)

  ## Init
i   <-1

t.<-c(0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,1,1.0,1.0,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,2,1.0,1.0,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,1,0.8,1.88,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,1,0.6,0.60,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,2,2025,1,1.0,1.0,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,2,2025,2,1.0,1.0,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,2,2025,1,0.8,1.88,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,2,2025,1,0.6,0.60,
      0.25,0.5, 0.5, 0.5, 0.75,0.75,0.75,0.75,3,2025,2,1.0,1.0,
      0.25,0.5, 0.5, 0.5, 0.75,0.75,0.75,0.75,4,2025,2,1.0,1.0,
      0.75,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,1,1.0,1.0,
      0.25,0.5, 0.5, 0.5, 0.75,0.75,0.75,0.75,1,2025,1,1.0,1.0,
      0.25,0.5, 0.25,0.5, 0.75,0.75,0.75,0.75,1,2025,2,1.0,1.0,
      0.25,0.5, 0.75,0.5, 0.75,0.75,0.75,0.75,1,2025,2,1.0,1.0,
      0.25,0.5, 0.5, 0.25,0.75,0.75,0.75,0.75,1,2025,2,1.0,1.0,
      0.25,0.5, 0.5, 0.75,0.75,0.75,0.75,0.75,1,2025,2,1.0,1.0)

runs<-data.frame(t(array(t.,dim=c(13,12),dimnames=list(c("F1","F2","F3","F4","SR1","SR2","SR3","SR4","Shift","yrShift","fShift","srShift.a","srShift.b"),NULL))))

my.dir<-"c:/papers/HerringLinking/"
#my.dir<-"z:/pc2900/Papers/HerringLinking/"
source(paste(my.dir,"R/InitRun.R",             sep=""))
load(  paste(my.dir,"data/conditioning.RData", sep=""))
dbFile<-paste(my.dir,"/db/Final2.dbf",sep="")
  
## Run #########################################################################
## get rid of old stuff
#file.remove(dbFile)
## initialise new database
SQLite(max.con = 16, fetch.default.rec = 500, force.reload = FALSE, shared.cache=FALSE)
## connect
TSCon <-dbConnect(dbDriver("SQLite"), dbname=dbFile)
#windows(120,20)

## runs
for (scenario in 1:4)
  {
  scen<-paste("run",scenario,"_",sep="")
  RunFunc(runs[scenario,],projPeriod,lapply(simPop,iter,501:1000),srs,iter(srDeviates,501:1000),lapply(sel,iter,501:1000),iter(cpueDeviates,501:1000),avail)
  scen<-paste("diff",scenario,"_",sep="")
  RunFunc(runs[scenario,],projPeriod,lapply(simPop,iter,501:1000),srs,iter(srDeviates,501:1000),lapply(sel,iter,501:1000),iter(cpueDeviates,501:1000),avail,diffussion)
  }
  
## close and clean up databases
dbListTables(TSCon)
dbDisconnect(TSCon)
file.info(dbFile)
  
