#testing fgindex
library(gam)

source("C:/data/nbmp2015/R/rcode/func_fgindex.R")
set.seed(5710)

load("C:/data/nbmp2015/r/rcode/hib_daub.rda",verbose=TRUE) 
attach(hib_daub)
#get smaller version to save time in development
#hib_daub=data.frame(hib_daub[1:1000,])
str(hib_daub)


fgdata(daub~site+year+months+tpc1)

asgenstat=fgindex(splinedf=6, nboot=40,nmonitor=5,baseyear=1999,title="as genstat",
               print=c("model","summary"))

fgplot(asgenstat, ymin=0)


