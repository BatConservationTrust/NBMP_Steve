# 19/9/16 GAM regional comparison for daubs
options(width=180)

#load functions and data
library(gam)

rm(list=ls())  #ensure memory clear
source("C:/data/nbmp2015/R/rcode/func_fgindex.R")
load("hib.rda")  #much faster than using xlsx file below
#need to restart R before re-running this line or get java error
#hib=read.xlsx("C:/data/nbmp2015/r/hib/hib2015.xlsx",sheetIndex=1)
                 
#********************************************************************************
#set all variables etc
#set below to UK values and overwrite for Wales, Scotland etc if needed
#********************************************************************************
nr=400   #n randomisations
nrchk=100 #check point where abandon if clearly not sig
baseyear=1999
splinedf=6  #16/11/13 go to 6 as default for all hibs (7 is correct for 90-13, 5 for 97 on)
firstyear=1998
minyears=2  #min number of years of valid data to include site
#set tcountries to those to be included
tcountries=c("England","Wales","Scotland")
monthstouse=c("Dec","Jan","Feb","Mar")
set.seed(157235) #so can reproduce exactly if required
nmon=50  #monitoring
tsp="daub"
#call count yvar to make easier to change species.  Use eval to reduce potential
#for error since only type species name once above
eval(parse(text=paste("hib$yvar=hib$",tsp)))
#give full model as character string (converted to full formula later once all vars defined)
charmodel="yvar~site+year+month"
#********************************************************************************
#end of definition section
#********************************************************************************

#with(hib,fullmodel=as.formula(charmodel))
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset


#********************************************************************************
#remove countries not required
#********************************************************************************
goodcountry=hib$country %in% tcountries
if (sum(goodcountry)==0) {
  rm(hib) #so can't miss warning!
  twarn=paste0("tcountry set to ",tcountry,". No data in subset")
  warning(twarn)  }
hib=hib[goodcountry,]  #subset to required country
table(hib$country,useNA="ifany")
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#build up subset for analysis
#********************************************************************************
goodyear=as.numeric(as.character(hib$year))>=firstyear
goodcount=!is.na(hib$yvar)
goodmonth=hib$month %in% monthstouse
#note na.rm is needed or gives missing sum if any count missing
sitewith=tapply(hib$yvar,hib$site,sum,na.rm=TRUE)>0
goodsite=sitewith[hib$site]
tapply(goodmonth, hib$month, sum)

hib=hib[goodyear&goodcount&goodmonth&goodsite,]
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#further reduce subset, taking just sites with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(hib$year,hib$site,tcount)
table(nyears)
hib=hib[(nyears[hib$site]>=minyears),]
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#no covariates model 
#********************************************************************************
with(hib,fgdata(yvar~site+year,distribution = "poisson"))
#any unused levels will cause fault, so drop them
hib$country=droplevels(hib$country)
reg=fgregion(print=c("summary"),region=hib$country,splinedf=splinedf,
             nmonitor =nmon,nrand=nr,ncheck=nrchk)

#********************************************************************************
#add code to compare separate curves for each region, with bootstrapping
#to highlight differences
#********************************************************************************


