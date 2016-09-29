# 7/8/16 GAM for daubs
options(width=180)

#load functions and data
library(gam)
library(xlsx)
rm(list=ls())  #ensure memory clear
source("C:/data/nbmp2015/R/rcode/func_fgindex.R")
load("hib.rda")  #much faster than using xlsx file below
#need to restart R before re-running this line or get java error
#hib=read.xlsx("C:/data/nbmp2015/r/hib/hib2015.xlsx",sheetIndex=1)
                 
#********************************************************************************
#set all variables etc
#set below to UK values and overwrite for Wales, Scotland etc if needed
#********************************************************************************
nboot1=10
nboot3=400   #n bootstraps for analysis used for trend
baseyear=1999
splinedf=6  #16/11/13 go to 6 as default for all hibs (7 is correct for 90-13, 5 for 97 on)
firstyear=1990
minyears=2  #min number of years of valid data to include site
tcountry='GB' #set to UK, England, Wales, Scotland
monthstouse=c("Dec","Jan","Feb","Mar")
set.seed(157235) #so can reproduce exactly if required
set.seed(45843) #crashes with value above
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

#filename for xlsx output
xlsfile=paste0("C:/data/nbmp2015/R/composite/xlsx/",substring(tcountry,1,3),
               "hib_",tsp,".xlsx")
xlsfile
#with(hib,fullmodel=as.formula(charmodel))
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#country information, then take appropriate subset
#might be worth putting all this in spreadsheet in future so don't need individual
#command file for each species
#********************************************************************************
if (tcountry=='Scotland'){
  firstyear=2005 #can go back to 1998 but v few so s.e.s unlikely to be sensible
  baseyear=2011}
if (tcountry=='Wales'){
  firstyear=1990}
if (tcountry=='England') {
 firstyear=1993}
#set tcountries variable - same for wales, scot, england, but list for GB, UK
tcountries=tcountry
nboot2=nboot3  #analysis 2 (covariates, no weights) used for Eng, Wales & Scotland
if (tcountry=='UK'){
  tcountries=c("England","Wales","Scotland","NI")
  nboot2=nboot1}
if (tcountry=='GB'){
  tcountries=c("England","Wales","Scotland")
  nboot2=nboot1}
goodcountry=hib$country %in% tcountries
#following ensures consistent with Genstat which included missing countries in GB
if (tcountry=="GB") goodcountry=goodcountry | is.na(hib$country)
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
analysis1=fgindex(splinedf=splinedf,nboot=nboot1,nmonitor=nmon,baseyear=baseyear,plot=TRUE,
        title="No covariates")
write.xlsx(analysis1,file=xlsfile,sheet='Analysis 1',row.names = FALSE)

#********************************************************************************
#With covariates model 
#********************************************************************************
with(hib,fgdata(as.formula(charmodel),distribution = "poisson"))
analysis2=fgindex(splinedf=splinedf,nboot=nboot2,nmonitor=nmon,baseyear=baseyear,plot=TRUE,
                  title="With covariates")
write.xlsx(analysis2,file=xlsfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)

#********************************************************************************
#Calculate waits and fit model With covariates & weighted 
#********************************************************************************
areas=read.xlsx("C:/data/nbmp2013/katepaper/areas for weighting.xlsx",sheetIndex=1,
                rowIndex = 7:12)
colabs=levels(hib$country)
areas=areas[areas[,1]%in% colabs,]
areas
#check order the same
if (mean(areas[1]==colabs)<1) rm(areas)
totobs=table(hib$country)
wt=areas$sqkmex/as.numeric(totobs)
data.frame(totobs,wt,areas$sqkmex)  #data.frame prints nicely in columns
wts=wt[hib$country] #expand to one value per survey

analysis3=fgindex(splinedf=splinedf,nboot=nboot3,nmonitor=nmon,baseyear=baseyear,plot=TRUE,
                  title="Weighted covariates",weight=wts)
write.xlsx(analysis3,file=xlsfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)

