#4/1/21 saving spot data
#11/8/17 getting waterways survey data into R
#based on similar Genstat program
#program order changed to get variables in sensible order
#datafile similar to that produced from Genstat in 2016 except:
# 1. time now called duration
# 2. sitename variable, as well as site
# Imputed counts not yet done - not sure needed now main emphasis is on binomial
# LCM data also outstanding
# 18/12/19 codes from wwMasterSite see 17/1/19 email from philip
# 11/2/20 add info on self-selected sites see 11/2/20 email from philip
timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

#library(xlsx)  #12/6/19 failed with 64 bit R
library(openxlsx)
library(RODBC)
library(lubridate)  #used for date functions - base versions limited
library(StreamMetabolism)  # for sunset times
library(rgdal)   #OS grid refs in function below
source("C:/data/dog/rgamcode/gridref_code.R")

###########################################################################
#set up database connection
#nb if using 32 bit office need to use 32bit R for ODBc drivers to work and
#32 bit R not supported by early versions of R Studio 1.2
#see https://community.rstudio.com/t/unable-to-select-32-bit-r-version-in-rstudio-v1-2-1047-1-preview/16404
#for now (12/6/19 just run in batch using 32 bit version)
###########################################################################
#need to do this in two stages, as otherwise get fault
#alternatives for accdb or mdb - quote out one not needed
db<-file.path("C:/databases/batlinks_newdb.accdb") #connect database.
channel<-odbcConnectAccess2007(db) #'channel' to identify database 
#db<-file.path("C:/databases/batlinks_newdb.mdb") #connect database.
#channel<-odbcConnectAccess(db) #'channel' to identify database 

###########################################################################
#First get survey info inc site data, and observer info
###########################################################################
wwc=sqlQuery( channel ,'SELECT * FROM Query_WaterwayCount order by ID')
cat(paste0("\n",nrow(wwc)," rows of data imported from Query_WaterwayCount"))

###########################################################################
# start building final 'daub' dataframe used for analysis
# Dates are imported as POSIX structure, stored as number of seconds since
# 1/1/1970, but switch to simpler date structure
###########################################################################
daub=data.frame(site=factor(wwc$Code),
               # as.character needed or get errors
               date=as.Date(as.character(wwc$CountDate)))
daub$year=factor(year(daub$date))
daub$dayno=as.numeric(daub$date)-as.numeric(as.Date(ISOdate(daub$year,1,1)))+1
table(daub$year)
daub$period=factor(wwc$DatePeriod)

# crashed if select *
sitdf=sqlQuery( channel ,'select ID,Code,SiteName,GridReference,CountryID,RegionID,
             WaterwayTypeID,Name,Mcode,orig from Query_WaterwaySite order by Code')
cat(paste0("\n",nrow(sitdf)," rows of data imported from Query_WaterwaySite"))

#check for orphans
wslev=unique(wwc$WaterwaySiteID)
print(mean(sitdf$ID%in%wslev))  # prop sites with surveys
print(mean(wslev%in%sitdf$ID))  # prop survey sites with site data
if (sum(!wslev%in%sitdf$ID)) {
  cat('\n Fault: surveys present with no site details \n')
  print(wwc[!wwc$WaterwaySiteID%in%sitdf$ID,c("ID","Code","CountDate")])
  rm(wwc)  #to ensure can't continue processing
}
if (sum(!sitdf$ID%in%wslev)) {
  cat('\n Warning: sites present with no surveys \n')
  print(sitdf[!sitdf$ID%in%wslev,c("ID","Code","SiteName","GridReference")])
  sitdf=sitdf[sitdf$ID%in%wslev,]  #remove details for those not used
}
daub$sitename=factor(paste(wwc$Code,wwc$SiteName))

statpo=sqlQuery(channel,'select * from dbo_WaterwayCountStatus order by ID')
daub$status=factor(wwc$WaterwayCountStatusID,levels=statpo$ID,labels=statpo$Description)
daub$status=droplevels(daub$status)

###########################################################################
# Geography.  Do everything at site level first, then to survey level
###########################################################################
ctrypo=sqlQuery(channel, query='select * from [dbo_Country] order by ID')
ctrypo$Name=as.character(ctrypo$Name) # imported as factor
# make name shorter
ctrypo$Name[ctrypo$Name=='British Crown Dependencies']='Crown Dependency'
sitdf$Country=factor(sitdf$CountryID,levels=ctrypo$ID,labels=ctrypo$Name)
sitdf$Country=droplevels(sitdf$Country)
# now convert site level to survey level
daub$country=factor(sitdf$Country[daub$site],levels=levels(sitdf$Country))
data.frame(table(daub$country))

regpo=sqlQuery(channel, query='select * from [dbo_Region] order by ID')
# merge 14 London with 13 SE
sitdf$Region=factor(sitdf$RegionID-(sitdf$RegionID==14),levels=regpo$ID,labels=regpo$Name)
sitdf$Region=droplevels(sitdf$Region)
daub$region=factor(sitdf$Region[daub$site],levels=levels(sitdf$Region))
tabcoreg=data.frame(table(daub$country,daub$region))
print(tabcoreg[tabcoreg$Freq>0,])

#version with NE combined with Y & H
reglab2=gsub('North East','North East/Y&H',levels(daub$region))
daub$region2=factor(gsub('North East','North East/Y&H',as.character(daub$region)),
                   levels=reglab2)
daub$region2[daub$region2=='Yorkshire and The Humber']='North East/Y&H'
daub$region2=droplevels(daub$region2)
#data.frame(table(daub$region2,daub$region))

daub$waterway=factor(sitdf$Name[daub$site])
daub$wwtype=factor(sitdf$WaterwayTypeID[daub$site],levels=1:2,labels=
                     c('Stream/River','Canal'))
#table(daub$waterway,daub$wwtype)

# grid references
daub$gridref=sitdf$GridReference[daub$site]  #survey level
#take out channel islands as not on standard grid
gridref2=daub$gridref
gridref2[daub$country=='Crown Dependency']=NA
# 2/1/18 now have a few from RoI (well maybe not as grid ref in Irish sea!)
osi=daub$region=='Northern Ireland'|daub$region=='Eire'
coord=gr2xy(gridref2,ireland=osi,date = daub$date)
daub$east=coord$easting/1000
daub$north=coord$northing/1000
daub$lat=coord$latitude
daub$long=coord$longitude
plot(daub$lat~daub$long,pch=c(1:11)[daub$region2],col=c(1:11)[daub$region2])

# 11/2/20 self-selected
daub$selfselect=sitdf$orig[daub$site]==0

############################################################################
#“Code” in “dbo_WaterwayMasterSite”. These are the site codes from the EA River 
# Habitat Survey site list, so if it has a code then it’s from this list, if 
# it doesn’t then it’s a self-selected site or created by BCT staff. 
############################################################################
daub$EAcode=sitdf$Mcode[daub$site]

############################################################################
# calculate time taken and check equal to imported Duration
# initially used parse_date_time gives POSIXct format, but now instead record as prop
# of day, as in Genstat 
###########################################################################
# daub$end=parse_date_time(wwc$EndTime,"T*")/(24*60)
daub$dist=wwc$Distance
daub$start=(as.numeric(substr(wwc$StartTime,1,2))+
             as.numeric(substr(wwc$StartTime,4,5))/60)/24
daub$end=(as.numeric(substr(wwc$EndTime,1,2))+
           as.numeric(substr(wwc$EndTime,4,5))/60)/24
# switch those with end before start, unless end is after midnight
# switch contains row nos of those to switch
switch=which((daub$start>daub$end) & (daub$end>0.75))
print(daub[switch,c("site","date","start","end")])
start=daub$start[switch]
daub$start[switch]=daub$end[switch]
daub$end[switch]=start
daub$settime=coord$set  #as prop of day
#1/4/19 NEED TO RECALC DURATION FOR SWITCHED VALUES
daub$duration=as.integer(wwc$Duration)
chktime=((daub$end-daub$start)+(daub$start>daub$end))*24*60
summary(round(chktime)==daub$duration)  #should be true or NA
# get in nice printing format by creating daubpr dataframe with text times
ttime=function(time){ #takes fraction of day & converts to hrs:mins
  hrs=floor(time*24)
  mins=as.character(round((time-hrs/24)*24*60)+100)
  ttime=paste0(hrs,':',substr(mins,2,3))
}
daubpr=daub[,c("site","date","duration")]
daubpr$start=ttime(daub$start)
daubpr$end=ttime(daub$end)
if (sum(round(chktime)!=daub$duration,na.rm=TRUE)>0) {
  wh=which(round(chktime)!=daub$duration&!is.na(daub$duration))
  print(daubpr[wh,])
}
print(daubpr[daub$duration>180&!(is.na(daub$duration)),])
daub$minsafter=(daub$start-daub$settime)*24*60

###########################################################################
# Met data
###########################################################################
daub$temp=wwc$Temperature
daub$temp[daub$temp<0]=NA   #-9.9 used for mvs in some cases

rnpo=sqlQuery( channel , 'select * from dbo_Rainfall order by ID')
daub$rain=factor(wwc$RainfallID,levels=rnpo$ID,labels=rnpo$Name)
daub$rain=droplevels(daub$rain)
data.frame(table(daub$rain))

cldpo=sqlQuery( channel , 'select * from dbo_CloudCover order by ID')
daub$cloud=factor(wwc$CloudCoverID,levels=cldpo$ID,labels=cldpo$Name)
daub$cloud=droplevels(daub$cloud)
data.frame(table(daub$cloud))

wndpo=sqlQuery( channel , 'select * from dbo_WindStrength order by ID')
daub$wind=factor(wwc$WindStrengthID,levels=wndpo$ID,labels=wndpo$Name)
daub$wind=droplevels(daub$wind)
data.frame(table(daub$wind))

###########################################################################
# Waterway characteristics
###########################################################################
daub$treeshelter=factor(wwc$WaterwayTreeShelterID)
data.frame(table(daub$treeshelter))

daub$smoothwater=factor(wwc$WaterwaySmoothnessID)
data.frame(table(daub$smoothwater))
# 4/1/19 clear renamed clearview to avoid clash of names with cloudy=clear/low
daub$clearview=wwc$SpotsClearView

daub$width=wwc$WaterWidth

###########################################################################
# Land classes - data was in old access database which can't easily be read
# into R direct.  Use genstat file C:\data\nbmp2016\nsp\_lcord.gen to
# convert to xlsx
###########################################################################
# first read mapping of grid square to land classes etc
#lc1po=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheet = 'lc1po',colIndex = 1:4)
# above fails as too large so use csv
lc1po=read.csv('C:/databases/nbmp_new/lcdata.csv')
lc1po$NE=lc1po$NORTH*1000+lc1po$EAST
lc1po=lc1po[order(lc1po$NE),]
ne=factor(floor(daub$north)*1000+floor(daub$east),levels=lc1po$NE)
ne[daub$country=="Northern Ireland"]=NA

lpo=read.xlsx('C:/databases/nbmp_new/lcdata2.xlsx',sheet = 'lpo')
daub$landclass=factor(lc1po$OldLCCode[ne],levels=lpo$Land_class,labels=as.character(lpo$Description))
#table(daub$landclass)

hspo=read.xlsx('C:/databases/nbmp_new/lcdata2.xlsx',sheet = 'hspo')
daub$strata2=factor(lc1po$Zone[ne],levels=hspo$Strata_id,labels=as.character(hspo$Description))
daub$strata2=droplevels(daub$strata2)
data.frame(summary(daub$strata2))

###########################################################################
# site.year information
###########################################################################
expo=sqlQuery( channel , 'SELECT * FROM dbo_VolunteerExperience')
#4/2/16 old database had 'no data' to avoid mvs so recreate
expo=rbind(expo,data.frame(ID=4,Description='No data',SortOrder=40))
daub$experience=factor(wwc$VolunteerExperienceID,levels=expo$ID,labels=expo$Description)
daub$experience[is.na(daub$experience)]='No data'
daub$experience=droplevels(daub$experience)
summary(daub$experience)

vspo=sqlQuery( channel , 'SELECT * FROM dbo_VolunteerSkill')
vspo=rbind(vspo,data.frame(ID=4,Description='No data',SortOrder=40))
daub$idskills=factor(wwc$VolunteerSkillID,levels=vspo$ID,labels=vspo$Description)
daub$idskills[is.na(daub$idskills)]='No data'
daub$idskills=droplevels(daub$idskills)
summary(daub$idskills)

#nb xls spreadsheet used to allow list of groups
#1/11/16 switch to dataset provided by Philip, which contains categories and characteristics in
#same sheet
detpo=read.xlsx('../composite/input/Detector details and categories.xlsx',
                sheet = 1)
# remove unused detectors 
detpo=detpo[detpo$Detector_ID%in%unique(wwc$DetectorID),]
print(mean(unique(wwc$DetectorID)%in%detpo$Detector_ID))
daub$detector=factor(wwc$DetectorID,levels=detpo$Detector_ID)
levels(daub$detector)=detpo$Name  #for some reason won't work if use labels in above

#add unknowns to detector categories
levels(detpo$Microphone_category)=c(levels(detpo$Microphone_category),"Unknown")
detpo$Microphone_category[is.na(detpo$Microphone_category)]="Unknown"
levels(detpo$Peak_Sensitivity_Category)=c(levels(detpo$Peak_Sensitivity_Category),"Unknown")
detpo$Peak_Sensitivity_Category[is.na(detpo$Peak_Sensitivity_Category)]="Unknown"
daub$microphone=factor(detpo$Microphone_category[daub$detector],levels=levels(detpo$Microphone_category))
daub$sensitivity=factor(detpo$Peak_Sensitivity_Category[daub$detector],
                       levels=levels(detpo$Peak_Sensitivity_Category))
table(daub$microphone,daub$sensitivity)
#and grouping
detgrpo=read.xlsx('../composite/input/Detector details and categories.xlsx',
                sheet = 2)
daub$detgrp=factor(detpo$New.group[daub$detector],levels=detgrpo$Group.ID,
                   labels=detgrpo$Description)
daub$detgrp[is.na(daub$detgrp)]='Unknown/none' #added in case not all set
tabdd=data.frame(table(daub$detector,daub$detgrp,useNA = 'always'))
print(tabdd[tabdd$Freq>0,])

daub$observer=as.factor(wwc$ObserverID)

###########################################################################
# spot level counts
###########################################################################
spotpo=sqlQuery(channel,'SELECT * FROM Query_WaterwaySpotCount order by WaterwayCountID')
cat(paste0("\n",nrow(spotpo)," rows of data imported from Query_WaterwaySpotCount"))

table(spotpo$SpotNotSurveyed,useNA='always')
spotspp=c('DaubCount', 'DaubUnsureCount')
# in case some spots not surveyed contain 0s -  reset to NA
spotpo[spotpo$SpotNotSurveyed==1,spotspp]=NA
# check
nmiss=rowSums(is.na(spotpo[,spotspp]))
table(nmiss,useNA='always')
#look at where some but not all missing (10/8/17 none present)
print(spotpo[nmiss!=0&nmiss!=2,])

# now match to fc$ID so can get sum for each species on each survey
# first need to check have same ids
spotpo$WaterwayCountID=factor(spotpo$WaterwayCountID)
wcidlev=levels(spotpo$WaterwayCountID)

if (sum(!wcidlev%in%wwc$ID,na.rm=TRUE)>0){
  cat("\n*** Surveys present in spot count table which are not in survey dataset ***\n")
  print(spotpo[!spotpo$WaterwayCountID%in%wwc$ID,])
  rm(spotpo)
}
if (sum(!wwc$ID%in%wcidlev,na.rm=TRUE)>0){
  cat("\n*** Some surveys have no corresponding entries in spot count table ***\n")
  print(wwc[!wwc$ID%in%wcidlev,c("ID","Code","CountDate")])
  rm(spotpo)
}

#now form tables of sums of counts and add to dataframe nsp
#first check that tabulating factor matches fc$ID so that get data in correct order
if (sum(wwc$ID!=wcidlev,na.rm=TRUE)>0) {
  cat("\n*** Order of surveys doesn't match ***\n")
  print(data.frame(fc$ID,fcidlev)[fc$ID!=fcidlev,])
  rm(spotpo)
}
daub$nspot=as.integer(tapply(spotpo$SpotNotSurveyed==0,spotpo$WaterwayCountID,sum,
                             na.rm=TRUE))
daub$Count=as.integer(tapply(spotpo$DaubCount,spotpo$WaterwayCountID,sum,na.rm=TRUE))
daub$Countuns=as.integer(tapply(spotpo$DaubUnsureCount,spotpo$WaterwayCountID,sum,
                                na.rm=TRUE))
daub$Countall=daub$Count+daub$Countuns
daub$nposure=as.integer(tapply(spotpo$DaubCount>0,spotpo$WaterwayCountID,sum,na.rm=TRUE))
daub$nposall=as.integer(tapply((spotpo$DaubCount+spotpo$DaubUnsureCount)>0,
                               spotpo$WaterwayCountID,sum,na.rm=TRUE))
# versions with 48 max per spot
sum(is.na(spotpo$DaubCount)); sum(spotpo$DaubCount,na.rm=TRUE);
     sum(spotpo$DaubCount>48,na.rm=TRUE)
spotpo$DaubCount[spotpo$DaubCount>48]=48
sum(is.na(spotpo$DaubCount)); sum(spotpo$DaubCount,na.rm=TRUE)
all48=spotpo$DaubCount+spotpo$DaubUnsureCount
sum(is.na(all48)); sum(all48,na.rm=TRUE)
all48[all48>48]=48
sum(is.na(all48)); sum(all48,na.rm=TRUE)
daub$Count48=as.integer(tapply(spotpo$DaubCount,spotpo$WaterwayCountID,sum,na.rm=TRUE))
daub$Countall48=as.integer(tapply(all48,spotpo$WaterwayCountID,sum,na.rm=TRUE))

# replace 0s with NA where no spots
daubspp=c("Count","Countuns","Countall","Count48","Countall48","nposure","nposall")
daub[daub$nspot==0,daubspp]=NA

#2/2/16 check where status says incomplete but all spots present, - set to missing
error=daub$nspot==10 & daub$status=='Incomplete route or count'
daub[error,c("site","date","nspot","nposall","Count")]
daub[error,daubspp]=NA

###########################################################################
# 5/1/21 export spot data to allow trend modelling at the level for ons work
#add wcid to dataframe daub so can link
###########################################################################
daub$wcid=wwc$ID
#date needs to be in ordinary date format
spotpo$CountDate=as.Date(as.character(spotpo$CountDate))
#Easier if FCID not factor
spotpo$WaterwayCountID=as.integer(as.character(spotpo$WaterwayCountID))

# 9/1/20 remove from daub & spotpo where status implies no data (see email 
# from philip 8/1/21)
badstatus=daub$status=="Site not surveyed this year"
badcount=daub[badstatus,"wcid"]
badspot=which(spotpo$WaterwayCountID%in%badcount)
sum(is.na(spotpo[,4]))
spotpo[badspot,(4:5)]=NA
sum(is.na(spotpo[,4]))
daub[badstatus,daubspp]=NA

save(spotpo,file="daub_spot.RData",version=2)

###########################################################################
# number of previous surveys for each observer
###########################################################################
yrlev=as.numeric(levels(daub$year))
nyr=nlevels(daub$year)
nv=nrow(daub)
one=rep(1,nv)
daub$previous=0  #so zero if unknown, eg in first year
for (i in nyr:2) {
  thisyr=daub$year==yrlev[i]
  one[thisyr]=0
  prev=tapply(one,daub$observer,sum)[daub$observer]
  daub$previous[thisyr]=prev[thisyr]
}

#call file daub.RData to distinguish from daub.rda from genstat
summary(daub[,daubspp])
print(nrow(daub))
#12/6/19 version 3 currently not read by genstat
save(daub,file="daub.RData",version=2)

sessionInfo()
timestamp()