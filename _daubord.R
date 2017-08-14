#11/8/17 getting waterways survey data into R
#based on similar Genstat program
#program order changed to get variables in sensible order
#datafile similar to that produced from Genstat in 2016 except:
# 1. time now called duration
# 2. sitename variable, as well as site
# Imputed counts not yet done - not sure needed now main emphasis is on binomial
# LCM data also outstanding
timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

library(xlsx)
library(RODBC)
library(lubridate)  #used for date functions - base versions limited
library(StreamMetabolism)  # for sunset times
library(rgdal)   #OS grid refs in function below
source("C:/data/dog/rgamcode/gridref_code.R")

###########################################################################
#set up database connection
###########################################################################
#need to do this in two stages, as otherwise get fault
db<-file.path("C:/databases/batlinks_newdb.accdb") #connect database.
channel<-odbcConnectAccess2007(db) #'channel' to identify database 

###########################################################################
#First get survey info inc site data, and observer info
###########################################################################
wwc=sqlQuery( channel ,'SELECT * FROM Query_WaterwayCount order by ID')

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
                WaterwayTypeID,Name from Query_WaterwaySite order by Code')
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
data.frame(table(daub$country,daub$region))

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
coord=gr2xy(gridref2,ireland=daub$region=='Northern Ireland',date = daub$date)
daub$east=coord$easting/1000
daub$north=coord$northing/1000
daub$lat=coord$latitude
daub$long=coord$longitude
plot(daub$lat~daub$long)

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
daub$settime=coord$set  #as prop of day
daub$duration=as.integer(wwc$Duration)
chktime=((daub$end-daub$start)+(daub$start>daub$end))*24*60
summary(round(chktime)==daub$duration)  #should be true or NA
if (sum(round(chktime)!=daub$duration,na.rm=TRUE)>0) {
  print(daub[round(chktime)!=daub$duration&!is.na(daub$duration),
             c('site','date','start','end','duration')])
}
print(daub[daub$duration>180&!(is.na(daub$duration)),
          c("site","date","start","end","duration")])
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

daub$clear=wwc$SpotsClearView

daub$width=wwc$WaterWidth

###########################################################################
# Land classes - data was in old access database which can't easily be read
# into R direct.  Use genstat file C:\data\nbmp2016\nsp\_lcord.gen to
# convert to xlsx
###########################################################################
# first read mapping of grid square to land classes etc
#lc1po=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheetName = 'lc1po',colIndex = 1:4)
# above fails as too large so use csv
lc1po=read.csv('C:/databases/nbmp_new/lcdata.csv')
lc1po$NE=lc1po$NORTH*1000+lc1po$EAST
lc1po=lc1po[order(lc1po$NE),]
ne=factor(floor(daub$north)*1000+floor(daub$east),levels=lc1po$NE)
ne[daub$country=="Northern Ireland"]=NA

lpo=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheetName = 'lpo')
daub$landclass=factor(lc1po$OldLCCode[ne],levels=lpo$Land_class,labels=as.character(lpo$Description))
#table(daub$landclass)

hspo=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheetName = 'hspo')
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
detpo=read.xlsx('C:/databases/nbmp_new/Detector details and categories.xlsx',sheetIndex = 1)
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

daub$observer=wwc$ObserverID

###########################################################################
# spot level counts
###########################################################################
spotpo=sqlQuery(channel,'SELECT * FROM Query_WaterwaySpotCount order by WaterwayCountID')
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

# replace 0s with NA where no spots
daubspp=c("Count","Countuns","Countall","nposure","nposall")
daub[daub$nspot==0,daubspp]=NA

#2/2/16 check where status says incomplete but all spots present, - set to missing
error=daub$nspot==10 & daub$status=='Incomplete route or count'
daub[error,c("site","date","nspot","nposall","Count")]
daub[error,daubspp]=NA

#call file daub.RData to distinguish from daub.rda from genstat
summary(daub[,daubspp])
print(nrow(daub))
save(daub,file="daub.RData")

timestamp()