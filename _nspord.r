#21/7/17 getting field survey data into R
#based on similar Genstat program
#datafile similar to that produced from Genstat in 2016 except:
# 1. time now called duration
# 2. sitename variable, as well as site
# 3. Cn now nspot and nwalk
#
# For historic reasons Genstat version had various variants on count variables. Here just have:
# 1. sumSP - sum of counts over 12 spots/walks (less if some missing).
# 2. s24SP/S48SP - sum of counts as above, but 24 threshold for spots, 48 for walks
# 3. nwSP - sum of spots/walks with species present
# SP is pip45, pip55, pipuns, pipall, noc, ser, nsuns, leis
# To make equivalent to old Genstat count variables, remove where nspot/nwalk <12

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
#First get survey data using query to produce 1 row per survey
###########################################################################
fc=sqlQuery( channel ,'select * from Query_FieldCount order by ID')

###########################################################################
# start building final 'nsp' dataframe used for analysis
# Dates are imported as POSIX structure, stored as number of seconds since
# 1/1/1970, but switch to simpler date structure
###########################################################################
nsp=data.frame(site=factor(fc$Code),
               # as.character needed or get errors
               date=as.Date(as.character(fc$CountDate)))
nsp$month=factor(month(nsp$date))
nsp$year=factor(year(nsp$date))
nsp$dayno=as.numeric(nsp$date)-as.numeric(as.Date(ISOdate(nsp$year,1,1)))+1
table(nsp$year,nsp$month)

###########################################################################
# Met data
###########################################################################
nsp$temp=fc$Temperature
nsp$temp[nsp$temp<0]=NA   #-9.9 used for mvs in some cases

rnpo=sqlQuery( channel , 'select * from dbo_Rainfall order by ID')
nsp$rain=factor(fc$RainfallID,levels=rnpo$ID,labels=rnpo$Name)
nsp$rain=droplevels(nsp$rain)
data.frame(table(nsp$rain))

cldpo=sqlQuery( channel , 'select * from dbo_CloudCover order by ID')
nsp$cloud=factor(fc$CloudCoverID,levels=cldpo$ID,labels=cldpo$Name)
nsp$cloud=droplevels(nsp$cloud)
data.frame(table(nsp$cloud))

wndpo=sqlQuery( channel , 'select * from dbo_WindStrength order by ID')
nsp$wind=factor(fc$WindStrengthID,levels=wndpo$ID,labels=wndpo$Name)
nsp$wind=droplevels(nsp$wind)
data.frame(table(nsp$wind))

###########################################################################
# site data
###########################################################################
sitpo=sqlQuery(channel,'select * from Query_FieldSite order by Code')
#check for orphans
fslev=unique(fc$FieldSiteID)
print(mean(sitpo$ID%in%fslev))  # prop sites with surveys
print(mean(fslev%in%sitpo$ID))  # prop survey sites with site data
if (sum(!fslev%in%sitpo$ID)) {
  cat('\n Fault: surveys present with no site details \n')
  print(fc[!fc$FieldSiteID%in%sitpo$ID,c("ID","Code","CountDate")])
  rm(nsp)  #to ensure can't continue processing
}
if (sum(!sitpo$ID%in%fslev)) {
  cat('\n Warning: sites present with no surveys \n')
  print(sitpo[!sitpo$ID%in%fslev,c("ID","Code","Location","GridReference")])
  sitpo=sitpo[sitpo$ID%in%fslev,]  #remove details for those not used
}
nsp$sitename=factor(paste(fc$Code,fc$SiteName))

statpo=sqlQuery(channel,'select * from dbo_FieldCountStatus order by ID')
nsp$status=factor(fc$FieldCountStatusID,levels=statpo$ID,labels=statpo$Description)
nsp$status=droplevels(nsp$status)

###########################################################################
# Geography.  Do everything at site level first, then to survey level
###########################################################################
ctrypo=sqlQuery(channel, query='select * from [dbo_Country] order by ID')
ctrypo$Name=as.character(ctrypo$Name) # make name shorter
ctrypo$Name[ctrypo$Name=='British Crown Dependencies']='Crown Dependency'
sitpo$Country=factor(sitpo$CountryID,levels=ctrypo$ID,labels=ctrypo$Name)
sitpo$Country=droplevels(sitpo$Country)
# now convert site level to survey level
nsp$country=factor(sitpo$Country[nsp$site],levels=levels(sitpo$Country))
data.frame(table(nsp$country))

ctypo=sqlQuery(channel, query='select * from [dbo_County] order by ID')
sitpo$County=factor(sitpo$CountyID,levels=ctypo$ID,labels=ctypo$Name)
sitpo$County=droplevels(sitpo$County)
nsp$county=factor(sitpo$County[nsp$site],levels=levels(sitpo$County))
#data.frame(table(nsp$county))

regpo=sqlQuery(channel, query='select * from [dbo_Region] order by ID')
# merge 14 London with 13 SE
sitpo$Region=factor(sitpo$RegionID-(sitpo$RegionID==14),levels=regpo$ID,labels=regpo$Name)
sitpo$Region=droplevels(sitpo$Region)
nsp$region=factor(sitpo$Region[nsp$site],levels=levels(sitpo$Region))
data.frame(table(nsp$country,nsp$region))

#version with NE combined with Y & H
reglab2=gsub('North East','North East/Y&H',levels(nsp$region))
nsp$region2=factor(gsub('North East','North East/Y&H',as.character(nsp$region)),
                   levels=reglab2)
nsp$region2[nsp$region2=='Yorkshire and The Humber']='North East/Y&H'
nsp$region2=droplevels(nsp$region2)
#data.frame(table(nsp$region2,nsp$region))

# grid references
nsp$gridref=sitpo$GridReference[nsp$site]  #survey level
#take out channel islands as not on standard grid
gridref2=nsp$gridref
gridref2[nsp$country=='Crown Dependency']=NA
coord=gr2xy(gridref2,ireland=nsp$region=='Northern Ireland',date = nsp$date)
nsp$east=coord$easting/1000
nsp$north=coord$northing/1000
nsp$lat=coord$latitude
nsp$long=coord$longitude
plot(nsp$lat~nsp$long)

############################################################################
# calculate time taken and check equal to imported Duration
# initially used parse_date_time gives POSIXct format, but now instead record as prop
# of day, as in Genstat 
###########################################################################
# nsp$end=parse_date_time(fc$EndTime,"T*")/(24*60)
nsp$start=(as.numeric(substr(fc$StartTime,1,2))+
             as.numeric(substr(fc$StartTime,4,5))/60)/24
nsp$end=(as.numeric(substr(fc$EndTime,1,2))+
           as.numeric(substr(fc$EndTime,4,5))/60)/24
nsp$settime=coord$set  #as prop of day
nsp$duration=as.integer(fc$Duration)
chktime=((nsp$end-nsp$start)+(nsp$start>nsp$end))*24*60
summary(round(chktime)==nsp$duration)  #should be true or NA
print(nsp[nsp$duration>180&!(is.na(nsp$duration)),
          c("site","date","start","end","duration")])
nsp$minsafter=(nsp$start-nsp$settime)*24*60

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
ne=factor(floor(nsp$north)*1000+floor(nsp$east),levels=lc1po$NE)
ne[nsp$country=="Northern Ireland"]=NA

lpo=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheetName = 'lpo')
nsp$landclass=factor(lc1po$OldLCCode[ne],levels=lpo$Land_class,labels=as.character(lpo$Description))
#table(nsp$landclass)

hspo=read.xlsx('C:/databases/nbmp_new/lcdata.xlsx',sheetName = 'hspo')
nsp$strata2=factor(lc1po$Zone[ne],levels=hspo$Strata_id,labels=as.character(hspo$Description))
data.frame(summary(nsp$strata2))

###########################################################################
# site.year information
###########################################################################
expo=sqlQuery( channel , 'SELECT * FROM dbo_VolunteerExperience')
#4/2/16 old database had 'no data' to avoid mvs so recreate
expo=rbind(expo,data.frame(ID=4,Description='No data',SortOrder=40))
nsp$experience=factor(fc$VolunteerExperienceID,levels=expo$ID,labels=expo$Description)
nsp$experience[is.na(nsp$experience)]='No data'
nsp$experience=droplevels(nsp$experience)
summary(nsp$experience)

vspo=sqlQuery( channel , 'SELECT * FROM dbo_VolunteerSkill')
vspo=rbind(vspo,data.frame(ID=4,Description='No data',SortOrder=40))
nsp$idskills=factor(fc$VolunteerSkillID,levels=vspo$ID,labels=vspo$Description)
nsp$idskills[is.na(nsp$idskills)]='No data'
nsp$idskills=droplevels(nsp$idskills)
summary(nsp$idskills)

#nb xls spreadsheet used to allow list of groups
#1/11/16 switch to dataset provided by Philip, which contains categories and characteristics in
#same sheet
detpo=read.xlsx('C:/databases/nbmp_new/Detector details and categories.xlsx',sheetIndex = 1)
# remove unused detectors 
detpo=detpo[detpo$Detector_ID%in%unique(fc$DetectorID),]
nsp$detector=factor(fc$DetectorID,levels=detpo$Detector_ID,labels=detpo$Name)

#add unknowns to detector categories
levels(detpo$Microphone_category)=c(levels(detpo$Microphone_category),"Unknown")
detpo$Microphone_category[is.na(detpo$Microphone_category)]="Unknown"
levels(detpo$Peak_Sensitivity_Category)=c(levels(detpo$Peak_Sensitivity_Category),"Unknown")
detpo$Peak_Sensitivity_Category[is.na(detpo$Peak_Sensitivity_Category)]="Unknown"
nsp$microphone=factor(detpo$Microphone_category[nsp$detector],levels=levels(detpo$Microphone_category))
nsp$sensitivity=factor(detpo$Peak_Sensitivity_Category[nsp$detector],
                       levels=levels(detpo$Peak_Sensitivity_Category))
table(nsp$microphone,nsp$sensitivity)

nsp$observer=fc$ObserverID
nsp$dist=fc$Distance
nsp$period=factor(fc$DatePeriod)

###########################################################################
# spot level counts
###########################################################################
spotpo=sqlQuery(channel,'SELECT * FROM Query_FieldSpotCount order by FCID')
table(spotpo$SpotNotSurveyed,useNA='always')
spotspp=c('CommonPipCount', 'SopranoPipCount', 'PipUnsureCount')
# some spots not surveyed contain 0s -  reset to NA
spotpo[spotpo$SpotNotSurveyed==1,spotspp]=NA
# check
nmiss=rowSums(is.na(spotpo[,spotspp]))
table(nmiss,useNA='always')
#look at where some but not all missing
print(spotpo[nmiss!=0&nmiss!=3,])
# replacing NA at 120185 in 2000 by 0 seems sensible, but stop prog if more than that to check
if (sum(nmiss!=0&nmiss!=3,na.rm=TRUE)>1) {rm(spotpo)
  print("**** check NA values in spot counts ****")}
spotpo[nmiss!=0&nmiss!=3,spotspp]=0

# now match to fc$ID so can get sum for each species on each survey
# first need to check have same ids
spotpo$FCID=factor(spotpo$FCID)
fcidlev=levels(spotpo$FCID)

if (sum(!fcidlev%in%fc$ID,na.rm=TRUE)>0){
  cat("\n*** Surveys present in spot count table which are not in survey dataset ***\n")
  print(spotpo[!spotpo$FCID%in%fc$ID,])
  rm(spotpo)
}
if (sum(!fc$ID%in%fcidlev,na.rm=TRUE)>0){
  cat("\n*** Some surveys have no corresponding entries in spot count table ***\n")
  print(fc[!fc$ID%in%fcidlev,c("ID","Code","CountDate")])
  rm(spotpo)
}

#now form tables of sums of counts and add to dataframe nsp
#first check that tabulating factor matches fc$ID so that get data in correct order
if (sum(fc$ID!=fcidlev,na.rm=TRUE)>0) {
  cat("\n*** Order of surveys doesn't match ***\n")
  print(data.frame(fc$ID,fcidlev)[fc$ID!=fcidlev,])
  rm(spotpo)
}
nsp$nspot=as.integer(tapply(spotpo$SpotNotSurveyed==0,spotpo$FCID,sum,na.rm=TRUE))
nsp$sumpip45=as.integer(tapply(spotpo$CommonPipCount,spotpo$FCID,sum,na.rm=TRUE))
nsp$sumpip55=as.integer(tapply(spotpo$SopranoPipCount,spotpo$FCID,sum,na.rm=TRUE))
nsp$sumpipuns=as.integer(tapply(spotpo$PipUnsureCount,spotpo$FCID,sum,na.rm=TRUE))
nsp$sumpipall=nsp$sumpip45+nsp$sumpip55+nsp$sumpipuns
nsp[nsp$site==120010,c("site","date","nspot","sumpip45","sumpip55","sumpipuns")]

#n spots with versions
nsp$nwpip45=as.integer(tapply(spotpo$CommonPipCount>0,spotpo$FCID,sum,na.rm=TRUE))
nsp$nwpip55=as.integer(tapply(spotpo$SopranoPipCount>0,spotpo$FCID,sum,na.rm=TRUE))
nsp$nwpipuns=as.integer(tapply(spotpo$PipUnsureCount>0,spotpo$FCID,sum,na.rm=TRUE))
nsp[nsp$site==120010,c("site","date","nspot","nwpip45","nwpip55","nwpipuns")]

# 24 max version
summary(spotpo[,spotspp])
for (s in spotspp) spotpo[spotpo[,s]>24&!is.na(spotpo[,s]),s]=24
nsp$s24pip45=as.integer(tapply(spotpo$CommonPipCount,spotpo$FCID,sum,na.rm=TRUE))
nsp$s24pip55=as.integer(tapply(spotpo$SopranoPipCount,spotpo$FCID,sum,na.rm=TRUE))
nsp$s24pipuns=as.integer(tapply(spotpo$PipUnsureCount,spotpo$FCID,sum,na.rm=TRUE))

# replace 0s with NA where no spots
pipspp=c("sumpip45","sumpip55","sumpipuns","s24pip45","s24pip55","s24pipuns",
         "nwpip45","nwpip55","nwpipuns")
nsp[nsp$nspot==0,pipspp]=NA

###########################################################################
# walk level counts
###########################################################################
walkpo=sqlQuery(channel,'SELECT * FROM Query_FieldWalkCount order by FCID')
table(walkpo$WalkNotSurveyed,useNA='always')
walkspp=c('NoctuleCount', 'SerotineCount', 'NoctuleSerotineUnsureCount')
#          'LeislerCount','LeislerUnsureCount')  ignore Leislers for now
# get country factor as species recorded different in ireland
walkpo$FCID=factor(walkpo$FCID)
wcountry=nsp$country[walkpo$FCID]
table(wcountry)
# tried inserting mvs for NI here but also need to do after summing since R thinks
# sum of 12 NAs is zero!
data.frame(tapply(is.na(walkpo$LeislerCount),wcountry,mean,na.rm=TRUE))
data.frame(tapply(is.na(walkpo$NoctuleCount),wcountry,mean,na.rm=TRUE))
data.frame(tapply(is.na(walkpo$SerotineCount),wcountry,mean,na.rm=TRUE))
walkpo$NoctuleCount[wcountry=='Northern Ireland']=NA
walkpo$SerotineCount[wcountry=='Northern Ireland']=NA
walkpo$NoctuleSerotineUnsureCount[wcountry=='Northern Ireland']=NA
# some walks not surveyed contain 0s -  reset to NA
walkpo[walkpo$WalkNotSurveyed==1,walkspp]=NA
# check
nmiss=rowSums(is.na(walkpo[,walkspp]))
table(nmiss[wcountry!='Northern Ireland'],useNA='always')
#look at where some but not all missing
print(walkpo[nmiss!=0&nmiss!=3,c(1:8)])
# replacing NA at 120303 in 2006 by 0 seems sensible, but stop prog if more than that to check
if (sum(nmiss!=0&nmiss!=3,na.rm=TRUE)>1) {rm(walkpo)
  print("**** check NA values in walk counts ****")}
walkpo[nmiss!=0&nmiss!=3,walkspp]=0

# now match to fc$ID so can get sum for each species on each survey
# first need to check have same ids
fcidlev=levels(walkpo$FCID)

if (sum(!fcidlev%in%fc$ID,na.rm=TRUE)>0){
  cat("\n*** Surveys present in walk count table which are not in survey dataset ***\n")
  print(walkpo[!walkpo$FCID%in%fc$ID,])
  rm(walkpo)
}
if (sum(!fc$ID%in%fcidlev,na.rm=TRUE)>0){
  cat("\n*** Some surveys have no corresponding entries in walk count table ***\n")
  print(fc[!fc$ID%in%fcidlev,c("ID","Code","CountDate")])
  rm(walkpo)
}

#now form tables of sums of counts and add to dataframe nsp
#first check that tabulating factor matches fc$ID so that get data in correct order
if (sum(fc$ID!=fcidlev,na.rm=TRUE)>0) {
  cat("\n*** Order of surveys doesn't match ***\n")
  print(data.frame(fc$ID,fcidlev)[fc$ID!=fcidlev,])
  rm(walkpo)
}
nsp$nwalk=as.integer(tapply(walkpo$WalkNotSurveyed==0,walkpo$FCID,sum,na.rm=TRUE))
nsp$sumnoc=as.integer(tapply(walkpo$NoctuleCount,walkpo$FCID,sum,na.rm=TRUE))
nsp$sumser=as.integer(tapply(walkpo$SerotineCount,walkpo$FCID,sum,na.rm=TRUE))
nsp$sumnsuns=as.integer(tapply(walkpo$NoctuleSerotineUnsureCount,walkpo$FCID,sum,na.rm=TRUE))
nsp$sumnsall=nsp$sumnoc+nsp$sumser+nsp$sumnsuns

nsp$nwnoc=as.integer(tapply(walkpo$NoctuleCount>0,walkpo$FCID,sum,na.rm=TRUE))
nsp$nwser=as.integer(tapply(walkpo$SerotineCount>0,walkpo$FCID,sum,na.rm=TRUE))
nsp$nwnsuns=as.integer(tapply(walkpo$NoctuleSerotineUnsureCount>0,walkpo$FCID,sum,na.rm=TRUE))

# 48 max version
summary(walkpo[,walkspp])
for (s in walkspp) walkpo[walkpo[,s]>48&!is.na(walkpo[,s]),s]=48
nsp$s48noc=as.integer(tapply(walkpo$NoctuleCount,walkpo$FCID,sum,na.rm=TRUE))
nsp$s48ser=as.integer(tapply(walkpo$SerotineCount,walkpo$FCID,sum,na.rm=TRUE))
nsp$s48nsuns=as.integer(tapply(walkpo$NoctuleSerotineUnsureCount,walkpo$FCID,sum,na.rm=TRUE))

#take out NI for noctule & serotine
nsspp=c("sumnoc","sumser","sumnsuns","s48noc","s48ser","s48nsuns",
        "nwnoc","nwser","nwnsuns")
nsp[nsp$country=='Northern Ireland',nsspp]=NA
# and where no walks
nsp[nsp$nwalk==0,nsspp]=NA

#call file nsp.RData to distinguish from nsp.rda from genstat
save(nsp,file="nsp.RData")

timestamp()
