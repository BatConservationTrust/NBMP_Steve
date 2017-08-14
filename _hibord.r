#6/7/17 getting hibernation data into R
#based on similar Genstat program
#datafile similar to that produced from Genstat in 2016 except:
# 1. completion now renamed duration
# 2. groupings of temperatures etc not included

timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

library(xlsx)
library(RODBC)
library(lubridate)  #used for day numbers
library(rgdal)   #OS grid refs in function below
source("C:/data/dog/rgamcode/gridref_code.R")

###########################################################################
#set maximum year to exclude early records from following year
###########################################################################
maxyear=2016

###########################################################################
#set up database connection
###########################################################################
#need to do this in two stages, as otherwise get fault
db<-file.path("C:/databases/batlinks_newdb.accdb") #connect database.
channel<-odbcConnectAccess2007(db) #'channel' to identify database 

###########################################################################
#First get roost and survey data using query to produce 1 row per survey
###########################################################################
tsql='select * from query_hib_year order by Code,ConvertedDate'
hibyear <- sqlQuery( channel , tsql)

# check for missing IDs & remove
# these occur if entry in dbo_HibernationCountYear but not dbo_HibernationCount
# To investigate, paste list into _check_orphans.gen"
if (sum(is.na(hibyear$ID))>0) {
  print(hibyear[is.na(hibyear$ID),c("Code","CountYear")])
  hibyear=hibyear[!is.na(hibyear$ID),]
}

# and years beyond current processing year
# some differences between calculated year and CountYear, so use calculated, plus
# summer months handled differently
vy=year(hibyear$ConvertedDate)+(month(hibyear$ConvertedDate)>=6)
diffyr=abs(hibyear$CountYear-vy)>0
data.frame(tapply(diffyr,month(hibyear$ConvertedDate),sum))
if (sum(vy>maxyear)>0) {
  print(hibyear[vy>maxyear,c("Code","CountYear","ConvertedDate")])
  hibyear=hibyear[!vy>maxyear,]
}

# 7/7/17 only one Irish record so remove
table(hibyear$Projection)
hibyear=hibyear[!hibyear$Projection=="OSI",]

###########################################################################
# now start building final 'hib' dataframe used for analysis
# Dates are imported as POSIX structure, stored as number of seconds since
# 1/1/1970, but switch to simpler date structure
###########################################################################
# get correct site names - need separate version for text and number since R
# doesn't have levels/labels distinction
hibyear$site=factor(hibyear$Code)
hibyear$sitename=factor(paste(hibyear$Code,hibyear$SiteName))
hib=data.frame(site=factor(hibyear$Code),
               sitename=factor(paste(hibyear$Code,hibyear$SiteName)),
# as.character needed or get errors
               date=as.Date(as.character(hibyear$ConvertedDate)),
#               date=hibyear$ConvertedDate,  #POSIX
               period=factor(hibyear$DatePeriod),
               nsurv=hibyear$NumObservers,nlic=hibyear$NumLicensed)

# sort temperatures
hib$exttemp=hibyear$ExternalTemperature
hib$tempwarm=hibyear$InternalTemperatureWarmest
hib$tempcool=hibyear$InternalTemperatureCoolest
#9/1/16 switch so coolest really is coolest, and switch to appropriate names"
switch=which(hib$tempwarm<hib$tempcool)
summary(hib$tempwarm<hib$tempcool)
hib$tempwarm[switch]=hibyear$InternalTemperatureCoolest[switch]
hib$tempcool[switch]=hibyear$InternalTemperatureWarmest[switch]

#count status
tsql='select * from dbo_HibernationCountStatus'
statdf <- sqlQuery( channel , tsql)
#NB levels is format in hibyear df, labels is how should be labelled
hib$status=factor(hibyear$HibernationCountStatusID,levels=as.character(statdf$ID),
                    labels=statdf$Description)
data.frame(summary(hib$status))

# 7/7/17 change 'completion' to duration
hib$duration=hibyear$Duration
#hib$dayno=lubridate::yday(hib$date) #see below - calc from 1st Sept
hib$month=factor(month(hib$date))
hib$year=factor(year(hib$date)+(hib$month>=6))
hib$dayno=as.numeric(hib$date)-as.numeric(as.Date(ISOdate(hib$year-1,9,1)))
table(hib$year,hib$month)
hib$observer=hibyear$ObserverID

###########################################################################
# keep only species with at least 100 obs
# Barbastele added 2008 - not many sites but worth looking at
# Others with high nos tend to be unidentified etc
# 2010 try all pips 
###########################################################################
tsql='select * from dbo_Species order by ID'
spdf <- sqlQuery( channel , tsql)
spdf

tsql='select * from dbo_HibernationCountSpecies order by ID'
hcdf <- sqlQuery( channel , tsql)
summary(hcdf$HibernationCountID %in% hibyear$ID)
hcdf=hcdf[hcdf$HibernationCountID %in% hibyear$ID,]
data.frame(table(hcdf$SpeciesID))
# recode all whiskered/brandts to 4
hcdf$SpeciesID[hcdf$SpeciesID %in% c(6,8)]=4
# and all pips to 12
hcdf$SpeciesID[hcdf$SpeciesID %in% c(13,14,32)]=12
goodspp=c("whiskbrandt","daub","natt","blongear","gthorse","lshorse","barb","allpip")
goodsplev=c(4,7,9,15,16,17,2,12)
#get rid of other species
hcdf$HibernationCountID=factor(hcdf$HibernationCountID,levels=hibyear$ID)
totbats=tapply(hcdf$TotalCount,hcdf$HibernationCountID,sum)
hcdf=hcdf[hcdf$SpeciesID %in% goodsplev,]
data.frame(table(hcdf$SpeciesID))
hcdf$SpeciesID=factor(hcdf$SpeciesID,levels=goodsplev,labels=goodspp)
sptot=tapply(hcdf$TotalCount,list(hcdf$HibernationCountID,hcdf$SpeciesID),sum)
# replace NAs with 0s, since if not listed not present
sptot[is.na(sptot)]=0
totbats[is.na(totbats)]=0

# add to main dataframe
hib=cbind(hib,sptot,totbats)
goodntot=c(goodspp,"totbats")
print(data.frame(colSums(hib[,goodntot])))  #check col totals
print(sum(hib[,goodspp]))  #should be slightly less than above

###########################################################################
#converting gridref to east, north
# initially tried osg_parse from rnrfa, but v slow (code does one at time) &
# doesn't cope with NAs etc.  So use rgdal writing function gr2xy
###########################################################################
geog=gr2xy(hibyear$GridReference)
hib=cbind(hib,geog)

###########################################################################
#type etc
###########################################################################
tsql3='select * from dbo_StructureType order by ID'
tpo <- sqlQuery( channel , tsql3)
hib$type=factor(hibyear$StructureTypeID,levels=tpo$ID,labels=tpo$Description)
hib$type=droplevels(hib$type)
data.frame(table(hib$type))

tsql4='select * from dbo_HibernationSiteCrevices order by ID'
tcpo <- sqlQuery( channel , tsql4)
hib$complexity=factor(hibyear$HibernationSiteCrevicesID,levels=tcpo$ID,
                      labels=tcpo$Description)
hib$complexity=droplevels(hib$complexity)
data.frame(table(hib$complexity))

tsql5='select * from dbo_HibernationSiteDisturbance order by ID'
tdpo <- sqlQuery( channel , tsql5)
hib$disturbance=factor(hibyear$HibernationSiteDisturbanceID,levels=tdpo$ID,
                       labels=tdpo$Description)
hib$disturbance=droplevels(hib$disturbance)
data.frame(table(hib$disturbance))

tsql6='select * from dbo_HibernationSiteSize order by ID'
tspo <- sqlQuery( channel , tsql6)
hib$size=factor(hibyear$HibernationSiteSizeID,levels=tspo$ID,
                labels=tspo$Description)
hib$size=droplevels(hib$size)
data.frame(table(hib$size))

tsql6='select * from dbo_HibernationSiteSurveyability order by ID'
tsvpo <- sqlQuery( channel , tsql6)
hib$surveyab=factor(hibyear$HibernationSiteSurveyabilityID,levels=tsvpo$ID,
                    labels=tsvpo$Description)
hib$surveyab=droplevels(hib$surveyab)
data.frame(table(hib$surveyab))

tsql7='select * from dbo_County order by ID'
copo <- sqlQuery( channel , tsql7)
hib$county=factor(hibyear$CountyID,levels=copo$ID,labels=copo$Name)
hib$county=droplevels(hib$county)
data.frame(table(hib$county))

tsql8='select * from dbo_Region order by ID'
regpo <- sqlQuery( channel , tsql8)
hib$region=factor(copo$RegionID[hibyear$CountyID],levels=regpo$ID,labels=regpo$Name)
#combine London with SE
hib$region[hib$region=='London']='South East'
hib$region=droplevels(hib$region)
data.frame(table(hib$region))

#version with NE combined with Y & H
reglab2=gsub('North East','North East/Y&H',levels(hib$region))
hib$region2=factor(gsub('North East','North East/Y&H',as.character(hib$region)),
                   levels=reglab2)
hib$region2[hib$region2=='Yorkshire and The Humber']='North East/Y&H'
hib$region2=droplevels(hib$region2)
data.frame(table(hib$region2))

tsql9='select * from dbo_Country order by ID'
ctrypo <- sqlQuery( channel , tsql9)
#get rid of unused regions
regpo=regpo[regpo$Name%in%levels(hib$region),]
hib$country=factor(regpo$CountryID[hib$region],levels=ctrypo$ID,labels=ctrypo$Name)
hib$country=droplevels(hib$country)
table(hib$region2,hib$country)

###########################################################################
# check status of codes, ensuring zero where implies no bats and missing
# where appropriate
# actually looks correct anyway, but this guards against future errors.
###########################################################################
nsurvey=table(hib$status)
Totbats=tapply(hib$totbats,hib$status,sum)
data.frame(nsurvey,Totbats)
#initially set to missing unless valid count
hib[hib$status!='Bats present and counts carried out',goodntot]=NA
#then add zeros where none present
nobats=c('No bats recorded on this survey date',
         'Site destroyed','Site changed, now unsuitable for bats')
hib[hib$status%in%nobats,goodntot]=0
nsurvey=table(hib$status)
Totbats=tapply(hib$totbats,hib$status,sum)
data.frame(nsurvey,Totbats)

###########################################################################
# remove tom mcowat sites
###########################################################################
mcowat=read.xlsx("c:/data/nbmp2016/hib/Sites excluded from analysis 2016.xlsx",
                 sheetIndex=1)
hib=hib[!hib$site%in%mcowat[,1],]

###########################################################################
# weighting by GB lowland area
# nb individual species weights are formed in gam program but this calculates
# an all-species weight hib$allweight, which is used if species weights too
# extreme.
###########################################################################
ipo=read.xlsx("C:/databases/nbmp_new/areas for weighting.xlsx",sheetIndex=1,
          startRow=6)
#just take eng, wales and scotland
ipo=ipo[ipo[[1]]%in%levels(hib$country),]
totobs=as.numeric(table(hib$country))  #number surveys in each country
wt=ipo$sqkmex/totobs #wt for each country
weight=wt[hib$country]  #one row for each survey
hib$allweight=weight/mean(weight,na.rm=TRUE)
summary(hib$allweight)
tapply(hib$allweight,hib$country,mean)

#call file hib.RData to distinguish from hib.rda from genstat
save(hib,file="hib.RData")

timestamp()