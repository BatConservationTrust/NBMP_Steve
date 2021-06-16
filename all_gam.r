# 25/1/17 GAM for all species for roost counts using settings spreadsheet
timestamp()
R.Version()$version.string
if(interactive()) options(width=80) else options(width=120) 
rm(list=ls())  #ensure memory clear

#load functions and data
library(gam)
predict.gam=predict.Gam #27/12/19 name of predict function has changed
#library(xlsx)  #12/6/19 failed with 64 bit R 
library(openxlsx)
library(fastglm)

source("C:/data/dog/rgamcode/func_fgindex.R")
#source("C:/data/dog/rgamcode/func_fgfastindex.R")
load("roost.rdata")  #much faster than using xlsx file below
#need to restart R before re-running this line or may get java error
settings=read.xlsx("../composite/input/gamsettings.xlsx",sheet = "roost")
                 
#********************************************************************************
#set key variables below to give correct row of settings dataframe
#Note that 'test' changes various options to allow quick and easy checking of
#program, and plots graph interactively
#textra="standard" in most cases, but allows versions with diff base year etc.
#********************************************************************************
tcountry='Scotland' #set to GB, UK, England, Wales, Scotland
tsp="blongear"  #abbreviated name eg pip45 used to name files, col C of gamsettings
textra="standard"
test=NA  # 0 to run normally, 1 for reduced bootstraps & interactive graphs, NA automatic
if (is.na(test)){  #if NA test=1 if interactive, test=0 if batch
  if(interactive()) test=1 else test=0 
}

row=which(settings$species==tsp&settings$country==tcountry&settings$extra==textra)
if (length(row)!=1) {
  print(row)
  rm(rc)
  stop("Fault: exactly one row must match variables above.")
}
settings=settings[row,]
#********************************************************************************
#other settings. May need to change seed for some species/countries
#19/1/19 fault in _roostord.r (& .gen) meant invalid stop reasons not picked
#up from 2013.  Now corrected but little or no difference to trends due to
#taking maximum
#jan 2020 eps is epsilon (convergence criterion) for gam fits of bootstrapped 
#data.  Default used by package gam is unnecessarily small. Optimal value seems
#to be around 1e-06 based on field survey work, but can go bigger if needed
#
#fast is fastglm method, -1 means don't use it, 3 generally best, NA means auto
#- method 3 unless covfile set
#********************************************************************************
fast=NA    
eps=1e-06
set.seed(157235) #so can reproduce exactly if required
#set.seed(65733) #18/1/21 need to try various for natterers for gb wted
minyears=2
nmon=50  #monitoring
xstop=c('Weather deteriorating','survey aborted')
if (test){
  settings$nb0=10;settings$nb1=10;settings$nb2=10;settings$nb3=10;nmon=1
  if (file.exists("del.xlsx")) file.remove("del.xlsx")
}
if (is.na(fast)) fast=-1+4*(is.na(settings$covfile))

#********************************************************************************
#end of definition section
#********************************************************************************
settings
#couldn't get this one to converge so use special prog natt_gam.r
if (tcountry=='England'&tsp=="natterers"){rm(settings)}

#filename for xlsx output, use del.xlsx for testing
#xlsfile is for output
#bootfile is for bootstrap estimates (used for annual comparisons)
xlsfile=paste0("../composite/r_xlsx/",settings$filename,".xlsx")
bootfile=paste0("../composite/r_boot/",settings$filename,".xlsx")
if (test) {
  xlsfile="del.xlsx";bootfile="del2.xlsx"
}else {
  grfile=paste0("all_gam.",tsp,".",substr(tcountry,1,3),".png")
  grrow=2+(settings$nb3>0)
  png(grfile,height=grrow*300,width=600)
  par(mfrow=c(grrow,2))  #2 x 2 matrix
}
xlsfile
bootfile
wbx=createWorkbook()
for (tsh in c('Analysis 0','Analysis 1','Analysis 2','Analysis 3','slow 2')) addWorksheet(wbx,tsh)
wbb=createWorkbook()
for (tsh in c('Analysis 0','Analysis 1','Analysis 2','Analysis 3')) addWorksheet(wbb,tsh)

length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#country information, then take appropriate subset
#********************************************************************************
#set tcountries variable - same for wales, scot, england, but list for GB, UK
tcountries=tcountry
if (tcountry=='UK') tcountries=c("England","Wales","Scotland","NI")
if (tcountry=='GB') tcountries=c("England","Wales","Scotland")
if (tcountry=='EW') tcountries=c("England","Wales")

goodcountry=rc$country %in% tcountries
#following ensures consistent with Genstat which included missing countries in GB
if (tcountry=="GB") goodcountry=goodcountry | is.na(rc$country)
if (sum(goodcountry)==0) {
  rm(rc) #so can't miss warning!
  twarn=paste0("tcountry set to ",tcountry,". No data in subset")
  warning(twarn)  }
rc=rc[goodcountry,]  #subset to required country
table(rc$country,useNA="ifany")
length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#roost data differently organised to other surveys, with species in different
#rows of batcount variable, rather than different columns.
#therefore first get rows for species of interest
#********************************************************************************
table(rc$species,useNA = "always")
rc=rc[rc$species==as.character(settings$spname) & (is.na(rc$species)==0),]
rc=droplevels(rc)
table(rc$species,rc$country,useNA = "always")
length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#build up subset for analysis
#genstat dataset only has dates in correct range, but include code for future,
#remembering that some records have unknown date - R behaves very oddly if NA
#present in subset
#********************************************************************************
goodyear=as.numeric(as.character(rc$year))>=settings$firstyear
goodcount=!is.na(rc$batcount)
#nb days currently subsetted out in _roostord.r so this not really needed
gooddays=(rc$dayno>=153 & rc$dayno<=181)
gooddays[is.na(rc$date)]=TRUE  # unknown dates allowed
goodstop=!(rc$stop %in% xstop)
#next row shouldn't be needed, but avoids any problem
goodstop[is.na(rc$stop)]=TRUE   #dummy records have unknown stop

summary(data.frame(goodyear,goodcount,gooddays,goodstop))
rc=rc[goodyear&goodcount&gooddays&goodstop,]
rc=droplevels(rc)
length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset
#check everything correct
table(list(rc$stop,rc$species),useNA = 'al')
summary(rc$dayno);table(is.na(rc$dayno)) #NAs should be present

#********************************************************************************
#further reduce subset, taking just roosts with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(rc$year,rc$site,tcount)
table(nyears)
rc=rc[(nyears[rc$site]>=minyears),]
rc=droplevels(rc)
length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#take out where never occurs, although should be very rare
#********************************************************************************
sumy=tapply(rc$batcount,rc$site,sum)
cat("\nProportion roosts where never present",round(mean(sumy==0),3),"\n")
rc=rc[(sumy[rc$site]>0),]
rc=droplevels(rc)
#length slightly different to genstat version since genstat takes out zero roosts at start
length(rc$batcount); sum(rc$batcount,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#now form data frame of variables for maximum count
#********************************************************************************
tmax=tapply(rc$batcount,list(rc$year,rc$site),max,na.rm=TRUE)
maxdf=data.frame(max=c(tmax))
maxdf$year=factor(rep(rownames(tmax),times=ncol(tmax)))
maxdf$site=factor(rep(colnames(tmax),each=nrow(tmax)))
maxdf$nobs=c(table(rc$year,rc$site))
#for dummy records and dummy dates always have single line of data but doesn't
#imply single count.  See email from Philip 29/11/11 for confirmation of this for
#NBCS records which have dummy dates
nbmpsurv=rc$datesource=="nbmp_survey"
#22/2/18 NAs occur due to missing stopreason in database but should be NBMP
nbmpsurv[is.na(nbmpsurv)]=TRUE
maxdf$NBMPsurvey=c(tapply((nbmpsurv),list(rc$year,rc$site),
                          mean,na.rm=TRUE))
maxdf$single=(maxdf$nobs==1) & (maxdf$NBMPsurvey==TRUE)
#need country info for weights
vco=tapply(as.numeric(rc$country),rc$site,max)
colev=levels(rc$country)
maxdf$country=factor(colev[vco[maxdf$site]],levels=colev)
maxdf=maxdf[!is.na(maxdf$max),]

#********************************************************************************
#check spline df, plotting graph only when testing
#********************************************************************************
with(maxdf,fgdata(max~site+year,distribution = "poisson"))
chkdf=fgindex(splinedf=3:10,nboot=0,baseyear=settings$baseyear,plot=test, fastglm = fast)

#********************************************************************************
#no covariates model using raw counts (roost dataframe, Analysis 0) 
#********************************************************************************
with(rc,fgdata(batcount~site+year,distribution = "poisson"))
title0=paste0(tcountry," ",tsp,": raw counts, no covs")
analysis0=fgindex(splinedf=settings$splinedf,nboot=settings$nb0,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,title=title0,
                  bsave=TRUE,glmboot = settings$bun0, fastglm = fast)
#write.xlsx(analysis0,file=xlsfile,sheet='Analysis 0',row.names = FALSE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 0',row.names = FALSE)
writeData(wbx, sheet = 'Analysis 0', x = analysis0)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 0', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE) 

#********************************************************************************
#no covariates model with maximum count (maxdf dataframe)
#********************************************************************************
with(maxdf,fgdata(max~site+year,distribution = "poisson"))
title1=paste0(tcountry," ",tsp,": max, no covariates")
analysis1=fgindex(splinedf=settings$splinedf,nboot=settings$nb1,nmonitor=nmon,
    baseyear=settings$baseyear,plot=TRUE,bepsilon = eps,
        title=title1,bsave=TRUE,glmboot = settings$bun1,fastglm=fast)
#write.xlsx(analysis1,file=xlsfile,sheet='Analysis 1',row.names = FALSE,append=TRUE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 1',row.names = FALSE,append=TRUE)
writeData(wbx, sheet = 'Analysis 1', x = analysis1)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 1', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#********************************************************************************
# 15/1/21 code for calculating offset to use fixed values of GB covariate 
# estimates.  From hib version
# relying on offset argument for fgdata and the gam/glm commands in fgindex
# doesn't work (no idea why), instead need to add to formula
#********************************************************************************
if (!is.na(settings$covfile)){
  print('calculating offset')
  covfile=as.character(settings$covfile)
  #covfile/gbfile is name of results file with GB covariate estimates
  gbfile=paste0("../composite/r_xlsx/",covfile,".xlsx")
  gbcov=read.xlsx(gbfile,sheet = 'Analysis 2')
  ft=sub('site+year+','',settings$model,fixed=TRUE)
  dm=model.matrix(as.formula(ft),data=maxdf)
  dmcolnames=colnames(dm)
  cov_offset=rep(0,nrow(maxdf))
  for (i in 2:ncol(dm)) {
    wh=which(gbcov$cov_name==dmcolnames[i])
    cov_offset=cov_offset+dm[,i]*gbcov$cov_est[wh]
  }
  modeloffset=paste0('max~site+year','+offset(cov_offset)')
  with(maxdf,fgdata(as.formula(modeloffset),distribution="poisson"))
  title2=paste0(tsp,": covariate estimates from ", covfile)
} else {  #standard model without offset
  settings$model=as.character(settings$model)  #goes to factor by default
  title2=paste0(tsp,": max with covariates")
  with(maxdf,fgdata(as.formula(settings$model),distribution = "poisson"))
}

#********************************************************************************
#With covariates model 
#********************************************************************************
analysis2=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,bepsilon = eps,
                  title=title2,bsave=TRUE,glmboot = settings$bun2, fastglm = fast)
#write.xlsx(analysis2,file=xlsfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
writeData(wbx, sheet = 'Analysis 2', x = analysis2)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 2', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#********************************************************************************
# 15/1/21 testing of fastglm, refit with fastglm=FALSE 
# with big datasets, just test point estimates, but with smaller compare properly
#********************************************************************************
if (nlevels(rc$site)>199) settings$nb2=0
analysis2slow=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,bepsilon = eps,
                  title=title2,bsave=TRUE,glmboot = settings$bun2)
writeData(wbx, sheet = 'slow 2', x = analysis2slow)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
diffind=abs(analysis2$Unsmoothed-analysis2slow$Unsmoothed)
summary(diffind)
if (max(diffind)>0){  #print all the time for checking, change to 0.001 if not needed
  printnice(data.frame(Year=analysis2$Year,unsmfast=analysis2$Unsmoothed,
                       unsmslow=analysis2slow$Unsmoothed,seunsmfast=analysis2$se_unsm,
                       seunsmslow=analysis2slow$se_unsm))
}
#********************************************************************************
#Calculate weights and fit model With covariates & weighted 
#********************************************************************************
if (length(tcountries)>1 & (settings$nb3>0)){
  areas=read.xlsx("C:/data/nbmp_old/nbmp2013/katepaper/areas for weighting.xlsx",sheet=1,
                  rows = 7:12)
  colabs=levels(rc$country)
  areas=areas[areas[,1]%in% colabs,]
  areas
#check order the same
  if (mean(areas[1]==colabs)<1) rm(areas)
  totobs=table(rc$country)
  wt=areas$sqkmex/as.numeric(totobs)
  wts=wt[maxdf$country] #expand to one value per survey
  #18/1/21 standardise wts as in old genstat to make easier to compare
  wts=wts/(mean(wts))
  stdwt=tapply(wts,maxdf$country,mean)
  printnice(data.frame(totobs,wt,stdwt,areas$sqkmex))  #data.frame prints nicely in columns
  if (tsp=="natterers"){
    nattdat=cbind(maxdf,wts)
    save(nattdat,file="natt.RData",version=2)  #save data to try in genstat
    wts=sqrt(wts)
    cat(paste("\n\nWARNING: weights reduced as too extreme for Natterers"))
    modwt=tapply(wts,maxdf$country,mean)
    printnice(data.frame(totobs,wt,stdwt,modwt,areas$sqkmex))  #data.frame prints nicely in columns
    }

  title3=paste0(tsp,": max, weighted covariates")
  analysis3=fgindex(splinedf=settings$splinedf,nboot=settings$nb3,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,bepsilon = eps,
                  title=title3,weight=wts,bsave=TRUE,glmboot = settings$bun3, fastglm = fast)
  #write.xlsx(analysis3,file=xlsfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)
  #write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)
  writeData(wbx, sheet = 'Analysis 3', x = analysis3)
  saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
  writeData(wbb, sheet = 'Analysis 3', x = fgglob.bsave )
  saveWorkbook(wbb, file = bootfile,overwrite = TRUE)
} #end of weighted gam if
dev.off()
timestamp() 

