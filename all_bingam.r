# 14/10/16 GAM for all species for field survey using settings spreadsheet
# 23/12/18 from 2017 exu version
# 26/1/21 fastglm added - with pip55 england actually slower with no covs, but much quicker with covs
# Initially used method=2 as for water
timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

#load functions and data
library(gam)
predict.gam=predict.Gam #27/12/19 name of predict function has changed
#library(xlsx)  #12/6/19 failed with 64 bit R
library(openxlsx)
library(fastglm)

source("C:/data/dog/rgamcode/func_fgindex.R")
load("nsp.rdata")  #much faster than using xlsx file below

settings=read.xlsx("../composite/input/gamsettings.xlsx",sheet = 'field')

#********************************************************************************
#set key variables below to give correct row of settings dataframe
#textra="standard" in most cases, but allows versions with diff base year etc.
#********************************************************************************
tcountry='Scotland' #set to GB, England, Wales, Scotland
tsp="pip45"
textra="standard"
row=which(settings$species==tsp&settings$country==tcountry&settings$extra==textra)
if (length(row)!=1) {
  print(row)
  rm(nsp)
  stop("Fault: exactly one row must match variables above.")
}
settings=settings[row,]
#settings$varname=as.character(settings$varname)
#********************************************************************************
#other settings.
#Note that 'test' changes various options to allow quick and easy checking of
#program
#jan 2020 eps is epsilon (convergence criterion) for gam fits of bootstrapped 
#data.  Default used by package gam is unnecessarily small. Optimal value seems
#to vary with dataset.  For pip45 scotland 1e-06 is optimal (see all_bingamtest.r)
#but for GB 1e-05 markedly quicker. For Noctule only 10% saving with 1e-05 so
#stick with 6.
#********************************************************************************
fast=-1     #fastglm method, -1 means don't use it, 3 generally best
eps=1e-06
test=NA  # 0 to run normally, 1 for reduced bootstraps & interactive graphs, NA automatic
if (is.na(test)) if(interactive()) test=1 else test=0
set.seed(157235) #so can reproduce exactly if required
#set.seed(6743) #so can reproduce exactly if required
minyears=2
nmon=50  #monitoring
nmindet=30  #threshold for min number of surveys to include detector
#use lower value if using GB covariates
if (!is.na(settings$covfile)) {nmindet=10}
exdet=c('Unknown')   #detectors to exclude
if (test){
  settings$nb1=10;settings$nb2=10;settings$nb3=10;nmon=1
  if (file.exists("del.xlsx")) file.remove("del.xlsx")
}
#call response variable tvar to make easier
tvar=as.character(settings$varname)
if (substring(tsp,1,3)=="pip") nvar='nspot' else nvar='nwalk'
#********************************************************************************
#end of definition section
#********************************************************************************
print(t(settings))  #transposing makes clearer

#filename for xlsx output, use del.xlsx for testing
xlsfile=paste0("../composite/r_xlsx/",settings$filename,".xlsx")
bootfile=paste0("../composite/r_boot/",settings$filename,".xlsx")
if (test) {  #added from hib 27/1/17
  xlsfile="del.xlsx";bootfile="del2.xlsx"
  grfile="not needed"
}else { #grfile name corrected 8/3/19
  grfile=paste0("all_bingam.",tsp,".",substr(tcountry,1,3),".png")
  grrow=2+(settings$nb3>0)
  png(grfile,height=grrow*300,width=480)
  par(mfrow=c(grrow,1))  #2 or 3 graphs above each other
}
data.frame(xlsfile,grfile,bootfile)
wbx=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3','slow 2')) addWorksheet(wbx,tsh)
wbb=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3')) addWorksheet(wbb,tsh)

#with(nsp,fullmodel=as.formula(charmodel))
length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)

#********************************************************************************
#country information, then take appropriate subset
#********************************************************************************
#set tcountries variable - same for wales, scot, england, but list for GB, UK
tcountries=tcountry
if (tcountry=='UK') tcountries=c("England","Wales","Scotland","NI")
if (tcountry=='GB') tcountries=c("England","Wales","Scotland")
if (tcountry=='EW') tcountries=c("England","Wales")

goodcountry=nsp$country %in% tcountries
#following ensures consistent with Genstat which included missing countries in GB
if (tcountry=="GB") goodcountry=goodcountry | is.na(nsp$country)
if (sum(goodcountry)==0) {
  rm(nsp) #so can't miss warning!
  twarn=paste0("tcountry set to ",tcountry,". No data in subset")
  warning(twarn)  }
nsp=nsp[goodcountry,]  #subset to required country
table(nsp$country,useNA="ifany")
length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)

#********************************************************************************
#build up subset for analysis
#********************************************************************************
goodyear=as.numeric(as.character(nsp$year))>=settings$firstyear
goodcount=(!is.na(nsp[[tvar]])) & nsp[[nvar]]>=6
gooddays=(nsp$dayno>=180 & nsp$dayno<=215)
#note for binomial analysis include all zero sites

nsp=nsp[goodyear&goodcount&gooddays,]
nsp=droplevels(nsp)  #removed unused levels
length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)

#********************************************************************************
#take out detectors
#********************************************************************************
detdf=data.frame(table(nsp$detector))
detdf
remrow=which(detdf$Freq<nmindet)
exdet2=append(exdet,as.character(detdf[remrow,1]))
exdet2
nsp=nsp[!nsp$detector%in%exdet2,]
nsp=droplevels(nsp)  #removed unused levels
data.frame(table(nsp$detector))

#********************************************************************************
#for serotine take out regions were doesn't occur regularly (at least 6 surveys)
#********************************************************************************
if (tsp=='serotine'){
  nreg=tapply(as.numeric(nsp[[tvar]]>0),nsp$region,sum)
  reglabs=levels(nsp$region)
  goodreg=nsp$region %in% reglabs[nreg>5]
  nsp=nsp[goodreg,]
  aftersubset=tapply(as.numeric(nsp[[tvar]]>0),nsp$region,sum)
  nsp=droplevels(nsp)  #removed unused levels
  length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)
  print(data.frame(nreg,aftersubset,row.names=reglabs))
  length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)
}

#********************************************************************************
#further reduce subset, taking just sites with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(nsp$year,nsp$site,tcount)
table(nyears)
nsp=nsp[(nyears[nsp$site]>=minyears),]
length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)

#********************************************************************************
#years (often 2001) with v few surveys can cause problems, so remove if indicated
#in settings sheet
#********************************************************************************
if (settings$removeyear>0){
  nsites=tapply(nsp$site,nsp$year,tcount)
  if (any(nsites<settings$removeyear)){
    cat("\nWARNING: some years removed due to less than",settings$removeyear,"sites:\n") 
    print(nsites[nsites<settings$removeyear])}
  nsp=nsp[(nsites[nsp$year]>=settings$removeyear),]
}
nsp=droplevels(nsp)
length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)

#********************************************************************************
#check spline df, plotting graph only when testing
#********************************************************************************
modelnocovs=paste0('cbind(',tvar,',nspot-',tvar,')~site+year')
with(nsp,fgdata(as.formula(modelnocovs),distribution = "binomial"))
chkdf=fgindex(splinedf=4:10,nboot=0,baseyear=settings$baseyear,plot=test, fastglm = fast)

#********************************************************************************
#no covariates model 
#********************************************************************************
title1=paste0(tsp,": no covariates ",tcountry)
with(nsp,fgdata(as.formula(modelnocovs),distribution = "binomial"))
analysis1=fgindex(splinedf=settings$splinedf,nboot=settings$nb1,nmonitor=nmon,
              baseyear=settings$baseyear,plot=TRUE,
      title=title1,bsave=TRUE,glmboot = settings$bun1,bepsilon = eps, fastglm = fast)
#write.xlsx(analysis1,file=xlsfile,sheet='Analysis 1',row.names = FALSE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 2',row.names = FALSE)
writeData(wbx, sheet = 'Analysis 1', x = analysis1)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 1', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#********************************************************************************
#With covariates model 
# nb define logical covariates as factors as otherwise get baffling error 
# messages when predicting
#********************************************************************************
nsp$below10C=as.factor(nsp$temp<=10)
nsp$vwind=c(1,1,2,3)[as.numeric(nsp$wind)]
nsp$Temp=cut(nsp$temp,c(0,10,15,20,99))
tapply(nsp$temp,nsp$Temp,min,na.rm=TRUE);tapply(nsp$temp,nsp$Temp,max,na.rm=TRUE)
nsp$minsafter[nsp$minsafter<=-30|nsp$minsafter>70]=NA
nsp$minsafterle20=as.factor(nsp$minsafter<=20)
#NAs problematic for wales
nsp$minsafterle20[is.na(nsp$minsafterle20)]="FALSE"
nsp$senslowmid=as.factor(nsp$sensitivity %in% c('Low-mid','Mid'))
nsp$senslow=as.factor(nsp$sensitivity=="Low")
table(nsp$sensitivity,nsp$senslow)

#replace idskills for 1998 with 1999 for observer where available
#must be a neater way of doing this!
#replace idskills for 1998 with 1999 for observer where available
# replace NA by "No data" otherwise missing vals in idpoor
nsp$idskills[is.na(nsp$idskills)]="No data"
obsskill99=tapply(as.numeric(nsp$idskills[nsp$year=="1999"]),
                  nsp$observer[nsp$year=="1999"],min,na.rm=TRUE)
skillrep=obsskill99[nsp$observer]
idlev=levels(nsp$idskills)
skillrep=factor(idlev[skillrep],levels=idlev)
nsp$idskills[nsp$year=="1998"]=skillrep[nsp$year=="1998"]
table(nsp$idskills[nsp$year==1998],useNA = "always")
nsp$idskillsnmv=nsp$idskills
nsp$idskillsnmv[nsp$idskillsnmv=='Unknown']='No data'
nsp$idskillsnmv[is.na(nsp$idskillsnmv)]='No data'
nsp$idskillsnmv=factor(nsp$idskillsnmv,levels=
                         c('Poor','OK','Good','Very Good','No data'))
table(nsp$idskillsnmv,useNA = "always")

nsp$breezy=as.factor(nsp$wind=='Breezy')
nsp$badweather=((nsp$temp<=10)+(nsp$wind=='Breezy')+(nsp$rain=='Showers'))
nsp$Badweather=as.factor(nsp$badweather)
nsp$idpoor=as.factor(nsp$idskills=='Poor')
nsp$idok=as.factor(nsp$idskills=='OK')
nsp$idpoorok=as.factor(nsp$idskills=='Poor'|nsp$idskills=='OK')
nsp$idpoor2=as.factor((nsp$idskills=='Poor')&(nsp$previous<=2))
nsp$anyrain=as.factor(nsp$rain %in% c("Light Rain","Showers"))
nsp$earlyclear=as.factor((nsp$cloud %in% c("Clear/low","Patchy"))&
                           (nsp$minsafterle20==TRUE))
#calc vcloud=mvrep(mvins(cloud;cloud.in.'Unknown');mean(cloud))
nsp$logduration=log10(nsp$duration)
nsp$logduration[nsp$duration<30 | nsp$duration>=240]=NA
nsp$blogduration=nsp$logduration
nsp$blogduration[nsp$duration<60]=log10(60)
nsp$blogduration[nsp$duration>120]=log10(120)
summary(nsp$blogduration)
nsp$week=as.factor(cut(nsp$dayno,c(0,187,194,201,208,365)))

#********************************************************************************
#now fit covariate model
#********************************************************************************
settings$model=as.character(settings$model)  #goes to factor by default

title2=paste0(tsp,": With covariates ",tcountry)
with(nsp,fgdata(as.formula(settings$model),distribution = "binomial"))
#missing values can cause problems, so exclude them
nmiss=apply(is.na(fgglob.data), 1, sum)  #n missing vals in each row
if (sum(nmiss,na.rm=TRUE)){
  cat('\n**** Missing values present in data matrix ****\n')
  print(table(nmiss,useNA = "always"))  # shows number rows with missing
  print(as.data.frame(apply(is.na(fgglob.data), 2, sum))) # shows vars with mvs
  nsp=nsp[nmiss==0,]
  nsp=droplevels(nsp)  #removed unused levels
  length(nsp[[tvar]]); sum(!is.na(nsp[[tvar]])); sum(nsp[[tvar]],na.rm=TRUE)
  with(nsp,fgdata(as.formula(settings$model),distribution = "binomial"))
}
#********************************************************************************
# 16/1/20 code for calculating offset to use fixed values of GB covariate 
# estimates
# Taken from ww version (which comes from hib) but use covariate estimates
# saved from fgglob.covest, as too many estimates to all fit on main sheet
# Note that with binomial distribution estimates slightly different to real
# fit with these covariates, since predicted at mean of offset.
#********************************************************************************
if (!is.na(settings$covfile)){
  print('calculating offset')
  covfile=as.character(settings$covfile)
  #covfile/gbfile is name of results file with GB covariate estimates
  gbfile=paste0("../composite/r_xlsx/",covfile,".xlsx")
  gbcov=read.xlsx(gbfile,sheet = 'cov')
  ft=sub('site+year+','',settings$model,fixed=TRUE)
  dm=model.matrix(as.formula(ft),data=nsp)
  dmcolnames=colnames(dm)
  cov_offset=rep(0,nrow(nsp))
  covval=rep(NA,ncol(dm)-1)
  for (i in 2:ncol(dm)) {
    wh=which(gbcov$cov_name==dmcolnames[i])
    cov_offset=cov_offset+dm[,i]*gbcov$cov_est[wh]
    covval[i]=gbcov$cov_est[wh]
  }
  print(data.frame(dmcolnames,covval))
  modeloffset=paste0(modelnocovs,'+offset(cov_offset)')
  with(nsp,fgdata(as.formula(modeloffset),distribution="binomial"))
  title2=paste0(tsp," ",tcountry,": covariate estimates from ", covfile)
}
analysis2=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
   print=c("summary","model","covariates"),baseyear=settings$baseyear,plot=TRUE,
   title=title2,bsave=TRUE,glmboot = settings$bun2,bepsilon = eps, fastglm = fast)
#write.xlsx(analysis2,file=xlsfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
writeData(wbx, sheet = 'Analysis 2', x = analysis2)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 2', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

if (nlevels(nsp$site)>199) settings$nb2=0
if (fast>=0) {    #no point if fastglm not used
  analysis2slow=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                        print=c("summary","model","covariates"),baseyear=settings$baseyear,
                        plot=TRUE,title=title2,bsave=TRUE,glmboot = settings$bun2)
  writeData(wbx, sheet = 'slow 2', x = analysis2slow)
  saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
  
  diffind=abs(analysis2$Unsmoothed-analysis2slow$Unsmoothed)
  summary(diffind)
  printnice(data.frame(Year=analysis2$Year,unsmfast=analysis2$Unsmoothed,
                       unsmslow=analysis2slow$Unsmoothed,seunsmfast=analysis2$se_unsm,
                       seunsmslow=analysis2slow$se_unsm))
} #FAST=TRUE if
#********************************************************************************
#Calculate weights and fit model With covariates & weighted 
#********************************************************************************
if (length(tcountries)>1 & (settings$nb3>0)){
  areas=read.xlsx("C:/data/nbmp_old/nbmp2013/katepaper/areas for weighting.xlsx",sheet=1,
                  rows = 7:12)
  colabs=levels(nsp$country)
  areas=areas[areas[,1]%in% colabs,]
  print(areas)
  #check order the same
  if (mean(areas[1]==colabs)<1) rm(areas)
  totobs=table(nsp$country)
  wt=areas$sqkmex/as.numeric(totobs)
  print(data.frame(totobs,wt,areas$sqkmex))  #data.frame prints nicely in columns
  wts=wt[nsp$country] #expand to one value per survey
  #18/1/21 standardise wts as in old genstat to make easier to compare
  wts=wts/(mean(wts))
  stdwt=tapply(wts,nsp$country,mean)
  printnice(data.frame(totobs,wt,stdwt,areas$sqkmex))  #data.frame prints nicely in columns
  
  title3=paste0(tsp,": Weighted covariates")
  analysis3=fgindex(splinedf=settings$splinedf,nboot=settings$nb3,nmonitor=nmon,
                    baseyear=settings$baseyear,plot=TRUE,title=title3,
                    weight=wts,bsave=TRUE,glmboot = settings$bun3,bepsilon = eps, fastglm = fast)
 # write.xlsx(analysis3,file=xlsfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)
 # write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 3',row.names = FALSE)
  writeData(wbx, sheet = 'Analysis 3', x = analysis3)
  saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
  writeData(wbb, sheet = 'Analysis 3', x = fgglob.bsave )
  saveWorkbook(wbb, file = bootfile,overwrite = TRUE)
} #end of weighted gam if
#put covariate estimates in separate sheet as versions above are incomplete
#where more covariate estimates than years.
#this will save those from analysis 3 or analysis 2 if no analysis 3
if (is.na(settings$covfile)){
  addWorksheet(wbx,'cov')
  writeData(wbx, sheet = 'cov', x = fgglob.covest)
  saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
}
dev.off()
sessionInfo()
timestamp()

