# 13/1/21 from nbmp2019/R/hib/hib_gam_covest - version allowing use of covariate estimates for GB
# 13/10/16 GAM for all species using settings spreadsheet 
# many improvements made december 2016, relating mainly to output/input rather than model
# 28/12/19 allows use of gb covariate estimates - confirmed worked by refitting
# gb models with offset
timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

#load functions and data
library(gam)
predict.gam=predict.Gam #27/12/19 name of predict function has changed
#library(xlsx)  #12/1/21 switch to openxls as failed with 64 bit R
library(openxlsx)
library(MASS)  #for ginv()
library(fastglm)

source("C:/data/dog/rgamcode/func_fgindex.R")
load("hib.rdata")  #much faster than using xlsx file below

# 19/12/17 temperature imputation still done in genstat due to issues producing fitted
# vals from complex reml models using lmer.
# 5/2/19 now includes Hadley central england temp data for 1-3 days before survey, plus
# himp versions imputed using HadCET database
# 26/12/19 now switching to modelled 1km square Hadley Centre data, as matched
# by BCT
load("tempimp.rdata")
#check match and delete if incorrect
if (sum(hib$site!=tempimp$sitech)){
  print(summary(hib$site==tempimp$sitech))
  rm(tempimp)}
if (sum(hib$id!=tempimp$impid)){
  print(summary(hib$date==tempimp$datech))
  rm(tempimp)}
hib=cbind(hib,tempimp)

#need to restart R before re-running this line or may get java error
settings=read.xlsx("../composite/input/gamsettings.xlsx",sheet = 'hib')
                 
#********************************************************************************
#set key variables below to give correct row of settings dataframe
#Note that 'test' changes various options to allow quick and easy checking of
#program, and plots graph interactively
#textra="standard" in most cases, but allows versions with diff base year etc.
#********************************************************************************
tcountry='Scotland' #set to UK, GB, England, Wales, Scotland
tsp="daub" 
textra="standard"
test=NA  # 0 to run normally, 1 for reduced bootstraps & interactive graphs, NA automatic
if (is.na(test)){  #if NA test=1 if interactive, test=0 if batch
  if(interactive()) test=1 else test=0
}

row=which(settings$species==tsp&settings$country==tcountry&settings$extra==textra)
if (length(row)!=1) {
  print(row)
  rm(hib)
  stop("Fault: exactly one row must match variables above.")
}
settings=settings[row,]
#********************************************************************************
#other settings.
#fast is fastglm method, -1 means don't use it, 3 generally best, NA means auto
#- method 3 unless covfile set
#********************************************************************************
fast=NA    
monthstouse=c(12,1:3)
set.seed(157235) #so can reproduce exactly if required
set.seed(6385)  #ble, daub for scotland
minyears=2
nmon=50  #monitoring
if (test){
  settings$nb1=10;settings$nb2=10;settings$nb3=10;nmon=5
  if (file.exists("del.xlsx")) file.remove("del.xlsx")
}
if (is.na(fast)) fast=-1+4*(is.na(settings$covfile))
#********************************************************************************
#end of definition section
#********************************************************************************
print(t(settings))  #transposing makes clearer

#filename for xlsx output, use del.xlsx for testing
#xlsfile is for output
#bootfile is for bootstrap estimates (used for annual comparisons)
xlsfile=paste0("../composite/r_xlsx/",settings$filename,".xlsx")
bootfile=paste0("../composite/r_boot/",settings$filename,".xlsx")
if (test) {
  xlsfile="del.xlsx";bootfile="del2.xlsx"
}else {
  grfile=paste0("hib_gam.",tsp,settings$extra,".",substr(tcountry,1,3),".png")
  grrow=2+(settings$nb3>0)
  png(grfile,height=grrow*300,width=480)
  par(mfrow=c(grrow,1))  #2 or 3 graphs above each other
}
xlsfile
bootfile
wbx=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3','slow 2')) addWorksheet(wbx,tsh)
wbb=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3')) addWorksheet(wbb,tsh)

#with(hib,fullmodel=as.formula(charmodel))
length(hib[[tsp]]); sum(hib[[tsp]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#country information, then take appropriate subset
#might be worth putting all this in spreadsheet in future so don't need individual
#command file for each species
#********************************************************************************
#set tcountries variable - same for wales, scot, england, but list for GB, UK
tcountries=tcountry
if (tcountry=='UK') tcountries=c("England","Wales","Scotland","NI")
if (tcountry=='GB') tcountries=c("England","Wales","Scotland")
if (tcountry=='EW') tcountries=c("England","Wales")

goodcountry=hib$country %in% tcountries
#following ensures consistent with Genstat which included missing countries in GB
if (tcountry=="GB") goodcountry=goodcountry | is.na(hib$country)
if (sum(goodcountry)==0) {
  rm(hib) #so can't miss warning!
  twarn=paste0("tcountry set to ",tcountry,". No data in subset")
  warning(twarn)  }
hib=hib[goodcountry,]  #subset to required country
table(hib$country,useNA="ifany")
length(hib[[tsp]]); sum(hib[[tsp]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#build up subset for analysis
#********************************************************************************
goodyear=as.numeric(as.character(hib$year))>=settings$firstyear
goodcount=!is.na(hib[[tsp]])
goodmonth=hib$month %in% monthstouse
#note na.rm is needed or gives missing sum if any count missing
# sitewith=tapply(hib[[tsp]],hib$site,sum,na.rm=TRUE)>0  done at end
# goodsite=sitewith[hib$site]
tapply(goodmonth, hib$month, sum)

hib=hib[goodyear&goodcount&goodmonth,]
length(hib[[tsp]]); sum(hib[[tsp]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#further reduce subset, taking just sites with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(hib$year,hib$site,tcount)
table(nyears)
hib=hib[(nyears[hib$site]>=minyears),]
hib=droplevels(hib)
length(hib[[tsp]]); sum(hib[[tsp]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#take out where never occurs
#********************************************************************************
sumy=tapply(hib[[tsp]],hib$site,sum)
cat("\nProportion sites where never present",round(mean(sumy==0),3),"\n")
hib=hib[(sumy[hib$site]>0),]
hib=droplevels(hib)
#length slightly different to genstat version since genstat takes out zero sites at start
length(hib[[tsp]]); sum(hib[[tsp]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#check spline df, plotting graph only when testing
#********************************************************************************
modelnocovs=paste0(tsp,'~site+year')
with(hib,fgdata(as.formula(modelnocovs),distribution = "poisson"))
chkdf=fgindex(splinedf=4:10,nboot=0,baseyear=settings$baseyear,plot=test, fastglm = fast)

#********************************************************************************
#no covariates model 
#********************************************************************************
with(hib,fgdata(as.formula(modelnocovs),distribution = "poisson"))
title1=paste0(tcountry," ",tsp,": No covariates")
analysis1=fgindex(splinedf=settings$splinedf,nboot=settings$nb1,nmonitor=nmon,
    baseyear=settings$baseyear,plot=TRUE,
        title=title1,bsave=TRUE,glmboot = settings$bun1, fastglm = fast)
#write.xlsx(analysis1,file=xlsfile,sheet='Analysis 1',row.names = FALSE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 1',row.names = FALSE)
writeData(wbx, sheet = 'Analysis 1', x = analysis1)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 1', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#********************************************************************************
#With covariates model 
#5/2/19 add polynomials of day
#********************************************************************************
hib$day1=hib$dayno
hib$day2=hib$dayno**2
hib$day3=hib$dayno**3
#interaction for lhs
hib$day1h13=hib$dayno*hib$haduk13
hib$day2h13=hib$dayno**2*hib$haduk13
hib$day3h13=hib$dayno**3*hib$haduk13
#11/2/19 and as factor
hib$hadgrp=cut(hib$haduk13,c(-99,5,7,999))
hib$haduk132=hib$haduk13**2
levels(hib$hadgrp)=c('<=5','5-7','>7')
hib$coolhimpgrp=cut(hib$coolhimp,c(-99,5,7,999))
levels(hib$coolhimpgrp)=c('<=5','5-7','>7')
mon=c(1,2,3,0)[hib$month]
mon=mon+(mon==3&hib$hadgrp=='5-7')+2*(mon==3&hib$hadgrp=='>7')
hib$monhadint=as.factor(mon)
hib$hotmarch=as.factor(hib$month==3&hib$hadgrp=='>7')
levels(hib$monhadint)=c('dec','jan','feb','mar cold','mar warm','mar hot')
mon=c(1,2,3,0)[hib$month]
hib$daygrp=as.factor(mon+(hib$month==3&hib$dayno>186))
levels(hib$daygrp)=c('dec','jan','feb','mar 1-15','mar 16-31')
#merge dec with jan for ble scotland
hib$daygrpnodec=as.factor(mon+(hib$month==3&hib$dayno>186)+(mon==0))
#trying grouping for lhs interaction - 27/12/19 no longer needed
#mon=c(0,1,2,3,6)[hib$daygrp]
#mon=mon+(mon>=3&hib$hadgrp=='5-7')+2*(mon>=3&hib$hadgrp=='>7')
#table(hib$hadgrp,mon)
    
settings$model=as.character(settings$model)  #goes to factor by default
title2=paste0(tsp," ",tcountry,": With covariates")
with(hib,fgdata(as.formula(settings$model),distribution = "poisson"))

#28/12/19 missing values can cause problems, so exclude them (from nbmp2008/field)
nmiss=apply(is.na(fgglob.data), 1, sum)  #n missing vals in each row
if (sum(nmiss,na.rm=TRUE)){
  cat('\n**** Missing values present in data matrix ****\n')
  cat('\nNumbers of rows with missing data:\n')
  print(table(nmiss,useNA = "always"))  # shows number rows with missing
  cat('\nNumbers of missing values for each variable:\n')
  print(apply(is.na(fgglob.data), 2, sum)) # and variables missing
  hib=hib[nmiss==0,]
  hib=droplevels(hib)  #removed unused levels
  cat('\n')
  with(hib,fgdata(as.formula(settings$model),distribution = "poisson"))
}

#********************************************************************************
# 28/12/19 code for calculating offset to use fixed values of GB covariate 
# estimates
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
  dm=model.matrix(as.formula(ft),data=hib)
  dmcolnames=colnames(dm)
  cov_offset=rep(0,nrow(hib))
  for (i in 2:ncol(dm)) {
    wh=which(gbcov$cov_name==dmcolnames[i])
    cov_offset=cov_offset+dm[,i]*gbcov$cov_est[wh]
  }
  modeloffset=paste0(modelnocovs,'+offset(cov_offset)')
  with(hib,fgdata(as.formula(modeloffset),distribution="poisson"))
  title2=paste0(tsp," ",tcountry,": cov estimates from ", covfile)
}

analysis2=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,
                  title=title2,bsave=TRUE,glmboot = settings$bun2, fastglm = fast) 
#write.xlsx(analysis2,file=xlsfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
writeData(wbx, sheet = 'Analysis 2', x = analysis2)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 2', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#28/12/19 was some code using fgwald here, but doesn't work properly - fgwald
#needs rewriting
#********************************************************************************
#compare fastglm with standard
#********************************************************************************
if (nlevels(hib$site)>199) settings$nb2=0
if (fast>=0) {    #no point if fastglm not used
zzz 24/2/21 if retain this make sure no plotting as overwrites one needed for checking
  analysis2slow=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                        baseyear=settings$baseyear,plot=TRUE,
                        title=title2,bsave=TRUE,glmboot = settings$bun2)
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
#Weights calculated below on basis of numbers of surveys contributing to trend
#Where this results in extreme weight distribution (i.e. one country > 10X another)
#consider switching to hib$allweight done on all surveys of any species
#********************************************************************************
if (length(tcountries)>1 & (settings$nb3>0)){
  areas=read.xlsx("C:/data/nbmp_old/nbmp2013/katepaper/areas for weighting.xlsx",
                  sheet=1, rows = 7:12)
  colabs=levels(hib$country)
  areas=areas[areas[,1]%in% colabs,]
  areas
#check order the same
  if (mean(areas[1]==colabs)<1) rm(areas)
  totobs=table(hib$country)
  wt=areas$sqkmex/as.numeric(totobs)
  print(data.frame(totobs,wt,areas$sqkmex))  #data.frame prints nicely in columns
  wts=wt[hib$country] #expand to one value per survey
  wts=wts/mean(wts)
  spwt=tapply(wts,hib$country,mean)
  allwt=tapply(hib$allweight,hib$country,mean)
  if (settings$wtmethod=='all') {
    wts=hib$allweight
    title3=paste0(tsp,": Weighted (all spp) covariates")
  } else {
    title3=paste0(tsp,": Weighted covariates")
  }  #wtmethod if
  wts_used=tapply(wts,hib$country,mean)
  printnice(data.frame(totobs,wt,areas$sqkmex,spwt,allwt,wts_used)) 
  
  analysis3=fgindex(splinedf=settings$splinedf,nboot=settings$nb3,nmonitor=nmon,
                  baseyear=settings$baseyear,plot=TRUE,
                  title=title3,weight=wts,bsave=TRUE,glmboot = settings$bun3)
  #write.xlsx(analysis3,file=xlsfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)
  #write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 3',row.names = FALSE,append=TRUE)
  writeData(wbx, sheet = 'Analysis 3', x = analysis3)
  saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
  writeData(wbb, sheet = 'Analysis 3', x = fgglob.bsave )
  saveWorkbook(wbb, file = bootfile,overwrite = TRUE)
} #end of weighted gam if
dev.off()

sessionInfo()
timestamp()


