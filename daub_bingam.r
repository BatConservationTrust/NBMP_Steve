#30/12/20 openxlsx used instead of xlsx as latter gave problems with 64 bit R, but too complex really
# 20/1/17 GAM for daubs using settings spreadsheet
timestamp()
R.Version()$version.string
options(width=180)
rm(list=ls())  #ensure memory clear

#load functions and data
library(gam)
predict.gam=predict.Gam #27/12/19 name of predict function has changed
#library(xlsx)
library(openxlsx)
library(fastglm)

source("C:/data/dog/rgamcode/func_fgindex.R")
load("daub.rdata")  

settings=read.xlsx("../composite/input/gamsettings.xlsx",sheet = 'ww')
 
#********************************************************************************
#set key variables below to give correct row of settings dataframe
#textra="standard" in most cases, but allows versions with diff base year etc.
#********************************************************************************
tcountry='NI' #set to UK, England, Wales, Scotland
tsp="Daubenton's"
textra="standard"
row=which(settings$species==tsp&settings$country==tcountry&settings$extra==textra)
if (length(row)!=1) {
  print(row)
  rm(daub)
  stop("Fault: exactly one row must match variables above.")
}
settings=settings[row,]
#settings$varname=as.character(settings$varname)
#********************************************************************************
#other settings.
#Note that 'test' changes various options to allow quick and easy checking of
#program
#jan 2020 eps is epsilon (convergence criterion) for gam fits of bootstrapped 
#data.  Default used by package gam is unnecessarily small and it can be 
#increased to 0.0001 without any impact on results
#********************************************************************************
eps=0.00001
fast=FALSE     #whether to use fastglm
test=NA  # 0 to run normally, 1 for reduced bootstraps & interactive graphs, NA automatic
if (is.na(test)) if(interactive()) test=1 else test=0
set.seed(157235) #so can reproduce exactly if required
minyears=2
#nmindet=30  #threshold for min number of surveys to include detector
exdet=c('Unknown')   #detectors to exclude
#nododgytimes=1  #1 to remove very long/short, early/late surveys
nmon=50  #monitoring
if (test){
  settings$nb1=10;settings$nb2=10;settings$nb3=10;nmon=1
  if (file.exists("del.xlsx")) file.remove("del.xlsx")
}
#call response variable tvar to make easier
tvar=as.character(settings$varname)
#********************************************************************************
#end of definition section
#********************************************************************************
print(t(settings))  #transposing makes clearer

#********************************************************************************
#filename for xlsx output, use del.xlsx for testing
#********************************************************************************
xlsfile=paste0("../composite/r_xlsx/",settings$filename,".xlsx")
bootfile=paste0("../composite/r_boot/",settings$filename,".xlsx")
if (test) {  #added from hib 27/1/17
  xlsfile="del.xlsx";bootfile="del2.xlsx"
  grfile="not needed"
}else {
  grfile=paste0("daub_bingam.",textra,".",substr(tcountry,1,3),".png")
  grrow=2+(settings$nb3>0)
  png(grfile,height=grrow*300,width=480)
  par(mfrow=c(grrow,1))  #2 or 3 graphs above each other
}
data.frame(xlsfile,grfile,bootfile)
wbx=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3','slow 2')) addWorksheet(wbx,tsh)
wbb=createWorkbook()
for (tsh in c('Analysis 1','Analysis 2','Analysis 3')) addWorksheet(wbb,tsh)

#with(daub,fullmodel=as.formula(charmodel))
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#country information, then take appropriate subset
#might be worth putting all this in spreadsheet in future so don't need individual
#command file for each species
#********************************************************************************
#set tcountries variable - same for wales, scot, england, but list for GB, UK
tcountries=tcountry
if (tcountry=='UK') tcountries=c("England","Wales","Scotland","Northern Ireland")
if (tcountry=='GB') tcountries=c("England","Wales","Scotland")
if (tcountry=='EW') tcountries=c("England","Wales")
if (tcountry=='NI') tcountries=c("Northern Ireland")

goodcountry=daub$country %in% tcountries
#following ensures consistent with Genstat which included missing countries in GB
if (tcountry=="GB") goodcountry=goodcountry | is.na(daub$country)
if (sum(goodcountry)==0) {
  rm(daub) #so can't miss warning!
  twarn=paste0("tcountry set to ",tcountry,". No data in subset")
  warning(twarn)  }
daub=daub[goodcountry,]  #subset to required country
table(daub$country,useNA="ifany")
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#build up subset for analysis
#********************************************************************************
goodyear=as.numeric(as.character(daub$year))>=settings$firstyear
goodcount=(!is.na(daub[[tvar]])) & daub$nspot>=5
gooddays=(daub$dayno>=205 & daub$dayno<=250)
#note for binomial analysis include all zero sites

daub=daub[goodyear&goodcount&gooddays,]
daub=droplevels(daub)  #removed unused levels
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#5/1/19 code added to remove dodgy times if required
#see dbbingamsures.gen - excluding gives slight precision gain in weighted analysis
goodduration=daub$duration>50 & daub$duration<180
goodduration[is.na(daub$duration)]=FALSE
goodstart=daub$minsafter>20 & daub$minsafter<90
goodstart[is.na(daub$minsafter)]=FALSE
table(goodduration,useNA='always')
table(goodstart,useNA='always')
if (settings$nododgytimes){
  daub=daub[goodduration&goodstart,]
  daub=droplevels(daub)  #removed unused levels
}
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#take out detectors
#********************************************************************************
if (settings$nmindet>0){
  detdf=data.frame(table(daub$detector))
  print(detdf)
  remrow=which(detdf$Freq<settings$nmindet)
  exdet2=append(exdet,as.character(detdf[remrow,1]))
  exdet2
  daub=daub[!daub$detector%in%exdet2,]
  daub=droplevels(daub)  #removed unused levels
  print(data.frame(table(daub$detector)))
}
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#further reduce subset, taking just sites with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(daub$year,daub$site,tcount)
table(nyears)
daub=daub[(nyears[daub$site]>=minyears),]
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#years (often 2001) with v few surveys can cause problems, so remove if indicated
#in settings sheet
#********************************************************************************
if (settings$removeyear>0){
  nsites=tapply(daub$site,daub$year,tcount)
  if (any(nsites<settings$removeyear)){
    cat("\nWARNING: some years removed due to less than",settings$removeyear,"sites:\n") 
    print(nsites[nsites<settings$removeyear])}
  daub=daub[(nsites[daub$year]>=settings$removeyear),]
}
daub=droplevels(daub)
length(daub[[tvar]]); sum(daub[[tvar]],na.rm=TRUE) #for checking correct subset

#********************************************************************************
#check spline df, plotting graph only when testing
#********************************************************************************
modelnocovs=paste0('cbind(',tvar,',nspot-',tvar,')~site+year')
with(daub,fgdata(as.formula(modelnocovs),distribution = "binomial"))
chkdf=fgindex(splinedf=4:10,nboot=0,baseyear=settings$baseyear,plot=test, fastglm = fast)

#********************************************************************************
#no covariates model 
#********************************************************************************
title1=paste0(tsp,": no covariates ",tcountry)
with(daub,fgdata(as.formula(modelnocovs),distribution = "binomial"))
analysis1=fgindex(splinedf=settings$splinedf,nboot=settings$nb1,nmonitor=nmon,
              baseyear=settings$baseyear,plot=TRUE, title=title1,bsave=TRUE,
              glmboot = settings$bun1,bepsilon = eps, fastglm = fast)
#write.xlsx(analysis1,file=xlsfile,sheet='Analysis 1',row.names = FALSE)
#write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 1',row.names = FALSE)
writeData(wbx, sheet = 'Analysis 1', x = analysis1)
saveWorkbook(wbx, file = xlsfile,overwrite = TRUE)
writeData(wbb, sheet = 'Analysis 1', x = fgglob.bsave )
saveWorkbook(wbb, file = bootfile,overwrite = TRUE)

#********************************************************************************
#With covariates model 
# nb define logical covariates as factors as otherwise get baffling error 
# messages when predicting
#********************************************************************************
daub$below10C=as.factor(daub$temp<=10)
daub$bminsafter=daub$minsafter
daub$bminsafter[daub$minsafter>50]=50

#replace idskills for 1998 with 1999 for observer where available
# replace NA by "No data" otherwise missing vals in idpoor
daub$idskills[is.na(daub$idskills)]="No data"
obsskill99=tapply(as.numeric(daub$idskills[daub$year=="1999"]),
                  daub$observer[daub$year=="1999"],min,na.rm=TRUE)
skillrep=obsskill99[daub$observer]
idlev=levels(daub$idskills)
skillrep=factor(idlev[skillrep],levels=idlev)
daub$idskills[daub$year=="1998"]=skillrep[daub$year=="1998"]
table(daub$idskills[daub$year==1998],useNA = "always")
table(daub$idskills,useNA = "always")

daub$breezy=as.factor(daub$wind %in% c('Breezy','Light/ Moderate'))
table(daub$wind,daub$breezy)
daub$anyrain=as.factor(daub$rain %in% c("Light Rain","Showers"))

# replace missing values in clear with site mean, or grand mean if no data for site
table(daub$clearview,useNA = "always")
vclear=round(tapply(daub$clearview,daub$site,mean,na.rm=TRUE))
vclear[is.na(vclear)]=round(mean(vclear,na.rm=TRUE))
siteclear=vclear[daub$site]
daub$clearview[is.na(daub$clearview)]=siteclear[is.na(daub$clearview)]
table(daub$clearview,useNA = "always")

#replace idskills for 1998 with 1999 for observer where available
#must be a neater way of doing this!
#replace idskills for 1998 with 1999 for observer where available
# replace NA by "No data" otherwise missing vals in idpoor
daub$idskills[is.na(daub$idskills)]="No data"
obsskill99=tapply(as.numeric(daub$idskills[daub$year=="1999"]),
                  daub$observer[daub$year=="1999"],min,na.rm=TRUE)
skillrep=obsskill99[daub$observer]
idlev=levels(daub$idskills)
skillrep=factor(idlev[skillrep],levels=idlev)
daub$idskills[daub$year=="1998"]=skillrep[daub$year=="1998"]
table(daub$idskills[daub$year==1998],useNA = "always")
daub$idskillsnmv=daub$idskills
daub$idskillsnmv[daub$idskillsnmv=='Unknown']='No data'
daub$idskillsnmv[is.na(daub$idskillsnmv)]='No data'
daub$idskillsnmv=factor(daub$idskillsnmv,levels=
                         c('Poor','OK','Good','Very Good','No data'))
table(daub$idskillsnmv,useNA = "always")

daub$idpoor=as.factor(daub$idskillsnmv=="Poor")
daub$idok=as.factor(daub$idskillsnmv=="OK")
daub$idpoorclear=(daub$idskillsnmv=="Poor")*daub$clearview
daub$idpoorok=(daub$idskillsnmv=="Poor") | (daub$idskillsnmv=="OK")

settings$model=as.character(settings$model)  #goes to factor by default
title2=paste0(tsp,": With covariates ",tcountry)
with(daub,fgdata(as.formula(settings$model),distribution = "binomial"))
#missing values can cause problems, so exclude them
nmiss=apply(is.na(fgglob.data), 1, sum)  #n missing vals in each row
if (sum(nmiss,na.rm=TRUE)){
  cat('\n**** Missing values present in data matrix ****\n')
  print(table(nmiss,useNA = "always"))  # shows number rows with missing
  apply(is.na(fgglob.data), 2, sum) # shows variables with missing
  daub=daub[nmiss==0,]
  daub=droplevels(daub)  #removed unused levels
  print(length(daub$yvar)); print(sum(daub$yvar,na.rm=TRUE))
  with(daub,fgdata(as.formula(settings$model),distribution = "binomial"))
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
  dm=model.matrix(as.formula(ft),data=daub)
  dmcolnames=colnames(dm)
  cov_offset=rep(0,nrow(daub))
  covval=rep(NA,ncol(dm)-1)
  for (i in 2:ncol(dm)) {
    wh=which(gbcov$cov_name==dmcolnames[i])
    cov_offset=cov_offset+dm[,i]*gbcov$cov_est[wh]
    covval[i]=gbcov$cov_est[wh]
  }
  print(data.frame(dmcolnames,covval))
  modeloffset=paste0(modelnocovs,'+offset(cov_offset)')
  with(daub,fgdata(as.formula(modeloffset),distribution="binomial"),offset=cov_offset)
  title2=paste0(tsp,": cov estimates from ", covfile," ",tcountry)
}
analysis2=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
            print=c("summary","model","covariates"),baseyear=settings$baseyear,
  plot=TRUE,title=title2,bsave=TRUE,glmboot = settings$bun2,bepsilon = eps, fastglm = fast)
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
if (nlevels(daub$site)>199) settings$nb2=0
if (fast==TRUE) {    #no point if fastglm not used
  analysis2slow=fgindex(splinedf=settings$splinedf,nboot=settings$nb2,nmonitor=nmon,
                  print=c("summary","model","covariates"),baseyear=settings$baseyear,
                  plot=TRUE,title=title2,bsave=TRUE,glmboot = settings$bun2)
  #write.xlsx(analysis2,file=xlsfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
  #write.xlsx(fgglob.bsave,file=bootfile,sheet='Analysis 2',row.names = FALSE,append=TRUE)
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
  colabs=levels(daub$country)
  #23/3/19 labelling of NI different
  areas[,1]=sub('NI','Northern Ireland',areas[,1])
  areas=areas[areas[,1]%in% colabs,]
  print(areas)
  #check order the same
  if (mean(areas[1]==colabs)<1) rm(areas)
  totobs=table(daub$country)
  wt=areas$sqkmex/as.numeric(totobs)
  print(data.frame(totobs,wt,areas$sqkmex))  #data.frame prints nicely in columns
  wts=wt[daub$country] #expand to one value per survey
  
  title3=paste0(tsp,": Weighted covariates ",tcountry)
  analysis3=fgindex(splinedf=settings$splinedf,nboot=settings$nb3,nmonitor=nmon,
                    baseyear=settings$baseyear,plot=TRUE,title=title3,weight=wts,
                    bsave=TRUE,glmboot = settings$bun3,bepsilon = eps, fastglm = fast)
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

