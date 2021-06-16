fgglob.version="v1.0 29/9/16 fgregion added"
fgglob.version="v1.1	03/10/2016 glmboot option added to fgindex"
fgglob.version="v1.2	18/10/2016 binomial added"
fgglob.version="v1.3	22/12/2016	bsave replaces bfile"
fgglob.version="v1.31	23/12/2016	covariate estimates saved and printed"
fgglob.version="v1.4	29/12/2016	sorting order of output variables"
fgglob.version="v1.41	31/12/2016	correcting bug in printing covariates"
fgglob.version="v1.5	9/1/2017 predicting over all sites, plus correcting bugs"
fgglob.version="v1.51	14/1/2017 print coefs with model & correct bug"
fgglob.version="v1.52	15/1/2017 fgregion changes, inc pvalue arg & return"
fgglob.version="v1.53	30/05/2017Sorts bug re tcount & message for invalid base year"
fgglob.version="v1.54 5/6/2017 Bugs in fgdata converting site, year to factors"
fgglob.version="v1.6 6/6/2017 Fixes bugs in prediction (for binomial)"
fgglob.version="v1.61 19/12/2017"  #r2 for multiple splinedf + minor bugs
fgglob.version="v1.62 21/12/2017"  #saves bootstrapped covariate estimates
fgglob.version="v1.64 23/2/2018"  #saves bootstrapped covariate estimates
fgglob.version="v1.7 27/2/2018"  #fgglob.covest
fgglob.version="v1.8 6/1/2020"  #bepsilon
fgglob.version="v1.81 6/1/2020"  #labelling of bootpred
fgglob.version="v1.9 15/1/21"   #fastglm option for fgindex
fgglob.version="v1.91 20/1/21"  #fastglm with binomial
fgglob.version="v1.92 27/1/21"  #fastglm with choice of method
#**************************************************************************** 
# Mode function (capital M because mode is built in function)
# from http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
#**************************************************************************** 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#**************************************************************************** 
# fgdata function: defines model for fgindex & fgregion
# Sets up following global variables for use in subsequent analyses
# fgglob.data main data frame
# fgglob.distribution for distribution
# fgglob.tvar strings showing info about analysis
# fgglob.formula formula for saturated model
#
# Note that formula should be the formula for a saturated model, with the sites
# factor after the ~ followed by the year factor e.g.
# count ~ site + year + covariates 
#NB NOT count ~ site + s(year, df)
#**************************************************************************** 
fgdata<-function(formula,
                 distribution=c("poisson","binomial"),
                 offset){
  vertext=fgglob.version
  #print(vertext,quote=FALSE)
  #match.arg gets full name if abbreviated in function call
  distribution=match.arg(distribution)
  if (!inherits(formula, "formula"))
    stop("Fault: argument formula needs to be an object of class formula.")
    fgglob.data<<-droplevels(data.frame(get_all_vars(formula)))
    fgnames=names(fgglob.data)

    #next 2 lines no longer needed as use match.arg
 #   if ((distribution!="poisson") & (distribution!="negbin") & (distribution!="binomial"))
 #     stop("Fault: Distribution must be set to 'poisson', 'negbin' or 'binomial'")
    if (distribution=="negbin") 
      stop("Fault: negative binomial not yet allowed") 
    #need to do checks and preferably allow scalar
    fgglob.distribution<<-distribution
    #position of site variable in dataframe
    spos=2+(distribution=="binomial")
    
    #ensure that y variate numeric - integers give problem for fastglm - NO NOT NEEDED
 #   fgglob.data[[1]]<<-as.numeric(fgglob.data[[1]])
 #   if (distribution=="binomial") fgglob.data[[2]]<<-as.numeric(fgglob.data[[2]])

    #check sites to make sure factor
    if(is.factor(fgglob.data[[spos]])==FALSE) 
      fgglob.data[[spos]]<<-factor(fgglob.data[[spos]])
    #check years to make sure factor
    if(is.factor(fgglob.data[[spos+1]])==FALSE) 
         fgglob.data[[spos+1]]<<-factor(fgglob.data[[spos+1]])
    yrlev=as.numeric(levels(fgglob.data[[spos+1]]))
    nyr=length(yrlev)
    
    nv=length(fgglob.data[[1]])
    #get numeric version of year for gams - unless go via text get ordinals
    tyears<-as(fgglob.data[[spos+1]],"character")
    fgglob.data$vyears<<-as(tyears,"numeric")
    #set weights to 1 for now, overwritten by fgindex if necessary
    fgglob.data$weights<<-rep(1,nv)
    if (missing(offset)) offset=rep(0,nv)
    fgglob.data$offset<<-offset

    #create dummy dataframe for predictions
    #v1.5 5/1/17 switch to predicting all site x year combinations, so can average
    #on natural scale - makes big difference for binomial models
    ns=nlevels(fgglob.data[[spos]])
    sitelev=levels(fgglob.data[[spos]])
    #next statement sets up df with correct number rows and correct column names
    #some columns (e.g. y variable) not needed but simpler to take all
    fgglob.preddf<<-fgglob.data[rep(1,nyr*ns),]
    fgglob.preddf[[spos]]<<-factor(rep(sitelev,each=nyr),levels=sitelev)
    fgglob.preddf[[spos+1]]<<-factor(rep(yrlev,ns),levels=yrlev)
    fgglob.preddf$vyears<<-rep(yrlev,ns) 
    ncpreddf=ncol(fgglob.preddf)
    lastcovpos=ncpreddf-3 #vyears,weights, offset at end
    firstcovpos=spos+3
    if (lastcovpos>=firstcovpos){
      for (i in firstcovpos:lastcovpos) {
        if (is.factor(fgglob.data[[i]])) fgglob.preddf[[i]]<<-Mode(fgglob.preddf[[i]])
        else fgglob.preddf[[i]]<<-mean(fgglob.data[[i]],na.rm=TRUE)
      }
    }
    
    tvar<-character(length=5)
    tvar[1]<-paste("Response variable:",fgnames[1])
    if (distribution=="binomial") paste("Binomial totals:",fgnames[2])
    else tvar[2]<-"Binomial totals:   not applicable"
    tvar[3]<-paste("Sites:",fgnames[spos])
    tvar[4]<-paste("Years:",fgnames[spos+1])
    nrow=nrow(fgglob.data)
    tvar[5]<-paste("Number of rows:",nrow)
    fgglob.tvar<<-tvar
    cat(paste(fgglob.tvar,"\n"))
    fgglob.formula<<-formula
}
#****************************************************************************
# fgdata end of function
#**************************************************************************** 

#****************************************************************************
#****************************************************************************
# fgindex function for fitting GAM models
# First use fgdata to define model
#****************************************************************************
#****************************************************************************
fgindex<-function(print="summary",
                  weights,       # default weights all equally
                  splinedf=-99,  # default 0.3*nyears
                  confidence=95,
                  nboot=0,       # number of bootstrap samples
                  nmonitor=0,    # monitoring progress e.g. 5 to show every 5th
                  baseyear=-99,  # -99 sets to second year of series
                  plot=FALSE,
                  title="NULL",
                  bsave=FALSE,   #Replaces bfile v1.3 22/12/16
                  glmboot=FALSE,
                  bepsilon=1e-04,
                  fastglm=-1)  { #-1 for standard glm, 0-3 for fastglm methods
  vertext="v0.2 18/6/16 with dataframe set up by fginit"
  vertext="v0.3 23/6/16 with dataframe for results"
  vertext="v0.4 6/7/16 switching to formula in fgdata"
  vertext="v0.5 14/7/16 (first working version), dataframe construction in fgdata"
  vertext="v0.6 4/8/16 multiple spline df"
  vertext="v0.7 12/09/2016 droplevels() in fgdata"
  vertext="v0.8 19/09/2016 Sorted bug if splinedf unset"
  vertext=fgglob.version
  #print(vertext,quote=FALSE)
  
  #check print arguments
  print=match.arg(print,choices=c("summary", "model","bootstrap","covariates","none"),
                  several.ok = TRUE)
  
  #21/1/21 offset not yet working with fastglm
  if (fastglm>=0 & regexpr('offset',deparse1(fgglob.formula)[1],ignore.case = TRUE)>1)
    stop("Fault: use fastglm=-1 if offset in formula")
  
  #****************************************************************************
  # Basic structures needed later
  #****************************************************************************
  spos=2+(fgglob.distribution=="binomial")
  #spos is position of site within fgglob.data, year is spos+1
  yrlev=as.numeric(levels(fgglob.data[[spos+1]]))
  nyr=length(yrlev)
  nsites=length(levels(fgglob.data[[spos]]))
  #Descs needs to have at least 12 values, if less years just take first nyr values
  desclength=max(nyr,12)
  Descs=character(length=desclength)
  Descs[1]=fgglob.tvar[1]
  Descs[2]=paste("Distribution:     ",fgglob.distribution)
  #  tform=paste0(deparse(fgglob.formula),collapse ='')
  tform=Reduce(paste, deparse(fgglob.formula,width.cutoff=500))
  #  tform=deparse(fgglob.formula,width.cutoff=500)
  #  tform=gsub(" ",tform,replacement="")
  Descs[3]=paste("Model:            ",tform)
  #splinedf
  ndf=length(splinedf)
  #fewster default for spline df
  if ((ndf==1)&(splinedf[1]<1)) {
    tsplsource=" (default)"
    splinedf=floor(0.3*nyr)}
  else tsplsource=" (supplied)"
  #code below allows for multiple splinedf
  Descs[4]=paste("Spline d.f.:      ",paste(splinedf,collapse=","),tsplsource)
  Descs[5]=paste("Time started:     ",format(Sys.time(), "%d %b %Y %X %Z"))
  Descs[6]=paste("N boot:           ",nboot)
  Descs[7]=paste("Confidence level: ",confidence)
  if (missing(title)) title=paste("GAM analysis ",format(Sys.time(), "%d %b %Y %X %Z"))
  Descs[8]=paste(" Title:            ",title)
  Descs[9]="GAM method:        full"  #genstat annual method not currently implemented
  if (fastglm>=0) Descs[9]=paste(Descs[9],' (fastglm method ',fastglm,")")
  Descs[10]=fgglob.tvar[2]
  Descs[11]=paste("Total sites:      ",format(nlevels(fgglob.data[[spos]])))
  Descs[12]=paste("Prog version:     ",vertext)
  fgresults=data.frame(Desc=as.character(Descs[1:nyr]),Year=yrlev)
  fgresults$Desc=as.character(fgresults$Desc)
  #,Ncount,Nsite,Rawmean,
  #                    se_mean,Smoothed,se_sm,low_sm,up_sm,Unsmoothed,se_unsm,
  #                   low_unsm,up_unsm,diff,low_diff,up_diff,chg,low_chg,up_chg)
  
  #fgresults$Desc=Descs[1:nyr]
  fgresults$Nobs=tapply(is.na(fgglob.data[[1]])==FALSE,fgglob.data[[spos+1]],sum)
  #following corrected 12/9/16 to refer to fgglob.data not names
  tcount=function(x) {length(unique(x))}
  fgresults$Nsite=tapply(fgglob.data[[spos]],fgglob.data[[spos+1]],tcount)
  fgresults$Mean=tapply(fgglob.data[[1]],fgglob.data[[spos+1]],mean)
  fgresults$SE_mean=tapply(fgglob.data[[1]],fgglob.data[[spos+1]],sd)/sqrt(fgresults$Nobs)
  
  nv=length(fgglob.data[[1]])
  ns=nlevels(fgglob.data[[spos]])
  
  if (nmonitor==0) nmonitor=99999  #0 gives fault with %%, 26/9/16 changed from -1
  
  #confidence limits
  slow=(1-confidence/100)/2
  sup=1-slow
  #if baseyear not set, take second year in series
  if (baseyear==-99)
    baseyear<-levels(fgglob.data[[spos+1]])[2]
  fgglob.baseyear=baseyear
  
  #create dataframe of variables for model
  if (!missing(weights)) fgglob.data$weights<<-weights
  
  # create formula for modela
  fgnames=names(fgglob.data)
  tspline=paste("s(vyears, ",format(splinedf[1]),")")
  tform=sub(fgnames[spos+1],tspline,tform)
  fform=as.formula(tform)
  fformglm=fgglob.formula
  
  if (sum(is.na(fgglob.data))>0){
    cat("\nWARNING: Missing values present in data matrix.  This can ")
    cat("sometimes cause problems. \nIf program crashes try removing NAs.\n") }
  
  #next few lines moved here 19/12/17 as now needed for multiple splinedf section
  if (sum(yrlev==fgglob.baseyear)==0){
    #exiting function using stop if no data for baseyear
    print(table(fgglob.data[[spos+1]]))
    cat(paste0("\n***** FAULT: no data for specified base year of ",
               fgglob.baseyear,". *****\n\n"))
    return(fgresults)
    stop("No data for specified base year",call.=FALSE)}
  baseyrno=which(yrlev==fgglob.baseyear)
  
  #****************************************************************************
  #fit models 
  #if multiple df just fit saturated (glmfit) model
  #
  #****************************************************************************
  if (ndf==1) {
    gamfit <- gam(fform, family = fgglob.distribution,data=fgglob.data,
                  offset=offset,weights=weights)  }
    
    if (fastglm>=0) {
      tformglm=as.character(fgglob.formula)
      txform=paste(tformglm[1],tformglm[3])
 #     print(txform)
      xmat=model.matrix(as.formula(txform),data=fgglob.data)
      if (fgglob.distribution=="binomial") {
        ymat=as.matrix(cbind(fgglob.data[[1]],(fgglob.data[[2]]-fgglob.data[[1]])),
                       nrow=nrow(fgglob.data),ncol=2)
      } else {
        ymat=as.matrix(fgglob.data[[1]],nrow=nrow(fgglob.data),ncol=1)
      }
      glmfit=fastglm(xmat,ymat,family=fgglob.distribution,method=fastglm,
                     offset=fgglob.data$offset,weights=fgglob.data$weights) 
      predmat=model.matrix(as.formula(txform),data=fgglob.preddf)
      predglm <- tapply(predict(glmfit,type="response",newdata=predmat,se.fit=FALSE),
                        fgglob.preddf[[spos+1]],mean)
      #print(tail(coefficients(glmfit)))
    } else {  #standard glm

        glmfit=glm(fformglm,family=fgglob.distribution,data=fgglob.data,
                   offset=offset,weights=weights) 
        predglm <- tapply(predict.glm(glmfit,type="response",newdata=fgglob.preddf),
                          fgglob.preddf[[spos+1]],mean)
        #print(tail(coefficients(glmfit)))
    }  #glm method if
  
  #****************************************************************************
  #  deal with multiple splines - fit spline models to saturated estimates
  #****************************************************************************
  if (ndf>1){
   # predglm <- tapply(predict.glm(glmfit,type="response",newdata=fgglob.preddf),
    #                  fgglob.preddf[[spos+1]],mean)
    aic=c(rep(NA,ndf));r2=c(rep(NA,ndf));adjr2=c(rep(NA,ndf))
    #structures to export results
    #will fail if n splinedf>n years, but this would be silly anyway
    fgresults$df=c(splinedf,rep(NA,(desclength-ndf)))
    fgresults$aic=c(rep(NA,desclength))
    fgresults$r2=c(rep(NA,desclength))
    #loop through fitting spline with required df to pred vals from saturated model
    for (i in 1:ndf) {
      fitspl=gam(predglm~s(yrlev,splinedf[i]))
      summ=summary(fitspl)
      aic[i]=summ$aic
      r2[i]=1-summ$deviance/summ$null.deviance
      fgresults$aic[i]=aic[i]
      fgresults$r2[i]=r2[i]
      fgresults[,i+9]=predict(fitspl)
      fgresults[,i+9]=fgresults[,i+9]/fgresults[baseyrno,i+9]*100
    } #end of for loop
    if (any(print=="summary")){
      cat("\n\n**** Comparing GAMs with different df ****\n\n")
      cat(paste(fgresults[c(8,1,10,2,3,4,9,7,6,11,12,13),1],"\n"))
      cat("\n\n**** AIC & Rsq for different d.f. ****\n\n")
      #print as data frame to get in columns
      print(as.data.frame(list(splinedf,aic,r2),col.names=c("d.f.","AIC","R-sq")))
    }#printing if
    if (plot==TRUE){fgplot(fgresults)  } 
    #exiting function using stop.  This is not ideal since throws a fault, but
    #saves the complication of putting subsequent programming in an if statement
    return(fgresults)
    stop("No bootstrapping with multiple spline d.f.",call.=FALSE)
  }  #end of multiple splinedf if()
  ppredgam=predict.gam(gamfit,type="response",newdata=fgglob.preddf)
  predgam <- tapply(ppredgam,
                    fgglob.preddf[[spos+1]],mean)
  
  fgresults$Smoothed=predgam/predgam[baseyrno]*100
  fgresults$se_sm=c(rep(NA,nyr))
  fgresults$low_sm=c(rep(NA,nyr))
  fgresults$up_sm=c(rep(NA,nyr))
  
  #predglm <- tapply(predict(glmfit,type="response",newdata=fgglob.preddf),
   #                 fgglob.preddf[[spos+1]],mean)
  fgresults$Unsmoothed=predglm/predglm[baseyrno]*100
  fgresults$se_unsm=c(rep(NA,nyr))
  fgresults$low_unsm=c(rep(NA,nyr))
  fgresults$up_unsm=c(rep(NA,nyr))
  
  #23/12/16 get covariate estimates
  est=coefficients(gamfit)
  estlab=names(est)
  # calc starting position for covariates in vector of estimates constant+(nsites-1)+linyr
  covstart=nsites+2  
  nvest=length(est)
  ncov=length(covstart:nvest)
  setcov=covstart<=nvest
  #print(data.frame(nyr,nsites,covstart,nvest,ncov,setcov))
  if (setcov) covest=est[covstart:nvest]
  
  ask="\n*******************************************************\n\n"
  if (any(print=="model")){
    cat("\n**************** GAM Model ****************\n\n")
    #   print(summary(fgglob.data))
    print(summary(gamfit))
    #   print(data.frame(yrlev,predgam))
    cat("\n**************** Coefficients from GAM Model ****************\n\n")
    nsp1=ns+1
    print(data.frame(Estimate=est[c(1,nsp1:nvest)]))
    cat("\nExcludes site estimates\n")
    cat(ask)
  }
  #****************************************************************************
  # On to bootstrapping if required, but first define remaining cols of fgresults
  # so that they exist and don't cause fault in fgprint
  # nb don't change order of definition or mucks up results output to xlsx
  #****************************************************************************
  fgresults$diff=c(rep(NA,nyr))
  fgresults$low_dif=c(rep(NA,nyr))
  fgresults$up_dif=c(rep(NA,nyr))
  fgresults$chg=c(rep(NA,nyr))
  fgresults$low_chg=c(rep(NA,nyr))
  fgresults$up_chg=c(rep(NA,nyr))
  fgresults$sig_chg=c(rep(NA,nyr))
  fgresults$sig_dif=c(rep(NA,nyr)) # out of natural order as not present in Genstat
  fgresults$cov_name=c(rep("",nyr))
  fgresults$cov_est=c(rep(NA,nyr))
  fgresults$low_cov=c(rep(NA,nyr))
  fgresults$up_cov=c(rep(NA,nyr))
  if (nboot<1){
    if (any(print=="summary")){
      fgprint(fgresults,round=TRUE)
      cat(ask)
    }
    #    if (plot==TRUE){fgplot(fgresults)  } gives fault
    return(fgresults)
    stop()  }
  
  #structures for bootstrapping
  siteunits=split(1:nv,fgglob.data[[spos]]) #unitnos for each site
  nrowsites=c(table(fgglob.data[[spos]]))
  
  #need to change sites in preddf so uses new site numbers 1:nsites not original levels
  local.preddf=fgglob.preddf
  local.preddf[[spos]]=factor(rep(1:nsites,each=nyr),levels=c(1:nsites))
  Xyrlev=paste0("X",yrlev)
  bootpred=matrix(nrow=nboot,ncol=nyr,dimnames=list(1:nboot,Xyrlev))
  xyrlev=paste0("x",yrlev)
  bglmpred=matrix(nrow=nboot,ncol=nyr,dimnames=list(1:nboot,xyrlev))
  if (setcov) bcovest=matrix(nrow=nboot,ncol=length(covstart:nvest))
  
  time0=Sys.time()
  for (i in 1:nboot) {
    bootsites=sample(1:nsites,nsites,replace=TRUE) #sites in bootstrap sample
    bootunits=unlist(siteunits[bootsites]) #append the different sets of units
    bsites=rep(1:nsites,times=nrowsites[bootsites])
    newdf=fgglob.data[bootunits,]
    #overwrite sites with new site nos, rather than site data has come from
    newdf[[spos]]=as.factor(bsites)
    
    bootgam <- gam(fform, family = fgglob.distribution,data=newdf,
                   offset=offset,weights=weights,epsilon=bepsilon)  
    
    bootpred[i,] <- tapply(predict.gam(bootgam,type="response",newdata=local.preddf),
                           local.preddf[[spos+1]],mean)
    if (setcov) {
      best=coefficients(bootgam) #changed est to best 14/1/17
      bcovest[i,]=best[covstart:nvest]
    }
    if (glmboot){
      if (fastglm>=0) {
        tformglm=as.character(fgglob.formula)
        txform=paste(tformglm[1],tformglm[3])
        xmat=model.matrix(as.formula(txform),data=newdf)
        if (fgglob.distribution=="binomial") {
          ymat=as.matrix(cbind(newdf[[1]],(newdf[[2]]-newdf[[1]])),
                         nrow=nrow(newdf),ncol=2)
        } else {
          ymat=as.matrix(newdf[[1]],nrow=nrow(newdf),ncol=1)
        }
        bootglm <- fastglm(xmat, ymat, family = fgglob.distribution,
                     offset=newdf$offset,weights=newdf$weights, method=fastglm)  
        predmat=model.matrix(as.formula(txform),data=local.preddf)
        bglmpred[i,] <- tapply(predict(bootglm,type="response",newdata=predmat, se.fit=FALSE),
                             local.preddf[[spos+1]],mean)
      } else {
        bootglm <- glm(fformglm, family = fgglob.distribution,data=newdf,
                       offset=offset,weights=weights)  
        bglmpred[i,] <- tapply(predict.glm(bootglm,type="response",newdata=local.preddf),
                               local.preddf[[spos+1]],mean)
        
      }
    }
    #output time etc every nmonitor bootstraps, plus details if print="bootstrap"
    if ((i%%nmonitor)==0) {
      ttime=format(Sys.time(), "%d %b %Y %X %Z")
      cat(paste("loop:",i," ",ttime,"\n"))
      if (any(print=="bootstrap")){
        cat("\nBootstrap sites:",bootsites)
        print(summary(bootgam))
        cat(ask)}}
  } 
  #end nboot loop
  timetaken=Sys.time()-time0
  fgresults[13,1]=paste("Time 100 loops:   ",
                        format(timetaken*100/nboot,nsmall=1))
  #express relative to baseyear
  bootpred=bootpred/bootpred[,baseyrno]*100
  if(nboot>5) fgresults$se_sm=apply(bootpred,2,sd)
  if ((nboot<99)&(nboot>5)){
    fgresults$low_sm=fgresults$Smoothed-2*fgresults$se_sm
    fgresults$up_sm=fgresults$Smoothed+2*fgresults$se_sm
    lmeth="limits from bootstrap s.e."
  } 
  if (nboot>98){
    fgresults$low_sm=apply(bootpred,2,quantile,prob=slow)
    fgresults$up_sm=apply(bootpred,2,quantile,prob=sup)
    lmeth="percentile limits"
  }
  #and similarly for the unsmoothed version from the GLM
  if (glmboot){
    bglmpred=bglmpred/bglmpred[,baseyrno]*100
    if(nboot>5) fgresults$se_unsm=apply(bglmpred,2,sd)
    if ((nboot<99)&(nboot>5)){
      fgresults$low_unsm=fgresults$Unsmoothed-2*fgresults$se_unsm
      fgresults$up_unsm=fgresults$Unsmoothed+2*fgresults$se_unsm
    } 
    if (nboot>98){
      fgresults$low_unsm=apply(bglmpred,2,quantile,prob=slow)
      fgresults$up_unsm=apply(bglmpred,2,quantile,prob=sup)
    }
  } #glmboot
  #****************************************************************************
  # Change and difference
  #****************************************************************************
  #temporarily write bootpred as global to help programming
  bootpred<<-bootpred
  nyr<<-nyr
  
  for (i in 2:nyr) {
    fgresults$diff[i]=fgresults$Smoothed[i]-fgresults$Smoothed[i-1]
    dif=bootpred[,i]-bootpred[,(i-1)]
    fgresults$low_dif[i]=quantile(dif,probs=slow)
    fgresults$up_dif[i]=quantile(dif,probs=sup)
    
    fromend=min(i-1,nyr-i)
    #now approximations to second derivative from Fewster et al p1977
    if (fromend==1){
      deriv2=bootpred[,i+1]-2*bootpred[,i]+bootpred[,i-1]
      fgresults$chg[i]=fgresults$Smoothed[i+1]-2*fgresults$Smoothed[i]+fgresults$Smoothed[i-1]}
    if (fromend==2){
      deriv2=(-bootpred[,i+2]+16*bootpred[,i+1]-30*bootpred[,i]+16*bootpred[,i-1]-bootpred[,i-2])/12
      fgresults$chg[i]=(-fgresults$Smoothed[i+2]+16*fgresults$Smoothed[i+1]-30*fgresults$Smoothed[i]+
                          16*fgresults$Smoothed[i-1]-fgresults$Smoothed[i-2])/12}
    if (fromend>2){
      deriv2=(2*bootpred[,i+3]-27*bootpred[,i+2]+270*bootpred[,i+1]-490*bootpred[,i]+
                270*bootpred[,i-1]-27*bootpred[,i-2]+2*bootpred[,i-3])/180
      fgresults$chg[i]=(2*fgresults$Smoothed[i+3]-27*fgresults$Smoothed[i+2]+270*fgresults$Smoothed[i+1]-
                          490*fgresults$Smoothed[i]+270*fgresults$Smoothed[i-1]-27*fgresults$Smoothed[i-2]+
                          2*fgresults$Smoothed[i-3])/180}
    if(fromend>0) {
      fgresults$low_chg[i]=quantile(deriv2,probs=slow)
      fgresults$up_chg[i]=quantile(deriv2,probs=sup) }
    deriv2=c(rep(NA,40)) #so values not carried forward to next loop
  }
  fgresults$sig_dif=ifelse(fgresults$low_dif>0 | fgresults$up_dif<0,"sig","NS")
  fgresults$sig_chg=ifelse(fgresults$low_chg>0 | fgresults$up_chg<0,"sig","NS")
  
  #covariates
  #26/2/18 allow for possibility that more covs than can fit in fgresults
  if (setcov){
    fgglob.covest<<-data.frame(cov_name=estlab[covstart:nvest])
    fgglob.covest$cov_est<<-est[covstart:nvest]
    # 27/2/18 na.rm=TRUE as otherwise get issues with rare categories
    ncmiss=sum(is.na(bcovest),na.rm=TRUE)
    if (ncmiss>0) {
      cat(paste0("\nWARNING: ",ncmiss,
                 " missing values in covariate bootstrap estimates.\n "))
    }
    fgglob.covest$low_cov<<-apply(bcovest,2,quantile,prob=slow,na.rm=TRUE)
    fgglob.covest$up_cov<<-apply(bcovest,2,quantile,prob=sup,na.rm=TRUE)
    covend=nvest
    nend=ncov
    if (ncov>nyr) {
      cat("\nWARNING: Too many covariates to save in results structure.\n ")
      cat("Full list of covariates stored in fgglob.covest.\n")
      nend=nyr
    }
    #cov_name kept defaulting to numeric
    covcols=c('cov_est','low_cov','up_cov')
    fgresults$cov_name=c(rep("",nyr))
    fgresults[1:nend,'cov_name']=as.character(fgglob.covest[1:nend,'cov_name'])
    fgresults[1:nend,covcols]=fgglob.covest[1:nend,covcols]
    colnames(bcovest)=estlab[covstart:nvest]
    #21/12/17 write to global so can use in Wald tests etc, 26/2/18 rename
    fgglob.bcovest<<-bcovest
  }
  
  if (any(print=="summary")){
    fgprint(fgresults,round=TRUE)
    cat(ask)
  }
  if (plot==TRUE){fgplot(fgresults)  } 
  
  if (bsave) {
    #    if (file.exists(bfile)) file.remove(bfile)  #14/10/16
    if (glmboot) {
      fgglob.bsave<<-cbind(bootpred,bglmpred)
    } else fgglob.bsave<<-bootpred
  }
  
  return(fgresults) } #end of fgindex

#****************************************************************************
#****************************************************************************
# fgprint
#****************************************************************************
#**************************************************************************** 
fgprint<-function(fgresults,round=TRUE){
  printdf=fgresults  #take copy so don't round supplied structure
  if (round==TRUE) {
    # if round numeric structure to NA goes to all NAs
    dec=c(NA,0,0,0,1,2,2,3,2,2 ,2,3,2,2 ,2,2,2 ,2,2,2,NA,NA,NA)
    nc=length(dec)
    for (i in 1:nc){
      if (is.numeric(printdf[[i]])) printdf[i]=round(printdf[i],dec[i])
    } #for loop
  } #if for rounding
  ask="\n****************************************************************\n\n"
  cat("\n\n**** GAM Bootstrapping Results ****\n\n")
  #use which to remove empty rows in Desc - or instead specify directly
  #cat(paste(printdf[which(printdf$Desc!=""),1],"\n"))
  cat(paste(printdf[c(8,1,10,2,3,4,9,7,6,11,12,13),1],"\n"))
  
  cat("\n\n**** GAM Smoothed & Unsmoothed Index ****\n\n")
  print(printdf[,c(2,3,4,5,7,8,9,10,11,12)])
  
  cat("\n\n**** Differences and change points ****\n\n")
  print(printdf[,c(2,15,16,17,22,18,19,20,21)])
  
  covpres=!is.na(printdf[24])
  if (sum(covpres,na.rm=TRUE)>=1){
    cat("\n\n**** Covariates & bootstrap confidence intervals ****\n\n")
    print(printdf[covpres,23:26])   
  }
 
} #end of fgprint

#****************************************************************************
#****************************************************************************
# fgplot
#**************************************************************************** 
#****************************************************************************
fgplot=function(fgresults,
                xmin=min(fgresults$Year),
                xmax=max(fgresults$Year),
                ymin,
                ymax){
  #extract title itself
  mtitle=substring(fgresults$Desc[8],20,99)
  
#***************************************************************************
#need either to plot full results or, for multiple dfs, just the aics
  if (names(fgresults)[7]=="df"){
#split graphics window to plot one above another
#margins must be set or get vast amounts of space
# 19/12/17 keep as one graph at a time or too squashed
#   parold=par(mfrow=c(2, 1),mar=c(3,4,2,2),oma=c(0,0,0,0))
   
#    ymin=min(fgresults$aic,na.rm=TRUE)
#    ymax=max(fgresults$aic,na.rm=TRUE)
    splinedf=fgresults$df[which(is.na(fgresults$df)==0)]
    ndf=length(splinedf)
    plot(fgresults$df,fgresults$r2,ylab="Rsq",xlab="Spline d.f.",main=
           "Comparing R-squared of spline models")
    
#and plot fits
    nc=ncol(fgresults)
    matplot(fgresults$Year,fgresults[,10:nc],type="l",lty=1,col=c(10:nc),
            ylab="index", xlab=" ",main="Comparing fits")
    legend("right",paste(format(splinedf)," d.f."),cex=0.6,lty=1,lwd=1,
           col=c(10:nc))
#    ymin=min(unlist(fgresults[,9:nc]),na.rm=TRUE)
#    ymax=max(unlist(fgresults[,9:nc]),na.rm=TRUE)
#    plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="Index",xlab=" ",main=mtitle)
#    for (i in 1:ndf) {points(fgresults$Year,fgresults[,i+8],type="l",col=i)}
#    par(parold)
  }
  else {
  #set up axes and then add lines/points one at a time
    #14/1/17 switch to ymin=0 unless specified
  if (missing(ymin)) ymin=0  #min(fgresults$low_sm,na.rm=TRUE)
  ymax=max(fgresults$up_sm,na.rm=TRUE)
  plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="Index",xlab=" ",main=mtitle)
  points(c(xmin-1,xmax+1),c(100,100),type="l")
  points(fgresults$Year,fgresults$Smoothed,type="l")
  #lty=3 gives dotted line type for conf interval
  points(fgresults$Year,fgresults$low_sm,type="l",lty=3)
  points(fgresults$Year,fgresults$up_sm,type="l",lty=3)
  points(fgresults$Year,fgresults$Unsmoothed,type="p",pch=16,col="green")
  # for sig differences plot between two years, using linear interpolated y value
  nrm1=nrow(fgresults)-1
  smdiff=(fgresults$Smoothed+fgresults$Smoothed[c(1,1:nrm1)])/2
  ydiff=ifelse(fgresults$sig_dif=="sig",smdiff,-100)
  xdiff=fgresults$Year-0.5
  points(xdiff,ydiff,type="p",pch=17,col="red")
  #changepoints
  ychange=ifelse(fgresults$sig_chg=="sig",fgresults$Smoothed,-100)
  points(fgresults$Year,ychange,type="p",pch=15,col="red")
  mtext(c("red squares are sig change points, where slope changes",
          "red triangles are where difference between years is sig"),
        side=1,line=c(2,3),cex=0.6)
  } #end of plotting for full results
} #end of fgplot

#****************************************************************************
#****************************************************************************
# fgregion
#****************************************************************************
#****************************************************************************
fgregion<-function(print="summary",
                   weights,       # default weights all equally
                   splinedf=-99,  # default 0.3*nyears
                   pvalue=0.05,   #pvalue (For early exit 14/1/17)
                   nrand=100,       # number of randomisations
                   ncheck=0,        # number at which to check for early halt
                   nmonitor=0,    # monitoring progress e.g. 5 to show every 5th
                   baseyear=-99,  # -99 sets to second year of series
                   region,        # regions to check for differences
                   plot=FALSE,
                   title="NULL")  {
  #edited from v0.8 of fgindex
  vertext=fgglob.version
  #print(vertext,quote=FALSE)
  
  #check print arguments
  print=match.arg(print,choices=c("summary", "model","rand","none","yeartab"),
                  several.ok = TRUE)
  
  #****************************************************************************
  # Basic structures needed later
  #****************************************************************************
  spos=2+(fgglob.distribution=="binomial")
  #spos is position of site within fgglob.data, year is spos+1
  yrlev=as.numeric(levels(fgglob.data[[spos+1]]))
  nyr=length(yrlev)
  #Descs needs to have at least 12 values, if less years just take first nyr values
  desclength=max(nyr,12)
  Descs=character(length=desclength)
  Descs[1]=fgglob.tvar[1]
  Descs[2]=paste("Distribution:     ",fgglob.distribution)
#  tform=deparse(fgglob.formula)
  tform=Reduce(paste, deparse(fgglob.formula,width.cutoff=500))
  Descs[3]=paste("Model:            ",tform)
  #splinedf
  ndf=length(splinedf)
  #fewster default for spline df
  if ((ndf==1)&(splinedf[1]<1)) {
    tsplsource=" (default)"
    splinedf=floor(0.3*nyr)}
  else tsplsource=" (supplied)"
  #code below allows for multiple splinedf
  Descs[4]=paste("Spline d.f.:      ",paste(splinedf,collapse=","),tsplsource)
  Descs[5]=paste("Time started:     ",format(Sys.time(), "%d %b %Y %X %Z"))
  Descs[7]=paste("N initial rand:   ",nrand)
  Descs[6]=paste("Check after rand: ",ncheck)
  if (missing(title)) title=paste("GAM regions ",format(Sys.time(), "%d %b %Y %X %Z"))
  Descs[8]=paste(" Title:            ",title)
  Descs[9]="GAM method:        full"  #genstat annual method not currently implemented
  Descs[11]=paste("Total sites:      ",format(nlevels(fgglob.data[[spos]])))
  Descs[12]=paste("Prog version:     ",vertext)
  #,Ncount,Nsite,Rawmean,
  #                    se_mean,Smoothed,se_sm,low_sm,up_sm,Unsmoothed,se_unsm,
  #                   low_unsm,up_unsm,diff,low_diff,up_diff,chg,low_chg,up_chg)
  
  nv=length(fgglob.data[[1]])
  
  if (nmonitor==0) nmonitor=99999  #0 gives fault with %%, 26/9/16 changed from -1
  
  #if baseyear not set, take second year in series
  if (baseyear==-99)
    baseyear<-levels(fgglob.data[[spos+1]])[2]
  fgglob.baseyear=baseyear
  
  #create dataframe of variables for model
  if (!missing(weights)) fgglob.data$weights<<-weights
  
  # create formula for model
  region=droplevels(region)
  fgnames=names(fgglob.data)
  tspline=paste("s(vyears, ",format(splinedf[1]),")")
  tform=sub(fgnames[spos+1],tspline,tform)
  fform=as.formula(tform)
  fformreg=as.formula(paste0(tform,"+region"))
  
  yeartab=table(fgglob.data[[spos+1]],region)
  
  #****************************************************************************
  #structures for rand test
  #****************************************************************************
  reglevs=levels(region)
  nreg=length(reglevs)
  sitereg=tapply(as.numeric(region),fgglob.data[spos],min)
  # 14/1/17use na.rm=TRUE below in case missing values or single obs
  if (sum(tapply(as.numeric(region),fgglob.data[spos],var),na.rm=TRUE)>0){
    stop("Each site can only be in one region.")  }
  devsep=c(rep(NA,nrand));dfsep=c(rep(NA,nrand))
  devcom=c(rep(NA,nrand));dfcom=c(rep(NA,nrand))
  ask = "\n*******************************************************\n\n"
  time0=Sys.time()
  for (j in 1:nrand) {
    if ((j%%nmonitor)==0) {
      ttime=format(Sys.time(), "%d %b %Y %X %Z")
      cat(paste("********** loop:",j," ",ttime," **********\n"))}
    #****************************************************************************
    #fit overall model
    #****************************************************************************
   # if (fgglob.distribution == "poisson") {
      gamfit <-gam(fformreg,family = fgglob.distribution,
                   data = fgglob.data,offset = offset,weights = weights)  
    
    #print summary of model if requested and if either first loop (real data) or nmon
    if (any(print == "model") & (((j %% nmonitor) == 0) | (j == 1))) {
      cat("\n**************** GAM Model ",j,"****************\n\n")
      print(summary(gamfit))
      cat(ask)    }
    devcom[j] = deviance(gamfit)
    dfcom[j] = df.residual(gamfit)
    
    #define these in loop so values not carried from last loop
    rdevsep=c(rep(NA,nreg));rdfsep=c(rep(NA,nreg)) #regional dev/df at each rand
    for (i in 1:nreg) {
      subfit = gam(fform,family = fgglob.distribution,
                   data = subset(fgglob.data, (region == reglevs[i])),
                   offset = offset,weights = weights) 
      #print deviance of model if requested and if either first loop (real data) or nmon
      
      rdevsep[i] = deviance(subfit)
      rdfsep[i] = df.residual(subfit)
      if (any(print == "model") & (((j %% nmonitor) == 0) | (j == 1))) {
        cat("\n**************** ",reglevs[i]," ****************\n\n")
        cat("\nseparate dev:   ",rdevsep[i],"\nseparate df:    ",rdfsep[i])
        cat(ask)
        print(rdfsep)}
    } #regions loop i
    devsep[j]=sum(rdevsep)
    dfsep[j]=sum(rdfsep)
    if (j==ncheck){  #check for early exit if cant be sig in full set
      ngt=sum((devcom-devsep)>=(devcom[1]-devsep[1]),na.rm=TRUE)
      if ((ngt-1)>(nrand*pvalue)) { #14/1/17 -1 to allow for real model
        nrand=j
        break}
      else{   #added 14/1/17
        cat("\n**** Running further randomisations ****\n")
        cat("\nncheck = ",ncheck,", ngt  ",ngt,"\n\n")}
    }  #ncheck if
    #********************************************************************************
    # finally randomise regions to sites for next loop
    #********************************************************************************
    rsitereg=sample(sitereg) #sample with no args set randomises data
    region=reglevs[rsitereg[fgglob.data[[spos]]]]
  } #nrand loop j
  timetaken=Sys.time()-time0
  devdif=devcom-devsep
  dfdif=dfcom-dfsep
  dfint=dfcom-dfsep  #df for interaction
  #asymptotic deviance test
  preal=pchisq(devdif[1],dfdif[1],lower.tail = FALSE)
  prand=mean(devdif>=devdif[1],na.rm = TRUE)
  ngteq=sum(devdif>=devdif[1],na.rm = TRUE)
  nless=sum(devdif<devdif[1],na.rm = TRUE)
  pvalues=structure(c(preal,prand),names=c("P dev","P rand"))
  Descs[10]=paste("N final rand:     ",nrand)  #nb recalculated if quit early
  
  if (any(print == "summary")) {
    cat("\n\n**** GAM Regional Analysis ****\n\n")
    #use which to remove empty rows in Desc - or instead specify directly
    #cat(paste(printdf[which(printdf$Desc!=""),1],"\n"))
    Descs[13]=paste("Time 100 loops:   ",
          format(timetaken*100/nrand,nsmall=1))
    printdf=data.frame(Descs)
    cat(paste(printdf[c(8,1,2,3,4,9,7,6,10,11,12,13),1],"\n"))
  }
  if (any(print == "yeartab")) print(yeartab)
  
  if (any(print == "summary")) {
    cat("\n**************** Deviances for observed data ****************\n")
    cat("\ncommon dev:     ",devcom[1],"\ncommon df:      ",dfcom[1])
    cat("\nseparate dev:   ",devsep[1],"\nseparate df:    ",dfsep[1])
    cat("\ndifference dev:   ",devdif[1],"\ndifference df:    ",dfint[1])
    cat("\n\n**************** Deviance test for differences in trend ****************\n")
    cat("\nAsymptotic P-value: ",round(preal,3))
    cat("\nRandomisation P-value: ",round(prand,3),"\n")
    quants=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)
    devquants=quantile(devdif,probs=quants,na.rm = TRUE)
    print(data.frame(quants,devquants))
  }
  #check for fails where df for interaction wrong
  #note that df not always exact integer so dont just use !=
  nfail=sum((abs(dfint-dfint[1]))>0.1,na.rm = TRUE)
  pfail=nfail/nrand
  if (pfail>0.05){
    print(table(round(dfint,2)))
    cat("\nRandomisations fail: ",nfail," due to incorrect df")
    stop("Too many randomisations fail.")  }
  if (nfail>0){
    cat("\nWARNING: ",nfail,"randomisations fail due to incorrect df") 
    print(table(round(dfint,2)))}
  #nvalues added 15/1/17
  nvalues=structure(c(nrand,ngteq,nless,nfail),
                    names=c("n rand","n>=real","n<real","n fail"))
  results=setNames(list(devcom,devsep,dfcom,dfsep,Descs,pvalues,nvalues),
                  c("devcom","devsep","dfcom","dfsep","Descs","pvalues","nvalues"))
  return(results)
} #end of fgregion

#****************************************************************************
#****************************************************************************
# fgwald
#****************************************************************************
#**************************************************************************** 
fgwald<-function(fgresults,covlabs,print=TRUE){
  wh=which(fgresults$cov_name%in%covlabs)
  if (length(wh)==0){
    #exiting function using stop if no estimates found
    stop("covlabs do not match covariate names",call.=FALSE)}
  
  if (length(wh)!=length(covlabs)){
    #exiting function using stop if no estimates found
    cat(paste0("\nWarning: following covlabs do not match covariate names:\n\n"))
    print(covlabs[!covlabs%in%fgresults$cov_name])
    #remove non-matching
    covlabs=covlabs[covlabs%in%fgresults$cov_name]
  }
  
  est=as.matrix(fgresults[wh,'cov_est'])
  df=length(est)
  vc=var(fgglob.covest[,covlabs])
  wald=round(as.numeric(crossprod(est,ginv(vc))%*%est),2)
  Pvalue=pchisq(wald,df, lower=FALSE)
  if (Pvalue<0.001) {Pvalue='<0.001'} else {
    Pvalue=round(Pvalue,3)
  }
  
  if (print==TRUE) {
  cat("\n**** Wald test for covariates in FGINDEX model ****\n\n")
  print(fgresults[fgresults$cov_name%in%covlabs,c('cov_name','cov_est')])
  cat(paste0("\nChi-squared = ",wald,' with ',df, ' d.f., P = ',Pvalue,"\n\n"))
  }

} #end of fgwald

#****************************************************************************
#****************************************************************************
# fgdiffplot
#****************************************************************************
#**************************************************************************** 
fgdiffplot<-function(bdata1,bdata2,tcountry1,tcountry2,title,print=TRUE,
                     plot=c('Indices','Diffdiff')){
  par(xpd=TRUE)
  #check print arguments
  plot=match.arg(plot,choices=c('Indices','Difference','Diffdiff'),
                  several.ok = TRUE)
  nyr1=ncol(bdata1)/2
  tyrs1=colnames(bdata1[,1:nyr1])
  nb1=nrow(bdata1)
  nyr2=ncol(bdata2)/2
  tyrs2=colnames(bdata2[,1:nyr2])
  nb2=nrow(bdata2)
  
  #********************************************************************************
  #check same years and nboot
  #********************************************************************************
  if (nb1!=nb2) {
    cat("\nN bootstraps differ\n")
    print(data.frame(nb1,nb2))
    stop("N bootstraps differ",call.=FALSE)
  }
  if (nyr1!=nyr2) {
    cat("\nWarning: number of years differ\n")
    print(data.frame(nyr1,nyr2))
    tyrs1=tyrs1[tyrs1%in%tyrs2]
    tyrs2=tyrs2[tyrs2%in%tyrs1]
    cat("\nFollowing years retained\n")
    print(tyrs2)  }
  
  #********************************************************************************
  # get differences
  #********************************************************************************
  bdiff=bdata1[,tyrs1]-bdata2[,tyrs1]  #nb tyrs1 & tyrs2 should be the same
  dmean=colMeans(bdiff)
  mean1=colMeans(bdata1[,tyrs1])
  mean2=colMeans(bdata2[,tyrs1])
  dlow=apply(bdiff,2,quantile,prob=slow)
  dup=apply(bdiff,2,quantile,prob=sup)
  yrs=as.numeric(substr(tyrs1,2,5))
  xmin=min(yrs);xmax=max(yrs)
  ymin=0;ymax=max(100,mean1,mean2)
  
  if ('Indices'%in%plot) {
    main=paste0(title,': Indices')
    plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="Index",xlab=" ",main=main)
    points(c(xmin-1,xmax+1),c(100,100),type="l")
    points(yrs,mean1,type="l",col='red')
    points(yrs,mean2,type="l",col='green')
    legend('topleft', c(tcountry1,tcountry2) , 
         lty=1, col=c('red', 'green'), bty='n', cex=.75)
  } # plot if
  #********************************************************************************
  # plot differences between region
  #********************************************************************************
  if ('Difference'%in%plot) {
    ymin=min(dlow);ymax=max(dup)
    main=paste0(title,': Difference between indices')
    plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="Difference",xlab=" ",
       main=main)
    points(c(xmin-1,xmax+1),c(0,0),type="l")
    points(yrs,dmean,type="l",col='black')
    points(yrs,dlow,type="l",col='black',lty=3)
    points(yrs,dup,type="l",col='black',lty=3)
    legend('topleft', c(tcountry1,tcountry2,'difference') , 
         lty=1, col=c('red', 'green','black'), bty='n', cex=.75)
    #lty=3 gives dotted line type for conf interval
  } #plot if
  #********************************************************************************
  # get differences in first temporal differences
  #********************************************************************************
  nyrs=length(yrs)
  ddmean=rep(NA,nyrs);dmean1=rep(NA,nyrs);dmean2=rep(NA,nyrs)
  ddlow=rep(NA,nyrs);ddup=rep(NA,nyrs)
  for (i in 2:nyrs) {
    im1=i-1
    #nb work with % change
    diff1=(bdata1[,i]-bdata1[,im1])/bdata1[,im1]*100
    diff2=(bdata2[,i]-bdata2[,im1])/bdata2[,im1]*100
    ddiff=diff1-diff2
    ddmean[i]=mean(ddiff)  #mean difference in %chg
    dmean1[i]=mean(diff1)   #mean %chg
    dmean2[i]=mean(diff2)
    ddlow[i]=quantile(ddiff,prob=slow)
    ddup[i]=quantile(ddiff,prob=sup)
  }
  data.frame(dmean1,dmean2,ddmean,ddlow,ddup)
  
  ymin=min(ddlow,dmean1,dmean2,na.rm=TRUE)
  ymax=max(ddup,dmean1,dmean2,na.rm=TRUE)
  main=paste0(title,': difference in % change')
  if ('Diffdiff'%in%plot) {
    plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="% Change/Difference",xlab=" ",
       main=main)
    points(c(xmin-1,xmax+1),c(0,0),type="l")
    points(yrs,dmean1,type='l',col='red')
    points(yrs,dmean2,type='l',col='green')
    points(yrs,ddmean,type="l",col='black',lwd=3)
    points(yrs,ddlow,type="l",col='black',lty=3,lwd=3)
    points(yrs,ddup,type="l",col='black',lty=3,lwd=3)
    legend('bottom', inset=c(0,-0.4),horiz=TRUE, c(tcountry1,tcountry2,'difference') , 
         lty=1, col=c('red', 'green','black'), bty='n', cex=.75,lwd=1,1,3)
  } #end of plot if
} #end of fgdiffplot
printnice<-function(x){
  cat('\n\n')
  print(x)
  cat('\n')
}
