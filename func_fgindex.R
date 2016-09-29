fggob.version="v1.0 29/9/16 fgregion added"
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
  vertext=fggob.version
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
    if (distribution=="binomial") 
      stop("Fault: binomial not yet allowed") 
    #need to do checks and preferably allow scalar
    fgglob.distribution<<-distribution
    #position of site variable in dataframe
    spos=2+(distribution=="binomial")
    
    #check sites to make sure factor
    if(is.factor(fgglob.data[[spos]])==FALSE) 
      fgglob.data[[spos]]<-factor(fgglob.data[[spos]])
    #check years to make sure factor
    if(is.factor(fgglob.data[[spos+1]])==FALSE) 
         fgglob.data[[spos+1]]<-factor(fgglob.data[[spos+1]])
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
    fgglob.preddf<<-data.frame(fgglob.data[1:nyr,])
    fgglob.preddf[]<<-fgglob.data[1,]
    fgglob.preddf$vyears<<-yrlev
    fgglob.preddf$year<<-as.factor(yrlev)
    
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

# fgindex function for fitting GAM models
# First use fgdata to define model

#****************************************************************************
fgindex<-function(print="summary",
                  weights,       # default weights all equally
                  splinedf=-99,  # default 0.3*nyears
                  confidence=95,
                  nboot=0,       # number of bootstrap samples
                  nmonitor=0,    # monitoring progress e.g. 5 to show every 5th
                  baseyear=-99,  # -99 sets to second year of series
                  plot=FALSE,
                  title="NULL")  {
vertext="v0.2 18/6/16 with dataframe set up by fginit"
vertext="v0.3 23/6/16 with dataframe for results"
vertext="v0.4 6/7/16 switching to formula in fgdata"
vertext="v0.5 14/7/16 (first working version), dataframe construction in fgdata"
vertext="v0.6 4/8/16 multiple spline df"
vertext="v0.7 12/09/2016 droplevels() in fgdata"
vertext="v0.8 19/09/2016 Sorted bug if splinedf unset"
vertext=fggob.version
#print(vertext,quote=FALSE)

#check print arguments
print=match.arg(print,choices=c("summary", "model","bootstrap","none"),
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
  tform=deparse(fgglob.formula)
  Descs[3]=paste("Model:    x       ",tform)
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
  Descs[10]=fgglob.tvar[2]
  Descs[11]=paste("Total sites:      ",format(nlevels(fgglob.data[[spos]])))
  Descs[12]=paste("Prog version:      ",vertext)
  fgresults=data.frame(Desc=as.character(Descs[1:nyr]),Year=yrlev)
  fgresults$Desc=as.character(fgresults$Desc)
  #,Ncount,Nsite,Rawmean,
   #                    se_mean,Smoothed,se_sm,low_sm,up_sm,Unsmoothed,se_unsm,
    #                   low_unsm,up_unsm,diff,low_diff,up_diff,chg,low_chg,up_chg)

  #fgresults$Desc=Descs[1:nyr]
  fgresults$Nobs=tapply(is.na(fgglob.data[[1]])==FALSE,fgglob.data[[spos+1]],sum)
  fgresults$Mean=tapply(fgglob.data[[1]],fgglob.data[[spos+1]],mean)
  fgresults$SE_mean=tapply(fgglob.data[[1]],fgglob.data[[spos+1]],sd)/sqrt(fgresults$Nobs)
  tcount=function(x) {length(unique(x))}
  #following corrected 12/9/16 to refer to fgglob.data not names
  fgresults$Nsite=tapply(fgglob.data[[spos]],fgglob.data[[spos+1]],tcount)

  nv=length(fgglob.data[[1]])
  
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

# create formulaa for modela
  fgnames=names(fgglob.data)
  tspline=paste("s(vyears, ",format(splinedf[1]),")")
  tform=sub(fgnames[spos+1],tspline,tform)
  fform=as.formula(tform)
  fformglm=fgglob.formula
  
#****************************************************************************
#fit models 
#if multiple df just fit saturated (glmfit) model
#****************************************************************************
  if (fgglob.distribution=="poisson"){
    if (ndf==1) {
     gamfit <- gam(fform, family = poisson(link = log),data=fgglob.data,
                offset=offset,weights=weights)  }
  glmfit=glm(fformglm,family=quasipoisson(link="log"),data=fgglob.data,
             offset=offset,weights=weights)}
  
#****************************************************************************
#  deal with multiple splines - fit spline models to saturated estimates
#****************************************************************************
  if (ndf>1){
    predglm <- predict.glm(glmfit,type="response",newdata=fgglob.preddf)
    aic=c(rep(NA,ndf))
    #structures to export results
    #will fail if n splinedf>n years, but this would be silly anyway
    fgresults$df=c(splinedf,rep(NA,(desclength-ndf)))
    fgresults$aic=c(rep(NA,desclength))
    #loop through fitting spline with required df to pred vals from saturated model
    for (i in 1:ndf) {
      fitspl=gam(predglm~s(fgglob.preddf$vyears,splinedf[i]))
      aic[i]=summary(fitspl)$aic
      fgresults$aic[i]=aic[i]
      fgresults[,i+8]=predict(fitspl)
    } #end of for loop
    if (any(print=="summary")){
      cat("\n\n**** Comparing GAMs with different df ****\n\n")
      cat(paste(fgresults[c(8,1,10,2,3,4,9,7,6,11,12,13),1],"\n"))
      cat("\n\n**** AIC for different d.f. ****\n\n")
    #print as data frame to get in columns
      print(as.data.frame(list(splinedf,aic),col.names=c("d.f.","AIC")))
      }#printing if

#exiting function using stop.  This is not ideal since throws a fault, but
#saves the complication of putting subsequent programming in an if statement
    return(fgresults)
    stop("No bootstrapping with multiple spline d.f.",call.=FALSE)
  }  #end of multiple splinedf if()

  predgam <- predict.gam(gamfit,type="response",newdata=fgglob.preddf)
  baseyrno=which(yrlev==fgglob.baseyear)
  fgresults$Smoothed=predgam/predgam[baseyrno]*100
  fgresults$se_sm=c(rep(NA,nyr))
  fgresults$low_sm=c(rep(NA,nyr))
  fgresults$up_sm=c(rep(NA,nyr))

  predglm <- predict.glm(glmfit,type="response",newdata=fgglob.preddf)
  fgresults$Unsmoothed=predglm/predglm[baseyrno]*100
  fgresults$se_unsm=c(rep(NA,nyr))
  fgresults$low_unsm=c(rep(NA,nyr))
  fgresults$up_unsm=c(rep(NA,nyr))
  
  ask="\n*******************************************************\n\n"
  if (any(print=="model")){
    cat("\n**************** GAM Model ****************\n\n")
 #   print(summary(fgglob.data))
    print(summary(gamfit))
 #   print(data.frame(yrlev,predgam))
    cat(ask)
  }
  #****************************************************************************
  # On to bootstrapping if required
  #****************************************************************************
  if (nboot<1){
    return(fgresults)
    stop()  }

#structures for bootstrapping
  siteunits=split(1:nv,fgglob.data[[spos]]) #unitnos for each site
  nrowsites=c(table(fgglob.data[[spos]]))
  nsites=length(levels(fgglob.data[[spos]]))
  #need to change sites in preddf so uses new site numbers 1:nsites not original levels
  fgglob.preddf[[spos]]=factor(rep(1,nyr),levels=c(1:nsites)) 
  bootpred=matrix(nrow=nboot,ncol=nyr)
  
  time0=Sys.time()
  for (i in 1:nboot) {
    bootsites=sample(1:nsites,nsites,replace=TRUE) #sites in bootstrap sample
    bootunits=unlist(siteunits[bootsites]) #append the different sets of units
    bsites=rep(1:nsites,times=nrowsites[bootsites])
    newdf=fgglob.data[bootunits,]
#overwrite sites with new site nos, rather than site data has come from
    newdf[[spos]]=as.factor(bsites)

    bootgam <- gam(fform, family = poisson(link = log),data=newdf,
                offset=offset,weights=weights)  

    bootpred[i,] <- predict.gam(bootgam,type="response",newdata=fgglob.preddf)
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
  #****************************************************************************
  # Change and difference
  #****************************************************************************
  fgresults$low_dif=c(rep(NA,nyr))
  fgresults$up_dif=c(rep(NA,nyr))
  fgresults$diff=c(rep(NA,nyr))
  fgresults$sig_dif=c(rep(NA,nyr))
  fgresults$chg=c(rep(NA,nyr))
  fgresults$low_chg=c(rep(NA,nyr))
  fgresults$up_chg=c(rep(NA,nyr))
  fgresults$sig_chg=c(rep(NA,nyr))
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
  #temporarily write bootpred as global to help programming
  bootpred<<-bootpred
  nyr<<-nyr
  
  if (any(print=="summary")){
    fgprint(fgresults,round=TRUE)
    cat(ask)
}
  if (plot==TRUE){fgplot(fgresults)  } 
  
  return(fgresults) }

#****************************************************************************
# fgprint
#**************************************************************************** 
fgprint<-function(fgresults,round=TRUE){
  printdf=fgresults
  if (round==TRUE) {
    dec=c(NA,0,0,1,2,0,2,3,2,2 ,2,3,2,2 ,2,2,2,NA ,2,2,2,NA)
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
  print(printdf[,c(2,3,6,7,8,9,10,11,12)])
  
  cat("\n\n**** Differences and change points ****\n\n")
  print(printdf[,c(2,17,15,16,18,19,20,21,22)])
}

#****************************************************************************
# fgplot
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
   parold=par(mfrow=c(2, 1),mar=c(3,4,2,2),oma=c(0,0,0,0))
   
#    ymin=min(fgresults$aic,na.rm=TRUE)
#    ymax=max(fgresults$aic,na.rm=TRUE)
    splinedf=fgresults$df[which(is.na(fgresults$df)==0)]
    ndf=length(splinedf)
    plot(fgresults$df,fgresults$aic,ylab="AIC",xlab="Spline d.f.",main=
           "Comparing AIC of spline models")
    
#and plot fits
    nc=ncol(fgresults)
    matplot(fgresults$Year,fgresults[,9:nc],type="l",lty=1,col=c(9:nc),ylab="index",
            xlab=" ",main="Comparing fits")
    legend("right",paste(format(splinedf)," d.f."),cex=0.6,lty=1,lwd=1,col=c(9:nc))
#    ymin=min(unlist(fgresults[,9:nc]),na.rm=TRUE)
#    ymax=max(unlist(fgresults[,9:nc]),na.rm=TRUE)
#    plot(c(xmin,xmax),c(ymin,ymax), type="n",ylab="Index",xlab=" ",main=mtitle)
#    for (i in 1:ndf) {points(fgresults$Year,fgresults[,i+8],type="l",col=i)}
    par(parold)
  }
  else {
  #set up axes and then add lines/points one at a time
  ymin=min(fgresults$low_sm,na.rm=TRUE)
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
          "red triangles are where difference between years is sig"),side=1,line=c(2,3))
  } #end of plotting for full results
}
fgregion<-function(print="summary",
                   weights,       # default weights all equally
                   splinedf=-99,  # default 0.3*nyears
                   confidence=95,
                   nrand=100,       # number of randomisations
                   ncheck=0,        # number at which to check for early halt
                   nmonitor=0,    # monitoring progress e.g. 5 to show every 5th
                   baseyear=-99,  # -99 sets to second year of series
                   region,        # regions to check for differences
                   plot=FALSE,
                   title="NULL")  {
  #edited from v0.8 of fgindex
  vertext=fggob.version
  #print(vertext,quote=FALSE)
  
  #check print arguments
  print=match.arg(print,choices=c("summary", "model","rand","none"),
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
  tform=deparse(fgglob.formula)
  Descs[3]=paste("Model:    x       ",tform)
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
  Descs[10]=fgglob.tvar[2]
  Descs[11]=paste("Total sites:      ",format(nlevels(fgglob.data[[spos]])))
  Descs[12]=paste("Prog version:     ",vertext)
  #,Ncount,Nsite,Rawmean,
  #                    se_mean,Smoothed,se_sm,low_sm,up_sm,Unsmoothed,se_unsm,
  #                   low_unsm,up_unsm,diff,low_diff,up_diff,chg,low_chg,up_chg)
  
  nv=length(fgglob.data[[1]])
  
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
  if (sum(tapply(as.numeric(region),fgglob.data[spos],var))>0){
    stop("Each site can only be in one region.")  }
  devsep=c(rep(NA,nrand));dfsep=c(rep(NA,nrand))
  devcom=c(rep(NA,nrand));dfcom=c(rep(NA,nrand))
  ask = "\n*******************************************************\n\n"
  
  for (j in 1:nrand) {
    if ((j%%nmonitor)==0) {
      ttime=format(Sys.time(), "%d %b %Y %X %Z")
      cat(paste("********** loop:",j," ",ttime," **********\n"))}
    #****************************************************************************
    #fit overall model
    #****************************************************************************
    if (fgglob.distribution == "poisson") {
      gamfit <-gam(fformreg,family = poisson(link = log),
                   data = fgglob.data,offset = offset,weights = weights)  }
    
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
      subfit = gam(fform,family = poisson(link = log),
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
      if (ngt>(nrand*0.05)) {
        nrand=j
        break}
    }  #ncheck if
    #********************************************************************************
    # finally randomise regions to sites for next loop
    #********************************************************************************
    rsitereg=sample(sitereg) #sample with no args set randomises data
    region=reglevs[rsitereg[fgglob.data[[spos]]]]
  } #nrand loop j
  devdif=devcom-devsep
  dfdif=dfcom-dfsep
  dfint=dfcom-dfsep  #df for interaction
  #asymptotic deviance test
  preal=pchisq(devdif[1],dfdif[1],lower.tail = FALSE)
  prand=mean(devdif>=devdif[1],na.rm = TRUE)
  plist=list(preal,prand)
  
  cat("\n\n**** GAM Regional Analysis ****\n\n")
  #use which to remove empty rows in Desc - or instead specify directly
  #cat(paste(printdf[which(printdf$Desc!=""),1],"\n"))
  printdf=data.frame(Descs)
  cat(paste(printdf[c(8,1,10,2,3,4,9,7,6,11,12,13),1],"\n"))
  print(yeartab)
  
  cat("\n**************** Deviances for observed data ****************\n")
  cat("\ncommon dev:     ",devcom[1],"\ncommon df:      ",dfcom[1])
  cat("\nseparate dev:   ",devsep[1],"\nseparate df:    ",dfsep[1])
  cat("\ndifference dev:   ",devdif[1],"\ndifference df:    ",dfint[1])
  cat("\n\n**************** Deviance test for interation ****************\n")
  cat("\nAsymptotic P-value: ",round(preal,3))
  cat("\nRandomisation P-value: ",round(prand,3),"\n")
  quants=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)
  devquants=quantile(devdif,probs=quants,na.rm = TRUE)
  print(data.frame(quants,devquants))
  
  #check for fails where df for interaction wrong
  #note that df not always exact integer so dont just use !=
  nfail=sum((abs(dfint-dfint[1]))>0.1,na.rm = TRUE)
  pfail=nfail/nrand
  if (pfail>0.05){
    table(dfint)
    cat("\nRandomisations fail: ",nfail," due to incorrect df")
    stop("Too many randomisations fail.")  }
  if (nfail>0){
    table(dfint)
    cat("\nWARNING: ",nfail,"randomisations fail due to incorrect df") }
  
  return(list(devcom,devsep,dfcom,dfsep,Descs,plist))
}

