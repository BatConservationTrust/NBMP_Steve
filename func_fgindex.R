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
  vertext="v0.5 14/7/16 (first working version), dataframe construction in fgdata"
  #print(vertext,quote=FALSE)
  #match.arg gets full name if abbreviated in function call
  distribution=match.arg(distribution)
  if (!inherits(formula, "formula"))
    stop("Fault: argument formula needs to be an object of class formula.")
    fgglob.data<<-data.frame(get_all_vars(formula))
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
#print(vertext,quote=FALSE)

#check print arguments
print=match.arg(print,choices=c("summary", "model","bootstrap","none"),
                several.ok = TRUE)

#****************************************************************************
# Basic structures needed later
#****************************************************************************
  spos=2+(fgglob.distribution=="binomial")
  yrlev=as.numeric(levels(fgglob.data[[spos+1]]))
  nyr=length(yrlev)
  #Descs needs to have at least 12 values, if less years just take first nyr values
  desclength=max(nyr,12)
  Descs=character(length=desclength)
  Descs[1]=fgglob.tvar[1]
  Descs[2]=paste("Distribution:     ",fgglob.distribution)
  tform=deparse(fgglob.formula)
  Descs[3]=paste("Model:    x       ",tform)
  #code below allows for multiple splinedf
  Descs[4]=paste("Spline d.f.:      ",paste(splinedf,collapse=","))
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
  fgresults$Nsite=tapply(site,year,tcount)
#fewster default for spline df
  if (splinedf<0) splinedf=floor(0.3*nyr)
  nv=length(fgglob.data[[1]])
  
  if (nmonitor==0) nmonitor=-1  #0 gives fault with %%

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
  tspline=paste("s(vyears, ",format(splinedf),")")
  tform=sub(fgnames[spos+1],tspline,tform)
  fform=as.formula(tform)
  fformglm=fgglob.formula
  
#fit models 
  if (fgglob.distribution=="poisson"){
  gamfit <- gam(fform, family = poisson(link = log),data=fgglob.data,
                offset=offset,weights=weights)  
  glmfit=glm(fformglm,family=quasipoisson(link="log"),data=fgglob.data,
             offset=offset,weights=weights)}

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
    stop  }

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
   # cat("\n**************** GAM Model Summary ****************\n\n")
    fgprint(fgresults,round=TRUE)
    cat(ask)
}
  
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
                ymin=min(fgresults$low_sm),
                ymax=max(fgresults$up_sm)){
  #extract title itself
  mtitle=substring(fgresults$Desc[8],20,99)
  #set up axes and then add lines/points one at a time
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
}
