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
  vertext="v0.8 19/09/2016 Sorted bug if splinedf unset" 
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
  Descs[12]=paste("Prog version:      ",vertext)
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
  