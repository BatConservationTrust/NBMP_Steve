# 16/9/16 mixed mode for daubs based on daubremlcov.gen
timestamp()
options(width=180)

#load functions and data
library(lme4)
#need lmerTest to get p-values for fixed terms
library(lmerTest)
rm(list=ls())  #ensure memory clear
source("C:/data/nbmp2015/R/rcode/func_fgindex.R")
load("hib.rda")  

#********************************************************************************
#set all variables etc
#set below to UK values and overwrite for Wales, Scotland etc if needed
#********************************************************************************
firstyear=1990
minyears=2  #min number of years of valid data to include site
monthstouse=c("Dec","Jan","Feb","Mar")
set.seed(157235) #so can reproduce exactly if required
tsp="daub"
#call count yvar to make easier to change species.  Use eval to reduce potential
#for error since only type species name once above
eval(parse(text=paste("hib$yvar=hib$",tsp)))
eval(parse(text=paste0("hib$ylog=log10(hib$",tsp,"+1)")))
#give full model as character string (converted to full formula later once all vars defined)
charmodel="yvar~site+year+month"

length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#build up subset for analysis
#********************************************************************************
goodyear=as.numeric(as.character(hib$year))>=firstyear
goodcount=!is.na(hib$yvar)
goodmonth=hib$month %in% monthstouse
#note na.rm is needed or gives missing sum if any count missing
sitewith=tapply(hib$yvar,hib$site,sum,na.rm=TRUE)>0
goodsite=sitewith[hib$site]
tapply(goodmonth, hib$month, sum)

hib=hib[goodyear&goodcount&goodmonth&goodsite,]
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset

#********************************************************************************
#further reduce subset, taking just sites with at least minyears of data
#********************************************************************************
tcount=function(x) {length(unique(x))}
nyears=tapply(hib$year,hib$site,tcount)
table(nyears)
hib=hib[(nyears[hib$site]>=minyears),]
length(hib$yvar); sum(hib$yvar,na.rm=TRUE) #for checking correct subset
hib=droplevels(hib)

#********************************************************************************
#fit reml model and check residuals
#********************************************************************************
model=lmer(ylog~year+(1|site/year)+(1|observer),data=hib)
summary(model)
plot(model)  #look for obvious high residuals

#survey level residuals
res=residuals(model)
p99=quantile(abs(res),0.99)
printdf=with(hib,data.frame(site,date,yvar,res))
yrlev=levels(hib$year)
#get list of all sites with residuals above 99th percentile in current year (assume
#those from earlier years already checked)
res_sites=unique(hib$site[(abs(res)>=p99)&(hib$year==max(yrlev))])
nres_sites=length(res_sites)
#now print all data for these sites so can assess if needs checking
for (i in 1:nres_sites){
  cat(paste0("site: ",res_sites[i]),"\n")
  print(printdf[printdf$site==res_sites[i],])}

#now same for site.year residuals
tef=ranef(model)
#this produces a vector with the site.year effects, taken from tef
vef=tef[[1]][[1]]
nsites=length(levels(hib$site))
nyr=length(yrlev)
dd=data.frame(hib[1:(nyr*nsites),1])
#year and site names are in the row labels
synames=rownames(tef[[1]])
Site=as.factor(substring(synames,6,99))
Year=as.factor(as.numeric(substring(synames,1,4)))
mean(levels(Site)==levels(hib$site))  #check levels are in correct order
p99=quantile(abs(vef),0.99)
res_sites=unique(Site[(abs(vef)>=p99)&(Year==max(yrlev))])
nres_sites=length(res_sites)
for (i in 1:nres_sites){
  cat(paste0("site: ",res_sites[i]),"\n")
  print(printdf[printdf$site==res_sites[i],])}

head(synames)
typeof(tef[[1]][[1]])

#********************************************************************************
#now covariate model, with year now as random term only
#********************************************************************************
covform="ylog~(1|year)+(1|site/year)+(1|observer)+month+type+region2"
covmodel=lmer(covform,data=hib)
summary(covmodel)
anova(covmodel)

#********************************************************************************
#following lines provide a quick check on significance of adding extra terms to
#model.  Interactions, quadratics, etc needed to be done individually, if needed.
#********************************************************************************
tvars=c("period","month","week","dayno","exttemp","tempcool","tempwarm","tpc1","tpc2",
        "tpc3","region2")
nvars=length(tvars)
#create structures to save results
Df1=c(rep(NA,nvars)); Df2=c(rep(NA,nvars)); F=c(rep(NA,nvars)); P=c(rep(NA,nvars))
fstats=data.frame(tvars,Df1,Df2,F,P)
for (i in 1:nvars) {
  covformplus=as.formula(paste0(covform,"+",tvars[i]))
#  print(covformplus)
  covmodelplus=lmer(covformplus,data=hib)
  cmp=anova(covmodelplus)
  rnames=row.names(cmp)
  fstats[i,c(2:5)]=cmp[rnames==tvars[i],c(3:6)]
}
#round to avoid excess decimals and exponential notation
fstats$P=round(fstats$P,4);fstats$F=round(fstats$F,2);fstats$Df2=round(fstats$Df2,1)
fstats
