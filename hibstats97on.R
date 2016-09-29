# 21/7/16 stats on hibernation data
options(width=180)
load("hib.rda")
monthstouse=c("Dec","Jan","Feb","Mar")

#dataset doesn't contain total numbers
hib$totbats=hib$whiskbrandt+hib$daub+hib$natt+hib$blongear+
    hib$gthorse+hib$lshorse+hib$barb+hib$allpip

#get initial stats on full dataset
tapply(hib$totbats,hib$year,sum,na.rm=TRUE)
tapply(!is.na(hib$totbats),hib$year,sum)

#**************************************************************************** 
#correlation matrix.  Use with() to avoid typing data frame name for each variable
#**************************************************************************** 
with(hib,round(cor(data.frame(whiskbrandt,daub,natt,blongear,gthorse,lshorse,barb,
                       allpip,totbats),use="complete"),3))

#**************************************************************************** 
#apply subsets for years, missing counts and months
#**************************************************************************** 
goodyear=as.numeric(as.character(hib$year))>1996
goodcount=!is.na(hib$totbats)

goodmonth=hib$month %in% monthstouse
tapply(goodmonth, hib$month, sum)

sum(hib$totbats,na.rm=TRUE)
hib=hib[goodyear&goodcount&goodmonth,]
table(hib$year)
sum(hib$totbats,na.rm=TRUE)

#**************************************************************************** 
#count number of unique years per site & plot in geog pos
#**************************************************************************** 
tcount=function(x) {length(unique(x))}
nyears=tapply(hib$year,hib$site,tcount)
grp=cut(nyears,c(0.5,2.5,5.5,10.5,100),lab=c("1 or 2 yrs","3-5 yrs","6-10 yrs",
                                             "11+ years"))
table(grp)
#would be good to add outline of coast, as in Genstat
plot(hib$east,hib$north,pch=unclass(grp))

#**************************************************************************** 
#now construct main table of numbers of each species
#**************************************************************************** 
ngoodsite=sum(nyears>0,na.rm=TRUE)
lgoodspp=with(hib,
          data.frame(whiskbrandt,daub,natt,blongear,gthorse,lshorse,barb,allpip))
tgoodspp=names(lgoodspp)
ngoodspp=length(tgoodspp)
title=paste("\nspecies     nwith  pwith (out of",ngoodsite,"sites)\n")
for (i in 1:ngoodspp) {
  if (i==1) cat(title)
  nbats=tapply(lgoodspp[,i],hib$site,sum)
  nwith=sum(nbats>0,na.rm=TRUE)
  pwith=nwith/ngoodsite
  cat(format(tgoodspp[i],width=12),format(nwith,width=4),round(pwith,4),"\n")
  }

#**************************************************************************** 
#types x years table - note this is number of surveys, not sites
#**************************************************************************** 
hib$year=droplevels(hib$year)  #removes unused year levels
tty=table(hib$type,hib$year)
(ttym=addmargins(tty))
round(prop.table(addmargins(tty,2),2),3)

#**************************************************************************** 
#number of years for each site (calculated above for plot)
#**************************************************************************** 
table(nyears)

#**************************************************************************** 
#number of sites per observer
#**************************************************************************** 
tnobs=tapply(hib$site,hib$observer,tcount)
table(tnobs)

#**************************************************************************** 
#number of sites per year
#**************************************************************************** 
tapply(hib$site,hib$year,tcount)

#**************************************************************************** 
#matrix of sites - this is important as it shows continuity in sites counted over
#time.
#R doesn't have symmetric matrix structure, but just use lower triangle of square
#array
#**************************************************************************** 
yrlev=levels(hib$year)
ny=length(yrlev)
ny2=ny*ny
sit=list() #list of sites with data in each year
nboth=array(c(rep(NA,ny*ny)),dim=c(ny,ny),dimnames=list(yrlev,yrlev))
for (i in 1:ny){
  sit[[i]]=unique(hib$site[hib$year==yrlev[i]])
  for (j in 1:ny)if (i>=j) nboth[i,j]=sum(sit[[i]] %in% sit[[j]])
}
nboth
