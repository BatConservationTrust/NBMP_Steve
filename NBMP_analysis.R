
indsp.func <- function(sp, dfvec = c(4, 7, 10, 15, 20, 33))
{
# indsp.func:
# Steve Langton August 2015 version
# Exactly as standard except for line calling GAM function, which is modified to use
# GAM package, not MGCV

# Calculates the data matrices of indices for different df.
# Arguments are the species data frame sp,  
# with column headers marked "site", "year", and "count", 
# and a list of required degrees of freedom (could be just one) for
# a selection of GAM fits. 
#
# All rows with NA's (missing values) must be removed from the data before
# applying this function.
# Result is a matrix with columns giving indices for each of the Nyears 
# years, with df as marked at the top of the column.
# 
# example: species data frame is called cb (for common bird):
# at command line:   indcb <- indsp.func(cb, c(4, 7, 10, 15, 20, 33))
#
# NOTE: SOME PARTS OF THIS FUNCTION MAY BE SPECIFIC TO THE GAM
# FITTED ON count ~ as.factor(site) + s(year).
# IF USING A DIFFERENT FORMULA, NEED TO MODIFY CODE: ESPECIALLY SEE
# NOTE MARKED *** BELOW.
#
	if(length(sp$count[is.na(sp$count)]) > 0) stop(
			"no missing data allowed")

# fit.func fits the GAM using the MGCV library and extracts the indices for a given 
# value of df:
# August 2015 modified to use GAM library rather than MGCV. Note df definition different.
	fit.func <- function(dfval)
	{
                gam.df <- gam(count~as.factor(site) + s(year, dfval), 
			family = poisson(link = log), data = sp)                  
                pred.df <- predict.gam(gam.df,type="terms")
                srow <- length(pred.df[1,])   
# *** For old MGCV versions, change line above to
#     srow <- length(pred.df[,1])
# srow is the row containing the smooth term in the predict object
#
# *** IF CHANGING THE ORDER OF VARIABLES IN FORMULA, OR FITTING GAM 
# WITH COVARIATES, NEED TO CHECK WHETHER IT IS STILL THE SAME ROW REQUIRED ***
#
                Nentries <- length(sp$count)
                reqd.entries <- (1:Nentries)[!duplicated(sp$year)] #
# reqd.entries gives the entries corresponding to the first occurrence of each
# year in the data frame sp.  The years are not in order, but given by 
# years.entries.  Both are then ordered before being returned to ind.df.  
                years.entries <- sp$year[!duplicated(sp$year)]
                years.pred <- pred.df[reqd.entries, srow]
# *** For old MGCV versions, change line above to
#     years.pred <- pred.df[srow, reqd.entries]
                years.ord <- sort(years.entries)
                pred.ord <- years.pred[order(years.entries)]                
                ind.df <- exp(pred.ord)/exp(pred.ord[1])        #

# indices are scaled around the base year: chosen here as year 1. 
		cat("df ", dfval, " complete \n")
		ind.df
	}
	ind.all <- sapply(dfvec, fit.func)
	dimnames(ind.all) <- list(character(0), paste("df", as.character(dfvec),
		sep = ""))
	ind.all
}

