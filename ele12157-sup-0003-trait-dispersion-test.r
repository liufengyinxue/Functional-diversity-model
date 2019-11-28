## Trait-based tests of community assembly in R

#The R functions posted here are free for you to download and use (under the terms of GNU GPL v2). The code offered is unsupported- I make no claims that they will work with your data or increase your general well being- but I am happy to answer questions. Please consider citing the related papers if you use them, modify them in your own research, or figure out how to solve a problem by looking at them.

## Simplification of functions used in:

#Kraft, N. J. B. and D. D. Ackerly. 2010. Functional trait and phylogenetic tests of community assembly across spatial scales in an Amazonian forest. Ecological Monographs 80:401-422.

#Kraft, N. J. B., R. Valencia, and D. Ackerly. 2008. Functional traits and niche-based tree community assembly in an Amazonian forest. Science 322:580-582.


# Fuctions come first in the file; example is at end. The key function is test_trait_data(); it requires all the other functions to run properly. 

## As a general warning, there is MINIMAL error/ safety checking in this code.  Use the example at the end of the file for formatting your data- I highly recommend that you wrap this code in something else to ensure that only complete datasets go into this function.  Also- Kurtosis is undefined for samples of less than 4, so the minimum community size for these tests is 4.  

## UPDATED June 5 2012 to correct an error in the calculation of SDNDr. 

#####################
### Functions #######
###########################################################



## BASIC TRAIT SPACE FUNCTIONS:

# Function to calculate the unbiased population estimate or the biased sample statistic of kurtosis. from http://www2.jura.uni-hamburg.de/instkrim/kriminologie/Mitarbeiter/Enzmann/Software/kurtosis.r Aug 2 2007


kurtosis=function(x, biased=F, na.rm=T)

{
 if (na.rm==T) x = x[!is.na(x)]
 n = length(x)
 if (n < 4)
 {
   if (na.rm==T)
   {
      cat('valid cases = ',n,'\nkurtosis is not defined for less than 4 valid cases!\n')
   }
   else
   {
      cat('cases = ',n,'\nkurtosis is not defined for less than 4 cases!\n')
   }
 }
 else
 {
   if (biased==T)
   {
      z = sqrt(n/(n-1))*scale(x)
      k = mean(z^4)-3
   }
   else
   {
      z = scale(x)
      k = sum(z^4)*n*(n+1)/((n-1)*(n-2)*(n-3))-3*(n-1)^2/((n-2)*(n-3))
   }
 k
 }
}

get_kurt=function(vect){
  return(kurtosis(vect))
	}

get_range=function(vect){
	r<-range(vect, na.rm=TRUE)
	return(r[2]-r[1])
	}
	
get_spacing=function(vect){
	vect<-vect[!is.na(vect)]
	vect[order(vect)]->v
	n<-length(v)
	spacing=NULL
	summary=NULL
	
	spacing[1]<-abs(v[1]-v[2])
	for (i in 2:(n-1)){
		first<-abs(v[i]-v[i-1])
		second<-abs(v[i]-v[i+1])
		spacing[i]<-min(first, second)
		}
	spacing[(length(spacing)+1)]<-abs(v[n]-v[n-1])
	summary$SDNN<-sd(spacing)
	summary$SDND<-sd(diff(sort(vect)))
	return(summary)
	}

get_var=function(vect){
	return(var(vect, na.rm=TRUE))
	}	
	
	
### Utility functions:
	
trim_pool_to_sample=function(restricted_sample,pool){
	sp<-unique(restricted_sample$sp)
	pool[(pool$sp %in% sp),]->np
	
	return(np)
	
	}


### Function to build a null distribution:

make_null=function(pool,richness,log=TRUE, reps=999, abweight=FALSE, abdata=NULL) {
	##abdata should be a sp x abund DF
	
	summary=NULL
	pool<-pool[!is.na(pool[2]),]
	if(log){
		pool[2]<-log10(pool[2])
		}
	
	if(abweight){
			merge(pool,abdata)->new
			new->pool
			}else{pool$abund<-1}
	

	
	for (i in 1:reps){
		sample(pool[,2], richness, prob=pool$abund)-> simcom
		##trait mean
		summary$mean[i]<-mean(simcom)
		##trait range
		summary$range[i]<-get_range(simcom)
		##spacing stats
		summary$var[i]<-get_var(simcom)
		space<-get_spacing(simcom)
		summary$SDNN[i]<-space$SDNN
		summary$SDNNr[i]<-space$SDNN/summary$range[i]
		summary$SDNDr[i]<-space$SDND/summary$range[i]
		summary$kurt[i]<-get_kurt(simcom)
		
				}
	summary<-as.data.frame(summary)
	return(summary)
	}
	

### Main function that puts it all together:

test_trait_data=function(sp_list, pool, log=TRUE, reps=999,abweight=FALSE, abdata=NULL, verbose=FALSE){
	
	
	community<-pool[pool$sp %in% sp_list,]
	
	dim(community)[1]->richness
	summary=NULL
	summary$test.richness<-richness
	summary$reps<-reps
	
	if(log){community[2]<-log10(community[2])}
	
	make_null(pool,richness, log=log, reps=reps, abweight=abweight, abdata=abdata)->nulldist
	summary$mean_rank<-sum(nulldist$mean<mean(community[,2]))
	summary$range_rank<-sum(nulldist$range<get_range(community[,2]))
	
	obspace<-get_spacing(community[,2])
	summary$SDNN_rank<-	sum(nulldist$SDNN<obspace$SDNN)
	summary$SDNNr_rank<-sum(nulldist$SDNNr<(obspace$SDNN/get_range(community[,2])))
	summary$SDNDr_rank<- sum(nulldist$SDNDr<(obspace$SDND/get_range(community[,2])))
	summary$kurt_rank<-	sum(nulldist$kurt<get_kurt(community[,2]))
	summary$var_rank<-	sum(nulldist$var<get_var(community[,2]))
	summary<-as.data.frame(summary)
	
	

	summary$obs_mean<-mean(community[,2])
	summary$null_mean_mean<-mean(nulldist$mean)
	summary$null_mean_sd<-sd(nulldist$mean)
		
	summary$obs_range<-get_range(community[,2])
	summary$null_mean_range<-mean(nulldist$range)
	summary$null_range_sd<-sd(nulldist$range)
		
	summary$obs_SDNN<-obspace$SDNN
	summary$null_mean_SDNN<-mean(nulldist$SDNN)
	summary$null_SDNN_sd<-sd(nulldist$SDNN)
		
	summary$obs_kurt<-get_kurt(community[,2])
	summary$null_mean_kurt<-mean(nulldist$kurt)
	summary$null_kurt_sd<-sd(nulldist$kurt)
		
	summary$obs_var<-get_var(community[,2])
	summary$null_mean_var<-mean(nulldist$var)
	summary$null_var_sd<-sd(nulldist$var)
		
	summary$obs_SDNNr<-summary$obs_SDNN/summary$obs_range
	summary$null_mean_SDNNr<-mean(nulldist$SDNNr)
	summary$null_SDNNr_sd<-sd(nulldist$SDNNr)
		
	summary$obs_SDNDr<-obspace$SDND/summary$obs_range
	summary$null_mean_SDNDr<-mean(nulldist$SDNDr)
	summary$null_SDNDr_sd<-sd(nulldist$SDNDr)
		
	## quantiles added for plotting:
	summary$null_SDNDr_025<-quantile(nulldist$SDNDr, c(.025))
	summary$null_SDNDr_975<-quantile(nulldist$SDNDr, c(.975))
		
	summary$mean_ES<-(summary$obs_mean-summary$null_mean_mean )/summary$null_mean_sd
	summary$range_ES<-(summary$obs_range-summary$null_mean_range )/summary$null_range_sd
	summary$var_ES<-(summary$obs_var-summary$null_mean_var )/summary$null_var_sd
	summary$SDNN_ES<-(summary$obs_SDNN-summary$null_mean_SDNN )/summary$null_SDNN_sd
	summary$SDNNr_ES<-(summary$obs_SDNNr-summary$null_mean_SDNNr )/summary$null_SDNNr_sd
	summary$SDNDr_ES<-(summary$obs_SDNDr-summary$null_mean_SDNDr )/summary$null_SDNDr_sd
	summary$kurt_ES<-(summary$obs_kurt-summary$null_mean_kurt )/summary$null_kurt_sd
	
	if(verbose){
		return(summary)
		}else{ 
			return(summary[,c("test.richness", "reps", "mean_rank", "mean_ES", "range_rank", "range_ES", "var_rank", "var_ES", "SDNN_rank", "SDNN_ES", "SDNNr_rank", "SDNNr_ES", "SDNDr_rank", "SDNDr_ES", "kurt_rank", "kurt_ES")])  
			
			}
	
}

## Description of arguments for test_trait_data()
#  sp_list  = vector of species names for local community
#  pool		= data frame for species pool, col 1 must be named "sp" and includes species names, used to match to sp_list; col 2 can be called anything and is a continuous trait
#  log		= Logical- should trait data be logged?
#  reps		= number of null communities to build- should be no less than 999
#  abweight = logical- should species be sampled in the null based on abundance? If true, requires a vector for abdata
#  abdata   = vector of the abundances (absolute or relative) of species in the pool, used to weight null model draws.  Must be in same order as species in pool dataframe.  Required if abweight is TRUE.
#  verbose  = Logical.  How much information do you want in the output?  If FALSE, summary includes the ranks of the observed in the null and effect sizes.  If TRUE, summary also includes the observed value, the mean of the null, and the sd of the null for each metric.

## Description of output:
# test.richness = number of species in the community.
# reps 			= number of null randomizations performed
# ..._rank		= rank of the observed metric (e.g. range_rank) in the null.  Used with 'reps' to calculate a p-value.
# ..._ES		= standard effect size for a metric (e.g. range_ES), calculated as (obs-expected)/ (sd of null)
# obs_...		= value of a metric for the real community (e.g. obs_range)
# null_mean_... = mean of the null distribution of a given metric (e.g. null_mean_range)
# null_..._sd	= standard deviation of the null distrbution for a given metric (e.g. null_range_sd)
