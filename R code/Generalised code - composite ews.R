require(reshape2)
require(parallel)

##########################################################################################
##Lets assume you have a 2 column matrix, the first is of times, the second of counts (do body size shit later)
##########################################################################################
##function to calculate CV


composite_ews<-function(dat, indicators, threshold, plotIt,...){
		
	CV <- function(x, na.rm){
		ave<-mean(x, na.rm=na.rm)
		dev<-sd(x, na.rm=na.rm)
	      	return((dev/ave))
	      }
	      
	##interpolate function
	interp<-function(days, obs){
		int.dat<-as.data.frame(approx(days, obs, n = length(obs), xout=seq(min(days), max(days), 1), method = "linear"))
		names(int.dat)<-c("time", "counts")
		return(int.dat)
	}
	
	##function for rolling mean
	rolling_mean <- function(x){
	     k = length(x);
	    result = rep(0, k);
	    for(i in 1 : k){
	        result[i] <- mean(x[1:i], na.rm=T);
	    }    
	   return(result);
	}
	##function for rolling sd
	rolling_sd <- function(x){
	    k = length(x);
	    result = rep(0, k);
	    for(i in 1 : k){
	        result[i] <- sd(x[1:i], na.rm=T);
	    }    
	    return(result);
	}

	
	############################################################
	##print some warnings
	if(length(dat[,1])<3){
		stop()
		print("Warning: Time series too short")}
	############################################################		
	
	dat<-data.frame("time"=dat[,1], "counts"=dat[,2])

	if(length(which(diff(dat$time)>1))>0){
		#print("Warning: some missing data. Count data interpolated where required")
	  dat<-interp(dat$time, dat$counts)}
    
    ##if there is no variation in the data (i.e. counts are all the same) then cut the data to where variation starts
    dat<-dat[min(which(diff(dat$counts)!=0)):(length(dat$counts)-min(which(diff(dat$counts)!=0))),]
  
  if(length(dat[,1]>2)){
  
	##blank objects to save results
	RES<-NULL
	roll.cv<-NULL
	roll.acf<-NULL
	roll.ar<-NULL
	roll.ar.rr<-NULL
	roll.return.rate<-NULL
	roll.density.ratio<-NULL

#	roll.body<-NULL
#	roll.body.sd<-NULL
	
	inds<-match.arg(indicators, choices=c("cv", "acf", "ar1", "dr", "rr"), several.ok=T)
	
	
	##looped to calculate the rolling change	
	for(i in 2:length(dat$time)){
		#i=3
		##subset the population of interest up until dat i
		dat.t<-subset(dat, time<=unique(sort(dat$time))[i])
		
		##calculate the CV at time t in the focal pop, relative to mean CV through time of that pop
		if(length(which(inds=="cv"))==1){
			roll.cv[[i]]<-CV(dat.t$counts, TRUE)
			cv<-(CV(dat.t$counts, TRUE)-(mean(roll.cv, na.rm=T)))/sd(roll.cv, TRUE)
		}
		## for autocorrelation
		if(length(which(inds=="acf"))==1){
			roll.acf[[i]]<-acf(dat.t$counts, lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
			acf<-(acf(dat.t$counts, lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]-mean(roll.acf, na.rm=T))/sd(roll.acf, TRUE)
		}
		
		## for ar1
		if(length(which(inds=="ar1"))==1){
			if(length(which(diff(dat.t$counts)!=0))>0){
			roll.ar[[i]]<-ar.ols(dat.t$counts, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]
			ar1<-(ar.ols(dat.t$counts, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]-mean(roll.ar, na.rm=T))/sd(roll.ar, TRUE)
			}else{ar1<-NA}
		}
		
		##for density ratio
		if(length(which(inds=="dr"))==1){
			if(length(which(diff(dat.t$counts)>0))>0){
				spectfft <- spec.ar(dat.t$counts, n.freq = length(dat.t$counts), plot = FALSE, order = 1)
				roll.density.ratio[[	i]] <- spectfft$spec[1]/spectfft$spec[length(dat.t$counts)]
				dr<-(roll.density.ratio[[i]]-(mean(roll.density.ratio, na.rm=T)))/sd(roll.density.ratio, TRUE)
				}else{dr<-NA}
		}
		
		##for return rate
		if(length(which(inds=="rr"))==1){
			if(length(which(diff(dat.t$counts)!=0))>0){
			roll.ar.rr[[i]]<-ar.ols(dat.t$counts, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]
			ar.t.rr<-(ar.ols(dat.t$counts, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]-mean(roll.ar.rr, na.rm=T))/sd(roll.ar.rr, TRUE)
			}else{ar.t.rr <-NA}
		    roll.return.rate[[i]]<-1/ar.t.rr
		    rr<-((1/ar.t.rr)-(mean(roll.return.rate, na.rm=T)))/sd(roll.return.rate, TRUE)
		}		
		    	    
		##save results
		RES[[i]]<-data.frame(dat.t[i,]
				,"cv"=if(length(which(inds=="cv"))==1){"cv"=cv}else{NA} 
				,"acf"=if(length(which(inds=="acf"))==1){"acf"=acf}else{NA}
				,"ar1"=if(length(which(inds=="ar1"))==1){"ar1"=ar1}else{NA}	
				,"dr"=if(length(which(inds=="dr"))==1){"dr"=dr}else{NA}
				,"rr"=if(length(which(inds=="rr"))==1){"rr"=rr}else{NA})
				
 		}}
 		results<-do.call("rbind", RES)
 		results$rr<--1*results$rr
 		output<-data.frame(results[,1:2], 
 				"metric.score"=rowSums(results[,3:length(results)], 
 				na.rm=T),
 				"metric.code"=paste(inds, collapse=" + "))
 		
 		output$rolling.mean<-rolling_mean(output$metric.score)
		output$rolling.sd<-rolling_sd(output$metric.score)
		output$threshold.crossed<-NA
		output$threshold.crossed[which(output$metric.score>(output$rolling.mean+(threshold*output$rolling.sd)))]<-1
		
		if(plotIt==T){
			dev.new()
       	 	par(mar = (c(1, 2, 0, 1) + 0.2), oma = c(4, 2, 3, 1))
			plot(output$time, output$rolling.mean, type="l", lwd=1, xlab="Time", ylab="metric score", 
			ylim=c(min(output[,c(3,5,6)], na.rm=T), max(output[,c(3,5,6)], na.rm=T)), lty="solid")
			lines(output$time, output$metric.score, type="l", lwd=2, col="skyblue")
			lines(output$time, output$rolling.mean+(output$rolling.sd*threshold), type="l", lwd=1, lty="dashed")
			lines(output$time, output$rolling.mean-(output$rolling.sd*threshold), type="l", lwd=1, lty="dashed")
		}
		
		return(output)
	}
	