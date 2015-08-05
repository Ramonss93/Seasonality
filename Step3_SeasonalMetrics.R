################################################################
#  Purpose: Extrcting LSP metrics							                 #
#  Project ID:  Shifting Seasons in a Changing World           #
#  Date: 11/11/2014                                            #
#  Author: Simeon Lisovski                                     #
################################################################


#--- initializing global variables ---#

library(zoo)
library(bbmle)
library(raster)
library(doParallel)
library(circular)


load('Data/evi_snow_M.RData')
	load("Data/period_years_M.RData")
  	period_years_M[(round(period_years_M*2,0)/2)==0 | (round(period_years_M*2,0)/2)>1] <- 0
	load("Data/power_years_M.RData")
load('Data/indices.RData')
	index.date <- indices$date

# --- day of the year for each evi raster
	# julian date
	day <- seq(as.POSIXct("1981-01-01"), as.POSIXct("2014-01-05"), by = "day")
		year <- as.numeric(format(day, "%Y"))
		week <- unsplit(lapply(split(year, f = cut(year, breaks = unique(year))), function(x) c(rep(1:53, each = 7))[1:length(x)]),
                              f = cut(year, breaks = unique(year)))

	mean.date <- aggregate(day, by = list(paste(year, week)), median)
  		first.week <- aggregate(data.frame(year, week), by = list(paste(year, week)), function(x) x[1])[,2:3]
     	first.week <- first.week[order(mean.date[,2]),]
      	mean.date <- mean.date[order(mean.date[,2]), 2]
		
	date.tab <- data.frame(year = as.numeric(format(mean.date, "%Y")),
						   week = first.week[,2],
						   mean.date = as.Date(mean.date),
						   jDay = as.numeric(format(mean.date, "%j")))
	date.tab <- merge(index.date, date.tab, by = c("year", "week"), sort = F, all.x = T)					   
						   

## --- aggregation (which raw pixel belong to which aggregated pixel)
	ind01 <- cbind(1:length(indices$wrld_land))
	     t1 <- indices$na_r
	     t1[] <- NA
	     t1[which(indices$wrld_land)] <- 1:sum(indices$wrld_land)
	     ind01 <- cbind(ind01, t1[])
	
	t2 <- indices$na_agg_r
	t2[] <- NA
	t2[which(indices$wrld_land_agg)] <- 1:sum(indices$wrld_land_agg)
	t3 <- resample(t2, indices$na_r, method = "ngb")
	
	t4 <- apply(period_years_M, 1, function(x) !all(is.na(x)))
	t2[which(indices$wrld_land_agg)] <- t4	
	t5 <- resample(t2, indices$na_r, method = "ngb")

	agg.ind <- cbind(ind01, t3[], t5[])
		agg.ind[agg.ind[,4]!=1] <- NA

#_________________________________________________________________________
#------------------------------------------------------------------------#
#--      SECTION 1: Curve fitting and LSP metric extraction    		   --# 
#------------------------------------------------------------------------#


###########################################################
#### Curve functions ######################################
###########################################################
	ads.curve <- function(t, parms) {
		t <- (t-mean(t))
		parms <- as.list(parms)
		(parms$c1) + 
			0.5*(parms$c2)*(tanh(parms$w1*(t - parms$mu)) - 
			tanh(parms$w2*(t - parms$v)))
	}


	ads.loglik <- function(c1, c2, w1, w2, mu, v) {
		fit <- ads.curve(1:nrow(y), parms = list(c1=c1, c2=c2, w1=w1, w2=w2, mu=mu, v=v))		
    	if(ncol(y)>1) {
    		-sum(apply(y, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log=TRUE)), na.rm = T)
    	} else {
    		-sum(dnorm(x = as.vector(y), mean = fit, sd = sd, log=TRUE), na.rm=T)
    		}
	}
		
	
	# fit cos-curve using least-square
	leastS.1 <- function(params, sd) {
			fit  <- params[1]*cos(pi*((1: length(x1))/(length(x1)/((length(x1)/52)*2))) + (pi+params[2]))
			if(is.matrix(Mx)) {
				-sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
				} else {
				-sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
				}
	}	
  	leastS.5 <- function(params, sd) {
    		fit  <- params[1]*cos(pi*((1: length(x1))/(length(x1)/((length(x1)/(52/2))*2))) + (pi+params[2]))
    		-sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }	
###########################################################
###########################################################


## parallel init ####### 
cl <- makeCluster(8)   #
registerDoParallel(cl) #
########################


min  <- matrix(ncol = 32, nrow = nrow(period_years_M))
max  <- min
sos1 <- min
eos1 <- min
sos2 <- min
eos2 <- min


err <- NA
k <- NA

ss <- apply(period_years_M, 1, function(x) {
	x <- round(x*2,0)/2
	x[x==0] <- NA
	if(sum(is.na(x))>10) NA else as.numeric(names(sort(table(x), decreasing = T))[1])
})

	sv <- which(!is.na(ss))[seq(1, sum(!is.na(ss)), length =200)]


for(seas in which(!is.na(ss) & ss==1)) {

cat(".")
if(seas%in%sv) {
	LSP.list <- list(min, max, sos1, eos1, sos2, eos2, err)
	save(LSP.list, file = "Data/LSP_loess23_0.5.RData")
}


if(!any(!is.na(min[seas,]))) {

tryCatch({

	x <- evi_snow_M[agg.ind[which(agg.ind[,3]==seas), 2],]
		 x[x<0] <- NA
		 
	dat <- data.frame(t = rep(1:ncol(x), nrow(x)), y = as.vector(t(x)))
	
  
  	###########################################################
	#### Define split columns and data suitability for fit ####
	###########################################################
	
	x0 <- loess(y ~ t, data = dat, span = 0.02)
	x1 <- predict(x0, newdata = data.frame(t = 1:ncol(x)))
	sd <- sd(abs(dat[!is.na(dat[,2]),2] - predict(x0)), na.rm = T)
	
	## normalize x
	Mx <- x1 - mean(x1, na.rm = T)

	fit0 <- optim(fn = leastS.1, par = c(a = 50, b = 0), sd = 0.001)
	if(fit0$par[1]==0){
		Mx <- x
		fit0 <- optim(fn = leastS.1, par = c(a = 50, b = 0), sd = 0.001)
		Mx <- x1 - mean(x1, na.rm = T)
	}		
	curve.1 <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/52)*2))) + 
                            (pi+fit0$par[2])) +  mean(x1, na.rm=T)
  
    spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0 | 
                   diff(curve.1[-length(curve.1)])>0 & diff(curve.1[-1])<0) +1
    
    out <- data.frame(t = 1:ncol(x), yday = date.tab$jDay, I = x1, cos.01 = curve.1)
    out$split1 <- ifelse(1:nrow(out)%in%spl, 1, 0) # minima/maxima
    out$split2 <- cut(out$t, breaks = c(0, spl, nrow(out)), labels = F) # segments
    
    ## indicate rise and set periods
    out$rise <- as.vector(unlist(apply(cbind(
      as.vector(unlist(lapply(split(out, f=out$split2), function(z) ifelse(coef(lm(z[,"cos.01"]~z[,"t"]))[2] > 0, TRUE, FALSE)))),
      unlist(lapply(split(out, f=out$split2), function(k) nrow(k)))), 1, function(r) rep(r[1], r[2]))))
    
    out$seas <- 1
    out <- cbind(out, t(x))


	## split 1a:
	spl1 <- which(out$split1==1 & out$rise==0)	
	  		tmp1 <- split(out, f = cut(1:nrow(out), breaks = c(0, spl1, nrow(out))))
		  	tmp4 <- split(out, f = cut(1:nrow(out), breaks = c(0, spl1, nrow(out))))					

			if(ss[seas]==1) {
			ytmp0 <- split(index.date, f = cut(1:nrow(out), breaks = c(0, spl1, nrow(out))))	
			ytmp1 <- matrix(unlist(lapply(ytmp0, function(x) cbind(as.numeric(sort(names(table(x[,1])), decreasing = F))[1], 
					 table(x$year)[which.max(table(x$year))], nrow(x)))),
				   			ncol = 3, byrow = T)
			
			
			year.ind <- unlist(lapply(split(data.frame(ytmp1[,2], 1:nrow(ytmp1)), f = cut(ytmp1[,1], breaks = 1981:2013)), 
						  function(x) ifelse(nrow(x)>1, x[which.max(x[,1]),2], x[,2])))	
					
					} else {
			spl1.temp <- which(diff(out$cos.01[-length(out$cos.01)])<0 & diff(out$cos.01[-1])>0) +1				
			ytmp0 <- split(index.date, f = cut(1:nrow(out), breaks = c(0, spl1.temp, nrow(out))))	
			ytmp1 <- matrix(unlist(lapply(ytmp0, function(x) cbind(as.numeric(sort(names(table(x[,1])), decreasing = F))[1], 
					 table(x$year)[which.max(table(x$year))], nrow(x)))),
				   			ncol = 3, byrow = T)
				
			year.ind <- unlist(lapply(split(data.frame(ytmp1[,2], 1:nrow(ytmp1)), f = cut(ytmp1[,1], breaks = 1981:2013)), 
						  function(x) ifelse(nrow(x)>1, x[which.max(x[,1]),2], x[,2])))	
			}		
  			
   			
   			tmp_loess <- foreach(f1=1:length(tmp1), .combine = rbind, .packages = "zoo") %dopar% {

   				z <- tmp1[[f1]]
					
					if((nrow(z)<(median(unlist(lapply(tmp1, nrow)))/3)*2) | 
			 		    sum(is.na(apply(z[,9:ncol(z)], 1, mean, na.rm = T)))>(nrow(z)/4) |
			 		    var(z$I, na.rm = T)<0.001) {
					
					rep(NA, 4) } else {
					
					m1 <- approxfun(x = seq(min(z$I, na.rm = T), max(z$I, na.rm = T), length = 100), y = seq(0, 1, length = 100))
					
					ind0 <- which(m1(z$I)>c(NA,  m1(z$I)[-nrow(z)]) &  m1(z$I)>c(m1(z$I)[-1], NA))
					ind0 <- ind0[which(z[ind0, "I"]>mean(z[,"I"], na.rm = T))]
						if(length(ind0)>1) {
							ind0 <- ind0[order(z[ind0, "I"])][1:2]						
							ind0 <- sort(ind0)
							z[ind0[1]:ind0[2], "I"] <- NA
							z[(ind0[1]-1):(ind0[2]+1), "I"] <- na.approx(z[(ind0[1]-1):(ind0[2]+1), "I"])
						}
	
					diff1 <- c(diff(m1(z$I)), NA)

					sos2 <- z[which.max(diff1), 1]
					eos2 <- z[which.min(diff1), 1]

	 		    	f.tm2 <- approxfun(x = z[,"t"], y = z[,"yday"])
					c(f.tm2(sos2), f.tm2(eos2), sos2, eos2)
					}
		   }		
						
	

		d_loess <- data.frame(x = 1:nrow(tmp_loess), y1 = tmp_loess[,1], y2 = tmp_loess[,2])
		
		circ.year <- approxfun(x = seq(1, 365, length = 360), y = 1:360)
    
    	t1 <- circular(d_loess[,2], type = "angles", units = "degrees")
    	t2 <- attr(median.circular(t1, na.rm = T), which = "medians")[1]
			if(t2>360) t2 <- 360
    	year.day  <- approxfun(x = 1:360, y = seq(1, 365, length = 360))
    
		out$sos_fix <- rep(year.day(t2), sum(ytmp1[,3]))
	
    
    	t1 <- circular(d_loess[,3], type = "angles", units = "degrees")
    	t2 <- attr(median.circular(t1, na.rm = T), which = "medians")[1]
		  	if(t2>360) t2 <- 360
	
		out$eos_fix <- rep(year.day(t2), sum(ytmp1[,3]))
		
	
	## split 1: negative ads curve
	tmp1 <- split(out, f = cut(1:nrow(out), breaks = c(0, spl1, nrow(out))))
			out$fit1 <- NA
			out$fit1_mu <- NA
			out$fit1_v <- NA
   	z.ind <- suppressWarnings(which(!is.na(as.numeric(names(out)))))		
   	
   			
		out[,c("fit1", "fit1_mu", "fit1_v")] <-  matrix(t(unlist(foreach(f1=1:length(tmp1), .packages = c("zoo", "bbmle")) %dopar% {
   			  
   			  z <- tmp1[[f1]]
					
			 		if((nrow(z) < 52) | 
			 		    sum(is.na(apply(z[, z.ind], 1, mean, na.rm = T)))>(nrow(z)/4) |
			 		    (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))<0.05) {
			 			
							t(matrix(nrow = nrow(z), ncol = 3))
							
			 			} else {
			 				
								
			 		    		
			 		    		min.t <- quantile(z[, z.ind], prob = c(0.05), na.rm = T) 
			 		    		max.t <- quantile(z[,z.ind], prob = c(0.95), na.rm = T)

							
								t.temp <- ((1:nrow(z))-mean(1:nrow(z)))
								mf <- approxfun(x = z$yday, y = t.temp)
								
									  sos0 <- unique(z$sos_fix)
									  	if(sos0>360) sos0 <- 360
									  
									  eos0 <- unique(z$eos_fix)
									  	if(eos0>360) eos0 <- 360

								mle <- suppressWarnings(mle2(ads.loglik, method = "L-BFGS-B", 
															 start=list(c1=quantile(z[,"I"], prob = 0.1, na.rm = T), 
																		c2=diff(quantile(z[,"I"], prob = c(0, 1), na.rm = T)), 
																		w1=0.3, w2=0.3, mu=mf(sos0), v=mf(eos0)),
															 upper = list(c1=quantile(z[,3], prob = 0.4, na.rm = T), 
															  	           c2=max.t-min.t, 
															  	           w1 = 2, w2 = 2, mu = mf(sos0)+3, v=mf(eos0)+3),
															 lower = list(c1=quantile(z[,3], prob = 0, na.rm = T)-0.01, 
															  			  c2=diff(quantile(z[,3], prob = c(0, 1), na.rm = T)), 
															  			  w1= 0.1, w2 = 0.1, mu = mf(sos0)-3, v=mf(eos0)-3),	
															 data = list(y=z[,z.ind], sd = 0.01)))


								t1 <- approxfun(y = z[,1], x = ((1:nrow(z))-mean(1:nrow(z))))
								
								mu <- rep(t1(coef(mle)[5]), nrow(z))
								v  <- rep(t1(coef(mle)[6]), nrow(z))

								t(matrix(c(ads.curve(seq(1, nrow(z), by = 1), coef(mle)), mu, v), ncol = 3))
						
						}
		  })), ncol = 3, byrow = T)	
	
	
	## split 2: negative ads curve
	spl2 <- which(out$split1==1 & out$rise==1)	
   			tmp2 <- split(out, f = cut(1:nrow(out), breaks = c(0, spl2, nrow(out))))
			
			out$fit2 <- NA
			out$fit2_mu <- NA
			out$fit2_v <- NA

			
   			out[,c("fit2", "fit2_mu", "fit2_v")] <-  matrix(t(unlist(foreach(f1=1:length(tmp2), .packages = c("zoo", "bbmle")) %dopar% {
   			  
   			  z <- tmp2[[f1]]
			 		
			 		if((nrow(z) < 52) | 
			 		    sum(is.na(apply(z[, z.ind], 1, mean, na.rm = T)))>(nrow(z)/4) |
			 		    (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))<0.05) {


			 		    	t(matrix(nrow = nrow(z), ncol = 3))
			 		    	
			 		    }  else {	
			 		    	
			 		      t.temp <- ((1:nrow(z))-mean(1:nrow(z)))
						  mf <- approxfun(x = z$yday, y = t.temp)
						
						  sos0 <- unique(z$sos_fix)
						  	if(length(sos0)>1) sos0 <- sos0[2]
						  	if(sos0>360) sos0 <- 360
							  
						  eos0 <- unique(z$eos_fix)[1]
						  	if(eos0>360) eos0 <- 360

								mle <- suppressWarnings(mle2(ads.loglik, method = "L-BFGS-B",  
															  start = list(c1=quantile(z[,3], prob = 0.9, na.rm = T), 
																		   c2=-diff(quantile(z[,3], prob = c(0, 1), na.rm = T)), 
																		   w1=0.3, w2=0.3, mu=mf(eos0), v=mf(sos0)),
															  upper = list(c1=quantile(z[,3], prob = 1, na.rm = T)+0.01, 
															  	           c2=-0.01, 
															  	           w1 = 2, w2 = 2, mu= mf(eos0)+3, v=mf(sos0)+3),
															  lower = list(c1=quantile(z[,3], prob = 0.6, na.rm = T), 
															  			   c2=-diff(quantile(z[,3], prob = c(0, 1), na.rm = T)), 
															  			   w1 = 0.1, w2 = 0.1, mu = mf(eos0)-3, v=mf(sos0)-3),	
															  data = list(y=z[,z.ind], sd = 0.01)))
								
								t1 <- approxfun(y = z[,1], x = ((1:nrow(z))-mean(1:nrow(z))))
								
								mu <- rep(t1(coef(mle)[5]), nrow(z))
								v  <- rep(t1(coef(mle)[6]), nrow(z))

								
								t(matrix(c(ads.curve(seq(1, nrow(z), by = 1), coef(mle)), mu, v), ncol = 3))
								
						}})), ncol = 3, byrow = T)

	
	## split 3: smooth transition between split 1 and split 2
	tmp3 <- split(out, f = cut(1:nrow(out), breaks = c(0, which(out$split1==1), nrow(out))))
	
	out$fit3 <- unlist(parLapply(cl, tmp3, function(z) {
	  
	  if(all(is.na(c(z$fit1, z$fit2)))) rep(NA, nrow(z)) else {
	    
	    if(sum(is.na(z$fit1))>(nrow(z)/4)) z$fit2 else {
	      if(sum(is.na(z$fit2))>(nrow(z)/4)) z$fit1 else { 
	        if(sum(z$rise)==0) {
	          ifelse(z$fit1>z$fit2, z$fit1-abs(z$fit1-z$fit2)*seq(0, 1, length=nrow(z)),
	                 z$fit1+abs(z$fit1-z$fit2)*seq(0, 1, length=nrow(z)))
	        } else {
	          ifelse(z$fit2>z$fit1, z$fit2-abs(z$fit2-z$fit1)*seq(0, 1, length=nrow(z)),
	                 z$fit2+abs(z$fit2-z$fit1)*seq(0, 1, length=nrow(z)))
	        }	  					  
	        
	      }}}
	}))


    
    ## extract LSP values
    tmp4 <- split(out, f = cut(1:nrow(out), breaks = c(0, spl1, nrow(out))))					
    tmp5 <- matrix(nrow = 6, ncol = length(tmp4))


		  for(phen in 1:length(tmp4)) {
          z <- tmp4[[phen]]

					
					if(sum(is.na(z$fit1)) > nrow(z)/2 |
					  (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))<0.05) {
			 		    	
			 		    	if(sum(is.na(apply(z[,z.ind], 1, median, na.rm = T)))<(nrow(z)/3) & 
			 		    	   (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))>0.05) {
			 		    		min.t <- quantile(z[,z.ind], prob = c(0.05), na.rm = T) 
			 		    		max.t <- quantile(z[,z.ind], prob = c(0.95), na.rm = T)
			 		    	} else {
			 		    		min.t <- NA 
			 		    		max.t <- NA
			 		    		}
			 		    	
			 		    	tmp5[,phen] <- c(min.t, max.t, rep(NA, 4))
			 		    	
			 		    } else { 
			 		    	
                			f.tm1 <- approxfun(x = 1:nrow(z), y = z[,"t"])
			 		    	f.tm2 <- approxfun(x = z[,"t"], y = z[,"yday"])
			 		    	
			 		    	# amplitude (max, min)
			 		    	if(sum(is.na(apply(z[,z.ind], 1, median, na.rm = T)))<(nrow(z)/3)) {
			 		    	min.t <- quantile(z[,z.ind], prob = c(0.05), na.rm = T) 
			 		    	max.t <- quantile(z[,z.ind], prob = c(0.95), na.rm = T)
			 		    	} else {
			 		    		min.t <- NA 
			 		    		max.t <- NA
			 		    	}
			 		    
			 		    tryCatch({
                 
							if(all(is.na(z$fit3[z$rise==1]))) {
								sos1.t <- NA
								sos2.t <- NA
							} else {
								sos1.t1 <- unique(z$fit1_mu[z$rise==1])
									if(length(sos1.t1)>1) sos1.t1 <- sos1.t1[1]
								sos1.t <- f.tm2(sos1.t1)
							}
							
							if(all(is.na(z$fit3[z$rise==0]))) {
								eos1.t <- NA
								eos2.t <- NA
							} else {
								eos1.t1 <- unique(z$fit1_v[z$rise==0])
									if(length(eos1.t1)>1) eos1.t1 <- eos1.t1[length(eos1.t1)]
								eos1.t <- f.tm2(eos1.t1)
							}
							
							if(!is.na(sos1.t)) {
								ind001 <- c(1:(which.min(abs(z[,1] - sos1.t1))+20))[which.max(
												z$fit3[c(1:(which.min(abs(z[,1] - sos1.t1))+20))])]
								tmp001.1 <- z[1:ind001,]
									f.tmp002 <- approxfun(x = tmp001.1$fit3, y = tmp001.1$t)
								sos2.t <-	f.tm2(f.tmp002(min(z$fit3[1:ind001], na.rm = T) + 
									((max(z$fit3[1:ind001], na.rm = T)-min(z$fit3[1:ind001], na.rm = T))*15)/100))
							}
							
							if(!is.na(eos1.t)) {
								ind001 <- c((which.min(abs(z[,1] - eos1.t1))-20):nrow(z))[which.max(
												z$fit3[(which.min(abs(z[,1] - eos1.t1))-20):nrow(z)])]
								tmp001.1 <- z[ind001:nrow(z),]
							 		f.tmp002 <- approxfun(x = tmp001.1$fit3, y = tmp001.1$t)
								eos2.t <-	f.tm2(f.tmp002(max(z$fit3, na.rm = T) - 
									((max(z$fit3, na.rm = T)-min(z$fit3[ind001:nrow(z)], na.rm = T))*15)/100))		
							  }
			 		    },  error = function(e) {
			 		        sos1.t <<- NA
               			    sos2.t <<- NA
                		    eos1.t <<- NA
                			eos2.t <<- NA})
			
				 tmp5[,phen] <- c(min.t, max.t, sos1.t,  sos2.t, eos1.t, eos2.t)
				 }	
		  	}
	
	min[seas,] <- tmp5[1,year.ind]
	max[seas,] <- tmp5[2,year.ind]
	sos1[seas,] <- tmp5[3,year.ind]	  	
	eos1[seas,] <- tmp5[5,year.ind]
	sos2[seas,] <- tmp5[4,year.ind]
	eos2[seas,] <- tmp5[6,year.ind]



}, error = function(e) {err <<- c(err, seas)})
}

} ## end function


LSP.list <- list(min, max, sos1, eos1, sos2, eos2, err)
# save(LSP.list, file = 'LSP_list.RData')


#_________________________________________________________________________
#------------------------------------------------------------------------#
#--      SECTION 2: Mean of SOS, EOS, Dur, Amp    		   			   --# 
#------------------------------------------------------------------------#


# load('Data/LSP_list.RData')

min <- LSP.list[[1]]
max <- LSP.list[[2]]
sos1 <- LSP.list[[3]]
eos1 <- LSP.list[[4]]
sos2 <- LSP.list[[5]]
eos2 <- LSP.list[[6]]


## Last 10 years (2004-2013)
   ind <- which(c(1982:2014)%in%c(2004:2013))
   
min.med <- indices$na_agg_r
	names(min.med) <- "min"
min.med[which(indices$wrld_land_agg)] <- apply(min[,ind], 1, median, na.rm = T)
# plot(min.med)

max.med <- indices$na_agg_r
	names(max.med) <- "max"
max.med[which(indices$wrld_land_agg)] <- apply(max[,ind], 1, median, na.rm = T)
# plot(max.med)

amp.med <- indices$na_agg_r
	names(amp.med) <- "amp"
amp.med[] <- max.med[]-min.med[]
# plot(amp.med)

	####################################################################
	### Circular statistic functions for SOS, EOS ######################
	####################################################################
	library(circular)
	
    circ.year <- approxfun(x = seq(1, 365, length = 360), y = 1:360)
    year.day  <- approxfun(x = 1:360, y = seq(1, 365, length = 360))
    
    
    circ.median <- function(x) {
    	if(sum(!is.na(x))>3) { 
    	x1 <- circ.year(x)
		cr <- circular(circ.year(x), type = "angles", units = "degrees")
    	year.day(attr(median(cr, na.rm = T), which = "medians")[1])
    	} else {
    		NA
    	}   	
    }
    
    
sos.med <- indices$na_agg_r	
	names(sos.med) <- "sos"
sos.med[which(indices$wrld_land_agg)] <- apply(sos2[,ind], 1, circ.median)

eos.med <- indices$na_agg_r	
	names(eos.med) <- "eos"
eos.med[which(indices$wrld_land_agg)] <- apply(eos2[,ind], 1, circ.median)

    
    
	duration <- function(x) {
  		if(any(is.na(x))) { NA }  else {
   		 	if(x[2]>x[1]) { 
     			 x[2]-x[1] } else {
        		(365-x[1]) + x[2]
      		}
  		}
	}


dur.med <- indices$na_agg_r
	names(dur.med) <- "dur"
dur.med[] <- apply(cbind(sos.med[], eos.med[]), 1, duration)


medianLSP_rs <- stack(min.med, max.med, amp.med,
					  sos.med, eos.med, dur.med)
# save(medianLSP_rs, file = "Data/medianLSP_rs.RData")



#_________________________________________________________________________
#------------------------------------------------------------------------#
#--      SECTION 3: Changes of SOS, EOS, Dur, Amp		   			   --# 
#------------------------------------------------------------------------#
	
	change_MMA <- function(y) {
  	  
  	  if(sum(is.na(y))>15) {
    	 rep(NA, 3) } else {		
			
			y <- y[-length(y)]
		  	x <- 1:length(y)
	  
	      	## remove outliers
        	qnt <- quantile(y, probs = c(0.25, 0.75), na.rm = T)
  			H <- 1.5*IQR(y, na.rm = T)
  			ind <- which(y < (qnt[1]-H) | y > (qnt[2]+H))
        	if(length(ind)>0) {
        		y <- y[-ind]
        		x <- x[-ind]
        	}
		mod <- lm(y~x)
        pr<- predict(mod)
        delta <- pr[length(pr)]-pr[1]
        p <- summary(mod)$coefficients[2,4]
        r <- summary(mod)$r.squared
    	c(delta, r, p)
		}
	}


min.delta <- indices$na_agg_r
min.delta.rs <- indices$na_agg_r
	names(min.delta) <- "min"
	names(min.delta.rs) <- "min"
min.delta[which(indices$wrld_land_agg)] <- t(apply(min, 1, change_MMA))[,1]
min.delta.rs[which(indices$wrld_land_agg)] <- t(apply(min, 1, change_MMA))[,2]


max.delta <- indices$na_agg_r
max.delta.rs <- indices$na_agg_r
	names(max.delta) <- "max"
	names(max.delta.rs) <- "max"
max.delta[which(indices$wrld_land_agg)] <- t(apply(max, 1, change_MMA))[,1]
max.delta.rs[which(indices$wrld_land_agg)] <- t(apply(max, 1, change_MMA))[,2]

	change_Amp <- function(y) {
	
		min1 <- y[1:30]
		max1 <- y[33:(length(y)-2)]
		
		amp <- apply(cbind(min1, max1), 1, function(z) z[2]-z[1])
		x <- 1:length(amp)
		
		if(sum(is.na(amp))>15) {
    	  rep(NA, 3) } else {
    	  	
			## remove outliers
        	qnt <- quantile(amp, probs = c(0.25, 0.75), na.rm = T)
  			H <- 1.5*IQR(amp, na.rm = T)
  			ind <- which(amp < (qnt[1]-H) | amp > (qnt[2]+H))
        	if(length(ind)>0) {
        		amp <- amp[-ind]
        		x <- x[-ind]
        	}
        
        mod <- lm(amp~x)
        pr<- predict(mod)
        delta <- pr[length(pr)]-pr[1]
        p <- summary(mod)$coefficients[2,4]
        r <- summary(mod)$r.squared
    	c(delta, r, p)
    	}
	}


amp.delta <- indices$na_agg_r
amp.delta.rs <- indices$na_agg_r
	names(amp.delta) <- "amp"
	names(amp.delta.rs) <- "amp"
amp.delta[which(indices$wrld_land_agg)] <- t(apply(cbind(min, max), 1, change_Amp))[,1]
amp.delta.rs[which(indices$wrld_land_agg)] <- t(apply(cbind(min, max), 1, change_Amp))[,2]


	
	# Changes in sos
	change_os <- function(y) {
  		
  		if(sum(is.na(y))>15) {
    	 rep(NA, 3) } else {
		 	
		 	y <- y[-length(y)]
		  	x <- 1:length(y)
		
			dist <- (min(y, na.rm = T)-y+180)%%360 -180
			dist <- - ifelse(dist>180, 360-dist, dist)
        
        	## remove outliers
        	qnt <- quantile(dist, probs = c(0.25, 0.75), na.rm = T)
  			H <- 1.5*IQR(dist, na.rm = T)
  			ind <- which(dist < (qnt[1]-H) | dist > (qnt[2]+H))
        	if(length(ind)>0) {
        		dist <- dist[-ind]
        		x <- x[-ind]
        	} 
        	
        mod <- lm(dist~x)
        pr<- predict(mod)
        delta <- pr[length(pr)]-pr[1]
        p <- summary(mod)$coefficients[2,4]
        r <- summary(mod)$r.squared
    	c(delta, r, p)
    	}
	}


### sos
sos.delta <- indices$na_agg_r
sos.delta.rs <- indices$na_agg_r
	names(sos.delta) <- "sos"
	names(sos.delta.rs) <- "sos"
	
sos.delta.m <- matrix(ncol = 3, nrow = nrow(sos2))	
	
	for(i in 1:nrow(sos2)) {
		sos.delta.m[i,] <- change_os(sos2[i,])
	}
	
sos.delta[which(indices$wrld_land_agg)] <- sos.delta.m[,1]
sos.delta.rs[which(indices$wrld_land_agg)] <- sos.delta.m[,2]
	sos.delta[] <- ifelse(sos.delta[]>75 | sos.delta[]< -75, NA, sos.delta[])
	sos.delta.rs[] <- ifelse(sos.delta[]>75 | sos.delta[]< -75, NA, sos.delta.rs[])

## eos
eos.delta <- indices$na_agg_r
eos.delta.rs <- indices$na_agg_r
	names(eos.delta) <- "eos"
	names(eos.delta.rs) <- "eos"

eos.delta.m <- matrix(ncol = 3, nrow = nrow(eos2))	
	
	for(i in 1:nrow(eos2)) {
		eos.delta.m[i,] <- change_os(eos2[i,])
	}

eos.delta[which(indices$wrld_land_agg)] <- eos.delta.m[,1]
eos.delta.rs[which(indices$wrld_land_agg)] <- eos.delta.m[,2]
	eos.delta[] <- ifelse(eos.delta[]>75 | eos.delta[]< -75, NA, eos.delta[])
	eos.delta.rs[] <- ifelse(eos.delta[]>75 | eos.delta[]< -75, NA, eos.delta.rs[])


## dur
dur.delta.m <- t(apply(cbind(sos2, eos2), 1, function(y) {
						
						s <- y[1:30]
						e <- y[33:(length(y)-2)]
						
						x <- 1:length(s)

						if(sum(is.na(s))>15 |sum(is.na(e))>15) {
    	 					rep(NA, 3) } else {
						
						dur1 <- apply(cbind(s, e), 1, duration)
						       
						       ## remove outliers
        						qnt <- quantile(dur1, probs = c(0.25, 0.75), na.rm = T)
  								H <- 1.5*IQR(dur1, na.rm = T)
  								ind <- which(dur1 < (qnt[1]-H) | dur1 > (qnt[2]+H))
        						if(length(ind)>0) {
        							dur1 <- dur1[-ind]
        							x <- x[-ind]
        						} 
						
						mod <- lm(dur1~x)
        				pr<- predict(mod)
        				delta <- pr[length(pr)]-pr[1]
        				p <- summary(mod)$coefficients[2,4]
        				r <- summary(mod)$r.squared
    					c(delta, r, p)
						}
					}))

dur.delta <- indices$na_agg_r
dur.delta.rs <- indices$na_agg_r
	names(dur.delta) <- "dur"
	names(dur.delta.rs) <- "dur"

dur.delta[which(indices$wrld_land_agg)] <- dur.delta.m[,1]
dur.delta.rs[which(indices$wrld_land_agg)] <- dur.delta.m[,2]

# plot(dur.delta)

deltaLSP_rs <- stack(min.delta, max.delta, amp.delta,
					 sos.delta, eos.delta, dur.delta)
# save(deltaLSP_rs, file = "Data/deltaLSP_rs.RData")

rsLSP_rs <- stack(min.delta.rs, max.delta.rs, amp.delta.rs,
					 sos.delta.rs, eos.delta.rs, dur.delta.rs)
# save(rsLSP_rs, file = "Data/rsLSP_rs.RData")