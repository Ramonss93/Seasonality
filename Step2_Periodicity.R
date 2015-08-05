#########################################################
#  Purpose: Cyclicity pattern in NDVI data              #
#  Project ID:  Shifting Seasons in a Changing World    #
#  Date: 17/09/2014                                     #
#  Author: Simeon Lisovski                              #
#########################################################


#--- initializing global variables ---#

library(biwavelet)
library(maptools)
	data(wrld_simpl)

load('Data/evi_snow_na_agg_M.RData')
load('Data/indices.RData')
  index.date <- indices$date

source('Map/map.R')

#_________________________________________________________________________
#------------------------------------------------------------------------#
#--  SECTION 1: Overall cyclicity - and from different time periods)   --# 
#------------------------------------------------------------------------#
	
per <- matrix(rep(1983:2012, each = 2), 
              ncol = 2, byrow = T)
       per <- cbind(per, c(3*c(1:nrow(per))-2), 3*c(1:nrow(per)))


## low NDVI value correction
evi_snow_na_agg_M[which(!indices[[6]]),] <- rep(NA, ncol(evi_snow_na_agg_M))

wavelet <- function(i) {

tryCatch({
		out1 <- rep(NA, nrow(per)*3)
		out2 <- rep(NA, nrow(per)*3)
		out3 <- rep(NA, nrow(per)*3)
			
		x <- evi_snow_na_agg_M[i,]
		if(!all(is.na(x))) {
	  	wt  <- wt(cbind(1:length(x),x))
	  	power  <- log2(wt$power.corr)
     
			time   <- wt$t 
			period <- wt$period/52

		for(j in 1:nrow(per)) {
		
	    tmp01.pow <- apply(power[ ,which(index.date[,1]%in%c(per[j,1]:per[j,2]))], 1, median, na.rm = T) 
		  		  
		  tmp01.sig  <- apply(wt$signif[ ,which(index.date[,1]%in%c(per[j,1]:per[j,2]))], 1, median, na.rm = T)
			# tmp01.pow[tmp01.sig<1] <- NA
			
      	ind <- c(1, which((!is.na(tmp01.pow[-length(tmp01.pow)]) &  is.na(tmp01.pow[-1])) |
                       ( is.na(tmp01.pow[-length(tmp01.pow)]) & !is.na(tmp01.pow[-1]))), length(tmp01.pow))
      
      	tmp02 <- split(data.frame(period, tmp01.pow, tmp01.sig), f = cut(1:length(tmp01.pow), breaks = unique(ind)))
      
        	pck.periods <- function(p) {
          	if(sum(!is.na(p[,2]))>1) {
            		ind01 <- c(ifelse((p[1, 2]>p[2, 2]), 1, NA), which((p[,2]>c(NA, p[-nrow(p),2]) & c(p[-1,2], NA)<p[,2])),
                       ifelse(p[nrow(p),2]>p[(nrow(p)-1),2], nrow(p), NA))
            	p[ind01[!is.na(ind01)],]
          	} else p[which.max(p[,2]),]
        	}
	    
      	tmp04 <- do.call(rbind, lapply(tmp02, pck.periods))[,1:3]

		
			if(!is.null(tmp04) & any(tmp04[,1]>0.85 & tmp04[,1]<1.15)) {
				ind02 <- which(tmp04[,1]>0.85 & tmp04[,1]<1.15)
					if(length(ind02)>1) ind02 <- ind02[which.min(abs(1 - tmp04[ind02, 1]))]
				
				out2[per[j, 3]+1] <- tmp04[ind02, 2]
				out1[per[j, 3]+1] <- tmp04[ind02, 1]
				out3[per[j, 3]+1] <- ifelse(tmp04[ind02, 3]>1, 1, 0)
				
				tmp04 <- tmp04[-ind02,] 
			}
		
			if(!is.null(tmp04) & nrow(tmp04)>0) {
				if(any(tmp04[,1]>0.35 & tmp04[,1]<0.65)) {
					ind03 <- which(tmp04[,1]>0.35 & tmp04[,1]<0.65)
						if(length(ind03)>1) ind03 <- ind03[which.min(abs(0.5 - tmp04[ind03, 1]))]
					
				out2[per[j, 3]] <- tmp04[ind03, 2]
				out1[per[j, 3]] <- tmp04[ind03, 1]
				out3[per[j, 3]] <- ifelse(tmp04[ind03, 3]>1, 1, 0)
				
				tmp04 <- tmp04[-ind03,] 	
				}
			
				if(nrow(tmp04)>0) {
					
					out2[per[j, 3]+2] <- tmp04[which.max(tmp04[,2]), 2]
					out1[per[j, 3]+2] <- tmp04[which.max(tmp04[,2]), 1]
					out3[per[j, 3]+2] <- ifelse(tmp04[which.max(tmp04[,2]), 3]>1, 1, 0)

				}
			}	
		} 
		
    }},
		error = function(x) err <<- c(err, i))
		c(out1, out2, out3)
}
		

err <- c()
all <- do.call("rbind", parallel:::mclapply(cbind(1:nrow(evi_snow_na_agg_M)), wavelet, mc.cores = 8))


period_M <- all[,c(1:90)]
power_M  <- all[,c(91:180)]
signif_M <- all[,c(181:ncol(all))]

# save(period_M, file = "Data/period_M.RData")
# save(power_M, file = "Data/power_M.RData")
# save(signif_M, file = "Data/signif_M.RData")
# load("Data/period_M.RData")
# load("Data/power_M.RData")
# load("Data/signif_M.RData")



#------------------------------------------------------------------------#
#--  SECTION 2: Matrix: Season per grid cell per year                  --# 
#------------------------------------------------------------------------#
	

  seas_year <- function(x) {
    if(any(!is.na(x))) {
      
      x1 <- x[1:90]
      x2 <- x[91:length(x)]
           
           period <- matrix(x1, ncol = 3, byrow = T) 
           power  <- matrix(x2, ncol = 3, byrow = T)
        
           period[apply(power, 1, function(x) all(is.na(x))), 1] <- 0
           power[apply(power, 1, function(x) all(is.na(x))), 1] <- 0
      
     c(apply(cbind(1:nrow(power), apply(power, 1, which.max)), 1, function(c) period[c[1],c[2]]),
       apply(cbind(1:nrow(power), apply(power, 1, which.max)), 1, function(c) power[c[1],c[2]]))
    } else rep(NA, 60) 
	}


all_years <- t(apply(cbind(period_M, power_M), 1, seas_year))

period_years_M <- all_years[,1:30]
power_years_M <- all_years[,31:ncol(all_years)]

# save(period_years_M, file = "Data/period_years_M.RData")
# save(power_years_M , file = "Data/power_years_M.RData")
# load("Data/period_years_M.RData")
# load("Data/power_years_M.RData")




#------------------------------------------------------------------------#
#--  SECTION 3: Change in Power                                        --# 
#------------------------------------------------------------------------#

power_change <- function(y) {
  if(any(!is.na(y))) {
	
	y <- y[1:28]
	  
    x <- 1:length(y)
    
    mod <- lm(y ~ x)

      p <- summary(mod)$coefficients[2, 4]
      r <- summary(mod)$r.squared
      delta <- predict(mod)[length(x)] - 
               predict(mod)[1]
    c(delta, r, p)
    } else rep(NA, 3)  
}


power_change_M <- t(apply(power_years_M, 1, power_change))
# save(power_change_M , file = "Data/power_change_M.RData")
# load("Data/power_change_M.RData")


#------------------------------------------------------------------------#
#--  SECTION 4: Change in Periodicity                                  --# 
#------------------------------------------------------------------------#
period_change_M <- matrix(ncol = 4, nrow = nrow(period_years_M))



period_change <- function(x) {
  if(all(is.na(x))) {
    out <- rep(NA, 4)
  } else {
    period <- x[1:30]
      period[is.na(period)] <- 0
      period <- round(period*2,0)/2
    power  <- x[31:length(x)]
      power <- abs(min(power))+power
	
	if(all(period==0)) {
		out <- rep(NA, 4)
	} else {
    a <- approxfun(x=seq(0, max(period, na.rm  =T), length = 10), y = seq(0, 1, length = 10))
    b <- approxfun(x=seq(0, 1, length = 10), y = seq(0, max(period, na.rm  =T), length = 10))
         period.a <- a(period)
    
    if(length(unique(period.a))==1) {
      out <- c(unique(period), unique(period), 1, sum(power)) 
    } else {
      dat <- data.frame(t = 1:length( period.a), p =  period.a)

      suppressWarnings(mod <- glm(p~t, data = dat, family = binomial, weights=power))
      # lines(dat[,1], predict(mod, type = "resp"), col = "red", lwd = 2)ß
      start <- round(b(as.numeric(predict(mod, type = "resp")[1]))*2, 0)/2
      end   <- round(b(as.numeric(predict(mod, type = "resp")[nrow(dat)]))*2, 0)/2
    
      # fit <- markovchainFit(period)$estimate@transitionMatrix
      out <- c(start, end, 1-(sum(abs(diff(period))>0)/length(period)), sum(power))
    }
    }
    
  }
out
}


for(i in 1:nrow(period_change_M)) {
	x <- cbind(period_years_M, power_years_M)[i,]
	period_change_M[i,] <- period_change(x)
}

period_change_M <- t(apply(cbind(period_years_M[1:5,], power_years_M[1:5,]), 1, period_change))
# save(period_change_M , file = "Data/period_change_M.RData")
# load("Data/period_change_M.RData")
