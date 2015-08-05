#########################################################
#  Purpose: Raw NDVI and snow cover data manipulation   #
#  Project ID:  Shifting Seasons in a Changing World    #
#  Date: 07/01/2015                                     #
#  Author: Simeon Lisovski                              #
#########################################################

## --- initializing global variables ---#

library(raster)
library(rgeos)
library(rgdal)
library(gdalUtils)
library(sp) 
library(maptools)
	data(wrld_simpl)

wrld  <- gIntersection(wrld_simpl, as(extent(-180, 180, 0, 90), 
					   "SpatialPolygons"), byid=T)
wrldP <- spTransform(wrld, CRS("+proj=laea +lat_0=90 +lon_0=0")) 

#_________________________________________________________________________
#------------------------------------------------------------------------#
#--      SECTION 1: NDVI Data      				       --# 
#------------------------------------------------------------------------#


## get file list and date index
fls <- list.files(pattern="*ND.hdf")
	years <- as.numeric(sapply(fls,substring,17,20))
	weeks <- as.numeric(sapply(fls,substring,22,23))
	dupl  <- !duplicated(paste(years, weeks, sep = "_"))
fls <- fls[dupl & years%in%c(1980:2013)]
	years <- as.numeric(sapply(fls,substring,17,20))
	weeks <- as.numeric(sapply(fls,substring,22,23))

## compile raster stack
evi_RS <- stack()

for(i in 1:length(fls)) {
	
	cat("evi manipulation: ", i, " of ", length(fls), "\r")
	
	tmp01 <- gdal_translate(fls[i], "tmp.tiff", 
	                        sd_index = 2, output_Raster = TRUE)[[1]]
	extent(tmp01) <- extent(-180, 180, -55.2, 75.324)
	proj4string(tmp01) <- proj4string(wrld_simpl)

	tmp02 <- suppressWarnings(projectRaster(tmp01, crs = CRS("+proj=moll +over")))	
  	names(tmp02) <- paste("date", years[i], weeks[i], sep = "_")
	
	evi_RS <- stack(evi_RS, tmp02)
}

# save(evi_RS, file = "evi_RS.RData")
# load("Data/evi_RS.RData")


#------------------------------------------------------------------------#
#--      SECTION 2: Snow Cover Data     			       --# 
#------------------------------------------------------------------------#

## get date index of evi raster stack and the index for corresponding snow cover data file
	t_seq <- seq(as.POSIXct("1980-01-01"), as.POSIXct("2014-01-01"), by = "day")
		  t_week <- as.numeric(format(t_seq, "%U"))
	t <- aggregate(t_seq, by = list(week = as.numeric(format(t_seq, "%U")), 
		       year = as.numeric(format(t_seq, "%Y"))), FUN = "mean")

evi_date <- t[match(substring(names(evi_RS), 6, 12), paste(t[,2], t[,1], sep = "_")),3]

## get projection details and one raster layer
	r01 <- evi_RS[[1]]
	r02 <- crop(r01, extent(-18125872, 18124628, 0, 8285575))
	r03 <- projectRaster(r02, crs = proj4string(wrldP))

nlayers_eviRS <- nlayers(evi_RS)
names_eviRS <- names(evi_RS)


## get file list and date index

fls <- list.files(pattern = "*.bin")	
	date1 <-  as.POSIXct(strptime(substring(fls, 21, 28), format = "%Y%m%d", tz = "GMT"))
	date2 <-  as.POSIXct(strptime(substring(fls, 30, 37), format = "%Y%m%d", tz = "GMT"))
date <- date1 + difftime(date1, date2)/2

## compile raster stack
snow_RS <- stack()

for(i in 1:nlayers_eviRS) {

	cat("snow manipulation: ", i, " of ", length(fls), "\r")
	
	ind01 <- which.min(abs(difftime(evi_date[i], date)))
	
	tmp01 <- readBin(fls[ind01], integer(), n = 720*720, size = 1, signed  = T)
	tmp02 <- matrix(tmp01, nrow = 720, ncol = 720, byrow = T)
	tmp03 <- raster(tmp02)
	  extent(tmp03) <- extent(wrldP)
  	  proj4string(tmp03) <- proj4string(wrldP)
	
	tmp04 <- resample(tmp03, r03)
	tmp05 <- projectRaster(tmp04, r01)
	tmp06 <- resample(tmp05, r01)
	
	tmp07 <- r01
	tmp07[] <- ifelse(!is.na(tmp06[]) & tmp06[]>0 & !is.na(tmp07[]) & tmp07[]>(-0.95), TRUE, FALSE)
	
	names(tmp07) <- names_eviRS[i]
	
	snow_RS <- stack(snow_RS, tmp07)
}

# save(snow_RS, file = "snow_RS.RData")
# load("Data/snow_RS.RData")

#-----------------------------------------------------------------------------------------------------------#
#--      SECTION 3: Matrices for NDVI and snow cover time series (rows = raster-cells, cols = time)       --# 
#-----------------------------------------------------------------------------------------------------------#

## index: raster layers
index.date <- data.frame(year = rep(1982:2013, each = 52), 
                        week = rep(1:52, length(1982:2013)))
evi.dat <- substring(names(evi_RS), 6, 12)
ind_evi <- apply(cbind(paste(index.date[,1], index.date[,2], sep = "_")), 1, 
                 FUN = function(x) {
                   out <- which(evi.dat==x[1])
                   ifelse(length(out)>0, out, NA)
                 })

## index: pixel on land
wrld_r <- rasterize(spTransform(wrld_simpl, CRS(proj4string(snow_RS[[1]]))), snow_RS[[1]])


## NDVI matrix
evi_M <- matrix(ncol = length(ind_evi), nrow = sum(!is.na(wrld_r[]))) 

for(i in 1:length(ind_evi)) {
  
  cat("evi matrix: ", i, " of ", length(ind_evi), "\r")
  
  if(!is.na(ind_evi[i])) {
    tmp01 <- evi_RS[[ind_evi[i]]]
    evi_M[,i] <- values(tmp01)[!is.na(wrld_r[])]
  }
}

save(evi_M, file = "evi_M.RData")


## Snow cover Matrix
snow_M <- matrix(ncol = length(ind_evi), nrow = sum(!is.na(wrld_r[]))) 

for(i in 1:length(ind_evi)) {
  
  cat("snow matrix: ", i, " of ", length(ind_evi), "\r")
  
  if(!is.na(ind_evi[i])) {
    tmp01 <- snow_RS[[ind_evi[i]]]
    snow_M[,i] <- values(tmp01)[!is.na(wrld_r[])]
  }
}

# save(snow_M, file = "snow_M.RData")
# load("Data/snow_M.RData")


#---------------------------------------------------------------------------------#
#--      SECTION 4: Snow cover correction and NA interpolation (missing vales)  --# 
#---------------------------------------------------------------------------------#

evi_snow_M <- evi_M
  evi_snow_M[!is.na(evi_snow_M) & evi_snow_M<(-0.25)] <- NA
  evi_snow_M[!is.na(snow_M) & snow_M==1] <- -1


## function: snow correction
snow_correction <- function(x) {
    
    ind01 <- x==(-1)
    ind02 <- ind01 & (c(x[-1], NA))==(-1) & !is.na(c(x[-1], NA))
    ind03 <- ind02 & (c(x[-c(1,2)], NA, NA))==(-1) & !is.na(c(x[-c(1,2)], NA, NA))
    
    ind04 <- ind01 & (c(NA, x[-length(x)]))==(-1) & !is.na(c(NA, x[-length(x)]))
    ind05 <- ind04 & (c(NA, NA, x[-c(length(x)-1,length(x))]))==(-1) & !is.na(c(NA, NA, x[-c(length(x)-1,length(x))]))
    
    tmp01 <- ifelse((ind03 | ind05), 0, x)
    ifelse(tmp01>(-1), tmp01, NA)
    
}


for(i in 1:nrow(evi_snow_M)) {
	cat(i, "    \r")
	x <- evi_snow_M[i,]
	if(any(!is.na(x) & x==(-1))) {
	evi_snow_M[i,] <- snow_correction(x)
	}
}

# save(evi_snow_M, file = "evi_snow_M.RData")
# load("Data/evi_snow_M.RData")



## function: NA interpolation
na.correction <- function(x) {
  
  nap <- sum(is.na(x))/length(x)
  qtl <- quantile(x, na.rm = T, probs = 0.95)
  
  if(nap<=0.8 & qtl>0.05) {
    
    f   <- approxfun(x = c(1:length(x))[!is.na(x)], y = x[!is.na(x)], rule = 2)
    ifelse(!is.na(x), x, f(c(1:length(x))))
    
  } else rep(NA, length(x))
}

evi_snow_na_M <- matrix(ncol = ncol(evi_snow_M), nrow = nrow(evi_snow_M))

for(i in 1:nrow(evi_snow_na_M)) {
	cat(i, "    \r")
	x <- evi_snow_M[i,]
	evi_snow_na_M[i,] <- na.correction(x)
}


# save(evi_snow_na_M, file = "evi_snow_na_M.RData")
# load("Data/evi_snow_na_M.RData")


#------------------------------------------------------------------------#
#--      SECTION 5: Spatial aggregation of snow and na corrected NDVI  --# 
#------------------------------------------------------------------------#

## downscale function
   downscale_evi <- function(x, ...) {
		if(sum(is.na(x))<=(length(x)/3)) median(x, na.rm =T) else NA
	}

## raster stack
evi_snow_na_agg_RS <- stack()

for(i in 1:ncol(evi_snow_na_M)) {
  
  cat(i, "\r")
  
  tmp2 <- wrld_r
  tmp2[] <- NA
  tmp2[][!is.na(wrld_r[])] <- evi_snow_na_M[,i]
  
  tmp3 <- aggregate(tmp2, fact = c(4, 4), downscale_evi)
  evi_snow_na_agg_RS <- stack(evi_snow_na_agg_RS, tmp3)
}

names(evi_snow_na_agg_RS) <- paste("date", index.date[,1], index.date[,2], sep = "_")

# save(evi_snow_na_agg_RS, file = "evi_snow_na_agg_RS.RData")
# load("Data/evi_snow_na_agg_RS.RData")


## time series matrix
wrld_r_agg <- aggregate(wrld_r, fact = c(4,4), mean)

evi_snow_na_agg_M <- matrix(ncol = nlayers(evi_snow_na_agg_RS), nrow = sum(!is.na(wrld_r_agg[])))

for(i in 1:ncol(evi_snow_na_agg_M)) {
  
  cat(i, "  \r")
  evi_snow_na_agg_M[,i] <- values(evi_snow_na_agg_RS[[i]])[!is.na(wrld_r_agg[])]
}

# save(evi_snow_na_agg_M, file = "evi_snow_na_agg_M.RData")
# load("Data/evi_snow_na_agg_M.RData")


#------------------------------------------------------------------------#
#--      SECTION 5: Indices for further processing		       --# 
#------------------------------------------------------------------------#

## pixel on land
wrld_land <- wrld_r
wrld_land[] <- ifelse(!is.na(wrld_land[]), TRUE, FALSE)

## aggregated pixel on land
wrld_land_agg <- aggregate(wrld_r, fact = c(4,4))
  wrld_land_agg[] <- ifelse(!is.na(wrld_land_agg[]), TRUE, FALSE)

## empty raster
na_r <- wrld_land
na_r[] <- NA

## empty aggregated raster
na_agg_r <- wrld_land_agg
na_agg_r[] <- NA

## pixel with vegetation (95% quantile of NDVI>0.25)
veg_agg <- apply(evi_snow_na_agg_M, 1, function(x) ifelse(quantile(x, probs = 0.95, na.rm = T)<0.25, FALSE, TRUE))


indices <- list(date = index.date, wrld_land = wrld_land[], wrld_land_agg = wrld_land_agg[], 
                na_r = na_r, na_agg_r = na_agg_r, veg_agg = veg_agg[])

# save(indices, file = "indices.RData")