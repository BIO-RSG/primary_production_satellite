# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-08-21

# Run primary production using the BIOv2 model.

# ******************************************************************************
#                   Required format of input Rdata file
#
# FOR IN SITU PRIMARY PRODUCTION:
#   tibble with the following columns:
#       event_id, sample_id, station, latitude, longitude, doy, year, par, depth_data, time_data
#           "doy" is the day of year
#           "par" is the surface par (currently we just use NASA OBPG's daily satellite PAR)
#           depth_data must be a list-col where each cell/element is a dataframe containing columns depth (a positive value), alphab, pmb, chl
#           time_data must be a list-col where each cell/element is a vector containing one column ("time"), the time(s) of day to calculate pp - EITHER ONE VALUE (sampling time), OR A LIST OF 23 VALUES (hours 1-23 in HHMM format)
#               "time" must be a 4-digit value: HHMM, AND MUST BE IN UTC
#
# FOR SATELLITE PRIMARY PRODUCTION:
#   tibble with the following columns:
#       bin, latitude, longitude, doy, year, par, depth_data, time_data
# ******************************************************************************

# UNITS:
#   par = einsteins/m^2/day (from NASA OBPG)
#   chl = mg/m^3
#   alphab: C/chl/sec/(watt/m^2) for PP at sampling time, C/chl/hr/(watt/m^2) for hourly PP
#   pmb = C/chl/sec for PP at sampling time, C/chl/hr for hourly PP
#   pp = mgC/sec for PP at sampling time, mgC/hr for hourly PP

# OUTPUT:
#   Rdata file containing the input dataframe, as well as columns containing the depth_data filled to the depths where PP is calculated, PP (a matrix where row=PP at depth, col=PP at time), and scalar irradiance E0z (also a matrix with row=E0z at depth, col=E0z at time)
#   csv file containing the input columns + total PP integrated over time and depth, + total PP integrated over time and depth in the top 30m of the water column


rm(list=ls())

#*******************************************************************************
# LIBRARIES

library(caTools)        # for integrating PP over depth and time using trapz()
library(data.table)     # for reading large csv input files with fast "fread" function
library(dplyr)          # for organizing input and output
library(oceancolouR)    # for various custom functions
library(pbapply)        # for progress bars
library(stringr)
library(Rcpp)
library(RcppArmadillo)
library(fst)
library(lubridate)
dformat <- "%Y%m%d"
source("scripts/model_BIOv2.R")
source("scripts/PAR_resolved.R")
sourceCpp("scripts/surface_irradiance.cpp")
sourceCpp("scripts/underwater_irradiance.cpp")
library(parallel)
num_cl <- NULL # set this to NULL to not use parallel processing


#*******************************************************************************
# VARIABLES THAT CAN BE CHANGED

# water column intervals to calculate PP (metres), and the min and max depth
# depth_data columns will be linearly interpolated to the depths in seq(Zmin,Zmax,by=dz)
dz <- 0.5
Zmin <- 0
Zmax <- 150

# Do a subsetted region? (good for testing)
subset <- FALSE
lat_bounds <- c(43,60)
lon_bounds <- c(-65,-43)

# calculating satellite primary production (or in situ)? (if satellite, you need column "bin" instead of columns "event_id", "sample_id", and "station")
satellite <- FALSE

# uniform chla and PI parameter profiles in water column? (this will use the values closest to the surface in your depth_data)
uniform_profile <- T

# Rdata file containing input data in the format specified at the top of this script, and output filename, EXCLUDING extension:

# IN SITU PP, replacing some values with satellite PP for "validation"
input_file <- "01_PP_input/in_situ/BIOv2/formatted_insitu_hourly_homogprofile_satChl_trophicregimePI_bySST_Q0.75_maxShelfTemp6_par1_UPDATE.Rdata"
output_file <- "02_PP_output/in_situ/BIOv2/formatted_insitu_hourly_satChl_trPIsst_Q0.75_maxtemp6_U_PP_BIOv2_par1_UPDATE"

# # SATELLITE PP:
# mydate <- "20180710"
# input_file <- paste0("01_PP_input/satellite/BIOv2/daily/modisaqua_PP_input_",mydate,".Rdata")
# output_file <- paste0("02_PP_output/satellite/BIOv2/daily/modisaqua_PP_input_",mydate,"_U")



#*******************************************************************************
# MAIN CODE

if (!is.null(num_cl)) {num_cl <- min(num_cl, detectCores()-1)}

dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)

input <- get(load(input_file))

if (subset) {
    input <- input %>%
        dplyr::filter(between(latitude, lat_bounds[1], lat_bounds[2]),
                      between(longitude, lon_bounds[1], lon_bounds[2]))
}

result <- pblapply(1:nrow(input), FUN=calculate_PP, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax, uniform_profile=uniform_profile, cl=num_cl)
# for (i in 1:nrow(input)) { # for diagnosing problem indices
#     test <- calculate_PP(i, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax, uniform_profile=uniform_profile)
# }

# add output to input dataframe, and integrate over time and depth to get summary values
output <- input %>%
    dplyr::mutate(depth_data_filled=lapply(result,"[[","depth_data_filled"),
                  PP=lapply(result,"[[","PP"),
                  E0z=lapply(result,"[[","E0z"))

cat("Input and output combined\n")

output$euphotic_depth <- pblapply(1:nrow(input), function(i) {
    E0z1 <- output$E0z[[i]]
    depth <- output$depth_data_filled[[i]]$depth
    tmp <- sapply(1:ncol(E0z1), function(j) depth[which.min(abs((E0z1[,j]-0.01*E0z1[1,j])))])
    final <- rep(NA,length(tmp))
    final[sapply(tmp,length)!=0] <- unlist(tmp[sapply(tmp,length)!=0])
    return(final)
}, cl=num_cl)
output$euphotic_depth_avg <- sapply(output$euphotic_depth, FUN=mean, na.rm=TRUE)

cat("Euphotic depth extracted\n")

output$PP = pblapply(1:nrow(input), function(i) {
    depths <- output$depth_data_filled[[i]]$depth
    euphotic_depths <- output$euphotic_depth[[i]]
    PP <- output$PP[[i]]
    tinds <- is.finite(euphotic_depths)
    for (j in which(tinds)) {
        PP[depths > euphotic_depths[j],j] <- NA
    }
    return(PP)
}, cl=num_cl)
output$pp_total = pbsapply(1:nrow(input), function(i) {
    depths <- output$depth_data_filled[[i]]$depth
    euphotic_depths <- output$euphotic_depth[[i]]
    PP <- output$PP[[i]]
    tinds <- is.finite(euphotic_depths)
    PP_pertime <- rep(NA,length(tinds))
    # for each time ind, integrate over depth
    for (j in which(tinds)) {
        zinds <- depths <= euphotic_depths[j]
        dsub <- depths[zinds]
        ppsub <- PP[zinds,j]
        PP_pertime[j] <- trapz(dsub,ppsub)
    }
    if (length(tinds)==1) {
        totalPP <- PP_pertime
    } else { # DAILY
        totalPP <- sum(PP_pertime,na.rm=TRUE)
    }
    return(totalPP)
}, cl=num_cl)

cat("PP extracted\n")

output$pp_percent_total_per_depth = pblapply(1:nrow(input), function(i) {
    depths <- output$depth_data_filled[[i]]$depth
    euphotic_depths <- output$euphotic_depth[[i]]
    PP <- output$PP[[i]]
    tinds <- is.finite(euphotic_depths)
    PP_percent_total <- matrix(nrow=nrow(PP),ncol=ncol(PP))
    # for each time ind, integrate over depth
    for (j in which(tinds)) {
        zinds <- depths <= euphotic_depths[j]
        dsub <- depths[zinds]
        ppsub <- PP[zinds,j]
        pptotal <- trapz(dsub,ppsub)
        final <- 100*c(0,sapply(2:length(dsub), function(k) trapz(dsub[1:k],ppsub[1:k])))/pptotal
        PP_percent_total[1:length(final),j] <- final
    }
    return(PP_percent_total)
}, cl=num_cl)

cat("PP % per depth extracted\n")

if (subset) {
    output_file <- paste0(output_file, "_", paste0(lat_bounds,collapse="-"),
                          "N_", paste0(sort(abs(lon_bounds)),collapse="-"),"W")
}
cat(paste0("Writing results to ",output_file,"\n"))

save(output, file=paste0(output_file,".Rdata"), compress=TRUE)

if (satellite) {
    fwrite(output %>% dplyr::select(bin,latitude,longitude,doy,year,par,euphotic_depth_avg,pp_total), file=paste0(output_file, ".csv"))
} else {
    fwrite(output %>% dplyr::select(event_id,sample_id,station,latitude,longitude,doy,year,par,euphotic_depth_avg,pp_total), file=paste0(output_file, ".csv"))
}
