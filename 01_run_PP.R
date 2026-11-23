# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-08-21

# Run primary production using the BIOv2 model.

# clear the environment
rm(list=ls())

#*******************************************************************************
# INPUT
#
# tibble in an Rdata file with the following columns:
#       latitude, longitude, doy, year, par, bathymetry, depth_data, time_data
#
# latitude and longitude = decimal degrees
# doy = day of year
# par = surface par, einsteins/m^2/day
# bathymetry = positive value, metres
# depth_data = list-col where each cell/element is a dataframe containing columns:
#     - depth (positive value, metres)
#     - alphab (C/chl/sec/(watt/m^2) for PP at sampling time, C/chl/hr/(watt/m^2) for hourly PP)
#     - pmb (C/chl/sec for PP at sampling time, C/chl/hr for hourly PP)
#     - chl (mg/m^3)
#     - optional (set calculate_iops=FALSE to use these):
#         - at (total absorption, 1/m)
#         - bbt (total backscattering, 1/m)
#         - aphy (phytoplankton absorption, 1/m)
#         ***for each depth/row, these must be matrices (nrow=1 and ncol=length(lambda)) with values at every wavelength lambda
#     ***THE DEPTHS IN DEPTH_DATA MUST BE THE SAME DEPTHS DEFINED BY seq(Zmin,Zmax,dz), and all data must be non-NA
# time_data = list-col where each record/row is a tibble containing column:
#     - time = the time(s) of day to calculate pp, must be a 4-digit value (HHMM) in UTC
#     ***EITHER ONE VALUE (sampling time), OR A LIST OF 23 VALUES (hours 1-23 in HHMM format)
#
#*******************************************************************************
# OUTPUT
#
# - csv file containing the input columns and:
#       - euphotic_depth_avg = average euphotic depth (metres)
#       - total PP integrated over time and depth (mgC/m^2 for PP at sampling time or mgC/m^2/day for hourly PP)
# - optional (set detailed=TRUE below): Rdata file containing the input columns and:
#       - depth_data filled to the depths where PP is calculated
#       - PP (matrix where row=PP at depth, col=PP at time, mgC/sec for PP at sampling time or mgC/hr for hourly PP)
#       - scalar irradiance E0z (matrix with row=E0z at depth, col=E0z at time)
#
#*******************************************************************************
# REQUIRED LIBRARIES, SCRIPTS, AND DATA FILES

rm(list=ls())
library(stringr)
library(lubridate)
library(caTools)        # for integrating PP over depth and time using trapz()
library(data.table)     # for reading large csv input files with fast "fread" function
library(dplyr)          # for organizing input and output
library(pbapply)        # for progress bars
library(oce)            # to calculate sun angle in PAR_resolved.R
library(Rcpp)           # to use chunks of C++ code
library(RcppArmadillo)  # to use chunks of C++ code
library(parallel)       # to use parallel processing (if your computer has enough cores)
pad0 <- function (s, len) stringr::str_pad(s, width = len, side = "left", pad = "0")
num_cl <- NULL # set this to NULL to not use parallel processing
# required C++ libraries: cmath, vector, and RcppArmadillo

source("scripts/model_BIOv2.R")
source("scripts/PAR_resolved.R")
sourceCpp("scripts/surface_irradiance.cpp")
sourceCpp("scripts/underwater_irradiance.cpp")
sourceCpp("scripts/math.cpp")
# scripts/convert_units.R # sourced from PAR_resolved.R
# scripts/insol.R # sourced from PAR_resolved.R

# Besides the input file in the format described above, you also need the following data files to run PP, which must be stored in the "data" subfolder of the current working directory:
#   data/pope97.dat # absorption of pure sea water according to Pope and Fry (1997)
#   data/gcirrad.txt # extraterrestrial irradiance and atmospheric absorption coefficients
#   data/LUT_PAR_lati_day_vis.Rdata # lookup table giving visibility for a given latitude and day to use in PAR_resolved


#*******************************************************************************
# VARIABLES THAT CAN BE CHANGED

# water column intervals to calculate PP (metres), and the min and max depth
# if Zmax>bathymetry, it will be reset to the bathymetry
# depth_data columns will be linearly interpolated to the depths in seq(Zmin,Zmax,by=dz)
dz <- 0.5
Zmin <- 0
Zmax <- 100

# if you also want resolved output (all depths and hours per day, stored in an Rdata file) set this to TRUE
detailed <- TRUE

# if you want IOPs to be calculated at each depth based on chl, set this to TRUE (otherwise you must provide columns at, bbt, and aphy within depth_data)
calculate_iops <- TRUE

# uniform profiles in water column? (this will use the values closest to the surface in your depth_data)
uniform_profile <- FALSE

# Rdata file containing input data in the format specified at the top of this script
input_file <- "PP_input_example_samplingtime.Rdata"

# output filename, EXCLUDING extension
output_file <- "PP_output_example_samplingtime"


#*******************************************************************************
# MAIN CODE

if (!is.null(num_cl)) {num_cl <- min(num_cl, detectCores()-1)}

dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)

input <- get(load(input_file))

result <- pblapply(1:nrow(input), FUN=calculate_PP, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax, uniform_profile=uniform_profile, calculate_iops=calculate_iops, cl=num_cl)
# for (i in 1:nrow(input)) { # for diagnosing problem indices
#     test <- calculate_PP(i, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax)
# }

# add output to input dataframe, and integrate over time and depth to get summary values
output <- input %>%
    dplyr::mutate(PP=lapply(result,"[[","PP"),
                  E0z=lapply(result,"[[","E0z"))
cat("Input and output combined\n")

output$euphotic_depth <- pblapply(1:nrow(input), function(i) {
    E0z1 <- output$E0z[[i]]
    depth <- output$depth_data[[i]]$depth
    tmp <- sapply(1:ncol(E0z1), function(j) depth[which.min(abs((E0z1[,j]-0.01*E0z1[1,j])))])
    final <- rep(NA,length(tmp))
    final[sapply(tmp,length)!=0] <- unlist(tmp[sapply(tmp,length)!=0])
    return(final)
}, cl=num_cl)
output$euphotic_depth_avg <- sapply(output$euphotic_depth, FUN=mean, na.rm=TRUE)
cat("Euphotic depth extracted\n")

output$PP = pblapply(1:nrow(input), function(i) {
    depths <- output$depth_data[[i]]$depth
    euphotic_depths <- output$euphotic_depth[[i]]
    PP <- output$PP[[i]]
    tinds <- is.finite(euphotic_depths)
    for (j in which(tinds)) {
        PP[depths > euphotic_depths[j],j] <- NA
    }
    return(PP)
}, cl=num_cl)
output$pp_total = pbsapply(1:nrow(input), function(i) {
    depths <- output$depth_data[[i]]$depth
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


cat(paste0("Writing results to ",output_file,"\n"))

if (detailed) {
    output$pp_percent_total_per_depth = pblapply(1:nrow(input), function(i) {
        depths <- output$depth_data[[i]]$depth
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
    save(output, file=paste0(output_file,".Rdata"), compress=TRUE)
}

fwrite(output %>% dplyr::select(latitude,longitude,doy,year,par,euphotic_depth_avg,pp_total), file=paste0(output_file, ".csv"))

