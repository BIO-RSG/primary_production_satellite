# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-08-21

# Run primary production using the BIOv2 model.
# THIS IS A SIMPLIFIED VERSION OF 01_RUN_PP.R, WHICH CALLS model_BIOv2_satellite_simple.R FOR SPEEDIER PROCESSING, WITHOUT SOME UNNECESSARY EXTRAS.
# - calculated at hourly intervals and integrated over day
# - uniform profile using satellite surface chla and oceanic regime PI parameters


# ******************************************************************************
#                   Required format of input Rdata file
#
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
# *fst file containing the input data, and output surface PP, total PP, and % PP down to different depths.
# fst is a format from the "fst" library in R, which reads/writes faster and write smaller files than csv. If you want to use csv instead, you can adjust the code near the bottom of this script.


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
source("scripts/model_BIOv2_satellite_simple.R")
source("scripts/PAR_resolved.R")
sourceCpp("scripts/surface_irradiance.cpp")
sourceCpp("scripts/underwater_irradiance.cpp")
library(parallel)
num_cl <- 10 # set this to NULL to not use parallel processing


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

mydates <- as_date(as_date("20140101",format="%Y%m%d"):as_date("20221231",format="%Y%m%d"))

for (mydate in mydates) {
    
    mydate <- as_date(mydate)
    
    input_file_dateless <- "01_PP_input/satellite/BIOv2/daily/modisaqua_PP_input_medianPI_dailyPAR_"
    output_file_dateless <- "02_PP_output/satellite/BIOv2/daily/modisaqua_PP_medianPI_dailyPAR_"
    
    input_file <- paste0(input_file_dateless,format(mydate,"%Y%m%d"),".Rdata")
    output_file <- paste0(output_file_dateless,format(mydate,"%Y%m%d"),"_U")
    
    if (!file.exists(input_file)) {
        cat("Input file missing:",input_file,"\n")
        next
    }
    
    
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
    
    if (nrow(input)==0) {
        cat("No valid data. Skipping...\n")
        next
    }
    
    result <- pblapply(1:nrow(input), FUN=calculate_PP, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax, cl=num_cl)
    # for (i in 1:nrow(input)) { # for diagnosing problem indices
    #     test <- calculate_PP(i, input=input, dz=dz, Zmin=Zmin, Zmax=Zmax)
    # }
    
    output <- dplyr::bind_cols(input %>% tidyr::unnest(depth_data) %>% dplyr::select(-depth,-time_data),do.call(dplyr::bind_rows,result))
    
    if (subset) {
        output_file <- paste0(output_file, "_", paste0(lat_bounds,collapse="-"),
                              "N_", paste0(sort(abs(lon_bounds)),collapse="-"),"W")
    }
    
    cat(paste0("Writing results to ",output_file,"\n"))
    # write results to csv and fst
    # fwrite(output, file=paste0(output_file, ".csv"))
    write_fst(output, path=paste0(output_file,".fst"), compress=100)
    
}
