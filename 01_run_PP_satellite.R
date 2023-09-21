cat("\014")     # Clear console
# rm(list=ls())   # Clear global environment

# Run primary production using the selected model and the input files that were
# created with the 01_format_PP_input.R script.

# CURRENTLY ONLY WORKS WITH BIO_V2
# other models use very similar formatting and will be incorporated later for easier use


#*******************************************************************************
# LIBRARIES, CUSTOM FUNCTIONS

library(ncdf4)          # to write output to netcdf files
library(data.table)     # for reading large csv input files with fast "fread" function
library(dplyr)          # for organizing input and output
library(oceancolouR)    # for various custom functions
library(pbapply)        # for progress bars
library(stringr)
library(Rcpp)
library(RcppArmadillo)
library(fst)
library(lubridate)
source("03b_BIOv2_model.R")
source("R/PAR_resolved.R")
sourceCpp("underwater_irradiance.cpp")
sourceCpp("surface_irradiance.cpp")
dformat <- "%Y%m%d"


#*******************************************************************************
# VARIABLES THAT CAN BE CHANGED

# Vector of date objects
# These will be converted to their equivalent 8day or monthly days, if those composites are selected.
d_list <- as_date("20150101",format=dformat):as_date("20150131",format=dformat)

# Length of composites used for PP calculation (daily, 8day, or monthly)
composite <- "monthly"

# model to use - BIOv2 only (need to do some reformatting before BIOv1 and Laval can be used with this)
model <- "BIOv2"

# Chla profile uniform or not?
uniform_profile <- TRUE

# If nonuniform, use Gaussian curve?
# (Must have properly named and formatted csv files containing parameters, created using format_input.R)
gauss <- TRUE



# Chla profile type in water column
# Options:
#       1. uniform (will use surface chla value)
#       2. gaussian (requires gaussian parameters in input table)
#       3. interpolated (requires table and chla and depth for basis of linear interpolation)
profile_type <- "interpolated"


# Use a single value for each PI parameter? (i.e. the same value for every pixel in the image)
use_constant_alphab <- TRUE
use_constant_pmb <- TRUE

# If using a constant alphab or pmb for the whole region, how should it be calculated?
#       OPTIONS: predefined (alphab=1, pmb=1.69), mean, or median
#       Note that if the image is subsetted below, the mean or median will be calculated using only the values inside those boundaries
alphab_method <- "median"
pmb_method <- "median"

# Output filename format = [nc_out]_NWA_[composite]_[date].nc
# [date] = YYYYMM or YYYYDDD for monthly or 8day respectively.
nc_out <- "BIOv2_"

# Do a subsetted region (for testing)
subset <- FALSE
lat_bounds <- c(43,60)
lon_bounds <- c(-65,-43)

# for nonuniform gaussian profiles, lowest allowed background chla (B0)
low_B0 <- 1e-5

# Number of clusters to use in parallel processing
# (2021-12-29 this is not used, instead we use Rcpp for faster processing)
num_cl <- 10 # detectCores()-1, or manually: 10 for hecla, 3 for Windows, 12 for Loki

# Input formatted csv files
input_path <- paste0("01_PP_input/", model, "/", composite, "/")

# Path for output files
output_path <- paste0("02_PP_output/", model, "/", composite, "/")



#*******************************************************************************
# MAIN CODE

# for writing to output
data("nwa_bins_4km")

# make sure the output folder exists
dir.create(output_path, showWarnings=FALSE, recursive=TRUE)

# get a full list of input files
files <- list.files(input_path)

# convert the selected dates to 8day or monthly, and subset to unique values
if (composite=="8day") {
    comp_nums <- week8(as_date(d_list))
    nondup_ind <- !duplicated(comp_nums)
    d_list <- d_list[nondup_ind]
    comp_nums <- comp_nums[nondup_ind]
    d_list <- week8_date(yearnum=year(as_date(d_list)), weeknum=comp_nums)
} else if (composite=="monthly") {
    comp_nums <- month(as_date(d_list))
    d_list <- d_list[!duplicated(comp_nums)]
}


# loop through dates
for (i in 1:length(d_list)) {
    
    date <- as_date(d_list[i])
    y <- year(date)
    m <- month(date)
    day <- day(date)
    doy <- yday(date)
    
    
    
    #*****************************************************************
    # LOAD INPUT

    
    
    
    #*****************************************************************
    # OPTION 1: USE NEW FORMATTED FILES
    datestr <- ifelse(composite=="monthly",
                      paste0(y, str_pad(m,width=2,side="left",pad="0")),
                      paste0(y, str_pad(doy,width=3,side="left",pad="0")))
    file <- files[grepl(datestr,files)]
    file_main <- paste0(input_path, file[which(grepl("satdata", file))])
    file_pi <- paste0(input_path, file[which(grepl("PIparams", file))])
    if (!file.exists(file_main) | !file.exists(file_pi)) {next}

    # variables needed for BIO_v2
    vars_to_select <- c("lat", "doy", "year", "par", "chl", "alphab", "pmb")

    # file_main = bin, lon, lat, bathy, chl, sst, par, tcc
    # file_pi = alphab and pmb (PI curve parameter)
    input <- bind_cols(fread(file_main), fread(file_pi)) %>% as.data.frame()

    # if non-uniform profile and gaussian curve, add BP parameter columns to input
    if (!uniform_profile & gauss) {
        file_bp <- paste0(input_path, file[which(grepl("BPparams", file))])
        if (!file.exists(file_bp)) {next}
        input <- bind_cols(input, fread(file_bp) %>% as.data.frame()) %>% dplyr::mutate(B0 = ifelse(B0 < low_B0, low_B0, B0))
        vars_to_select <- c(vars_to_select, "zm", "B0", "h", "sigma")
    }
    #*****************************************************************
    
    
    # #*******************************************************************
    # # OPTION 2: USE CARLA'S FORMATTED FILES FOR COMPARISON
    # datestr <- biov1_dates[i]
    # 
    # # temporary reading/formatting for Carla's files
    # input <- bind_cols(fread(paste0(input_path,datestr,"_LHTBI.csv")),
    #                    fread(paste0(input_path,datestr,"_pi-20141113.csv")),
    #                    fread(paste0(input_path,datestr,"_pp_fs.csv")),
    #                    fread(paste0(input_path,datestr,"_cp-20141117.csv"))) %>%
    #     dplyr::select(Name,Longitude,Latitude,height,chlor_a,sst,par,CLD_FR,
    #                   alphaB,PmB,zm,B0_sc,h_sc,sigma) %>%
    #     dplyr::mutate(B0_sc = ifelse(B0_sc < low_B0, low_B0, B0_sc))
    # colnames(input) <- c("bin", "lon", "lat", "bathy", "chl", "sst", "par",
    #                      "tcc", "alphab", "pmb", "zm", "B0", "h", "sigma")
    # 
    # vars_to_select <- c("lat", "doy", "year", "par", "chl",
    #                     "alphab", "pmb", "zm", "B0", "h", "sigma", "tcc")
    # 
    # # Carla's input files have erroneous data marked as -9000,-9001,-9002...
    # input[input < -500] <- NA
    # #*******************************************************************
    
    
    
    
    # subset (better for testing)
    if (subset) {
        input <- input %>%
            dplyr::filter(between(lat, lat_bounds[1], lat_bounds[2]),
                          between(lon, lon_bounds[1], lon_bounds[2]))
    }
    

    # if you're using one pair of PI parameters for every pixel, set them
    # all to either the mean value of the region to process, the median value,
    # or another predefined value (alphab=1 and PmB=1.69, from the Laval model)
    if (use_constant_alphab) {
        input$alpha <- ifelse(alphab_method=="median", median(input$alphab, na.rm=TRUE),
                              ifelse(alphab_method=="mean", mean(input$alphab, na.rm=TRUE),
                                     ifelse(alphab_method=="predefined", 1, NA)))
    }
    if (use_constant_pmb) {
        input$pmb <- ifelse(pmb_method=="median", median(input$pmb, na.rm=TRUE),
                            ifelse(pmb_method=="mean", mean(input$pmb, na.rm=TRUE),
                                   ifelse(pmb_method=="predefined", 1.69, NA)))
    }
    
    # grab only the columns needed for the BIO model
    input_subset <- input %>%
        dplyr::mutate(doy=doy, year=y) %>%
        dplyr::select(all_of(vars_to_select)) %>%
        as.data.frame()
    
    # make sure the necessary variables for the selected model are all valid
    valid_inds <- rowSums(is.finite(as.matrix(input_subset)))==ncol(input_subset)
    input <- input[valid_inds,]
    input_subset <- input_subset[valid_inds,]
    
    stop()
    
    #*******************************************************************
    # CALCULATE PP FOR EACH VALID INDEX
    
    if (uniform_profile) {
        result <- pbapply(input_subset, MARGIN=1, FUN=full_PP_calc)
    } else {
        result <- pbapply(input_subset, MARGIN=1, FUN=full_PP_calc_NU)
    }
    
    result <- as.numeric(result)
    
    # # find the index creating an error
    # for (i in 53246:55000) {
    #     print(i)
    #     print(as.numeric(full_PP_calc_NU(df=input_subset[i,])))
    # }
    
    
    #*******************************************************************
    # WRITE OUTPUT TO RDA AND NETCDF
    
    output_name <- paste0(output_path, nc_out, "_NWA_", composite, "_", datestr)
    if (subset) {
        output_name <- paste0(output_name, "_", paste0(lat_bounds,collapse="-"),"N_", paste0(sort(abs(lon_bounds)),collapse="-"),"W")
    }
    cat(paste0("Writing results to ",output_name,"\n"))
    
    # format output: remove invalid data, place on NWA grid
    output <- dplyr::left_join(data.frame(bin=nwa_bins_4km, stringsAsFactors = FALSE),
                               dplyr::bind_cols(input, data.frame(PP=result)) %>% dplyr::filter(PP >= 0),
                               by="bin")
    
    # write the rows with valid PP (in the subsetted region, if applicable),
    # and all the columns (all input and output) to Rdata file
    rda_output <- output %>% dplyr::filter(is.finite(PP))
    save(rda_output, file=paste0(output_name, ".rda"), compress=TRUE)
    
    # create netcdf dimensions
    dim_bin <- ncdim_def(name="binDataDim",
                         units = "",
                         vals = 1:nrow(output),
                         create_dimvar = FALSE)
    dim_time <- ncdim_def(name = "time",
                          units = "seconds since 1970-01-01 00:00:00.000 +0000",
                          vals = as.numeric(as_date(paste0(y,doy),format="%Y%j"))*24*60*60,
                          unlim = TRUE,
                          calendar = "standard")
    # create netcdf output variable
    output_var <- ncvar_def(name="PP",
                            units="mg C / m^2 / day",
                            dim=list(dim_bin, dim_time),
                            missval=NA,
                            longname="Primary Production using BIO model",
                            compression=5)
    # create new output netcdf
    ncout <- nc_create(paste0(output_name, ".nc"), output_var, force_v4=TRUE)
    # put variable in file
    ncvar_put(ncout, output_var, vals=output$PP)
    # close file
    nc_close(ncout)
    
}
