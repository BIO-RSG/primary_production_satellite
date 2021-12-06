# Stephanie.Clay@dfo-mpo.gc.ca
# 2021-12-03

library(dplyr)
library(data.table)
library(oceancolouR)
library(fst)
library(stringr)
library(lubridate)
source("01b_get_param.R") # to get profile and PI curve parameters
dformat <- "%Y%m%d"



# NOTE THAT WE CURRENTLY ONLY HAVE 8DAY OR MONTHLY TCC DOWNLOADED (not daily)
# THE NEW BIO MODEL DOESN'T USE IT THOUGH



#*******************************************************************************
# DESCRIPTION
#
# Create standard input for later use in BIO or Laval primary production scripts.
# NWA REGION ONLY.
#
# Output files (stored in ../data/formatted/):
#       1. sat_data_YYYYMM_formatted_composite_created[created date].csv
#               Columns: bin, lon, lat, bathy, chl, sst, par, tcc, yelsub
#       2. parameters_BP_YYYYMM_formatted_inteval_created[created date].csv
#               Columns: alphab, pmb
#       3. parameters_PI_YYYYMM_formatted_inteval_created[created date].csv
#               Columns: zm, B0, h, sigma
# 
# bin, lon, lat, bathy = bin number, latitude, longitude, bathymetry
# chl = (satellite) chlorophyll-a
# sst = (satellite) sea surface temperature
# par = (satellite) photosynthetically active radiation
# tcc = (satellite) total cloud cover
# alphab = (PI curve parameter) initial slope normalized to biomass concentration
# pmb = (PI curve parameter) assimiliation number, i.e. maximum photosynthetic rate, normalized to biomass concentration
# zm, B0, h, sigma = chlorophyll-a profile parameters


#*******************************************************************************
# VARIABLES TO CHANGE

# Vector of date objects
# These will be converted to their equivalent 8day or monthly days, if those composites are selected.
d_list <- as_date("20180701",format=dformat):as_date("20180731",format=dformat)

# Options: daily, 8day or monthly (i.e. get data for weekly PP processing, or monthly?)
composite <- "8day"

# VIIRS-SNPP or MODIS
sensor <- "VIIRS-SNPP"

input_path <- "/mnt/data3/claysa"
chl_path <- file.path(input_path, "VIIRS-SNPPfilled", "CHL_POLY4", "NWA", "8day", "annual_fst")
par_path <- file.path(input_path, "VIIRS-SNPP", "PAR", "PANCAN", "annual_fst")
sst_path <- file.path(input_path, "VIIRS-SNPP", "SST", "PANCAN", "annual_fst")
tcc_path <- NULL # set to NULL if it's not needed
# tcc_path <- "00b_formatted_TCC"

chl_filled <- TRUE

# number of clusters for parallel processing
num_cl <- 10

# valid data ranges
tcc_range <- c(0,1)
chl_range <- c(0.0095,64.0005)


#*******************************************************************************

# (nn_pre.f90: issue warnings for extreme values, but still use them.  A large
# number of warnings might indicate a data format problem.)

# get the region of a variable based on the number of rows in the annual fst file
which_region <- function(x) {
    ind <- nrow(x)/num_pix
    return(names(num_pix)[floor(ind)==ind & ind <= 366])
}

# get variables and format them
get_var <- function(file, composite, avg_cols) {
    var <- read_fst(file)
    reg <- which_region(var)
    var <- matrix(var$var, nrow=num_pix[[reg]])
    reg_bins <- get_bins(region=tolower(reg), variables="bin")
    bin_inds <- reg_bins %in% nwa_bins
    var <- var[bin_inds,]
    if (avg_cols) {
        var <- avg_columns(var, year=y, composite=composite)
    }
    return(var)
}

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

# Reference parameters for nonuniform biomass profile
bp_par <- read.csv("data/02a_LUT_cp-20141117.csv") %>%
    mutate(dates=as.Date(paste0(year,"-",mon,"-",day))) %>%
    filter(chlsur>=0,!is.na(chlsur),!is.na(tempsur))#,
           #between(sigma,0,250),between(zm,-100,200),between(rho,0,1))

# Reference parameters for PI curve
pi_par <- read.csv("data/02a_LUT_pi-20141113.csv") %>%
    mutate(dates = as.Date(paste0(year,"-",mon,"-",day)),
           dn = format(as.Date(dates,origin="1970-01-01"),"%j")) %>%
    filter(alpha >= 0, !is.na(pm), !is.na(alpha))#,
           #between(alpha,0,0.6),between(PmB,0,15))

# Path for output files
output_path <- file.path("01_PP_input", composite)
dir.create(output_path, showWarnings=FALSE, recursive=TRUE)

# End of the output name (after the date which is automatically used, before the extension .csv)
output_name_end <- paste0(composite, "_NWA_", sensor, "_created",Sys.Date())

num_pix <- sapply(num_pix, "[[", 1)
nwa_bins <- get_bins(region="nwa", variables="bin")


#*******************************************************************************

for (i in 1:length(d_list)) {
    
    date <- as_date(d_list[i])
    y <- year(date)
    m <- month(date)
    day <- day(date)
    doy <- yday(date)
    
    # if this date is for a new year, load and process a new file
    prev_year <- ifelse(i==1, 0, year(as_date(d_list[i-1])))
    if (i==1 | (i > 1 & y != prev_year)) {
        chl_file <- list.files(chl_path, pattern=as.character(y), full.names=TRUE)
        par_file <- list.files(par_path, pattern=as.character(y), full.names=TRUE)
        sst_file <- list.files(sst_path, pattern=as.character(y), full.names=TRUE)
        if (!file.exists(chl_file) | !file.exists(par_file) | !file.exists(sst_file)) {next}
        chl <- get_var(chl_file, composite, avg_cols=ifelse(composite != "daily" & !chl_filled, TRUE, FALSE))
        par <- get_var(par_file, composite, avg_cols=ifelse(composite != "daily", TRUE, FALSE))
        sst <- get_var(sst_file, composite, avg_cols=ifelse(composite != "daily", TRUE, FALSE))
    }
    
    date_ind <- ifelse(composite=="daily", doy,
                       ifelse(composite=="8day", week8(date),
                              ifelse(composite=="monthly", m)))
    if (date_ind > min(ncol(chl), ncol(sst), ncol(par))) {next}
    
    sat_vars <- bind_cols(get_bins(region="nwa", variables=c("bin", "longitude", "latitude", "bathymetry")),
                          chl=chl[,date_ind], sst=sst[,date_ind], par=par[,date_ind]) %>%
        dplyr::rename(lon=longitude,
                      lat=latitude,
                      bathy=bathymetry) %>%
        dplyr::mutate(bathy=abs(bathy)) %>%
        dplyr::filter(is.finite(chl) & is.finite(sst) & is.finite(par))
    
    if (is.null(tcc_path)) {
        sat_vars$tcc <- 0 # ignore tcc, set to clear sky
    } else {
        # get tcc filename and check if it exists
        # NOTE: tcc filenames use YYYYMMDD, but are not affected by leap years
        tcc_file <- ifelse(composite=="8day",
                           paste0(tcc_path, "/", composite, "/MYDAL2_E_CLD_FR_",y,"-",str_pad(m,width=2,side="left",pad="0"),"-",day,".csv"),
                           paste0(tcc_path, "/", composite, "/MYDAL2_M_CLD_FR_",y,"-",str_pad(m,width=2,side="left",pad="0"),".csv"))
        if (!file.exists(tcc_file)) {next}
        sat_vars <- left_join(x=sat_vars, y=fread(tcc_file), by="bin") %>%
            dplyr::filter(is.finite(tcc))
    }
    
    sat_vars <- sat_vars %>%
        dplyr::filter(between(tcc, tcc_range[1], tcc_range[2]) & between(chl, chl_range[1], chl_range[2]))
    
    
    #-------------------------------------------------------------------
    # GET PARAMETERS - "valid_chl" SHOULD BE LOGGED CHLA
    
    # This is using the current method coded in R, which is slightly
    # different from the one in the original complex Fortran BIO code.
    # How? We may never know
    
    cat("Computing biomass profile and PI curve parameters...\n")
    
    valid_sst <- sat_vars$sst
    valid_chl <- log(sat_vars$chl)
    
    # Biomass profile
    bp_param <- get_param(logchla=valid_chl,
                          sst=valid_sst,
                          month=m,
                          day=as.numeric(doy),
                          ref_tab=bp_par,
                          param_type="bp",
                          num_cl=min(parallel::detectCores()-1, num_cl))
    
    # PI curve
    pi_param <- get_param(logchla=valid_chl,
                          sst=valid_sst,
                          month=m,
                          day=as.numeric(doy),
                          ref_tab=pi_par,
                          param_type="pi",
                          num_cl=min(parallel::detectCores()-1, num_cl))
    
    
    #-------------------------------------------------------------------
    # FORMAT AND WRITE OUTPUT TO .CSV FILE
    
    datestr <- ifelse(composite=="monthly", paste0(y, str_pad(m,width=2,side="left",pad="0")), paste0(y, doy))
    
    cat(paste0("Writing output for ", datestr, "...\n"))
    
    # Write satellite variables, biomass profile parameters, and PI curve parameters
    fwrite(sat_vars, file = paste0(output_path, "/", datestr, "_satdata_", output_name_end, ".csv"))
    fwrite(bp_param[,c(3,4,1,2)], file = paste0(output_path, "/", datestr, "_BPparams_", output_name_end, ".csv"))
    fwrite(pi_param[,c(2,1)], file = paste0(output_path, "/", datestr, "_PIparams_", output_name_end, ".csv"))

}

