# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-08-21

# format in situ data and write to Rdata file
# required variables:
#    event_id, sample_id, station, latitude, longitude, doy, year, par, depth_data, time_data
#       "doy" is the day of year
#       "par" is the surface par (currently we just take this from the satellite value)
#       depth_data must be a list-col where each cell/element is a dataframe containing columns depth, alphab, pmb, chl
#       time_data must be a list-col where each cell/element is a vector containing one column ("time"), the time(s) of day to calculate pp
#           "time" must be a 4-digit value: HHMM


# 2023-09-11
# CURRENTLY THIS TAKES THE SHALLOWEST CHL AND THE AVG OF THE PI PARAMETERS IN THE WATER COLUMN


rm(list=ls())
library(tibble)
library(dplyr)
library(lubridate)
library(oceancolouR)
library(ncdf4)
library(geodist)

biochem_file <- "../LabSea2022Bloom/in_situ_data/final_dataset/LabSea2022project_biochem_and_PIparams_withMLD_onlyAR7W_updatedSept623.csv"

par_files <- list.files("01_PP_input/in_situ", pattern=".PAR.", full.names=TRUE)
# collect satellite par pixels within this distance, and average them
max_par_dist1 <- 5000
max_par_dist2 <- 10000 # if no valid par within range, expand search

max_depth <- 200

output_file <- "01_PP_input/in_situ/BIOv2/formatted_insitu_atsamplingtime_homogprofile.Rdata"
# output_file <- "01_PP_input/in_situ/BIOv2/formatted_insitu_hourly_homogprofile.Rdata"

# units are per hour and will be multiplied by this
# if calculating PP at sampling time, want per second (pi_conv=1/3600)
pi_conv <- 1/3600
# pi_conv <- 1

# vector of times (HHMM) when PP should be calculated
# set to NULL to use the sampling time from the biochem file
times <- NULL
# times <- as.numeric(paste0(1:23,"00"))


#*******************************************************************************
# READ AND FORMAT INPUT DATA

# # check which records use beta in the PI model
# test <- read.csv(biochem_file) %>% dplyr::distinct(STATION,EVENT_ID,SAMPLE_ID,YEAR,MONTH,DAY,LATDEC,LONGDEC,DEPTH,MLD,PI_model,Pmax,alpha,beta)
# write.csv(test,file=gsub(".Rdata","_PIbetarecords.csv",output_file),row.names=FALSE)

df <- read.csv(biochem_file) %>%
    dplyr::filter(DEPTH <= max_depth) %>%
    dplyr::mutate(date=as_date(paste0(pad0(YEAR,4),pad0(MONTH,2),pad0(DAY,2)),format="%Y%m%d")) %>%
    dplyr::mutate(doy=yday(date)) %>%
    dplyr::rename(event_id=EVENT_ID,
                  sample_id=SAMPLE_ID,
                  station=STATION,
                  latitude=LATDEC,
                  longitude=LONGDEC,
                  year=YEAR,
                  depth=DEPTH,
                  alphab=alpha,
                  pmb=Pmax,
                  chl=Turner_CHLA,
                  time=TIME) %>%
    dplyr::select(event_id,sample_id,station,latitude,longitude,doy,year,depth,alphab,pmb,chl,time) %>%
    tidyr::pivot_longer(cols=c(chl,alphab,pmb)) %>%
    tidyr::drop_na(value) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider() %>%
    # make sure PI parameters have proper units
    dplyr::mutate(alphab=alphab*pi_conv,
                  pmb=pmb*pi_conv) %>%
    dplyr::arrange(year,doy,station,event_id,depth) %>%
    dplyr::group_by(year,doy,station,event_id) %>%
    # only use sample_id, lat, lon, and time from the surface
    dplyr::summarize(sample_id=first(sample_id),
                  latitude=first(latitude),
                  longitude=first(longitude),
                  time=first(time),
                  depth=first(depth),
                  alphab=mean(alphab,na.rm=TRUE),
                  pmb=mean(pmb,na.rm=TRUE),
                  chl=first(chl)) %>%
    dplyr::ungroup() %>%
    tidyr::nest(depth_data=c(depth,alphab,pmb,chl))

if (is.null(times)) {
    df <- df %>% tidyr::nest(time_data=c(time))
} else {
    df <- df %>% dplyr::mutate(time_data=list(tibble(time=times)))
}

df <- lapply(1:nrow(df), function(i) {
    tmp <- df[i,]
    tmpdd <- tmp$depth_data[[1]]
    if (any(is.finite(tmpdd$alphab)) & any(is.finite(tmpdd$chl))) return(tmp)
    return(NULL)
}) %>% do.call(what=dplyr::bind_rows)


#*******************************************************************************
# GET SATELLITE SURFACE PAR FOR EACH RECORD AND SAVE TO OUTPUT

final <- lapply(1:nrow(df), function(i) {
    tmp <- df[i,]
    # get par from satellite file
    date <- as_date(paste0(pad0(tmp$year,4),pad0(tmp$doy,3)),format="%Y%j")
    par_file <- par_files[grepl(format(date,"%Y%m%d"),par_files)]
    nc <- nc_open(par_file)
    par_dat <- data.frame(bin=ncvar_get(nc,"bin"), par=ncvar_get(nc,"par"))
    nc_close(nc)
    # get satellite par
    coords <- tmp %>% dplyr::select(longitude,latitude)
    par <- dplyr::left_join(get_bins(region="nwa",variables=c("bin","latitude","longitude")), par_dat, by="bin") %>%
        dplyr::filter(between(longitude,coords$longitude-0.5,coords$longitude+0.5) &
                          between(latitude,coords$latitude-0.5,coords$latitude+0.5)) %>%
        dplyr::mutate(dist=as.numeric(geodist::geodist(x=data.frame(longitude,latitude),y=coords))) %>%
        tidyr::drop_na(par)
    tmp_par <- par %>% dplyr::filter(dist<=max_par_dist1)
    if (nrow(tmp_par)==0) {tmp_par <- par %>% dplyr::filter(dist<=max_par_dist2)}
    if (nrow(tmp_par)==0) return(NULL)
    tmp <- tmp %>% dplyr::mutate(par=mean(tmp_par$par,na.rm=TRUE))
    return(tmp %>% dplyr::select(event_id,sample_id,station,latitude,longitude,doy,year,par,depth_data,time_data))
}) %>% do.call(what=dplyr::bind_rows)

save(final, file=output_file, compress=TRUE)

