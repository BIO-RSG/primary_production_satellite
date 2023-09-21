source("R/zenith_BIO.R")
source("R/zenith_LAVAL.R")
source("R/daylat.R")
library(oce)
library(ggplot2)

latipxl = 50.
daypxl = 173
yearpxl = 2018

xhr = 1:23

dl_vars = daylat(day=daypxl,lat=latipxl) # Jan1 = Day 0
delta = dl_vars$delta
phi = dl_vars$phi
biov1_zenith = unlist(sapply(1:23, function(i) zenith_bio(TIME=xhr[i], LONG=0, delta=delta, phi=phi))) * 180/pi

laval_zenith = unlist(sapply(1:23, function(i) zenith_laval(doy=daypxl, hour=xhr[i], xcoor=0, ycoor=latipxl))[1,])

datepxl = str_replace(as.character(doyday(yearpxl,daypxl) + hours(xhr)),"AST","UTC")
biov2_zenith = 90 - oce::sunAngle(datepxl,longitude = 0, latitude = latipxl)$altitude

#*******************************************************************************

all = data.frame(time = rep(xhr, 3),
                 model = c(rep("biov1", 23), rep("laval", 23), rep("biov2", 23)),
                 zenith_angle = c(biov1_zenith, laval_zenith, biov2_zenith),
                 stringsAsFactors = FALSE)

p <- ggplot(all, aes(x=time, y=zenith_angle, color=model)) +
    geom_point(alpha=0.6, size=4) +
    theme_bw()

print(p)
