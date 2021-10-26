
#*******************************************************************************
# Load libraries, custom source code, and define variables

library(insol)
library(lubridate)
library(oce)
library(stringr)
library(fields)
source("R/convert_units.R")

# Gcmod.r is the Gregg and Carder model, that is used to retrieve the spectral
# irradiance just under the sea surface: E_d(0-, lambda) in Einstiein.m-2
# We first use the visibility (in km) as a parameter to scale above-water Satellite
# PAR to the above-water modeled PAR, once the visibility is found, we derive E_d(0-,lambda)
# source("R/Gcmod.r")
# Note: instead of sourcing this file, the functions have been converted to C++ for
# faster processing (see surface_irradiance.cpp), and the atmospheric data for
# Gregg-Carder is loaded here:
# Opens and reads light data for the radiative transfer model of Gregg and Carder 1990.
# This is the data of Table 1 in the L&O paper, with values added to extend the model to 800 nm.
data_atm <- read.table("data/gcirrad.txt",skip=4,header=TRUE)
inds <- 51:351 # table goes from 350-801nm, but we only want 400-700nm (indices 51-351)
lam <- data_atm[inds,1]
# convert from W/cm2/um to W/m2/nm (*10000 for cm2-->m2, /1000 for um-->nm, *10000/1000 = *10)
# Note: do NOT convert the units of the other input here (this is how the input
# should be used according to the Gregg-Carder 1990 paper)
Fobar <- data_atm[inds,2] * 10
oza <- data_atm[inds,3]
ag <- data_atm[inds,4]
aw <- data_atm[inds,5]

# Here We download clear sky irradiance, it will be used to derive Ed(lambda,0-)
# as we sill scale the visibility to match model and satellite PAR, the variable is PARday
load("data/LUT_PAR_lati_day_vis.Rdata")
# latibin = seq(35,80,0.25) # x axis for PARday 
# daybin = 1:365 # y axis for PARday 
visibility = 10^seq(-1.6,1.8,30) # z axis for PARday 

# We define lambda for the visible spectrum (i.e. PAR)
lambda = 400:700
# Atmospheric parameters
rh = 0.8  # humidity default value
am = 1 # aerosol model (1-10), 1 is default value (marine model)
wv = 1.5 # precipitable water
wsm = 0.0 # Average wind speed over 24 hours, default value (m/s)
ws = 0.0 # instant wind speed (m/s)
# vis = 40  # visibility in km; enter 15 for default - this is calculated instead using PARday to get pxlvis in PAR_resolved()
ro3 = 100 # ozone in Dobson Units, or enter -99 for climatological ozone
# values (computed from the day and Earth position)
iblw = 1 # = 0 to compute above-surface values, or 1 to compute below-surface values
p0 = 29.92 # reference atmospheric pressure at sea level
pres = 29.9 # Real time atmospheric pressure
nw = 1.34 # refractive index of seawater


# EXAMPLE AT END OF SCRIPT



#*******************************************************************************
# FUNCTION

PAR_resolved <- function(latipxl, daypxl, yearpxl, SatPAR) {#, Cloud) {
    
    # computer sun angle for day
    xhr = 1:23
    datepxl = doyday(yearpxl,daypxl) + hours(xhr)
    julianpxl = JD(doyday(yearpxl,daypxl))
    tz(datepxl) <- "UTC" # Here using timezone UTC
    
    # Here we convert PAR from Einstein.m-2.d-1 to mW.cm-2.um-1
    # see seawifs tech report vol 22 page 49 
    # https://oceancolor.gsfc.nasa.gov/atbd/par/seawifs_par_wfigs.pdf
    SatPAR = SatPAR/1.193 
    
    # Here we compute the sun zenith angle in degree, which is 90 - elevation
    zendR = 90 - oce::sunAngle(datepxl,longitude = 0, latitude = latipxl)$altitude
    zendw = asin( 1/nw * sin(zendR * pi/180))
    
    # We compute the ozone content, it can be derived from climatology or use a real time value
    sco3 = comp_sco3_c(latipxl, 0, daypxl, ro3) # Hoz in Gregg-Carder 1990 (ozone scale height)
    Fo = solar_irr_c(Fobar, daypxl) # Extraterrestrial irradiance in W/m2/nm
    
    # Find the visibility *pixvil that corresponds to  SatPAR for latitude and day of year:
    latidx_match = seq(35.0, 80.0, 0.25) # lats of the PARday file, from 35 to 80 in 0.25 degree increments
    if (length(latidx_match) != nrow(PARday)) {
        print("Error with PARday file and line 92")
        break
    } else if (latipxl > 80 | latipxl < 35) {
        return(list(Eqdifw=rep(NA,length(lam)), Eqdirw=rep(NA,length(lam)),
                    zendR=zendR, zendw=zendw))
    }
    latidx = which.min(abs(latidx_match - latipxl))
    # BEFORE: latidx = round((latipxl-35)/0.25) # remove 35 deg because the look-up-table starts at 35 degs
    
    PARvis = PARday[latidx,daypxl,]
    pxlvis = approx(PARvis, visibility, xout = SatPAR, rule=c(2,2))$y
    
    # get planar total, direct, and diffuse irradiance in w/m2/nm, using gregg-carder
    # si <- surface_irr(zendR, lam, oza, ag, aw, sco3, p0, wv, rh, am, wsm, ws, pxlvis, Fo)
    si <- surface_irr_c(zendR, lam, oza, ag, aw, sco3, p0, wv, rh, am, wsm, ws, pxlvis, Fo)
    # above surface
    Eddir <- si$Eddir
    Eddif <- si$Eddif
    # below surface
    Eqdirw <- si$Eqdirw
    Eqdifw <- si$Eqdifw
    
    # # THIS IS THE EXACT PROCESS USED IN BIOV1 for the same thing
    # Ed <- Watt_per_m2_to_Mol_per_m2_perhour(Ed, lambda) # change units to mol per m^2 per hour
    # Ed <- rowSums(Ed, na.rm=TRUE) # for each hour, integrate across wavebands
    # Ed <- sum(Ed, na.rm=TRUE) # integrate over the day to get to einsteins/m2/day
    # ScalingPARfactor = SatPAR * 1.193 / Ed # create scaling factor
    Ed <- (Eddir + Eddif)/cos(zendR*pi/180) # scalar direct and diffuse, above water surface
    ScalingPARfactor = SatPAR * 1.193 / sum(Watt_per_m2_to_Mol_per_m2_perhour(Ed, lambda), na.rm=TRUE)
    
    # UNDERWATER DIRECT AND DIFFUSE --- SCALE BY SATELLITE PAR
    Eqdirw = Eqdirw * ScalingPARfactor
    Eqdifw = Eqdifw * ScalingPARfactor
    
    # UNDERWATER DIRECT AND DIFFUSE --- PLANAR TO SCALAR
    # We divide the direct irradiance by the in-water sun zenith angle to convert from planar to scalar
    Eqdirw = Eqdirw/cos(zendw)
    # We divide the diffuse irradiance by 0.833 to convert from planar to scalar, as in Platt and Sathyendranath 1997, page 2624 last paragraph
    Eqdifw = Eqdifw/0.83
    
    return(list(Eqdifw=Eqdifw, Eqdirw=Eqdirw, zendR=zendR, zendw=zendw))
    
}


# # returns total, direct, and diffuse irradiance in w/m2/nm
# surface_irr <- function(zendR, lam, oza, ag, aw, sco3, p0, wv, rh, am, wsm, ws, pxlvis, Fo) {
# 
#     rad = 180/pi
# 
#     # Ed(0+,lambda)
#     Ed = matrix(0,23,301)
#     Eddir = matrix(0,23,301)
#     Eddif = matrix(0,23,301)
#
#     # Eq(0-)
#     Eqdirw = matrix(0,23,301) # direct
#     Eqdifw = matrix(0,23,301) # diffus
# 
#     for (i in 1:23)
#     {
# 
#         if (zendR[i] < 90.)
#         {
# 
#             # above water spectral irradiance
#             irr.gc = atmodd_v2(rod=0.0, ros=0.0, rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag,
#                                aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
#                                ws = ws, vis = pxlvis, Fo = Fo)
#             Ed[i,] = irr.gc$Ed
#             Eddir[i,] = irr.gc$Edir
#             Eddif[i,] = irr.gc$Eddif
# 
#             # Under water spectral irradiance
#             res.sfc <- sfcrfl(rad,zendR[i],ws)
#             irrm.gc = atmodd_v2(rod=res.sfc$rod, ros=res.sfc$ros, rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag,
#                                 aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
#                                 ws = ws, vis = pxlvis, Fo = Fo)
#             Eqdirw[i,] = irrm.gc$Edir
#             Eqdifw[i,] = irrm.gc$Edif
# 
#         }
#     }
# 
#     return(list(Ed=Ed, Eddir=Eddir, Eddif=Eddif, Eqdirw=Eqdirw, Eqdifw=Eqdifw))
# 
# }





# #*******************************************************************************
# # EXAMPLE ####
# 
# latipxl = 50 # decimal degrees
# daypxl = 173
# yearpxl = 2018
# SatPAR = 24 # einsteins/m^2/day
# 
# # load custom functions and necessary libraries
# library(oceancolouR)    # for various custom functions
# library(Rcpp)
# library(RcppArmadillo)
# library(lubridate)
# sourceCpp("surface_irradiance.cpp")
# # get surface irradiance and solar zenith angles
# result <- PAR_resolved(latipxl, daypxl, yearpxl, SatPAR)
# # plot of total planar irradiance just below the surface
# xhr = 1:23
# Z = seq(0,120,by=0.5)
# image.plot(xhr, lambda, result$Eqdirw + result$Eqdifw, xlab = "Time (h)", ylab = "wavelength (nm)",axes = F)
# axis(1)
# axis(2,at = seq(400,700,100),labels = seq(400,700,100))
# box()
