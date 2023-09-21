
#*******************************************************************************
# DATA ####

# We define lambda (in nm) for the visible spectrum (i.e. PAR)
lambda = 400:700

# Absorption:
# Pure seawater is taken from Pope and Fry 1997
awp <- read.table("data/pope97.dat",skip=6) # absorption of pure sea water according to Pope and Fry (1997)
asw <- approx(awp[,1],awp[,2],lambda)$y # asw is absorption by pure seawater at wavelengths of interest

# EXAMPLE AT END OF SCRIPT



#*******************************************************************************
# PP FUNCTION ####

# df                         # Dataframe containing columns "lat", "doy", "year", "par", "chl", "alphab", "pmb"
#                              (latitude, day of year, year, surface PAR, surface chla, surface PI parameters)
# dz                         # Depth intervals (metres)
# Z = seq(Zmin,Zmax,by=dz)   # Vector of depths (metres)
# profile_type               # How to calculate chla in water column (uniform, gaussian, or interpolated)
# Rsize = length(lambda)     # Number of wavelengths (400 to 700 nm, 5nm increments)
# nstep = length(Z)          # Length of depths vector, where dz=0.5m and Zmax=250m
# xhr = 1:23                 # hours

calculate_PP <- function(df, lambda=400:700, dz=0.5, Z=seq(0,120,by=dz), profile_type="gaussian",
                         Rsize=length(lambda), nstep=length(Z), xhr=1:23) {
    
    latipxl <- df$lat
    daypxl <- df$doy
    yearpxl <- df$year
    SatPAR <- df$par
    chlpix <- df$chl
    alphaB <- df$alphab
    PBm <- df$pmb
    
    
    #*********************************
    # GET SURFACE PAR
    
    presolved <- PAR_resolved(latipxl=latipxl, daypxl=daypxl, yearpxl=yearpxl, SatPAR=SatPAR)
    Eqdifw <- presolved$Eqdifw
    if (all(is.na(Eqdifw))) {return(NA)}
    Eqdirw <- presolved$Eqdirw
    zendR <- presolved$zendR
    zendw <- presolved$zendw
    
    
    #*********************************
    # GET CHLA IN THE WATER COLUMN
    
    if (profile_type=="uniform") {
        chlpix <- unique(chlpix)
        if (length(chlpix)>1) {
            stop("Uniform chla profile selected but nonuniform chla profile supplied.")
        }
    } else if (profile_type=="gauss") {
        z_m = df[8]
        B_0 = df[9]
        h = df[10]
        sigma = df[11]
        chlpix = shifted_gaussian(tv=Z, B0=B_0, h=h, sigma=sigma, tmax=z_m)
        if (length(chlpix)==0 | all(!is.finite(chlpix))) {
            stop("Gaussian chla profile selected but parameters missing.")
        }
    } else if (profile_type=="interpolated") {
        tmp <- chlpix
        if (nrow(tmp)==0) {
            return(NA)
        } else if (nrow(tmp) == 1) {
            chlpix <- rep(tmp$chl,length(Z))
            warning("Chla provided at only one depth, uniform profile returned.")
        } else {
            chlpix <- approx(x=tmp$depth, y=tmp$chl, xout=Z, rule=2)$y
        }
    }
    
    
    #*********************************
    # GET UNDERWATER IRRADIANCE AND PP, using C++
    
    # Some information on "constant" IOPs such as absorption and backscattering from pure seawater
    # Seawater backscattering 
    BW500 = 0.00288
    BW = BW500*(lambda/500)^(-4.3)
    
    # for phytoplankton absorption we used the model from Devred et al. 2006 !!! THIS NEED TO BE UPDATED
    # AC = pc1*(1-exp(-rate*chlz) + pc2*chlz with:
    # ATLANTIC ZONE
    lamc = seq(400,700,5)
    pc1 = c(0.0744, 0.0773, 0.0778, 0.0776, 0.0773, 0.0745, 0.0746, 0.0748, 0.0723, 0.0668, 0.0612, 0.0571, 0.0546, 0.0516, 0.0488, 0.0446, 0.0401, 0.036, 0.0322, 0.0283, 0.0247, 0.0217, 0.019, 0.017, 0.0151, 0.0134, 0.012, 0.0106, 0.0093, 0.0081, 0.007, 0.0061, 0.0056, 0.0056, 0.0064, 0.007, 0.0089, 0.0104, 0.0103, 0.0118, 0.0105, 0.0102, 0.0113, 0.0112, 0.0116, 0.0119, 0.0128, 0.0135, 0.0116, 0.0103, 0.0104, 0.012, 0.0168, 0.0224, 0.0277, 0.0296, 0.0269, 0.0207, 0.0135, 0.0082, 0.0072)
    pc2 = c(0.0104, 0.0105, 0.0112, 0.0115, 0.0114, 0.0117, 0.0119, 0.0124, 0.0125, 0.0121, 0.0114, 0.0109, 0.0106, 0.0103, 0.0096, 0.009, 0.0084, 0.0081, 0.0079, 0.0078, 0.0077, 0.0073, 0.007, 0.0066, 0.0063, 0.0061, 0.0058, 0.0056, 0.0053, 0.005, 0.0046, 0.0042, 0.0037, 0.0034, 0.0031, 0.003, 0.0029, 0.0028, 0.0028, 0.0025, 0.0024, 0.0025, 0.0027, 0.003, 0.0033, 0.0034, 0.0035, 0.0035, 0.0036, 0.0036, 0.0037, 0.0043, 0.0056, 0.0074, 0.0089, 0.0092, 0.0081, 0.0061, 0.0039, 0.0023, 0.0011)
    rate = c(0.5582, 0.5706, 0.611, 0.6465, 0.6668, 0.7038, 0.7335, 0.7657, 0.7975, 0.8253, 0.8571, 0.8993, 0.9557, 1.0389, 1.0929, 1.1765, 1.259, 1.3455, 1.3954, 1.4338, 1.4485, 1.4396, 1.434, 1.415, 1.4161, 1.399, 1.365, 1.3307, 1.3254, 1.293, 1.2783, 1.2116, 1.0704, 0.8557, 0.6616, 0.5854, 0.47, 0.429, 0.4673, 0.4162, 0.4494, 0.4422, 0.3926, 0.4014, 0.4103, 0.4147, 0.4125, 0.4179, 0.5266, 0.6462, 0.6428, 0.5651, 0.4828, 0.5287, 0.5984, 0.6433, 0.6482, 0.5811, 0.4876, 0.4428, 0.3376)
    
    # pc1f = approx(lamc,pc1,lambda)$y
    # pc2f = approx(lamc,pc2,lambda)$y
    # ratef = approx(lamc,rate,lambda)$y
    # this uses a C++ function, which is generally veeeeeeeeeery slightly faster than approx()
    # (but the speed increase adds up when processing millions of pixels)
    pc1f = higher_res(lamc,lambda,pc1)
    pc2f = higher_res(lamc,lambda,pc2)
    ratef = higher_res(lamc,lambda,rate)
    
    # the slope of the exponential decrease is from Bricaud et al. 1981 Limnology and Oceanography
    Sy = 0.014
    # fraction of ays443 to aph443
    fracys = 0.35
    
    day_hrs <- which(zendR <= 90)
    
    # Eqdirw and Eqdifw are scalar here, not planar
    # They are also in Watt/m2/nm rather than units with Einsteins
    # As long as the units are compatible with alphaB and PmB units, this works out
    # In BIOCHEM, alphaB units are mgC/mgChl/hr, and PmB units are C/Chl/hr/Watt/m^2
    if (profile_type=="uniform") {
        PP <- underwater_irradiance_c(Eqdirw, Eqdifw, chlpix, nstep, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendR, zendw, alphaB, PBm, BW, min(day_hrs)-1, max(day_hrs)-1, dz)
    } else {
        PP <- underwater_irradiance_NU_c(Eqdirw, Eqdifw, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendR, zendw, alphaB, PBm, BW, min(day_hrs)-1, max(day_hrs)-1, dz)
    }
    
    # # R version for uniform profile:
    # # Here is the table where we store the light at each depth (i.e., 50cm), hour and lambda
    # E0dir = array(NaN,c(23,nstep,Rsize))
    # E0dif = array(NaN,c(23,nstep,Rsize))
    # # Here we just initialise the light at the top of the water column from the atmospheric model
    # E0dir[,1,] = Eqdirw
    # E0dif[,1,] = Eqdifw
    # # need to double-check if this requires scalar or planar irradiance input!!
    # PP <- underwater_irradiance(E0dir, E0dif, day_hrs, nstep, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendw, alphaB, PBm, BW)
    
    # multiply by depth interval to get the final value in metres
    finalPP <- sum(PP)*dz
    
    return(finalPP)
    
}


# #*******************************************************************************
# # EXAMPLE ####
# 
# # units:
# #       lat = decimal degrees
# #       par = einsteins/m^2/day
# #       chl = mg/m^3
# #       alphab = C/chl/hr/watt/m^2
# #       pmb = mgC/mg chl/hr
# #       zm = meters
# #       B0 = mg/m^3
# df <- data.frame(matrix(c(39.02083,91,2017,38.13178,5.12604,0.122075,2.12867,
#                           6.00070,0.295,243.056,19.1074),
#                         nrow=1))
# colnames(df) <- c("lat","doy","year","par","chl","alphab","pmb","zm","B0","h","sigma")
# # load custom functions and necessary libraries
# library(oceancolouR)    # for various custom functions
# library(Rcpp)
# library(RcppArmadillo)
# library(lubridate)
# source("R/PAR_resolved.R")
# sourceCpp("underwater_irradiance.cpp")
# sourceCpp("surface_irradiance.cpp")
# # try with uniform profile
# full_PP_calc(df, profile_type="uniform")
# # try with nonuniform profile
# full_PP_calc(df, profile_type="gaussian")
