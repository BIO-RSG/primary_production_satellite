
#*******************************************************************************
# DATA ####

# We define lambda (in nm) for the visible spectrum (i.e. PAR)
lambda = 400:700

# Absorption:
# Pure seawater is taken from Pope and Fry 1997
awp <- read.table("data/pope97.dat",skip=6) # absorption of pure sea water according to Pope and Fry (1997)
asw <- approx(awp[,1],awp[,2],lambda)$y # asw is absorption by pure seawater at wavelengths of interest


#*******************************************************************************
# PP FUNCTION ####

# THIS RETURNS A MATRIX OF PP VALUES (ROW=DEPTH, COLUMN=TIME)

# i                 # row index of input dataframe
# input             # Dataframe containing the following columns:
#                       latitude, longitude, doy, year, par, bathymetry, depth_data, time_data
#                       (see 01_run_PP.R for details)
# dz                # Depth interval (in metres)
# Zmin, Zmax        # First and last depth over which PP will be calculated and integrated
# uniform_profile   # TRUE/FALSE, assume the water profile is homogeneous? (if TRUE, surface depth_data will be repeated at depth)

calculate_PP <- function(i, input, dz=0.5, Zmin=0, Zmax=150, lambda=400:700, uniform_profile=TRUE, calculate_iops=TRUE) {
    
    df <- input[i,]
    latipxl <- df$latitude
    daypxl <- df$doy
    yearpxl <- df$year
    SatPAR <- df$par
    dd <- df$depth_data[[1]] %>%
        dplyr::arrange(depth) %>% # make sure dataframe is arranged from shallow to deep values
        dplyr::relocate(depth) # make sure first column is depth
    time <- df$time_data[[1]]$time
    bathy <- df$bathymetry
    
    if (calculate_iops) {
        dd <- dd %>% dplyr::select(depth,chl,alphab,pmb)
    }
    
    dd <- dd %>% tidyr::drop_na()
    
    # make sure you don't calculate PP below the bottom depth
    if (bathy<Zmax) {
        Zmax <- bathy
        dd <- dd %>% dplyr::filter(depth<=bathy)
        if (nrow(dd)==0) return(NA)
    }
    
    Rsize <- length(lambda)
    Z <- seq(Zmin,Zmax,by=dz)
    nstep <- length(Z)
    
    if (!all(dd$depth[1:length(Z)]==Z)) stop("Required input data missing for some depths")
    
    chlpix <- dd$chl
    alphaB <- dd$alphab
    PBm <- dd$pmb
    if (!calculate_iops) {
        if (uniform_profile) {
            dd[,2:ncol(dd)] <- dd[rep(1,nrow(dd)),2:ncol(dd)]
        }
        at <- dd$at %>% do.call(what=rbind)
        bbt <- dd$bbt %>% do.call(what=rbind)
        aphy <- dd$aphy %>% do.call(what=rbind)
    }
    
    
    #*********************************
    # GET SURFACE PAR
    
    datepxl <- as_datetime(paste0(pad0(yearpxl,4),pad0(daypxl,3),pad0(time,4)),format="%Y%j%H%M")
    tz(datepxl) <- "UTC"
    presolved <- PAR_resolved(latipxl=latipxl, daypxl=daypxl, yearpxl=yearpxl, SatPAR=SatPAR, datepxl=datepxl)
    Eqdifw <- presolved$Eqdifw
    if (all(is.na(Eqdifw))) {return(NA)}
    Eqdirw <- presolved$Eqdirw
    zendR <- presolved$zendR
    zendw <- presolved$zendw
    
    
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
    
    # Grab a subset of the input during day hours so you don't have to calculate PP at night too
    day_times <- which(zendR <= 90)
    if (length(day_times)==0) {
        warning("None of the selected times are during daylight hours, returning NULL.")
        return(NULL)
    }
    tmp_Eqdirw <- Eqdirw[day_times,]
    tmp_Eqdifw <- Eqdifw[day_times,]
    tmp_zendR <- zendR[day_times]
    tmp_zendw <- zendw[day_times]
    if (length(day_times)==1) {
        tmp_Eqdirw <- matrix(tmp_Eqdirw, nrow=1)
        tmp_Eqdifw <- matrix(tmp_Eqdifw, nrow=1)
    }
    
    # Eqdirw and Eqdifw are scalar here, not planar
    # They are also in Watt/m2/nm rather than units with Einsteins
    # As long as the units are compatible with alphaB and PmB units, this works out
    # In BIOCHEM, alphaB units are C/Chl/hr, and PmB units are C/Chl/hr/Watt/m^2
    if (calculate_iops) {
        if (uniform_profile) {
            ui_res <- underwater_irradiance_c(tmp_Eqdirw, tmp_Eqdifw, chlpix[1], nstep, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, tmp_zendR, tmp_zendw, alphaB[1], PBm[1], BW, length(day_times), dz)
        } else {
            ui_res <- underwater_irradiance_NU_c(tmp_Eqdirw, tmp_Eqdifw, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, tmp_zendR, tmp_zendw, alphaB, PBm, BW, length(day_times), dz)
        }
    } else {
        ui_res <- underwater_irradiance_precalIOP_c(tmp_Eqdirw, tmp_Eqdifw, chlpix, lambda, tmp_zendR, tmp_zendw, alphaB, PBm, length(day_times), dz, at, bbt, aphy)
    }
    
    PP <- E0z <- matrix(nrow=length(Z),ncol=length(datepxl))
    PP[,day_times] <- t(ui_res$PP)
    E0z[,day_times] <- t(ui_res$E0z) # E0z is the sum of planar irradiance from 400-700nm at each depth and time interval
    
    return(list(depth_data_filled=dd,PP=PP,E0z=E0z))
    
}
