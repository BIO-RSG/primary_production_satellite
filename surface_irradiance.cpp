#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

// Calculate ozone content
// Output units: kilo Dobson units???
// [[Rcpp::export]]
double comp_sco3_c(double ylat,double xlon,int jday,double ro3) {
    // ylat = latitude (decimal degrees)
    // xlon = longitude (decimal degrees)
    // jday = day of year
    // ro3 = ozone (Dobson Units)
    double rlat = ylat*PI/180;
    double rlon = xlon*PI/180;
    double to3;
    if (ro3 > 0) {
        to3 = ro3;
    } else {
        to3 = 235.0 + (150.0+40.0*sin(0.9865*(jday-30.0)) + 20.0*sin(3.0*(rlon)))*sin(1.28*rlat)*sin(1.28*rlat);
    }
    double sco3 = to3*1e-3;
    return(sco3);
}


// Correct for Earth-Sun distance
// Output units: Watt/m2/nm
// [[Rcpp::export]]
NumericVector solar_irr_c(NumericVector Fobar, int jday) {
    // Fobar = extraterrestrial solar irradiance, Watt/m2/nm
    // jday = day of year
    if(jday >365) {jday=365.0;} // modify date if it is last day of a leap year
    NumericVector Fo = Fobar * pow((1.0+0.0167*cos(2*PI*(jday-3)/365.0)),2);
    return(Fo);
}


// Computes surface reflectance for direct (rod) and diffuse (ros)
// components separately, as a function of theta, wind speed or stress.
// Output units: None
// [[Rcpp::export]]
List sfcrfl_c(double theta, double ws) {
    
    // theta (solar zenith angle, degrees)
    // ws = instant wind speed (m/s)
    
    double rn = 1.341;      //index of refraction of pure seawater (unitless)
    double roair = 1.2E3;   //density of air (g/m3)
    
    double rof;
    double rosps;
    double rospd;
    double rtheta;
    double rthetar;
    double rmin;
    double rpls;
    double sinp;
    double tanp;
    double a;
    double b;

    //  Foam and diffuse reflectance
    if (ws > 4.0) {
        if (ws < 7.0) {
            rof = roair*(6.2E-4 + 1.56E-3/ws)*2.2E-5*ws*ws - 4.0E-4;
        } else {
            rof = (roair*(0.49E-3 + 0.065E-3*ws)*4.5E-5 - 4.0E-5)*ws*ws;
        }
        rosps = 0.057;
    } else {
        rof = 0.0;
        rosps = 0.066;
    }

    //  Direct
    //  Fresnel reflectance for theta < 40, ws < 2 m/s
    if (theta < 40.0 || ws < 2.0) {
        if (theta == 0.0) {
            rospd = 0.0211;
        } else {
            rtheta = theta*PI/180;
            rthetar = asin(sin(rtheta)/rn);
            rmin = rtheta - rthetar;
            rpls = rtheta + rthetar;
            sinp = (sin(rmin)*sin(rmin))/(sin(rpls)*sin(rpls));
            tanp = (tan(rmin)*tan(rmin))/(tan(rpls)*tan(rpls));
            rospd = 0.5*(sinp + tanp);
        }
    } else { // Empirical fit otherwise
        a = 0.0253;
        b = -7.14E-4*ws + 0.0618;
        rospd = a*exp(b*(theta-40.0));
    }

    //  Reflectance totals
    double rodaux = rospd + rof;
    double rosaux = rosps + rof;
    
    List Output;
    Output["rod"] = rodaux;
    Output["ros"] = rosaux;
    return(Output);
    
}


// Computes aerosol parameters according to a simplified version of the Navy marine aerosol model.
// Output units: None
// [[Rcpp::export]]
List navaer_c(double rh,double am,double wsm,double ws,double vis) {
    
    // rh = relative humidity (unitless)
    // am = aerosol model (air mass type, 1-10, unitless)
    // wsm = average wind speed over 24 hours (m/s)
    // ws = instant wind speed (m/s)
    // vis = visibility (km)
    
    NumericVector ro = {0.03,0.24,2.0}; // mode radius (microns)
    NumericVector r = {0.1,1.0,10.0}; // aerosol size distribution radii (microns)
    double rlam = 0.55; // waveband (microns)
    
    //  Relative humidity factor
    if (rh >= 100.0) rh = 99.9;
    double frh = pow(((2.0 - rh/100.0)/(6.0*(1.0-rh/100.0))),0.333);

    //  Size distribution amplitude components
    NumericVector a = {2000.0*am*am,5.866*(wsm-2.2),0.01527*(ws-2.2)*0.05};
    if (a[1] < 0.5) a[1] = 0.5;
    if (a[2] < 1.4E-5) a[2] = 1.4E-5;

    //  Compute size distribution at three selected radii according to Navy method
    NumericVector dndr(3);
    for (int i=0; i < 3; i++) {
        dndr = dndr + a[i]*exp((-1)*pow(log(r/frh/ro[i]),2))/frh;
    }
    
    //  Least squares approximation
    double gama = (sum(log10(r)*log10(dndr)))/(sum(pow(log10(r),2))); // Junge exponent
    double alphaux = -(gama+3.0); // Angstrom exponent

    //  Compute beta
    double betaux = (3.91/vis)*pow(rlam,alphaux);

    //  Compute asymmetry parameter -- a function of alpha
    double asympaux;
    if (alphaux > 1.2) {
        asympaux = 0.65;
    } else if (alphaux > 0.0) {
        asympaux = 0.82;
    } else {
        asympaux = -0.14167*alphaux + 0.82;
    }
    
    //  Single scattering albedo at 550; function of RH
    double waux = (-0.0032*am + 0.972)*exp(3.06E-4*rh);
    
    List Output;
    Output["alpha"] = alphaux;
    Output["beta"] = betaux;
    Output["asymp"] = asympaux;
    Output["wa"] = waux;
    return(Output);

}


// Model for atmospheric transmittance of solar irradiance through a maritime atmosphere.
// Computes direct and diffuse separately. Includes water vapor and oxygen absorption.
// Output units:??????????????????????????????
// [[Rcpp::export]]
List atmodd_v2_c(double rod, double ros, NumericVector lam, double theta,
                 NumericVector oza, NumericVector ag, NumericVector aw,
                 double sco3, double p, double p0, double wv, double rh, double am,
                 double wsm, double ws, double vis, NumericVector Fo) {
    
    // rod = surface reflectance for direct irradiance (unitless)
    // ros = surface reflectance for diffuse irradiance (unitless)
    // lam = wavebands (nm)
    // theta = solar zenith angle (degrees)
    // oza = ozone absorption coefficient (1/cm)
    // ag = oxygen absorption coefficient (1/cm)
    // aw = water vapor absorption coefficient (1/cm)
    // sco3 = ozone scale height (cm)
    // p = atmospheric pressure at sea level (mb)
    // p0 = reference atmospheric pressure at sea level (mb)
    // wv = precipitable water (cm)
    // rh = relative humidity (unitless)
    // am = aerosol model (air mass type, 1-10, unitless)
    // wsm = average wind speed over 24 hours (m/s)
    // ws = instant wind speed (m/s)
    // vis = visibility (km)
    // Fo = extraterrestrial irradiance (Watt/m2/nm)
    
    //  Compute atmospheric path lengths (air mass); pressure-corrected
    //  Modified March, 1994 according to Kasten and Young 1989.
    double cosunz = cos(theta*PI/180);
    double rm = 1.0/(cosunz+0.50572*pow((96.07995-theta),-1.6364));
    double rmp = p/p0*rm; // M'theta in Gregg-Carder 1990
    double rmo = (1.0+22.0/6370.0)/pow((cosunz*cosunz+44.0/6370.0),0.5); // Moz in Gregg-Carder 1990

    //  Obtain aerosol parameters; simplified Navy aerosol model
    List res_aer = navaer_c(rh,am,wsm,ws,vis);
    double alpha = res_aer["alpha"];
    double beta = res_aer["beta"];
    double asymp = res_aer["asymp"];
    double wa = res_aer["wa"];
    double eta = -alpha;
    
    //   Forward scattering probability
    double alg = log(1.0-asymp);
    double afs = alg*(1.459+alg*(.1595+alg*.4129));
    double bfs = alg*(.0783+alg*(-.3824-alg*.5874));
    double Fa = 1.0 - 0.5*exp((afs+bfs*cosunz)*cosunz);

    // Compute spectral irradiance
    // Rayleigh, by Bird's method (note: lambda must be in units of nanometres,
    // and then it's converted to microns below)
    NumericVector rlam = lam*1.0E-3;
    NumericVector rtra = exp(-(1.0/(115.6406*pow(rlam,4) - 1.335*pow(rlam,2)))*rmp);
    //    Ozone
    NumericVector to = oza*sco3;   //optical thickness
    NumericVector otra = exp(-to*rmo);   //transmittance
    //   Aerosols
    NumericVector ta = beta*pow(rlam,eta);
    NumericVector atra = exp(-ta*rm);
    NumericVector taa = exp(-(1.0-wa)*ta*rm);
    NumericVector tas = exp(-wa*ta*rm);
    //   Oxygen/gases
    NumericVector gtra = exp((-1.41*ag*rmp)/(pow((1.0 + 118.3*ag*rmp),0.45)));
    //   Water Vapor
    NumericVector wtra = exp((-0.2385*aw*wv*rm)/(pow((1.0+20.07*aw*wv*rm),0.45)));

    //  Direct irradiance (multiplied by cos(SOLZ) so it's planar, not scalar)
    NumericVector Ediraux = Fo*cosunz*rtra*otra*atra*gtra*wtra*(1.0-rod);

    //   Diffuse irradiance (planar)
    NumericVector dray = Fo*cosunz*gtra*wtra*otra*taa*0.5*(1.0-pow(rtra,.95));
    NumericVector daer = Fo*cosunz*gtra*wtra*otra*pow(rtra,1.5)*taa*Fa*(1.0-tas);
    NumericVector Edifaux = (dray + daer)*(1.0-ros);
    
    //  Total irradiance
    NumericVector Edaux = Ediraux + Edifaux;
    
    List Output;
    Output["Edir"] = Ediraux;
    Output["Edif"] = Edifaux;
    Output["Ed"] = Edaux;
    return(Output);

}


// Calculates total, direct, and diffuse surface irradiance.
// Output units: Watt/m2/nm
// [[Rcpp::export]]
List surface_irr_c(NumericVector zendR, NumericVector lam,
                   NumericVector oza, NumericVector ag, NumericVector aw,
                   double sco3, double p0, double wv, double rh, double am,
                   double wsm, double ws, double pxlvis, NumericVector Fo) {
    
    // zendR = solar zenith angle (degrees)
    // lam = wavebands (nm)
    // oza = ozone absorption coefficient (1/cm)
    // ag = oxygen absorption coefficient (1/cm)
    // aw = water vapor absorption coefficient (1/cm)
    // sco3 = ozone scale height (cm)
    // p0 = reference atmospheric pressure at sea level (mb)
    // wv = precipitable water (cm)
    // rh = relative humidity (unitless)
    // am = aerosol model (air mass type, 1-10, unitless)
    // wsm = average wind speed over 24 hours (m/s)
    // ws = instant wind speed (m/s)
    // pxlvis = visibility (km)
    // Fo = extraterrestrial irradiance (Watt/m2/nm)
    
    // Ed(0+,lambda)
    NumericMatrix Ed(23,301);
    NumericMatrix Eddir(23,301);
    NumericMatrix Eddif(23,301);
    
    // Eq(0-)
    NumericMatrix Eqdirw(23,301); // direct
    NumericMatrix Eqdifw(23,301); // diffus
    
    List irr_gc;
    List irrm_gc;
    List res_sfc;
    double rod;
    double ros;
    
    NumericVector Eqdirw_res(23);
    NumericVector Eqdifw_res(23);
    NumericVector Ed_res(301);
    NumericVector Eddir_res(301);
    NumericVector Eddif_res(301);
    
    for (int i=0; i < 23; i++) {
        
        if (zendR[i] < 90.) {
            
            // above water spectral irradiance
            rod = 0.0;
            ros = 0.0;
            irr_gc = atmodd_v2_c(rod, ros, lam, zendR[i], oza, ag, aw, sco3,
                                 p0, p0, wv, rh, am, wsm, ws, pxlvis, Fo);
            Ed_res = irr_gc["Ed"];
            Ed(i, _) = Ed_res;
            Eddir_res = irr_gc["Edir"];
            Eddir(i, _) = Eddir_res;
            Eddif_res = irr_gc["Edif"];
            Eddif(i, _) = Eddif_res;
            
            // Under water spectral irradiance.
            res_sfc = sfcrfl_c(zendR[i],ws);
            rod = res_sfc["rod"];
            ros = res_sfc["ros"];
            irrm_gc = atmodd_v2_c(rod, ros, lam, zendR[i], oza, ag, aw, sco3,
                                  p0, p0, wv, rh, am, wsm, ws, pxlvis, Fo);
            Eqdirw_res = irrm_gc["Edir"];
            Eqdirw(i, _) = Eqdirw_res;
            Eqdifw_res = irrm_gc["Edif"];
            Eqdifw(i, _) = Eqdifw_res;

        }
    }
    
    List Output;
    Output["Ed"] = Ed;
    Output["Eddir"] = Eddir;
    Output["Eddif"] = Eddif;
    Output["Eqdirw"] = Eqdirw;
    Output["Eqdifw"] = Eqdifw;
    return(Output);
    
}

