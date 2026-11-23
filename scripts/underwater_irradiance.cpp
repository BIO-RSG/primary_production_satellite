#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <vector>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List underwater_irradiance_c(NumericMatrix Eqdirw, NumericMatrix Eqdifw, double chlz, int nstep,
                             NumericVector lambda, NumericVector pc1f, NumericVector pc2f,
                             NumericVector ratef, double fracys, double Sy, NumericVector asw,
                             NumericVector zendR, NumericVector zendw, double alphaB, double PBm,
                             NumericVector BW, double tsize, double dz) {
    
    double Rsize = lambda.size();
    NumericMatrix PP(tsize, nstep);
    NumericMatrix E0z(tsize, nstep); //scalar irradiance at depth z
    NumericVector Kddir(Rsize);
    NumericVector Kddif(Rsize);
    NumericVector tmp_Eqdirw(Rsize);
    NumericVector tmp_Eqdifw(Rsize);
    double PIz;
    double apbar;
    
    // To compute the underwater light field, we need to know total absorption and total backscattering
    // a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength
    
    // backscattering is modelled according to Loisel and Morel 1998
    double bbtilda = (0.78 + 0.42 * (-1) * log10(chlz)) * 0.01; // bbtilda et le rapport entre bb et b.
    if (bbtilda < 0.0005) {
        bbtilda = 0.0005;
    } else if (bbtilda > 0.01) {
        bbtilda = 0.01;
    }
    double BC660 = 0.407 * pow(chlz, 0.795); // ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)
    
    NumericVector BC = BC660*pow((660/lambda),(-log10(chlz))); // BC spectral dependence depends on chlz,
    // this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
    BC[BC < 0] = 0;
    NumericVector bbt = BC*bbtilda + BW*0.5; // Here, we get total backscattering.
    
    //phytoplankton absorption //## !!!!! THIS HAS TO BE UPDATED.
    NumericVector aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz;
    
    // Yellow substances and detritus absorption is derived from chl
    // Note that C++ is zero-indexed, so aphy[43] is the same as aphy[44] in R
    NumericVector ays = fracys * aphy[43] * exp(-Sy * (lambda -443));
    
    // Now we can compute total absorption
    NumericVector at = asw + ays + aphy;
    
    // Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
    apbar = mean(aphy);
    NumericVector alphaBc = alphaB * aphy / apbar;
    
    
    for (int it=0; it < tsize; it++) { // loop on the time of day (for daylight zendR)
        
        // get E at surface, to be reduced at each depth based on Kd
        tmp_Eqdirw = Eqdirw(it, _);
        tmp_Eqdifw = Eqdifw(it, _);
        
        // calculate PP for the first layer
        // PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_Eqdirw) + 1.20 * sum(alphaBc * tmp_Eqdifw);
        PIz = sum(alphaBc * (tmp_Eqdirw + tmp_Eqdifw));
        PP(it, 0) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));  // Caverhill and Platt 1991, equation 4
        E0z(it, 0) = sum(tmp_Eqdirw + tmp_Eqdifw);
        
        // calculate attenuation coefficient for each wavelength (for uniform profile, it's the same at every depth)
        Kddir = exp(-dz * (at + bbt) / cos(zendw[it]));
        Kddif = exp(-dz * (at + bbt) / 0.83);
        
        for (int iz=1; iz < nstep; iz++) { // loop on depth, we start at the second step as the first one is Eqdw (note: C++ uses zero-indexing)
            
            tmp_Eqdirw = tmp_Eqdirw * Kddir; 
            tmp_Eqdifw = tmp_Eqdifw * Kddif;
            
            // PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_Eqdirw) + 1.20 * sum(alphaBc * tmp_Eqdifw);
            PIz = sum(alphaBc * (tmp_Eqdirw + tmp_Eqdifw));
            PP(it, iz) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));  // Caverhill and Platt 1991, equation 4
            E0z(it, iz) = sum(tmp_Eqdirw + tmp_Eqdifw);
            
        } // end loop on depth
        
    } // end loop on time of day (xhr)
    
    
    List Output;
    Output["PP"] = PP;
    Output["E0z"] = E0z;
    return(Output);
    
}

// [[Rcpp::export]]
List underwater_irradiance_NU_c(NumericMatrix E0dir, NumericMatrix E0dif, NumericVector chlpix,
                                NumericVector lambda, NumericVector pc1f, NumericVector pc2f,
                                NumericVector ratef, double fracys, double Sy, NumericVector asw,
                                NumericVector zendR, NumericVector zendw, NumericVector alphaB, NumericVector PBm,
                                NumericVector BW, double tsize, double dz) {
    
    double nstep = chlpix.size();
    double Rsize = lambda.size();
    NumericMatrix PP(tsize, nstep);
    NumericMatrix E0z(tsize, nstep); //scalar irradiance at depth z
    NumericVector BC(Rsize);
    NumericVector bbt(Rsize);
    NumericVector aphy(Rsize);
    NumericVector ays(Rsize);
    NumericVector at(Rsize);
    NumericVector Kddif(Rsize);
    NumericVector alphaBc(Rsize);
    NumericVector E0dir_num(Rsize);
    double chlz;
    double alphaBz;
    double PBmz;
    double bbtilda;
    double BC660;
    double apbar;
    double PIz;
    NumericMatrix tmp_E0dir(tsize, Rsize);
    NumericMatrix tmp_E0dif(tsize, Rsize);
    
    for (int iz=1; iz < nstep; iz++) { // loop on depth, we start at the second step as the first one is Eqdw (note: C++ uses zero-indexing)
        
        chlz = chlpix[iz];
        alphaBz = alphaB[iz];
        PBmz = PBm[iz];
        
        // To compute the underwater light field, we need to know total absorption and total backscattering
        // a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength
        
        // backscattering is modelled according to Loisel and Morel 1998
        bbtilda = (0.78 + 0.42 * (-1) * log10(chlz)) * 0.01; // bbtilda et le rapport entre bb et b.
        if (bbtilda < 0.0005) {
            bbtilda = 0.0005;
        } else if (bbtilda > 0.01) {
            bbtilda = 0.01;
        }
        BC660 = 0.407 * pow(chlz, 0.795); // ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)
        
        BC = BC660*pow((660/lambda),(-log10(chlz))); // BC spectral dependence depends on chlz,
        // this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
        // also if chlz = 0 this will give an error
        BC[BC < 0] = 0;
        bbt = BC*bbtilda + BW*0.5; // Here, we get total backscattering.
        
        //phytoplankton absorption //## !!!!! THIS HAS TO BE UPDATED.
        aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz;
        
        // Yellow substances and detritus absorption is derived from chl
        // Note that C++ is zero-indexed, so aphy[43] is the same as aphy[44] in R
        ays = fracys * aphy[43] * exp(-Sy * (lambda -443));
        
        // Now we can compute total absorption
        at = asw + ays + aphy;
        
        // Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
        apbar = mean(aphy);
        alphaBc = alphaBz * aphy / apbar;
        
        Kddif = exp(-dz * (at + bbt) / 0.83);
        E0dir_num = -dz * (at + bbt);
        
        if (iz==1) {
            for (int it=0; it < tsize; it++) {
                // calculate PP for surface irradiance
                PIz = sum(alphaBc * (E0dir(it, _) + E0dif(it, _)));
                PP(it, 0) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2));
                E0z(it, 0) = sum(E0dir(it, _) + E0dif(it, _));
                // calculate attenuated light for second depth
                tmp_E0dir(it, _) = E0dir(it, _) * exp(E0dir_num / cos(zendw[it]));
                tmp_E0dif(it, _) = E0dif(it, _) * Kddif;
                // calculate PP for second depth
                PIz = sum(alphaBc * (tmp_E0dir(it, _) + tmp_E0dif(it, _)));
                PP(it, 1) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2)); // Caverhill and Platt 1991, equation 4
                E0z(it, 1) = sum(tmp_E0dir(it, _) + tmp_E0dif(it, _));
            }
        } else {
            for (int it=0; it < tsize; it++) {
                // calculate attenuated light for this depth
                tmp_E0dir(it, _) = tmp_E0dir(it, _) * exp(E0dir_num / cos(zendw[it]));
                tmp_E0dif(it, _) = tmp_E0dif(it, _) * Kddif;
                // calculate PP for this depth
                PIz = sum(alphaBc * (tmp_E0dir(it, _) + tmp_E0dif(it, _)));
                PP(it, iz) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2)); // Caverhill and Platt 1991, equation 4
                E0z(it, iz) = sum(tmp_E0dir(it, _) + tmp_E0dif(it, _));
            }
        }
        
    } // end loop on depth
    
    List Output;
    Output["PP"] = PP;
    Output["E0z"] = E0z;
    return(Output);
    
}

// [[Rcpp::export]]
List underwater_irradiance_precalIOP_c(NumericMatrix E0dir, NumericMatrix E0dif, NumericVector chlpix,
                                       NumericVector lambda, NumericVector zendR, NumericVector zendw,
                                       NumericVector alphaB, NumericVector PBm, double tsize, double dz,
                                       NumericMatrix at, NumericMatrix bbt, NumericMatrix aphy) {
    
    double nstep = chlpix.size();
    double Rsize = lambda.size();
    NumericMatrix PP(tsize, nstep);
    NumericMatrix E0z(tsize, nstep); //scalar irradiance at depth z
    NumericVector atz(Rsize);
    NumericVector bbtz(Rsize);
    NumericVector aphyz(Rsize);
    NumericVector Kddif(Rsize);
    NumericVector alphaBc(Rsize);
    NumericVector Knum(Rsize);
    double chlz;
    double alphaBz;
    double PBmz;
    double apbar;
    double PIz;
    NumericMatrix tmp_E0dir(tsize, Rsize);
    NumericMatrix tmp_E0dif(tsize, Rsize);
    
    for (int iz=1; iz < nstep; iz++) { // loop on depth, we start at the second step as the first one is Eqdw (note: C++ uses zero-indexing)
        
        // get chla, photosynthetic parameters, and IOPs at this depth for light and PP calculations
        chlz = chlpix[iz];
        alphaBz = alphaB[iz];
        PBmz = PBm[iz];
        atz = at(iz, _);
        bbtz = bbt(iz, _);
        aphyz = aphy(iz, _);
        
        // Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
        apbar = mean(aphyz);
        alphaBc = alphaBz * aphyz / apbar;
        
        Knum = -dz * (atz + bbtz);
        Kddif = exp(Knum / 0.83);
        
        if (iz==1) {
            for (int it=0; it < tsize; it++) {
                // calculate PP for surface irradiance
                PIz = sum(alphaBc * (E0dir(it, _) + E0dif(it, _)));
                PP(it, 0) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2));
                E0z(it, 0) = sum(E0dir(it, _) + E0dif(it, _));
                // calculate attenuated light for second depth
                tmp_E0dir(it, _) = E0dir(it, _) * exp(Knum / cos(zendw[it]));
                tmp_E0dif(it, _) = E0dif(it, _) * Kddif;
                // calculate PP for second depth
                PIz = sum(alphaBc * (tmp_E0dir(it, _) + tmp_E0dif(it, _)));
                PP(it, 1) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2)); // Caverhill and Platt 1991, equation 4
                E0z(it, 1) = sum(tmp_E0dir(it, _) + tmp_E0dif(it, _));
            }
        } else {
            for (int it=0; it < tsize; it++) {
                // calculate attenuated light for this depth
                tmp_E0dir(it, _) = tmp_E0dir(it, _) * exp(Knum / cos(zendw[it]));
                tmp_E0dif(it, _) = tmp_E0dif(it, _) * Kddif;
                // calculate PP for this depth
                PIz = sum(alphaBc * (tmp_E0dir(it, _) + tmp_E0dif(it, _)));
                PP(it, iz) = chlz * PIz/sqrt(1 + pow((PIz/PBmz),2)); // Caverhill and Platt 1991, equation 4
                E0z(it, iz) = sum(tmp_E0dir(it, _) + tmp_E0dif(it, _));
            }
        }
        
    } // end loop on depth
    
    List Output;
    Output["PP"] = PP;
    Output["E0z"] = E0z;
    return(Output);
    
}

