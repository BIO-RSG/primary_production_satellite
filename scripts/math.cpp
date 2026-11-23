#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <vector>
using namespace std;
using namespace Rcpp;

// Source: http://www.cplusplus.com/forum/general/216928/
// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
// boolean argument extrapolate determines behaviour beyond ends of array (if needed)
// [[Rcpp::export]]
double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )
{
    int size = xData.size();
    
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > xData[i+1] ) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    
    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient
    
    return yL + dydx * ( x - xL );                                              // linear interpolation
}

// [[Rcpp::export]]
vector<double> higher_res(vector<double> xold, vector<double> xnew, vector<double> yold) {
    vector<double> ynew;
    for ( double x : xnew ) {
        double ynew_tmp = interpolate( xold, yold, x, true );
        ynew.push_back( ynew_tmp );
    }
    return(ynew);
}