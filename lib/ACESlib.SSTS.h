
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.SSTS.a1.1.0</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - SSTS</ACESuserName>

//
// Contains functions used for forward and inverse tone scale 
//

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.CustomDefs.h"


// Textbook monomial to basis-function conversion matrix.
const array <array <double, 3>, 3> M1 = { {
  {  0.5, -1.0, 0.5 },
  { -1.0,  1.0, 0.5 },
  {  0.5,  0.0, 0.0 }
} };

struct TsPoint
{
    double x;        // ACES
    double y;        // luminance
    double slope;    // 
};

struct TsParams
{
    TsPoint Min;
    TsPoint Mid;
    TsPoint Max;
    double coefsLow[6];
    double coefsHigh[6];    
};



// TODO: Move all "magic numbers" (i.e. values in interpolation tables, etc.) to top 
// and define as constants

const double MIN_STOP_SDR = -6.5;
const double MAX_STOP_SDR = 6.5;

const double MIN_STOP_RRT = -15.0;
const double MAX_STOP_RRT = 18.0;

const double MIN_LUM_SDR = 0.02;
const double MAX_LUM_SDR = 48.0;

const double MIN_LUM_RRT = 0.0001;
const double MAX_LUM_RRT = 10000.0;


double lookup_ACESmin( double minLum )
{
    const array <array <double, 2>, 2> minTable = { { { log10f(MIN_LUM_RRT), MIN_STOP_RRT },
                                   { log10f(MIN_LUM_SDR), MIN_STOP_SDR } } };

    return 0.18*powf( 2.0, interpolate1D( minTable, log10f( minLum)));
}

double lookup_ACESmax( double maxLum )
{
    const array <array <double, 2>, 2>  maxTable = { { { log10f(MAX_LUM_SDR), MAX_STOP_SDR },
                                   { log10f(MAX_LUM_RRT), MAX_STOP_RRT } } };

    return 0.18*powf( 2.0, interpolate1D( maxTable, log10f( maxLum)));
}

array <double, 5> init_coefsLow(
    TsPoint TsPointLow,
    TsPoint TsPointMid
)
{
    array <double, 5> coefsLow;

    double knotIncLow = (log10f(TsPointMid.x) - log10f(TsPointLow.x)) / 3.0;
    // double halfKnotInc = (log10f(TsPointMid.x) - log10f(TsPointLow.x)) / 6.0;

    // Determine two lowest coefficients (straddling minPt)
    coefsLow[0] = (TsPointLow.slope * (log10f(TsPointLow.x)-0.5*knotIncLow)) + ( log10f(TsPointLow.y) - TsPointLow.slope * log10f(TsPointLow.x));
    coefsLow[1] = (TsPointLow.slope * (log10f(TsPointLow.x)+0.5*knotIncLow)) + ( log10f(TsPointLow.y) - TsPointLow.slope * log10f(TsPointLow.x));
    // NOTE: if slope=0, then the above becomes just 
        // coefsLow[0] = log10f(TsPointLow.y);
        // coefsLow[1] = log10f(TsPointLow.y);
    // leaving it as a variable for now in case we decide we need non-zero slope extensions

    // Determine two highest coefficients (straddling midPt)
    coefsLow[3] = (TsPointMid.slope * (log10f(TsPointMid.x)-0.5*knotIncLow)) + ( log10f(TsPointMid.y) - TsPointMid.slope * log10f(TsPointMid.x));
    coefsLow[4] = (TsPointMid.slope * (log10f(TsPointMid.x)+0.5*knotIncLow)) + ( log10f(TsPointMid.y) - TsPointMid.slope * log10f(TsPointMid.x));
    
    // Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    array <array <double, 2>, 2> bendsLow = { { {MIN_STOP_RRT, 0.18},
                             {MIN_STOP_SDR, 0.35} } };
    double pctLow = interpolate1D( bendsLow, log2f(TsPointLow.x/0.18));
    coefsLow[2] = log10f(TsPointLow.y) + pctLow*(log10f(TsPointMid.y)-log10f(TsPointLow.y));

    return coefsLow;
} 

array <double, 5> init_coefsHigh(
    TsPoint TsPointMid, 
    TsPoint TsPointMax
)
{
    array <double, 5> coefsHigh;

    double knotIncHigh = (log10f(TsPointMax.x) - log10f(TsPointMid.x)) / 3.0;
    // double halfKnotInc = (log10f(TsPointMax.x) - log10f(TsPointMid.x)) / 6.0;

    // Determine two lowest coefficients (straddling midPt)
    coefsHigh[0] = (TsPointMid.slope * (log10f(TsPointMid.x)-0.5*knotIncHigh)) + ( log10f(TsPointMid.y) - TsPointMid.slope * log10f(TsPointMid.x));
    coefsHigh[1] = (TsPointMid.slope * (log10f(TsPointMid.x)+0.5*knotIncHigh)) + ( log10f(TsPointMid.y) - TsPointMid.slope * log10f(TsPointMid.x));

    // Determine two highest coefficients (straddling maxPt)
    coefsHigh[3] = (TsPointMax.slope * (log10f(TsPointMax.x)-0.5*knotIncHigh)) + ( log10f(TsPointMax.y) - TsPointMax.slope * log10f(TsPointMax.x));
    coefsHigh[4] = (TsPointMax.slope * (log10f(TsPointMax.x)+0.5*knotIncHigh)) + ( log10f(TsPointMax.y) - TsPointMax.slope * log10f(TsPointMax.x));
    // NOTE: if slope=0, then the above becomes just
        // coefsHigh[0] = log10f(TsPointHigh.y);
        // coefsHigh[1] = log10f(TsPointHigh.y);
    // leaving it as a variable for now in case we decide we need non-zero slope extensions
    
    // Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    array <array <double, 2>, 2> bendsHigh = { { {MAX_STOP_SDR, 0.89},
                              {MAX_STOP_RRT, 0.90} } };
    double pctHigh = interpolate1D( bendsHigh, log2f(TsPointMax.x/0.18));
    coefsHigh[2] = log10f(TsPointMid.y) + pctHigh*(log10f(TsPointMax.y)-log10f(TsPointMid.y));
    
    return coefsHigh;
}


double shift( double in, double expfShift)
{
    return powf(2.0,(log2f(in)-expfShift));
}


TsParams init_TsParams(
    double minLum,
    double maxLum,
    double expfShift = 0
)
{
    TsPoint MIN_PT = { lookup_ACESmin(minLum), minLum, 0.0};
    TsPoint MID_PT = { 0.18, 4.8, 1.55};
    TsPoint MAX_PT = { lookup_ACESmax(maxLum), maxLum, 0.0};
    array <double, 5> cLow = init_coefsLow( MIN_PT, MID_PT);
    array <double, 5> cHigh = init_coefsHigh( MID_PT, MAX_PT);
    MIN_PT.x = shift(lookup_ACESmin(minLum),expfShift);
    MID_PT.x = shift(0.18,expfShift);
    MAX_PT.x = shift(lookup_ACESmax(maxLum),expfShift);

    TsParams P = {
        {MIN_PT.x, MIN_PT.y, MIN_PT.slope},
        {MID_PT.x, MID_PT.y, MID_PT.slope},
        {MAX_PT.x, MAX_PT.y, MAX_PT.slope},
        {cLow[0], cLow[1], cLow[2], cLow[3], cLow[4], cLow[4]},
        {cHigh[0], cHigh[1], cHigh[2], cHigh[3], cHigh[4], cHigh[4]}
    };
         
    return P;
}


double ssts
( 
    double x,
    TsParams C
)
{
    const int N_KNOTS_LOW = 4;
    const int N_KNOTS_HIGH = 4;

    // Check for negatives or zero before taking the log. If negative or zero,
    // set to HALF_MIN.
    double logx = log10f( max(x, HALF_MIN )); 

    double logy;

    if ( logx <= log10f(C.Min.x) ) { 

        logy = logx * C.Min.slope + ( log10f(C.Min.y) - C.Min.slope * log10f(C.Min.x) );

    } else if (( logx > log10f(C.Min.x) ) && ( logx < log10f(C.Mid.x) )) {

        double knot_coord = (N_KNOTS_LOW-1) * (logx-log10f(C.Min.x))/(log10f(C.Mid.x)-log10f(C.Min.x));
        int j = knot_coord;
        double t = knot_coord - j;

        array <double, 3> cf = { C.coefsLow[ j], C.coefsLow[ j + 1], C.coefsLow[ j + 2]};

        array <double, 3> monomials = { t * t, t, 1.0 };
        logy = dot_f3_f3( monomials, mult_f3_f33( cf, M1));

    } else if (( logx >= log10f(C.Mid.x) ) && ( logx < log10f(C.Max.x) )) {

        double knot_coord = (N_KNOTS_HIGH-1) * (logx-log10f(C.Mid.x))/(log10f(C.Max.x)-log10f(C.Mid.x));
        int j = knot_coord;
        double t = knot_coord - j;

        array <double, 3> cf = { C.coefsHigh[ j], C.coefsHigh[ j + 1], C.coefsHigh[ j + 2]};

        array <double, 3> monomials = { t * t, t, 1.0 };
        logy = dot_f3_f3( monomials, mult_f3_f33( cf, M1));

    } else { //if ( logIn >= log10f(C.Max.x) ) { 

        logy = logx * C.Max.slope + ( log10f(C.Max.y) - C.Max.slope * log10f(C.Max.x) );

    }

    return pow10f(logy);

}


double inv_ssts
( 
    double y,
    TsParams C
)
{  
    const int N_KNOTS_LOW = 4;
    const int N_KNOTS_HIGH = 4;

    const double KNOT_INC_LOW = (log10f(C.Mid.x) - log10f(C.Min.x)) / (N_KNOTS_LOW - 1.0);
    const double KNOT_INC_HIGH = (log10f(C.Max.x) - log10f(C.Mid.x)) / (N_KNOTS_HIGH - 1.0);

    // KNOT_Y is luminance of the spline at each knot
    double KNOT_Y_LOW[ N_KNOTS_LOW];
    for (int i = 0; i < N_KNOTS_LOW; i = i+1) {
    KNOT_Y_LOW[ i] = ( C.coefsLow[i] + C.coefsLow[i+1]) / 2.0;
    };

    double KNOT_Y_HIGH[ N_KNOTS_HIGH];
    for (int i = 0; i < N_KNOTS_HIGH; i = i+1) {
    KNOT_Y_HIGH[ i] = ( C.coefsHigh[i] + C.coefsHigh[i+1]) / 2.0;
    };

    double logy = log10f( max(y,1e-10));

    double logx;
    if (logy <= log10f(C.Min.y)) {

        logx = log10f(C.Min.x);

    } else if ( (logy > log10f(C.Min.y)) && (logy <= log10f(C.Mid.y)) ) {

        unsigned int j;
        array <double, 3> cf;
        if ( logy > KNOT_Y_LOW[ 0] && logy <= KNOT_Y_LOW[ 1]) {
            cf[ 0] = C.coefsLow[0];  cf[ 1] = C.coefsLow[1];  cf[ 2] = C.coefsLow[2];  j = 0;
        } else if ( logy > KNOT_Y_LOW[ 1] && logy <= KNOT_Y_LOW[ 2]) {
            cf[ 0] = C.coefsLow[1];  cf[ 1] = C.coefsLow[2];  cf[ 2] = C.coefsLow[3];  j = 1;
        } else if ( logy > KNOT_Y_LOW[ 2] && logy <= KNOT_Y_LOW[ 3]) {
            cf[ 0] = C.coefsLow[2];  cf[ 1] = C.coefsLow[3];  cf[ 2] = C.coefsLow[4];  j = 2;
        } 

        const array <double, 3> tmp = mult_f3_f33( cf, M1);

        double a = tmp[ 0];
        double b = tmp[ 1];
        double c = tmp[ 2];
        c = c - logy;

        const double d = sqrtf( b * b - 4. * a * c);

        const double t = ( 2. * c) / ( -d - b);

        logx = log10f(C.Min.x) + ( t + j) * KNOT_INC_LOW;

    } else if ( (logy > log10f(C.Mid.y)) && (logy < log10f(C.Max.y)) ) {

        unsigned int j;
        array <double, 3> cf;
        if ( logy >= KNOT_Y_HIGH[ 0] && logy <= KNOT_Y_HIGH[ 1]) {
            cf[ 0] = C.coefsHigh[0];  cf[ 1] = C.coefsHigh[1];  cf[ 2] = C.coefsHigh[2];  j = 0;
        } else if ( logy > KNOT_Y_HIGH[ 1] && logy <= KNOT_Y_HIGH[ 2]) {
            cf[ 0] = C.coefsHigh[1];  cf[ 1] = C.coefsHigh[2];  cf[ 2] = C.coefsHigh[3];  j = 1;
        } else if ( logy > KNOT_Y_HIGH[ 2] && logy <= KNOT_Y_HIGH[ 3]) {
            cf[ 0] = C.coefsHigh[2];  cf[ 1] = C.coefsHigh[3];  cf[ 2] = C.coefsHigh[4];  j = 2;
        } 

        const array <double, 3> tmp = mult_f3_f33( cf, M1);

        double a = tmp[ 0];
        double b = tmp[ 1];
        double c = tmp[ 2];
        c = c - logy;

        const double d = sqrtf( b * b - 4. * a * c);

        const double t = ( 2. * c) / ( -d - b);

        logx = log10f(C.Mid.x) + ( t + j) * KNOT_INC_HIGH;

    } else { //if ( logy >= log10f(C.Max.y) ) {

        logx = log10f(C.Max.x);

    }

    return pow10f( logx);

}


array <double, 3> ssts_f3
( 
    array <double, 3> &x,
    TsParams C
)
{
    array <double, 3> out;
    out[0] = ssts( x[0], C);
    out[1] = ssts( x[1], C);
    out[2] = ssts( x[2], C);

    return out;
}


array <double, 3> inv_ssts_f3
( 
    const array <double, 3> &x,
    TsParams C
)
{
    array <double, 3> out;
    out[0] = inv_ssts( x[0], C);
    out[1] = inv_ssts( x[1], C);
    out[2] = inv_ssts( x[2], C);

    return out;
}