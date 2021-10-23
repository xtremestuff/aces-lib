
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.ODT_Common.a1.1.0</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - ODT Common</ACESuserName>

//
// Contains functions and constants shared by forward and inverse ODT transforms 
//

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.CustomDefs.h"
#include "ACESlib.Transform_Common.h"


// Target white and black points for cinema system tonescale
const float CINEMA_WHITE = 48.0f;
const float CINEMA_BLACK = pow10f(log10f(0.02f)); // CINEMA_WHITE / 2400. 
    // CINEMA_BLACK is defined in this roundabout manner in order to be exactly equal to 
    // the result returned by the cinema 48-nit ODT tonescale.
    // Though the min point of the tonescale is designed to return 0.02, the tonescale is 
    // applied in log-log space, which loses precision on the antilog. The tonescale 
    // return value is passed into Y_2_linCV, where CINEMA_BLACK is subtracted. If 
    // CINEMA_BLACK is defined as simply 0.02, then the return value of this subfunction
    // is very, very small but not equal to 0, and attaining a CV of 0 is then impossible.
    // For all intents and purposes, CINEMA_BLACK=0.02.




// Gamma compensation factor
const float DIM_SURROUND_GAMMA = 0.9811f;

// Saturation compensation factor
const float ODT_SAT_FACTOR = 0.93f;
const array<array <float, 3>, 3> ODT_SAT_MAT = calc_sat_adjust_matrix( ODT_SAT_FACTOR, AP1_RGB2Y);




const array<array <float, 3>, 3> D60_2_D65_CAT = calculate_cat_matrix( AP0.white, REC709_PRI.white);




float Y_2_linCV( float Y, float Ymax, float Ymin) 
{
  return (Y - Ymin) / (Ymax - Ymin);
}

float linCV_2_Y( float linCV, float Ymax, float Ymin) 
{
  return linCV * (Ymax - Ymin) + Ymin;
}

array <float, 3> Y_2_linCV_f3(const array <float, 3> &Y, float Ymax, float Ymin)
{
    array <float, 3> linCV;
    linCV[0] = Y_2_linCV( Y[0], Ymax, Ymin);
    linCV[1] = Y_2_linCV( Y[1], Ymax, Ymin);
    linCV[2] = Y_2_linCV( Y[2], Ymax, Ymin);
    return linCV;
}

array <float, 3> linCV_2_Y_f3(array <float, 3> linCV, float Ymax, float Ymin)
{
    array <float, 3> Y;
    Y[0] = linCV_2_Y( linCV[0], Ymax, Ymin);
    Y[1] = linCV_2_Y( linCV[1], Ymax, Ymin);
    Y[2] = linCV_2_Y( linCV[2], Ymax, Ymin);
    return Y;
}

array <float, 3> darkSurround_to_dimSurround(array <float, 3> linearCV)
{
  array <float, 3> XYZ = mult_f3_f44( linearCV, AP1_2_XYZ_MAT);
  array <float, 3> xyY = XYZ_2_xyY(XYZ);
  xyY[2] = clamp( xyY[2], 0.0f, HALF_POS_INF);
  xyY[2] = powf( xyY[2], DIM_SURROUND_GAMMA);
  XYZ = xyY_2_XYZ(xyY);

  return mult_f3_f44( XYZ, XYZ_2_AP1_MAT);
}

array <float, 3> dimSurround_to_darkSurround(array <float, 3> linearCV)
{
  array <float, 3> XYZ = mult_f3_f44( linearCV, AP1_2_XYZ_MAT);

  array <float, 3> xyY = XYZ_2_xyY(XYZ);
  xyY[2] = clamp( xyY[2], 0.0f, HALF_POS_INF);
  xyY[2] = powf( xyY[2], 1.0f/DIM_SURROUND_GAMMA);
  XYZ = xyY_2_XYZ(xyY);

  return mult_f3_f44( XYZ, XYZ_2_AP1_MAT);
}




/* ---- Functions to compress highlights ---- */
// allow for simulated white points without clipping

float roll_white_fwd( 
    float in,      // color value to adjust (white scaled to around 1.0)
    float new_wht, // white adjustment (e.g. 0.9 for 10% darkening)
    float width    // adjusted width (e.g. 0.25 for top quarter of the tone scale)
  )
{
    const float x0 = -1.0f;
    const float x1 = x0 + width;
    const float y0 = -new_wht;
    const float y1 = x1;
    const float m1 = (x1 - x0);
    const float a = y0 - y1 + m1;
    const float b = 2.0f * ( y1 - y0) - m1;
    const float c = y0;
    const float t = (-in - x0) / (x1 - x0);
    float out = 0.0f;
    if ( t < 0.0f)
        out = -(t * b + c);
    else if ( t > 1.0f)
        out = in;
    else
        out = -(( t * a + b) * t + c);
    return out;
}

float roll_white_rev( 
    float in,      // color value to adjust (white scaled to around 1.0)
    float new_wht, // white adjustment (e.g. 0.9 for 10% darkening)
    float width    // adjusted width (e.g. 0.25 for top quarter of the tone scale)
  )
{
    const float x0 = -1.0f;
    const float x1 = x0 + width;
    const float y0 = -new_wht;
    const float y1 = x1;
    const float m1 = (x1 - x0);
    const float a = y0 - y1 + m1;
    const float b = 2.0f * ( y1 - y0) - m1;
    float c = y0;
    float out = 0.0f;
    if ( -in < y0)
        out = -x0;
    else if ( -in > y1)
        out = in;
    else {
        c = c + in;
        const float discrim = sqrtf( b * b - 4.0f * a * c);
        const float t = ( 2.0f * c) / ( -discrim - b);
        out = -(( t * ( x1 - x0)) + x0);
    }
    return out;
}