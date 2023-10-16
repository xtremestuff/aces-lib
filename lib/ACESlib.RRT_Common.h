
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.RRT_Common.a1.1.0</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - RRT Common</ACESuserName>

//
// Contains functions and constants shared by forward and inverse RRT transforms
//

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.CustomDefs.h"



// "Glow" module constants
const double RRT_GLOW_GAIN = 0.05;
const double RRT_GLOW_MID = 0.08;

// Red modifier constants
const double RRT_RED_SCALE = 0.82;
const double RRT_RED_PIVOT = 0.03;
const double RRT_RED_HUE = 0.0;
const double RRT_RED_WIDTH = 135.0;

// Desaturation contants
const double RRT_SAT_FACTOR = 0.96;
const array <array <double, 3>, 3> RRT_SAT_MAT = calc_sat_adjust_matrix( RRT_SAT_FACTOR, AP1_RGB2Y);




// ------- Glow module functions
double glow_fwd( double ycIn, double glowGainIn, double glowMid)
{
   double glowGainOut;

   if (ycIn <= 2.0/3.0 * glowMid) {
     glowGainOut = glowGainIn;
   } else if ( ycIn >= 2.0 * glowMid) {
     glowGainOut = 0.0;
   } else {
     glowGainOut = glowGainIn * (glowMid / ycIn - 1.0/2.0);
   }

   return glowGainOut;
}

double glow_inv( double ycOut, double glowGainIn, double glowMid)
{
    double glowGainOut;

    if (ycOut <= ((1.0 + glowGainIn) * 2.0/3.0 * glowMid)) {
    glowGainOut = -glowGainIn / (1.0 + glowGainIn);
    } else if ( ycOut >= (2.0 * glowMid)) {
    glowGainOut = 0.0;
    } else {
    glowGainOut = glowGainIn * (glowMid / ycOut - 1.0/2.0) / (glowGainIn / 2.0 - 1.0);
    }

    return glowGainOut;
}

double sigmoid_shaper( double x)
{
    // Sigmoid function in the range 0 to 1 spanning -2 to +2.

    double t = max( 1.0 - fabs( x / 2.0), 0.0);
    double y = 1.0 + sign(x) * (1.0 - t * t);

    return y / 2.0;
}


// ------- Red modifier functions
double cubic_basis_shaper
( 
  double x, 
  double w   // full base width of the shaper function (in degrees)
)
{
    array <array <double, 4>, 4> M = { { { -1.0 / 6.0,  3.0 / 6.0, -3.0 / 6.0,  1.0 / 6.0 },
                      {  3.0 / 6.0, -6.0 / 6.0,  3.0 / 6.0,  0.0 / 6.0 },
                      { -3.0 / 6.0,  0.0 / 6.0,  3.0 / 6.0,  0.0 / 6.0 },
                      {  1.0 / 6.0,  4.0 / 6.0,  1.0 / 6.0,  0.0 / 6.0 } } };
  
    array <double, 5> knots = { -w/2.0,
                     -w/4.0,
                     0.0,
                     w/4.0,
                     w/2.0 };
  
  double y = 0.0;
  if ((x > knots[0]) && (x < knots[4])) {  
    double knot_coord = (x - knots[0]) * 4.0/w;  
    int j = knot_coord;
    double t = knot_coord - j;
      
    array <double, 4> monomials = { t*t*t, t*t, t, 1.0 };

    // (if/else structure required for compatibility with CTL < v1.5.)
    if ( j == 3) {
      y = monomials[0] * M[0][0] + monomials[1] * M[1][0] + 
          monomials[2] * M[2][0] + monomials[3] * M[3][0];
    } else if ( j == 2) {
      y = monomials[0] * M[0][1] + monomials[1] * M[1][1] + 
          monomials[2] * M[2][1] + monomials[3] * M[3][1];
    } else if ( j == 1) {
      y = monomials[0] * M[0][2] + monomials[1] * M[1][2] + 
          monomials[2] * M[2][2] + monomials[3] * M[3][2];
    } else if ( j == 0) {
      y = monomials[0] * M[0][3] + monomials[1] * M[1][3] + 
          monomials[2] * M[2][3] + monomials[3] * M[3][3];
    } else {
      y = 0.0;
    }
  }
  
  return y * 3.0/2.0;
}

double center_hue( double hue, double centerH)
{
  double hueCentered = hue - centerH;
  if (hueCentered < -180.0) hueCentered = hueCentered + 360.0;
  else if (hueCentered > 180.0) hueCentered = hueCentered - 360.0;
  return hueCentered;
}

double uncenter_hue( double hueCentered, double centerH)
{
  double hue = hueCentered + centerH;
  if (hue < 0.0) hue = hue + 360.0;
  else if (hue > 360.0) hue = hue - 360.0;
  return hue;
}



array <double, 3> rrt_sweeteners( array <double, 3> in)
{
    array <double, 3> aces = in;
    
    // --- Glow module --- //
    double saturation = rgb_2_saturation( aces);
    double ycIn = rgb_2_yc( aces);
    double s = sigmoid_shaper( (saturation - 0.4) / 0.2);
    double addedGlow = 1.0 + glow_fwd( ycIn, RRT_GLOW_GAIN * s, RRT_GLOW_MID);

    aces = mult_f_f3( addedGlow, aces);

    // --- Red modifier --- //
    double hue = rgb_2_hue( aces);
    double centeredHue = center_hue( hue, RRT_RED_HUE);
    double hueWeight = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH);

    aces[0] = aces[0] + hueWeight * saturation * (RRT_RED_PIVOT - aces[0]) * (1.0 - RRT_RED_SCALE);

    // --- ACES to RGB rendering space --- //
    aces = clamp_f3( aces, 0.0, HALF_POS_INF);
    array <double, 3> rgbPre = mult_f3_f44( aces, AP0_2_AP1_MAT);
    rgbPre = clamp_f3( rgbPre, 0.0, HALF_MAX);
    
    // --- Global desaturation --- //
    rgbPre = mult_f3_f33( rgbPre, RRT_SAT_MAT);
    return rgbPre;
}

array <double, 3> inv_rrt_sweeteners(array <double, 3> in)
{
    array <double, 3> rgbPost = in;
    
    // --- Global desaturation --- //
    rgbPost = mult_f3_f33( rgbPost, invert_f33(RRT_SAT_MAT));

    rgbPost = clamp_f3( rgbPost, 0.0, HALF_MAX);

    // --- RGB rendering space to ACES --- //
    array <double, 3> aces = mult_f3_f44( rgbPost, AP1_2_AP0_MAT);

    aces = clamp_f3( aces, 0.0, HALF_MAX);

    // --- Red modifier --- //
    double hue = rgb_2_hue( aces);
    double centeredHue = center_hue( hue, RRT_RED_HUE);
    double hueWeight = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH);

    double minChan;
    if (centeredHue < 0.0) { // min_f3(aces) = aces[1] (i.e. magenta-red)
      minChan = aces[1];
    } else { // min_f3(aces) = aces[2] (i.e. yellow-red)
      minChan = aces[2];
    }

    double a = hueWeight * (1.0 - RRT_RED_SCALE) - 1.0;
    double b = aces[0] - hueWeight * (RRT_RED_PIVOT + minChan) * (1.0 - RRT_RED_SCALE);
    double c = hueWeight * RRT_RED_PIVOT * minChan * (1.0 - RRT_RED_SCALE);

    aces[0] = ( -b - sqrtf( b * b - 4.0 * a * c)) / ( 2.0 * a);

    // --- Glow module --- //
    double saturation = rgb_2_saturation( aces);
    double ycOut = rgb_2_yc( aces);
    double s = sigmoid_shaper( (saturation - 0.4) / 0.2);
    double reducedGlow = 1.0 + glow_inv( ycOut, RRT_GLOW_GAIN * s, RRT_GLOW_MID);

    aces = mult_f_f3( ( reducedGlow), aces);
    return aces;
}