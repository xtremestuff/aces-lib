
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
const float RRT_GLOW_GAIN = 0.05f;
const float RRT_GLOW_MID = 0.08f;

// Red modifier constants
const float RRT_RED_SCALE = 0.82f;
const float RRT_RED_PIVOT = 0.03f;
const float RRT_RED_HUE = 0.0f;
const float RRT_RED_WIDTH = 135.0f;

// Desaturation contants
const float RRT_SAT_FACTOR = 0.96f;
const array <array <float, 3>, 3> RRT_SAT_MAT = calc_sat_adjust_matrix( RRT_SAT_FACTOR, AP1_RGB2Y);




// ------- Glow module functions
float glow_fwd( float ycIn, float glowGainIn, float glowMid)
{
   float glowGainOut;

   if (ycIn <= 2.0f/3.0f * glowMid) {
     glowGainOut = glowGainIn;
   } else if ( ycIn >= 2.0f * glowMid) {
     glowGainOut = 0.0f;
   } else {
     glowGainOut = glowGainIn * (glowMid / ycIn - 1.0f/2.0f);
   }

   return glowGainOut;
}

float glow_inv( float ycOut, float glowGainIn, float glowMid)
{
    float glowGainOut;

    if (ycOut <= ((1.0f + glowGainIn) * 2.0f/3.0f * glowMid)) {
    glowGainOut = -glowGainIn / (1.0f + glowGainIn);
    } else if ( ycOut >= (2.0f * glowMid)) {
    glowGainOut = 0.0f;
    } else {
    glowGainOut = glowGainIn * (glowMid / ycOut - 1.0f/2.0f) / (glowGainIn / 2.0f - 1.0f);
    }

    return glowGainOut;
}

float sigmoid_shaper( float x)
{
    // Sigmoid function in the range 0 to 1 spanning -2 to +2.

    float t = max( 1.0f - fabs( x / 2.0f), 0.0f);
    float y = 1.0f + sign(x) * (1.0f - t * t);

    return y / 2.0f;
}


// ------- Red modifier functions
float cubic_basis_shaper
( 
  float x, 
  float w   // full base width of the shaper function (in degrees)
)
{
    array <array <float, 4>, 4> M = { { { -1.0f / 6.0f,  3.0f / 6.0f, -3.0f / 6.0f,  1.0f / 6.0f },
                      {  3.0f / 6.0f, -6.0f / 6.0f,  3.0f / 6.0f,  0.0f / 6.0f },
                      { -3.0f / 6.0f,  0.0f / 6.0f,  3.0f / 6.0f,  0.0f / 6.0f },
                      {  1.0f / 6.0f,  4.0f / 6.0f,  1.0f / 6.0f,  0.0f / 6.0f } } };
  
    array <float, 5> knots = { -w/2.0f,
                     -w/4.0f,
                     0.0f,
                     w/4.0f,
                     w/2.0f };
  
  float y = 0.0f;
  if ((x > knots[0]) && (x < knots[4])) {  
    float knot_coord = (x - knots[0]) * 4.0f/w;  
    int j = knot_coord;
    float t = knot_coord - j;
      
    array <float, 4> monomials = { t*t*t, t*t, t, 1.0f };

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
      y = 0.0f;
    }
  }
  
  return y * 3.0f/2.0f;
}

float center_hue( float hue, float centerH)
{
  float hueCentered = hue - centerH;
  if (hueCentered < -180.0f) hueCentered = hueCentered + 360.0f;
  else if (hueCentered > 180.0f) hueCentered = hueCentered - 360.0f;
  return hueCentered;
}

float uncenter_hue( float hueCentered, float centerH)
{
  float hue = hueCentered + centerH;
  if (hue < 0.0f) hue = hue + 360.0f;
  else if (hue > 360.0f) hue = hue - 360.0f;
  return hue;
}



array <float, 3> rrt_sweeteners( array <float, 3> in)
{
    array <float, 3> aces = in;
    
    // --- Glow module --- //
    float saturation = rgb_2_saturation( aces);
    float ycIn = rgb_2_yc( aces);
    float s = sigmoid_shaper( (saturation - 0.4f) / 0.2f);
    float addedGlow = 1.0f + glow_fwd( ycIn, RRT_GLOW_GAIN * s, RRT_GLOW_MID);

    aces = mult_f_f3( addedGlow, aces);

    // --- Red modifier --- //
    float hue = rgb_2_hue( aces);
    float centeredHue = center_hue( hue, RRT_RED_HUE);
    float hueWeight = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH);

    aces[0] = aces[0] + hueWeight * saturation * (RRT_RED_PIVOT - aces[0]) * (1.0f - RRT_RED_SCALE);

    // --- ACES to RGB rendering space --- //
    aces = clamp_f3( aces, 0.0f, HALF_POS_INF);
    array <float, 3> rgbPre = mult_f3_f44( aces, AP0_2_AP1_MAT);
    rgbPre = clamp_f3( rgbPre, 0.0f, HALF_MAX);
    
    // --- Global desaturation --- //
    rgbPre = mult_f3_f33( rgbPre, RRT_SAT_MAT);
    return rgbPre;
}

array <float, 3> inv_rrt_sweeteners(array <float, 3> in)
{
    array <float, 3> rgbPost = in;
    
    // --- Global desaturation --- //
    rgbPost = mult_f3_f33( rgbPost, invert_f33(RRT_SAT_MAT));

    rgbPost = clamp_f3( rgbPost, 0.0f, HALF_MAX);

    // --- RGB rendering space to ACES --- //
    array <float, 3> aces = mult_f3_f44( rgbPost, AP1_2_AP0_MAT);

    aces = clamp_f3( aces, 0.0f, HALF_MAX);

    // --- Red modifier --- //
    float hue = rgb_2_hue( aces);
    float centeredHue = center_hue( hue, RRT_RED_HUE);
    float hueWeight = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH);

    float minChan;
    if (centeredHue < 0.0f) { // min_f3(aces) = aces[1] (i.e. magenta-red)
      minChan = aces[1];
    } else { // min_f3(aces) = aces[2] (i.e. yellow-red)
      minChan = aces[2];
    }

    float a = hueWeight * (1.0f - RRT_RED_SCALE) - 1.0f;
    float b = aces[0] - hueWeight * (RRT_RED_PIVOT + minChan) * (1.0f - RRT_RED_SCALE);
    float c = hueWeight * RRT_RED_PIVOT * minChan * (1.0f - RRT_RED_SCALE);

    aces[0] = ( -b - sqrtf( b * b - 4.0f * a * c)) / ( 2.0f * a);

    // --- Glow module --- //
    float saturation = rgb_2_saturation( aces);
    float ycOut = rgb_2_yc( aces);
    float s = sigmoid_shaper( (saturation - 0.4f) / 0.2f);
    float reducedGlow = 1.0f + glow_inv( ycOut, RRT_GLOW_GAIN * s, RRT_GLOW_MID);

    aces = mult_f_f3( ( reducedGlow), aces);
    return aces;
}