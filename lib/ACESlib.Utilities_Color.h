
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.Utilities_Color.a1.1.0</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - Color Utilities</ACESuserName>

//
// Color related constants and functions
//

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.CustomDefs.h"
#include "ACESlib.Utilities.h"



/* ---- Chromaticities of some common primary sets ---- */

const Chromaticities AP0 = // ACES Primaries from SMPTE ST2065-1
{
  { 0.73470f,  0.26530f},
  { 0.00000f,  1.00000f},
  { 0.00010f, -0.07700f},
  { 0.32168f,  0.33767f}
};

const Chromaticities AP1 = // Working space and rendering primaries for ACES 1.0
{
  { 0.713f,    0.293f},
  { 0.165f,    0.830f},
  { 0.128f,    0.044f},
  { 0.32168f,  0.33767f}
};

const Chromaticities REC709_PRI =
{
  { 0.64000f,  0.33000f},
  { 0.30000f,  0.60000f},
  { 0.15000f,  0.06000f},
  { 0.31270f,  0.32900f}
};

const Chromaticities P3D60_PRI =
{
  { 0.68000f,  0.32000f},
  { 0.26500f,  0.69000f},
  { 0.15000f,  0.06000f},
  { 0.32168f,  0.33767f}
};

const Chromaticities P3D65_PRI =
{
  { 0.68000f,  0.32000f},
  { 0.26500f,  0.69000f},
  { 0.15000f,  0.06000f},
  { 0.31270f,  0.32900f}
};

const Chromaticities P3DCI_PRI =
{
  { 0.68000f,  0.32000f},
  { 0.26500f,  0.69000f},
  { 0.15000f,  0.06000f},
  { 0.31400f,  0.35100f}
};

const Chromaticities ARRI_ALEXA_WG_PRI =
{
  { 0.68400f,  0.31300f},
  { 0.22100f,  0.84800f},
  { 0.08610f, -0.10200f},
  { 0.31270f,  0.32900f}
};

const Chromaticities REC2020_PRI = 
{
  { 0.70800f,  0.29200f},
  { 0.17000f,  0.79700f},
  { 0.13100f,  0.04600f},
  { 0.31270f,  0.32900f}
};

const Chromaticities RIMMROMM_PRI = 
{
  { 0.7347f,  0.2653f},
  { 0.1596f,  0.8404f},
  { 0.0366f,  0.0001f},
  { 0.3457f,  0.3585f}
};

const Chromaticities SONY_SGAMUT3_PRI =
{
  { 0.730f,  0.280f},
  { 0.140f,  0.855f},
  { 0.100f, -0.050f},
  { 0.3127f,  0.3290f}
};

const Chromaticities SONY_SGAMUT3_CINE_PRI =
{
  { 0.766f,  0.275f},
  { 0.225f,  0.800f},
  { 0.089f, -0.087f},
  { 0.3127f,  0.3290f}
};

// Note: No official published primaries exist as of this day for the
// Sony VENICE SGamut3 and Sony VENICE SGamut3.Cine colorspaces. The primaries
// have thus been derived from the IDT matrices.
const Chromaticities SONY_VENICE_SGAMUT3_PRI =
{
  { 0.740464264304292f,  0.279364374750660f},
  { 0.089241145423286f,  0.893809528608105f},
  { 0.110488236673827f, -0.052579333080476f},
  { 0.312700000000000f,  0.329000000000000f}
};

const Chromaticities SONY_VENICE_SGAMUT3_CINE_PRI =
{
  { 0.775901871567345f,  0.274502392854799f},
  { 0.188682902773355f,  0.828684937020288f},
  { 0.101337382499301f, -0.089187517306263f},
  { 0.312700000000000f,  0.329000000000000f}
};

const Chromaticities CANON_CGAMUT_PRI =
{
  { 0.7400f,  0.2700f},
  { 0.1700f,  1.1400f},
  { 0.0800f, -0.1000f},
  { 0.3127f,  0.3290f}
};

const Chromaticities RED_WIDEGAMUTRGB_PRI =
{
  { 0.780308f,  0.304253f},
  { 0.121595f,  1.493994f},
  { 0.095612f, -0.084589f},
  { 0.3127f,  0.3290f}
};

const Chromaticities PANASONIC_VGAMUT_PRI =
{
  { 0.730f,  0.280f},
  { 0.165f,  0.840f},
  { 0.100f, -0.030f},
  { 0.3127f,  0.3290f}
};




/* ---- Conversion Functions ---- */
// Various transformations between color encodings and data representations
//

// Transformations between CIE XYZ tristimulus values and CIE x,y 
// chromaticity coordinates
array <float,3> XYZ_2_xyY( array<float,3> XYZ)
{  
  array<float,3> xyY;
  float divisor = (XYZ[0] + XYZ[1] + XYZ[2]);
  if (divisor == 0.0f) divisor = (float)1e-10;
  xyY[0] = XYZ[0] / divisor;
  xyY[1] = XYZ[1] / divisor;  
  xyY[2] = XYZ[1];
  
  return xyY;
}

array <float,3> xyY_2_XYZ( array<float,3> xyY)
{
  array <float,3> XYZ;
  XYZ[0] = xyY[0] * xyY[2] / max( xyY[1], (float)1e-10);
  XYZ[1] = xyY[2];  
  XYZ[2] = (1.0f - xyY[0] - xyY[1]) * xyY[2] / max( xyY[1], (float)1e-10);

  return XYZ;
}


// Transformations from RGB to other color representations
float rgb_2_hue( array<float,3> rgb) 
{
  // Returns a geometric hue angle in degrees (0-360) based on RGB values.
  // For neutral colors, hue is undefined and the function will return a quiet NaN value.
  float hue;
  if (rgb[0] == rgb[1] && rgb[1] == rgb[2]) {
    hue = (float)FLT_NAN; // RGB triplets where RGB are equal have an undefined hue
  } else {
    hue = (float)(180.0/M_PI) * atan2( sqrtf(3.0f)*(rgb[1]-rgb[2]), 2.0f*rgb[0]-rgb[1]-rgb[2]);
  }
    
  if (hue < 0.0f) hue = hue + 360.0f;

  return hue;
}

float rgb_2_yc( array <float,3> rgb, float ycRadiusWeight = 1.75f)
{
  // Converts RGB to a luminance proxy, here called YC
  // YC is ~ Y + K * Chroma
  // Constant YC is a cone-shaped surface in RGB space, with the tip on the 
  // neutral axis, towards white.
  // YC is normalized: RGB 1 1 1 maps to YC = 1
  //
  // ycRadiusWeight defaults to 1.75, although can be overridden in function 
  // call to rgb_2_yc
  // ycRadiusWeight = 1 -> YC for pure cyan, magenta, yellow == YC for neutral 
  // of same value
  // ycRadiusWeight = 2 -> YC for pure red, green, blue  == YC for  neutral of 
  // same value.

  float r = rgb[0]; 
  float g = rgb[1]; 
  float b = rgb[2];
  
  float chroma = sqrtf(b*(b-g)+g*(g-r)+r*(r-b));

  return ( b + g + r + ycRadiusWeight * chroma) / 3.0f;
}



/* ---- Chromatic Adaptation ---- */

const array<array <float, 3>, 3> CONE_RESP_MAT_BRADFORD = { {
  { 0.89510f, -0.75020f,  0.03890f},
  { 0.26640f,  1.71350f, -0.06850f},
  {-0.16140f,  0.03670f,  1.02960f}
} };

const array<array <float, 3>, 3> CONE_RESP_MAT_CAT02 = { {
  { 0.73280f, -0.70360f,  0.00300f},
  { 0.42960f,  1.69750f,  0.01360f},
  {-0.16240f,  0.00610f,  0.98340f}
} };

array<array <float, 3>, 3> calculate_cat_matrix
  (const array <float, 2> &src_xy,         // x,y chromaticity of source white
      const array <float, 2> &des_xy,         // x,y chromaticity of destination white
      const array<array <float, 3>, 3> coneRespMat = CONE_RESP_MAT_CAT02 //Changed Default to CAT02
  )
{
  //
  // Calculates and returns a 3x3 Von Kries chromatic adaptation transform 
  // from src_xy to des_xy using the cone response primaries defined 
  // by coneRespMat. By default, coneRespMat is set to CONE_RESP_MAT_CAT02. 
  // The default coneRespMat can be overridden at runtime. 
  //

  const array <float, 3> src_xyY = { src_xy[0], src_xy[1], 1. };
  const array <float, 3> des_xyY = { des_xy[0], des_xy[1], 1. };

  array <float, 3> src_XYZ = xyY_2_XYZ( src_xyY );
  array <float, 3> des_XYZ = xyY_2_XYZ( des_xyY );

  array <float, 3> src_coneResp = mult_f3_f33( src_XYZ, coneRespMat);
  array <float, 3> des_coneResp = mult_f3_f33( des_XYZ, coneRespMat);

  array<array <float, 3>, 3> vkMat = { {
      { des_coneResp[0] / src_coneResp[0], 0.0, 0.0 },
      { 0.0, des_coneResp[1] / src_coneResp[1], 0.0 },
      { 0.0, 0.0, des_coneResp[2] / src_coneResp[2] }
  } };

  array<array <float, 3>, 3> cat_matrix = mult_f33_f33( coneRespMat, mult_f33_f33( vkMat, invert_f33( coneRespMat ) ) );

  return cat_matrix;
}



array<array <float, 3>, 3> calculate_rgb_to_rgb_matrix
(Chromaticities SOURCE_PRIMARIES,
    Chromaticities DEST_PRIMARIES,
    const array<array <float, 3>, 3> coneRespMat = CONE_RESP_MAT_BRADFORD
)
{
    //
    // Calculates and returns a 3x3 RGB-to-RGB matrix from the source primaries to the 
    // destination primaries. The returned matrix is effectively a concatenation of a 
    // conversion of the source RGB values into CIE XYZ tristimulus values, conversion to
    // cone response values or other space in which reconciliation of the encoding white is 
    // done, a conversion back to CIE XYZ tristimulus values, and finally conversion from 
    // CIE XYZ tristimulus values to the destination RGB values.
    //
    // By default, coneRespMat is set to CONE_RESP_MAT_BRADFORD. 
    // The default coneRespMat can be overridden at runtime. 
    //

    const array<array <float, 4>, 4> RGBtoXYZ_44 = RGBtoXYZ(SOURCE_PRIMARIES, 1.0);
    const array<array <float, 3>, 3> RGBtoXYZ_MAT =
    { { {RGBtoXYZ_44[0][0], RGBtoXYZ_44[0][1], RGBtoXYZ_44[0][2]},
      {RGBtoXYZ_44[1][0], RGBtoXYZ_44[1][1], RGBtoXYZ_44[1][2]},
      {RGBtoXYZ_44[2][0], RGBtoXYZ_44[2][1], RGBtoXYZ_44[2][2]} } };

  // Chromatic adaptation from source white to destination white chromaticity
  // Bradford cone response matrix is the default method
  const array<array <float, 3>, 3> CAT = calculate_cat_matrix( SOURCE_PRIMARIES.white,
                                                DEST_PRIMARIES.white,
                                                coneRespMat );

  const array<array <float, 4>, 4> XYZtoRGB_44 = XYZtoRGB( DEST_PRIMARIES, 1.0);
  const array<array <float, 3>, 3> XYZtoRGB_MAT =
  { { {XYZtoRGB_44[0][0], XYZtoRGB_44[0][1], XYZtoRGB_44[0][2]},
      {XYZtoRGB_44[1][0], XYZtoRGB_44[1][1], XYZtoRGB_44[1][2]},
      {XYZtoRGB_44[2][0], XYZtoRGB_44[2][1], XYZtoRGB_44[2][2]}} };

return mult_f33_f33( RGBtoXYZ_MAT, mult_f33_f33( CAT, XYZtoRGB_MAT));
  return mult_f33_f33( RGBtoXYZ_MAT, mult_f33_f33( CAT, XYZtoRGB_MAT));
}



array<array <float, 3>, 3> calc_sat_adjust_matrix ( float sat, array<float, 3> rgb2Y)
{
  //
  // This function determines the terms for a 3x3 saturation matrix that is
  // based on the luminance of the input.
  //
  array<array <float, 3>, 3> M;
  M[0][0] = (1.0f - sat) * rgb2Y[0] + sat;
  M[1][0] = (1.0f - sat) * rgb2Y[0];
  M[2][0] = (1.0f - sat) * rgb2Y[0];
  
  M[0][1] = (1.0f - sat) * rgb2Y[1];
  M[1][1] = (1.0f - sat) * rgb2Y[1] + sat;
  M[2][1] = (1.0f - sat) * rgb2Y[1];
  
  M[0][2] = (1.0f - sat) * rgb2Y[2];
  M[1][2] = (1.0f - sat) * rgb2Y[2];
  M[2][2] = (1.0f - sat) * rgb2Y[2] + sat;

  M = transpose_f33(M);
  return M;
} 





/* ---- Signal encode/decode functions ---- */

float moncurve_f( float x, float gamma, float offs )
{
  // Forward monitor curve
  float y;
  const float fs = (( gamma - 1.0f) / offs) * powf( offs * gamma / ( ( gamma - 1.0f) * ( 1.0f + offs)), gamma);
  const float xb = offs / ( gamma - 1.0f);
  if ( x >= xb) 
    y = powf( ( x + offs) / ( 1.0f + offs), gamma);
  else
    y = x * fs;
  return y;
}

float moncurve_r( float y, float gamma, float offs )
{
  // Reverse monitor curve
  float x;
  const float yb = powf( offs * gamma / ( ( gamma - 1.0f) * ( 1.0f + offs)), gamma);
  const float rs = powf( ( gamma - 1.0f) / offs, gamma - 1.0f) * powf( ( 1.0f + offs) / gamma, gamma);
  if ( y >= yb) 
    x = ( 1.0f + offs) * powf( y, 1.0f / gamma) - offs;
  else
    x = y * rs;
  return x;
}

array <float, 3> moncurve_f_f3(array <float, 3> x, float gamma, float offs)
{
    array <float, 3> y;
    y[0] = moncurve_f( x[0], gamma, offs);
    y[1] = moncurve_f( x[1], gamma, offs);
    y[2] = moncurve_f( x[2], gamma, offs);
    return y;
}

array <float, 3> moncurve_r_f3(array <float, 3> y, float gamma, float offs)
{
    array <float, 3> x;
    x[0] = moncurve_r( y[0], gamma, offs);
    x[1] = moncurve_r( y[1], gamma, offs);
    x[2] = moncurve_r( y[2], gamma, offs);
    return x;
}

float bt1886_f( float V, float gamma, float Lw, float Lb)
{
  // The reference EOTF specified in Rec. ITU-R BT.1886
  // L = a(max[(V+b),0])^g
  float a = powf( powf( Lw, 1.0f/gamma) - powf( Lb, 1.0f/gamma), gamma);
  float b = powf( Lb, 1.0f/gamma) / ( powf( Lw, 1.0f/gamma) - powf( Lb, 1.0f/gamma));
  float L = a * powf( max( V + b, 0.0f), gamma);
  return L;
}

float bt1886_r( float L, float gamma, float Lw, float Lb)
{
  // The reference EOTF specified in Rec. ITU-R BT.1886
  // L = a(max[(V+b),0])^g
  float a = powf( powf( Lw, 1.0f/gamma) - powf( Lb, 1.0f/gamma), gamma);
  float b = powf( Lb, 1.0f/gamma) / ( powf( Lw, 1.0f/gamma) - powf( Lb, 1.0f/gamma));
  float V = powf( max( L / a, 0.0f), 1.0f/gamma) - b;
  return V;
}

array <float, 3> bt1886_f_f3(array <float, 3> V, float gamma, float Lw, float Lb)
{
    array <float, 3> L;
    L[0] = bt1886_f( V[0], gamma, Lw, Lb);
    L[1] = bt1886_f( V[1], gamma, Lw, Lb);
    L[2] = bt1886_f( V[2], gamma, Lw, Lb);
    return L;
}

array <float, 3> bt1886_r_f3(array <float, 3> L, float gamma, float Lw, float Lb)
{
    array <float, 3> V;
    V[0] = bt1886_r( L[0], gamma, Lw, Lb);
    V[1] = bt1886_r( L[1], gamma, Lw, Lb);
    V[2] = bt1886_r( L[2], gamma, Lw, Lb);
    return V;
}

// SMPTE Range vs Full Range scaling formulas
float smpteRange_to_fullRange( float in)
{
    const float REFBLACK = (  64.0f / 1023.0f);
    const float REFWHITE = ( 940.0f / 1023.0f);

    return (( in - REFBLACK) / ( REFWHITE - REFBLACK));
}

float fullRange_to_smpteRange( float in)
{
    const float REFBLACK = (  64.0f / 1023.0f);
    const float REFWHITE = ( 940.0f / 1023.0f);

    return ( in * ( REFWHITE - REFBLACK) + REFBLACK );
}

array <float, 3> smpteRange_to_fullRange_f3( array <float, 3> rgbIn)
{
    array <float, 3> rgbOut;
    rgbOut[0] = smpteRange_to_fullRange( rgbIn[0]);
    rgbOut[1] = smpteRange_to_fullRange( rgbIn[1]);
    rgbOut[2] = smpteRange_to_fullRange( rgbIn[2]);

    return rgbOut;
}

array <float, 3> fullRange_to_smpteRange_f3( array <float, 3> rgbIn)
{
    array <float, 3> rgbOut;
    rgbOut[0] = fullRange_to_smpteRange( rgbIn[0]);
    rgbOut[1] = fullRange_to_smpteRange( rgbIn[1]);
    rgbOut[2] = fullRange_to_smpteRange( rgbIn[2]);

    return rgbOut;
}


// SMPTE 431-2 defines the DCDM color encoding equations. 
// The equations for the decoding of the encoded color information are the 
// inverse of the encoding equations
// Note: Here the 4095 12-bit scalar is not used since the output of CTL is 0-1.
array <float, 3> dcdm_decode( array <float, 3> XYZp)
{
    array <float, 3> XYZ;
    XYZ[0] = (52.37f/48.0f) * powf( XYZp[0], 2.6f);  
    XYZ[1] = (52.37f/48.0f) * powf( XYZp[1], 2.6f);  
    XYZ[2] = (52.37f/48.0f) * powf( XYZp[2], 2.6f);  

    return XYZ;
}

array <float, 3> dcdm_encode(array <float, 3> XYZ)
{
    array <float, 3> XYZp;
    XYZp[0] = powf( (48.0f/52.37f) * XYZ[0], 1.0f/2.6f);
    XYZp[1] = powf( (48.0f/52.37f) * XYZ[1], 1.0f/2.6f);
    XYZp[2] = powf( (48.0f/52.37f) * XYZ[2], 1.0f/2.6f);

    return XYZp;
}



// Base functions from SMPTE ST 2084-2014

// Constants from SMPTE ST 2084-2014
const float pq_m1 = 0.1593017578125f; // ( 2610.0 / 4096.0 ) / 4.0;
const float pq_m2 = 78.84375f; // ( 2523.0 / 4096.0 ) * 128.0;
const float pq_c1 = 0.8359375f; // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
const float pq_c2 = 18.8515625f; // ( 2413.0 / 4096.0 ) * 32.0;
const float pq_c3 = 18.6875f; // ( 2392.0 / 4096.0 ) * 32.0;

const float pq_C = 10000.0f;

// Converts from the non-linear perceptually quantized space to linear cd/m^2
// Note that this is in float, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex 
// sections of SMPTE ST 2084-2014
float ST2084_2_Y( float N )
{
  // Note that this does NOT handle any of the signal range
  // considerations from 2084 - this assumes full range (0 - 1)
  float Np = powf( N, 1.0f / pq_m2 );
  float L = Np - pq_c1;
  if ( L < 0.0f )
    L = 0.0f;
  L = L / ( pq_c2 - pq_c3 * Np );
  L = powf( L, 1.0f / pq_m1 );
  return L * pq_C; // returns cd/m^2
}

// Converts from linear cd/m^2 to the non-linear perceptually quantized space
// Note that this is in float, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex 
// sections of SMPTE ST 2084-2014
float Y_2_ST2084( float C )
//pq_r
{
  // Note that this does NOT handle any of the signal range
  // considerations from 2084 - this returns full range (0 - 1)
  float L = C / pq_C;
  float Lm = powf( L, pq_m1 );
  float N = ( pq_c1 + pq_c2 * Lm ) / ( 1.0f + pq_c3 * Lm );
  N = powf( N, pq_m2 );
  return N;
}

array <float, 3> Y_2_ST2084_f3(array <float, 3> in)
{
  // converts from linear cd/m^2 to PQ code values
  
    array <float, 3> out;
  out[0] = Y_2_ST2084( in[0]);
  out[1] = Y_2_ST2084( in[1]);
  out[2] = Y_2_ST2084( in[2]);

  return out;
}

array <float, 3> ST2084_2_Y_f3(array <float, 3> in)
{
  // converts from PQ code values to linear cd/m^2
  
    array <float, 3> out;
  out[0] = ST2084_2_Y( in[0]);
  out[1] = ST2084_2_Y( in[1]);
  out[2] = ST2084_2_Y( in[2]);

  return out;
}


// Conversion of PQ signal to HLG, as detailed in Section 7 of ITU-R BT.2390-0
array <float, 3> ST2084_2_HLG_1000nits_f3(array <float, 3> PQ)
{
    // ST.2084 EOTF (non-linear PQ to display light)
    array <float, 3> displayLinear = ST2084_2_Y_f3( PQ);

    // HLG Inverse EOTF (i.e. HLG inverse OOTF followed by the HLG OETF)
    // HLG Inverse OOTF (display linear to scene linear)
    float Y_d = 0.2627f*displayLinear[0] + 0.6780f*displayLinear[1] + 0.0593f*displayLinear[2];
    const float L_w = 1000.0f;
    const float L_b = 0.0f;
    const float alpha = (L_w-L_b);
    const float beta = L_b;
    const float gamma = 1.2f;
    
    array <float, 3> sceneLinear;
    if (Y_d == 0.0f) { 
        /* This case is to protect against powf(0,-N)=Inf error. The ITU document
        does not offer a recommendation for this corner case. There may be a 
        better way to handle this, but for now, this works. 
        */ 
        sceneLinear[0] = 0.0f;
        sceneLinear[1] = 0.0f;
        sceneLinear[2] = 0.0f;        
    } else {
        sceneLinear[0] = powf( (Y_d-beta)/alpha, (1.0f-gamma)/gamma) * ((displayLinear[0]-beta)/alpha);
        sceneLinear[1] = powf( (Y_d-beta)/alpha, (1.0f-gamma)/gamma) * ((displayLinear[1]-beta)/alpha);
        sceneLinear[2] = powf( (Y_d-beta)/alpha, (1.0f-gamma)/gamma) * ((displayLinear[2]-beta)/alpha);
    }

    // HLG OETF (scene linear to non-linear signal value)
    const float a = 0.17883277f;
    const float b = 0.28466892f; // 1.-4.*a;
    const float c = 0.55991073f; // 0.5-a*log(4.*a);

    array <float, 3> HLG;
    if (sceneLinear[0] <= 1./12) {
        HLG[0] = sqrtf(3.0f*sceneLinear[0]);
    } else {
        HLG[0] = a*log(12.0f*sceneLinear[0]-b)+c;
    }
    if (sceneLinear[1] <= 1.0f/12.0f) {
        HLG[1] = sqrtf(3.0f*sceneLinear[1]);
    } else {
        HLG[1] = a*log(12.0f*sceneLinear[1]-b)+c;
    }
    if (sceneLinear[2] <= 1.0f/12.0f) {
        HLG[2] = sqrtf(3.0f*sceneLinear[2]);
    } else {
        HLG[2] = a*log(12.0f*sceneLinear[2]-b)+c;
    }

    return HLG;
}


// Conversion of HLG to PQ signal, as detailed in Section 7 of ITU-R BT.2390-0
array <float, 3> HLG_2_ST2084_1000nits_f3(array <float, 3> HLG)
{
    const float a = 0.17883277f;
    const float b = 0.28466892f; // 1.-4.*a;
    const float c = 0.55991073f; // 0.5-a*log(4.*a);

    const float L_w = 1000.0f;
    const float L_b = 0.0f;
    const float alpha = (L_w-L_b);
    const float beta = L_b;
    const float gamma = 1.2f;

    // HLG EOTF (non-linear signal value to display linear)
    // HLG to scene-linear
    array <float, 3> sceneLinear;
    if ( HLG[0] >= 0.0f && HLG[0] <= 0.5f) {
        sceneLinear[0] = powf(HLG[0],2.0f)/3.0f;
    } else {
        sceneLinear[0] = (expf((HLG[0]-c)/a)+b)/12.0f;
    }        
    if ( HLG[1] >= 0.0f && HLG[1] <= 0.5f) {
        sceneLinear[1] = powf(HLG[1],2.0f)/3.0f;
    } else {
        sceneLinear[1] = (expf((HLG[1]-c)/a)+b)/12.0f;
    }        
    if ( HLG[2] >= 0.0f && HLG[2] <= 0.5f) {
        sceneLinear[2] = powf(HLG[2],2.0f)/3.0f;
    } else {
        sceneLinear[2] = (expf((HLG[2]-c)/a)+b)/12.0f;
    }        
    
    float Y_s = 0.2627f*sceneLinear[0] + 0.6780f*sceneLinear[1] + 0.0593f*sceneLinear[2];

    // Scene-linear to display-linear
    array <float, 3> displayLinear;
    displayLinear[0] = alpha * powf( Y_s, gamma-1.0f) * sceneLinear[0] + beta;
    displayLinear[1] = alpha * powf( Y_s, gamma-1.0f) * sceneLinear[1] + beta;
    displayLinear[2] = alpha * powf( Y_s, gamma-1.0f) * sceneLinear[2] + beta;
        
    // ST.2084 Inverse EOTF
    array <float, 3> PQ = Y_2_ST2084_f3( displayLinear);

    return PQ;
}