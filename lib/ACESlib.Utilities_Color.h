
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
  { 0.73470,  0.26530},
  { 0.00000,  1.00000},
  { 0.00010, -0.07700},
  { 0.32168,  0.33767}
};

const Chromaticities AP1 = // Working space and rendering primaries for ACES 1.0
{
  { 0.713,    0.293},
  { 0.165,    0.830},
  { 0.128,    0.044},
  { 0.32168,  0.33767}
};

const Chromaticities REC709_PRI =
{
  { 0.64000,  0.33000},
  { 0.30000,  0.60000},
  { 0.15000,  0.06000},
  { 0.31270,  0.32900}
};

const Chromaticities P3D60_PRI =
{
  { 0.68000,  0.32000},
  { 0.26500,  0.69000},
  { 0.15000,  0.06000},
  { 0.32168,  0.33767}
};

const Chromaticities P3D65_PRI =
{
  { 0.68000,  0.32000},
  { 0.26500,  0.69000},
  { 0.15000,  0.06000},
  { 0.31270,  0.32900}
};

const Chromaticities P3DCI_PRI =
{
  { 0.68000,  0.32000},
  { 0.26500,  0.69000},
  { 0.15000,  0.06000},
  { 0.31400,  0.35100}
};

const Chromaticities ARRI_ALEXA_WG_PRI =
{
  { 0.68400,  0.31300},
  { 0.22100,  0.84800},
  { 0.08610, -0.10200},
  { 0.31270,  0.32900}
};

const Chromaticities REC2020_PRI = 
{
  { 0.70800,  0.29200},
  { 0.17000,  0.79700},
  { 0.13100,  0.04600},
  { 0.31270,  0.32900}
};

const Chromaticities RIMMROMM_PRI = 
{
  { 0.7347,  0.2653},
  { 0.1596,  0.8404},
  { 0.0366,  0.0001},
  { 0.3457,  0.3585}
};

const Chromaticities SONY_SGAMUT3_PRI =
{
  { 0.730,  0.280},
  { 0.140,  0.855},
  { 0.100, -0.050},
  { 0.3127,  0.3290}
};

const Chromaticities SONY_SGAMUT3_CINE_PRI =
{
  { 0.766,  0.275},
  { 0.225,  0.800},
  { 0.089, -0.087},
  { 0.3127,  0.3290}
};

// Note: No official published primaries exist as of this day for the
// Sony VENICE SGamut3 and Sony VENICE SGamut3.Cine colorspaces. The primaries
// have thus been derived from the IDT matrices.
const Chromaticities SONY_VENICE_SGAMUT3_PRI =
{
  { 0.740464264304292,  0.279364374750660},
  { 0.089241145423286,  0.893809528608105},
  { 0.110488236673827, -0.052579333080476},
  { 0.312700000000000,  0.329000000000000}
};

const Chromaticities SONY_VENICE_SGAMUT3_CINE_PRI =
{
  { 0.775901871567345,  0.274502392854799},
  { 0.188682902773355,  0.828684937020288},
  { 0.101337382499301, -0.089187517306263},
  { 0.312700000000000,  0.329000000000000}
};

const Chromaticities CANON_CGAMUT_PRI =
{
  { 0.7400,  0.2700},
  { 0.1700,  1.1400},
  { 0.0800, -0.1000},
  { 0.3127,  0.3290}
};

const Chromaticities RED_WIDEGAMUTRGB_PRI =
{
  { 0.780308,  0.304253},
  { 0.121595,  1.493994},
  { 0.095612, -0.084589},
  { 0.3127,  0.3290}
};

const Chromaticities PANASONIC_VGAMUT_PRI =
{
  { 0.730,  0.280},
  { 0.165,  0.840},
  { 0.100, -0.030},
  { 0.3127,  0.3290}
};




/* ---- Conversion Functions ---- */
// Various transformations between color encodings and data representations
//

// Transformations between CIE XYZ tristimulus values and CIE x,y 
// chromaticity coordinates
array <double,3> XYZ_2_xyY( array<double,3> XYZ)
{  
  array<double,3> xyY;
  double divisor = (XYZ[0] + XYZ[1] + XYZ[2]);
  if (divisor == 0.0) divisor = (double)1e-10;
  xyY[0] = XYZ[0] / divisor;
  xyY[1] = XYZ[1] / divisor;  
  xyY[2] = XYZ[1];
  
  return xyY;
}

array <double,3> xyY_2_XYZ( array<double,3> xyY)
{
  array <double,3> XYZ;
  XYZ[0] = xyY[0] * xyY[2] / max( xyY[1], (double)1e-10);
  XYZ[1] = xyY[2];  
  XYZ[2] = (1.0 - xyY[0] - xyY[1]) * xyY[2] / max( xyY[1], (double)1e-10);

  return XYZ;
}


// Transformations from RGB to other color representations
double rgb_2_hue( array<double,3> rgb) 
{
  // Returns a geometric hue angle in degrees (0-360) based on RGB values.
  // For neutral colors, hue is undefined and the function will return a quiet NaN value.
  double hue;
  if (rgb[0] == rgb[1] && rgb[1] == rgb[2]) {
    hue = (double)FLT_NAN; // RGB triplets where RGB are equal have an undefined hue
  } else {
    hue = (double)(180.0/M_PI) * atan2( sqrtf(3.0)*(rgb[1]-rgb[2]), 2.0*rgb[0]-rgb[1]-rgb[2]);
  }
    
  if (hue < 0.0) hue = hue + 360.0;

  return hue;
}

double rgb_2_yc( array <double,3> rgb, double ycRadiusWeight = 1.75)
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

  double r = rgb[0]; 
  double g = rgb[1]; 
  double b = rgb[2];
  
  double chroma = sqrtf(b*(b-g)+g*(g-r)+r*(r-b));

  return ( b + g + r + ycRadiusWeight * chroma) / 3.0;
}



/* ---- Chromatic Adaptation ---- */

const array<array <double, 3>, 3> CONE_RESP_MAT_BRADFORD = { {
  { 0.89510, -0.75020,  0.03890},
  { 0.26640,  1.71350, -0.06850},
  {-0.16140,  0.03670,  1.02960}
} };

const array<array <double, 3>, 3> CONE_RESP_MAT_CAT02 = { {
  { 0.73280, -0.70360,  0.00300},
  { 0.42960,  1.69750,  0.01360},
  {-0.16240,  0.00610,  0.98340}
} };

array<array <double, 3>, 3> calculate_cat_matrix
  (const array <double, 2> &src_xy,         // x,y chromaticity of source white
      const array <double, 2> &des_xy,         // x,y chromaticity of destination white
      const array<array <double, 3>, 3> coneRespMat = CONE_RESP_MAT_CAT02 //Changed Default to CAT02
  )
{
  //
  // Calculates and returns a 3x3 Von Kries chromatic adaptation transform 
  // from src_xy to des_xy using the cone response primaries defined 
  // by coneRespMat. By default, coneRespMat is set to CONE_RESP_MAT_CAT02. 
  // The default coneRespMat can be overridden at runtime. 
  //

  const array <double, 3> src_xyY = { src_xy[0], src_xy[1], 1. };
  const array <double, 3> des_xyY = { des_xy[0], des_xy[1], 1. };

  array <double, 3> src_XYZ = xyY_2_XYZ( src_xyY );
  array <double, 3> des_XYZ = xyY_2_XYZ( des_xyY );

  array <double, 3> src_coneResp = mult_f3_f33( src_XYZ, coneRespMat);
  array <double, 3> des_coneResp = mult_f3_f33( des_XYZ, coneRespMat);

  array<array <double, 3>, 3> vkMat = { {
      { des_coneResp[0] / src_coneResp[0], 0.0, 0.0 },
      { 0.0, des_coneResp[1] / src_coneResp[1], 0.0 },
      { 0.0, 0.0, des_coneResp[2] / src_coneResp[2] }
  } };

  array<array <double, 3>, 3> cat_matrix = mult_f33_f33( coneRespMat, mult_f33_f33( vkMat, invert_f33( coneRespMat ) ) );

  return cat_matrix;
}



array<array <double, 3>, 3> calculate_rgb_to_rgb_matrix
(Chromaticities SOURCE_PRIMARIES,
    Chromaticities DEST_PRIMARIES,
    const array<array <double, 3>, 3> coneRespMat = CONE_RESP_MAT_CAT02 //Changed Default to CAT02
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
    // By default, coneRespMat is set to CONE_RESP_MAT_CAT02. 
    // The default coneRespMat can be overridden at runtime. 
    //

    const array<array <double, 4>, 4> RGBtoXYZ_44 = RGBtoXYZ(SOURCE_PRIMARIES, 1.0);
    const array<array <double, 3>, 3> RGBtoXYZ_MAT =
    { { {RGBtoXYZ_44[0][0], RGBtoXYZ_44[0][1], RGBtoXYZ_44[0][2]},
      {RGBtoXYZ_44[1][0], RGBtoXYZ_44[1][1], RGBtoXYZ_44[1][2]},
      {RGBtoXYZ_44[2][0], RGBtoXYZ_44[2][1], RGBtoXYZ_44[2][2]} } };

  // Chromatic adaptation from source white to destination white chromaticity
  // Bradford cone response matrix is the default method
  const array<array <double, 3>, 3> CAT = calculate_cat_matrix( SOURCE_PRIMARIES.white,
                                                DEST_PRIMARIES.white,
                                                coneRespMat );

  const array<array <double, 4>, 4> XYZtoRGB_44 = XYZtoRGB( DEST_PRIMARIES, 1.0);
  const array<array <double, 3>, 3> XYZtoRGB_MAT =
  { { {XYZtoRGB_44[0][0], XYZtoRGB_44[0][1], XYZtoRGB_44[0][2]},
      {XYZtoRGB_44[1][0], XYZtoRGB_44[1][1], XYZtoRGB_44[1][2]},
      {XYZtoRGB_44[2][0], XYZtoRGB_44[2][1], XYZtoRGB_44[2][2]}} };

return mult_f33_f33( RGBtoXYZ_MAT, mult_f33_f33( CAT, XYZtoRGB_MAT));
  return mult_f33_f33( RGBtoXYZ_MAT, mult_f33_f33( CAT, XYZtoRGB_MAT));
}



array<array <double, 3>, 3> calc_sat_adjust_matrix ( double sat, array<double, 3> rgb2Y)
{
  //
  // This function determines the terms for a 3x3 saturation matrix that is
  // based on the luminance of the input.
  //
  array<array <double, 3>, 3> M;
  M[0][0] = (1.0 - sat) * rgb2Y[0] + sat;
  M[1][0] = (1.0 - sat) * rgb2Y[0];
  M[2][0] = (1.0 - sat) * rgb2Y[0];
  
  M[0][1] = (1.0 - sat) * rgb2Y[1];
  M[1][1] = (1.0 - sat) * rgb2Y[1] + sat;
  M[2][1] = (1.0 - sat) * rgb2Y[1];
  
  M[0][2] = (1.0 - sat) * rgb2Y[2];
  M[1][2] = (1.0 - sat) * rgb2Y[2];
  M[2][2] = (1.0 - sat) * rgb2Y[2] + sat;

  M = transpose_f33(M);
  return M;
} 





/* ---- Signal encode/decode functions ---- */

double moncurve_f( double x, double gamma, double offs )
{
  // Forward monitor curve
  double y;
  const double fs = (( gamma - 1.0) / offs) * powf( offs * gamma / ( ( gamma - 1.0) * ( 1.0 + offs)), gamma);
  const double xb = offs / ( gamma - 1.0);
  if ( x >= xb) 
    y = powf( ( x + offs) / ( 1.0 + offs), gamma);
  else
    y = x * fs;
  return y;
}

double moncurve_r( double y, double gamma, double offs )
{
  // Reverse monitor curve
  double x;
  const double yb = powf( offs * gamma / ( ( gamma - 1.0) * ( 1.0 + offs)), gamma);
  const double rs = powf( ( gamma - 1.0) / offs, gamma - 1.0) * powf( ( 1.0 + offs) / gamma, gamma);
  if ( y >= yb) 
    x = ( 1.0 + offs) * powf( y, 1.0 / gamma) - offs;
  else
    x = y * rs;
  return x;
}

array <double, 3> moncurve_f_f3(array <double, 3> x, double gamma, double offs)
{
    array <double, 3> y;
    y[0] = moncurve_f( x[0], gamma, offs);
    y[1] = moncurve_f( x[1], gamma, offs);
    y[2] = moncurve_f( x[2], gamma, offs);
    return y;
}

array <double, 3> moncurve_r_f3(array <double, 3> y, double gamma, double offs)
{
    array <double, 3> x;
    x[0] = moncurve_r( y[0], gamma, offs);
    x[1] = moncurve_r( y[1], gamma, offs);
    x[2] = moncurve_r( y[2], gamma, offs);
    return x;
}

double bt1886_f( double V, double gamma, double Lw, double Lb)
{
  // The reference EOTF specified in Rec. ITU-R BT.1886
  // L = a(max[(V+b),0])^g
  double a = powf( powf( Lw, 1.0/gamma) - powf( Lb, 1.0/gamma), gamma);
  double b = powf( Lb, 1.0/gamma) / ( powf( Lw, 1.0/gamma) - powf( Lb, 1.0/gamma));
  double L = a * powf( max( V + b, 0.0), gamma);
  return L;
}

double bt1886_r( double L, double gamma, double Lw, double Lb)
{
  // The reference EOTF specified in Rec. ITU-R BT.1886
  // L = a(max[(V+b),0])^g
  double a = powf( powf( Lw, 1.0/gamma) - powf( Lb, 1.0/gamma), gamma);
  double b = powf( Lb, 1.0/gamma) / ( powf( Lw, 1.0/gamma) - powf( Lb, 1.0/gamma));
  double V = powf( max( L / a, 0.0), 1.0/gamma) - b;
  return V;
}

array <double, 3> bt1886_f_f3(array <double, 3> V, double gamma, double Lw, double Lb)
{
    array <double, 3> L;
    L[0] = bt1886_f( V[0], gamma, Lw, Lb);
    L[1] = bt1886_f( V[1], gamma, Lw, Lb);
    L[2] = bt1886_f( V[2], gamma, Lw, Lb);
    return L;
}

array <double, 3> bt1886_r_f3(array <double, 3> L, double gamma, double Lw, double Lb)
{
    array <double, 3> V;
    V[0] = bt1886_r( L[0], gamma, Lw, Lb);
    V[1] = bt1886_r( L[1], gamma, Lw, Lb);
    V[2] = bt1886_r( L[2], gamma, Lw, Lb);
    return V;
}

// SMPTE Range vs Full Range scaling formulas
double smpteRange_to_fullRange( double in)
{
    const double REFBLACK = (  64.0 / 1023.0);
    const double REFWHITE = ( 940.0 / 1023.0);

    return (( in - REFBLACK) / ( REFWHITE - REFBLACK));
}

double fullRange_to_smpteRange( double in)
{
    const double REFBLACK = (  64.0 / 1023.0);
    const double REFWHITE = ( 940.0 / 1023.0);

    return ( in * ( REFWHITE - REFBLACK) + REFBLACK );
}

array <double, 3> smpteRange_to_fullRange_f3( array <double, 3> rgbIn)
{
    array <double, 3> rgbOut;
    rgbOut[0] = smpteRange_to_fullRange( rgbIn[0]);
    rgbOut[1] = smpteRange_to_fullRange( rgbIn[1]);
    rgbOut[2] = smpteRange_to_fullRange( rgbIn[2]);

    return rgbOut;
}

array <double, 3> fullRange_to_smpteRange_f3( array <double, 3> rgbIn)
{
    array <double, 3> rgbOut;
    rgbOut[0] = fullRange_to_smpteRange( rgbIn[0]);
    rgbOut[1] = fullRange_to_smpteRange( rgbIn[1]);
    rgbOut[2] = fullRange_to_smpteRange( rgbIn[2]);

    return rgbOut;
}


// SMPTE 431-2 defines the DCDM color encoding equations. 
// The equations for the decoding of the encoded color information are the 
// inverse of the encoding equations
// Note: Here the 4095 12-bit scalar is not used since the output of CTL is 0-1.
array <double, 3> dcdm_decode( array <double, 3> XYZp)
{
    array <double, 3> XYZ;
    XYZ[0] = (52.37/48.0) * powf( XYZp[0], 2.6);  
    XYZ[1] = (52.37/48.0) * powf( XYZp[1], 2.6);  
    XYZ[2] = (52.37/48.0) * powf( XYZp[2], 2.6);  

    return XYZ;
}

array <double, 3> dcdm_encode(array <double, 3> XYZ)
{
    array <double, 3> XYZp;
    XYZp[0] = powf( (48.0/52.37) * XYZ[0], 1.0/2.6);
    XYZp[1] = powf( (48.0/52.37) * XYZ[1], 1.0/2.6);
    XYZp[2] = powf( (48.0/52.37) * XYZ[2], 1.0/2.6);

    return XYZp;
}



// Base functions from SMPTE ST 2084-2014

// Constants from SMPTE ST 2084-2014
const double pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
const double pq_m2 = 78.84375; // ( 2523.0 / 4096.0 ) * 128.0;
const double pq_c1 = 0.8359375; // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
const double pq_c2 = 18.8515625; // ( 2413.0 / 4096.0 ) * 32.0;
const double pq_c3 = 18.6875; // ( 2392.0 / 4096.0 ) * 32.0;

const double pq_C = 10000.0;

// Converts from the non-linear perceptually quantized space to linear cd/m^2
// Note that this is in double, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex 
// sections of SMPTE ST 2084-2014
double ST2084_2_Y( double N )
{
  // Note that this does NOT handle any of the signal range
  // considerations from 2084 - this assumes full range (0 - 1)
  double Np = powf( N, 1.0 / pq_m2 );
  double L = Np - pq_c1;
  if ( L < 0.0 )
    L = 0.0;
  L = L / ( pq_c2 - pq_c3 * Np );
  L = powf( L, 1.0 / pq_m1 );
  return L * pq_C; // returns cd/m^2
}

// Converts from linear cd/m^2 to the non-linear perceptually quantized space
// Note that this is in double, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex 
// sections of SMPTE ST 2084-2014
double Y_2_ST2084( double C )
//pq_r
{
  // Note that this does NOT handle any of the signal range
  // considerations from 2084 - this returns full range (0 - 1)
  double L = C / pq_C;
  double Lm = powf( L, pq_m1 );
  double N = ( pq_c1 + pq_c2 * Lm ) / ( 1.0 + pq_c3 * Lm );
  N = powf( N, pq_m2 );
  return N;
}

array <double, 3> Y_2_ST2084_f3(array <double, 3> in)
{
  // converts from linear cd/m^2 to PQ code values
  
    array <double, 3> out;
  out[0] = Y_2_ST2084( in[0]);
  out[1] = Y_2_ST2084( in[1]);
  out[2] = Y_2_ST2084( in[2]);

  return out;
}

array <double, 3> ST2084_2_Y_f3(array <double, 3> in)
{
  // converts from PQ code values to linear cd/m^2
  
    array <double, 3> out;
  out[0] = ST2084_2_Y( in[0]);
  out[1] = ST2084_2_Y( in[1]);
  out[2] = ST2084_2_Y( in[2]);

  return out;
}


// Conversion of PQ signal to HLG, as detailed in Section 7 of ITU-R BT.2390-0
array <double, 3> ST2084_2_HLG_1000nits_f3(array <double, 3> PQ)
{
    // ST.2084 EOTF (non-linear PQ to display light)
    array <double, 3> displayLinear = ST2084_2_Y_f3( PQ);

    // HLG Inverse EOTF (i.e. HLG inverse OOTF followed by the HLG OETF)
    // HLG Inverse OOTF (display linear to scene linear)
    double Y_d = 0.2627*displayLinear[0] + 0.6780*displayLinear[1] + 0.0593*displayLinear[2];
    const double L_w = 1000.0;
    const double L_b = 0.0;
    const double alpha = (L_w-L_b);
    const double beta = L_b;
    const double gamma = 1.2;
    
    array <double, 3> sceneLinear;
    if (Y_d == 0.0) { 
        /* This case is to protect against powf(0,-N)=Inf error. The ITU document
        does not offer a recommendation for this corner case. There may be a 
        better way to handle this, but for now, this works. 
        */ 
        sceneLinear[0] = 0.0;
        sceneLinear[1] = 0.0;
        sceneLinear[2] = 0.0;        
    } else {
        sceneLinear[0] = powf( (Y_d-beta)/alpha, (1.0-gamma)/gamma) * ((displayLinear[0]-beta)/alpha);
        sceneLinear[1] = powf( (Y_d-beta)/alpha, (1.0-gamma)/gamma) * ((displayLinear[1]-beta)/alpha);
        sceneLinear[2] = powf( (Y_d-beta)/alpha, (1.0-gamma)/gamma) * ((displayLinear[2]-beta)/alpha);
    }

    // HLG OETF (scene linear to non-linear signal value)
    const double a = 0.17883277;
    const double b = 0.28466892; // 1.-4.*a;
    const double c = 0.55991073; // 0.5-a*log(4.*a);

    array <double, 3> HLG;
    if (sceneLinear[0] <= 1./12) {
        HLG[0] = sqrtf(3.0*sceneLinear[0]);
    } else {
        HLG[0] = a*log(12.0*sceneLinear[0]-b)+c;
    }
    if (sceneLinear[1] <= 1.0/12.0) {
        HLG[1] = sqrtf(3.0*sceneLinear[1]);
    } else {
        HLG[1] = a*log(12.0*sceneLinear[1]-b)+c;
    }
    if (sceneLinear[2] <= 1.0/12.0) {
        HLG[2] = sqrtf(3.0*sceneLinear[2]);
    } else {
        HLG[2] = a*log(12.0*sceneLinear[2]-b)+c;
    }

    return HLG;
}


// Conversion of HLG to PQ signal, as detailed in Section 7 of ITU-R BT.2390-0
array <double, 3> HLG_2_ST2084_1000nits_f3(array <double, 3> HLG)
{
    const double a = 0.17883277;
    const double b = 0.28466892; // 1.-4.*a;
    const double c = 0.55991073; // 0.5-a*log(4.*a);

    const double L_w = 1000.0;
    const double L_b = 0.0;
    const double alpha = (L_w-L_b);
    const double beta = L_b;
    const double gamma = 1.2;

    // HLG EOTF (non-linear signal value to display linear)
    // HLG to scene-linear
    array <double, 3> sceneLinear;
    if ( HLG[0] >= 0.0 && HLG[0] <= 0.5) {
        sceneLinear[0] = powf(HLG[0],2.0)/3.0;
    } else {
        sceneLinear[0] = (expf((HLG[0]-c)/a)+b)/12.0;
    }        
    if ( HLG[1] >= 0.0 && HLG[1] <= 0.5) {
        sceneLinear[1] = powf(HLG[1],2.0)/3.0;
    } else {
        sceneLinear[1] = (expf((HLG[1]-c)/a)+b)/12.0;
    }        
    if ( HLG[2] >= 0.0 && HLG[2] <= 0.5) {
        sceneLinear[2] = powf(HLG[2],2.0)/3.0;
    } else {
        sceneLinear[2] = (expf((HLG[2]-c)/a)+b)/12.0;
    }        
    
    double Y_s = 0.2627*sceneLinear[0] + 0.6780*sceneLinear[1] + 0.0593*sceneLinear[2];

    // Scene-linear to display-linear
    array <double, 3> displayLinear;
    displayLinear[0] = alpha * powf( Y_s, gamma-1.0) * sceneLinear[0] + beta;
    displayLinear[1] = alpha * powf( Y_s, gamma-1.0) * sceneLinear[1] + beta;
    displayLinear[2] = alpha * powf( Y_s, gamma-1.0) * sceneLinear[2] + beta;
        
    // ST.2084 Inverse EOTF
    array <double, 3> PQ = Y_2_ST2084_f3( displayLinear);

    return PQ;
}