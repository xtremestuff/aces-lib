
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.Tonescales.a1.0.3</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - Tonescales</ACESuserName>

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.Utilities.h"



// Textbook monomial to basis-function conversion matrix.
const array <array <float, 3>, 3> M = { {
  {  0.5f, -1.0f, 0.5f },
  { -1.0f,  1.0f, 0.5f },
  {  0.5f,  0.0f, 0.0f }
} };



struct SplineMapPoint
{
  float x;
  float y;
};

struct SegmentedSplineParams_c5
{
  array <float, 6> coefsLow;    // coefs for B-spline between minPoint and midPoint (units of log luminance)
  array <float, 6> coefsHigh;   // coefs for B-spline between midPoint and maxPoint (units of log luminance)
  SplineMapPoint minPoint; // {luminance, luminance} linear extension below this
  SplineMapPoint midPoint; // {luminance, luminance} 
  SplineMapPoint maxPoint; // {luminance, luminance} linear extension above this
  float slopeLow;       // log-log slope of low linear extension
  float slopeHigh;      // log-log slope of high linear extension
};

struct SegmentedSplineParams_c9
{
  array <float, 10> coefsLow;    // coefs for B-spline between minPoint and midPoint (units of log luminance)
  array <float, 10> coefsHigh;   // coefs for B-spline between midPoint and maxPoint (units of log luminance)
  SplineMapPoint minPoint; // {luminance, luminance} linear extension below this
  SplineMapPoint midPoint; // {luminance, luminance} 
  SplineMapPoint maxPoint; // {luminance, luminance} linear extension above this
  float slopeLow;       // log-log slope of low linear extension
  float slopeHigh;      // log-log slope of high linear extension
};


const SegmentedSplineParams_c5 RRT_PARAMS =
{
  // coefsLow[6]
  { -4.0000000000f, -4.0000000000f, -3.1573765773f, -0.4852499958f, 1.8477324706f, 1.8477324706f },
  // coefsHigh[6]
  { -0.7185482425f, 2.0810307172f, 3.6681241237f, 4.0000000000f, 4.0000000000f, 4.0000000000f },
  { 0.18f*powf(2.0f,-15.0f), 0.0001f},    // minPoint
  { 0.18f,                4.8f},    // midPoint  
  { 0.18f*powf(2.0f, 18.0f), 10000.0f},    // maxPoint
  0.0f,  // slopeLow
  0.0f   // slopeHigh
};


float segmented_spline_c5_fwd
  ( 
    float x,
    SegmentedSplineParams_c5 C = RRT_PARAMS
  )
{
  const int N_KNOTS_LOW = 4;
  const int N_KNOTS_HIGH = 4;

  // Check for negatives or zero before taking the log. If negative or zero,
  // set to HALF_MIN.
  float logx = log10f( max(x, HALF_MIN )); 

  float logy;

  if ( logx <= log10f(C.minPoint.x) ) { 

    logy = logx * C.slopeLow + ( log10f(C.minPoint.y) - C.slopeLow * log10f(C.minPoint.x) );

  } else if (( logx > log10f(C.minPoint.x) ) && ( logx < log10f(C.midPoint.x) )) {

    float knot_coord = (N_KNOTS_LOW-1) * (logx-log10f(C.minPoint.x))/(log10f(C.midPoint.x)-log10f(C.minPoint.x));
    int j = (int)knot_coord;
    float t = knot_coord - j;

    array <float, 3> cf = { C.coefsLow[ j], C.coefsLow[ j + 1], C.coefsLow[ j + 2]};
    // NOTE: If the running a version of CTL < 1.5, you may get an 
    // exception thrown error, usually accompanied by "Array index out of range" 
    // If you receive this error, it is recommended that you update to CTL v1.5, 
    // which contains a number of important bug fixes. Otherwise, you may try 
    // uncommenting the below, which is longer, but equivalent to, the above 
    // line of code.
    //
    // float cf[ 3];
    // if ( j <= 0) {
    //     cf[ 0] = C.coefsLow[0];  cf[ 1] = C.coefsLow[1];  cf[ 2] = C.coefsLow[2];
    // } else if ( j == 1) {
    //     cf[ 0] = C.coefsLow[1];  cf[ 1] = C.coefsLow[2];  cf[ 2] = C.coefsLow[3];
    // } else if ( j == 2) {
    //     cf[ 0] = C.coefsLow[2];  cf[ 1] = C.coefsLow[3];  cf[ 2] = C.coefsLow[4];
    // } else if ( j == 3) {
    //     cf[ 0] = C.coefsLow[3];  cf[ 1] = C.coefsLow[4];  cf[ 2] = C.coefsLow[5];
    // } else if ( j == 4) {
    //     cf[ 0] = C.coefsLow[4];  cf[ 1] = C.coefsLow[5];  cf[ 2] = C.coefsLow[6];
    // } else if ( j == 5) {
    //     cf[ 0] = C.coefsLow[5];  cf[ 1] = C.coefsLow[6];  cf[ 2] = C.coefsLow[7];
    // } else if ( j == 6) {
    //     cf[ 0] = C.coefsLow[6];  cf[ 1] = C.coefsLow[7];  cf[ 2] = C.coefsLow[8];
    // } 
    
    array <float, 3> monomials = { t * t, t, 1. };
    logy = dot_f3_f3( monomials, mult_f3_f33( cf, M));

  } else if (( logx >= log10f(C.midPoint.x) ) && ( logx < log10f(C.maxPoint.x) )) {

    float knot_coord = (N_KNOTS_HIGH-1) * (logx-log10f(C.midPoint.x))/(log10f(C.maxPoint.x)-log10f(C.midPoint.x));
    int j = knot_coord;
    float t = knot_coord - j;

    array <float, 3> cf = { C.coefsHigh[ j], C.coefsHigh[ j + 1], C.coefsHigh[ j + 2]};
    // NOTE: If the running a version of CTL < 1.5, you may get an 
    // exception thrown error, usually accompanied by "Array index out of range" 
    // If you receive this error, it is recommended that you update to CTL v1.5, 
    // which contains a number of important bug fixes. Otherwise, you may try 
    // uncommenting the below, which is longer, but equivalent to, the above 
    // line of code.
    //
    // float cf[ 3];
    // if ( j <= 0) {
    //     cf[ 0] = C.coefsHigh[0];  cf[ 1] = C.coefsHigh[1];  cf[ 2] = C.coefsHigh[2];
    // } else if ( j == 1) {
    //     cf[ 0] = C.coefsHigh[1];  cf[ 1] = C.coefsHigh[2];  cf[ 2] = C.coefsHigh[3];
    // } else if ( j == 2) {
    //     cf[ 0] = C.coefsHigh[2];  cf[ 1] = C.coefsHigh[3];  cf[ 2] = C.coefsHigh[4];
    // } else if ( j == 3) {
    //     cf[ 0] = C.coefsHigh[3];  cf[ 1] = C.coefsHigh[4];  cf[ 2] = C.coefsHigh[5];
    // } else if ( j == 4) {
    //     cf[ 0] = C.coefsHigh[4];  cf[ 1] = C.coefsHigh[5];  cf[ 2] = C.coefsHigh[6];
    // } else if ( j == 5) {
    //     cf[ 0] = C.coefsHigh[5];  cf[ 1] = C.coefsHigh[6];  cf[ 2] = C.coefsHigh[7];
    // } else if ( j == 6) {
    //     cf[ 0] = C.coefsHigh[6];  cf[ 1] = C.coefsHigh[7];  cf[ 2] = C.coefsHigh[8];
    // } 

    array <float, 3> monomials = { t * t, t, 1. };
    logy = dot_f3_f3( monomials, mult_f3_f33( cf, M));

  } else { //if ( logIn >= log10f(C.maxPoint.x) ) { 

    logy = logx * C.slopeHigh + ( log10f(C.maxPoint.y) - C.slopeHigh * log10f(C.maxPoint.x) );

  }

  return pow10f(logy);
  
}


float segmented_spline_c5_rev
  ( 
    float y,
    SegmentedSplineParams_c5 C = RRT_PARAMS
  )
{  
  const int N_KNOTS_LOW = 4;
  const int N_KNOTS_HIGH = 4;

  const float KNOT_INC_LOW = (log10f(C.midPoint.x) - log10f(C.minPoint.x)) / (N_KNOTS_LOW - 1.);
  const float KNOT_INC_HIGH = (log10f(C.maxPoint.x) - log10f(C.midPoint.x)) / (N_KNOTS_HIGH - 1.);
  
  // KNOT_Y is luminance of the spline at each knot
  float KNOT_Y_LOW[ N_KNOTS_LOW];
  for (int i = 0; i < N_KNOTS_LOW; i = i+1) {
    KNOT_Y_LOW[ i] = ( C.coefsLow[i] + C.coefsLow[i+1]) / 2.;
  };

  float KNOT_Y_HIGH[ N_KNOTS_HIGH];
  for (int i = 0; i < N_KNOTS_HIGH; i = i+1) {
    KNOT_Y_HIGH[ i] = ( C.coefsHigh[i] + C.coefsHigh[i+1]) / 2.;
  };

  float logy = log10f( max(y,1e-10));

  float logx;
  if (logy <= log10f(C.minPoint.y)) {

    logx = log10f(C.minPoint.x);

  } else if ( (logy > log10f(C.minPoint.y)) && (logy <= log10f(C.midPoint.y)) ) {

    unsigned int j;
    array <float, 3> cf;
    if ( logy > KNOT_Y_LOW[ 0] && logy <= KNOT_Y_LOW[ 1]) {
        cf[ 0] = C.coefsLow[0];  cf[ 1] = C.coefsLow[1];  cf[ 2] = C.coefsLow[2];  j = 0;
    } else if ( logy > KNOT_Y_LOW[ 1] && logy <= KNOT_Y_LOW[ 2]) {
        cf[ 0] = C.coefsLow[1];  cf[ 1] = C.coefsLow[2];  cf[ 2] = C.coefsLow[3];  j = 1;
    } else if ( logy > KNOT_Y_LOW[ 2] && logy <= KNOT_Y_LOW[ 3]) {
        cf[ 0] = C.coefsLow[2];  cf[ 1] = C.coefsLow[3];  cf[ 2] = C.coefsLow[4];  j = 2;
    } 
    
    const array <float, 3> tmp = mult_f3_f33( cf, M);

    float a = tmp[ 0];
    float b = tmp[ 1];
    float c = tmp[ 2];
    c = c - logy;

    const float d = sqrtf( b * b - 4.0f * a * c);

    const float t = ( 2.0f * c) / ( -d - b);

    logx = log10f(C.minPoint.x) + ( t + j) * KNOT_INC_LOW;

  } else if ( (logy > log10f(C.midPoint.y)) && (logy < log10f(C.maxPoint.y)) ) {

    unsigned int j;
    array <float, 3> cf;
    if ( logy > KNOT_Y_HIGH[ 0] && logy <= KNOT_Y_HIGH[ 1]) {
        cf[ 0] = C.coefsHigh[0];  cf[ 1] = C.coefsHigh[1];  cf[ 2] = C.coefsHigh[2];  j = 0;
    } else if ( logy > KNOT_Y_HIGH[ 1] && logy <= KNOT_Y_HIGH[ 2]) {
        cf[ 0] = C.coefsHigh[1];  cf[ 1] = C.coefsHigh[2];  cf[ 2] = C.coefsHigh[3];  j = 1;
    } else if ( logy > KNOT_Y_HIGH[ 2] && logy <= KNOT_Y_HIGH[ 3]) {
        cf[ 0] = C.coefsHigh[2];  cf[ 1] = C.coefsHigh[3];  cf[ 2] = C.coefsHigh[4];  j = 2;
    } 
    
    const array <float, 3> tmp = mult_f3_f33( cf, M);

    float a = tmp[ 0];
    float b = tmp[ 1];
    float c = tmp[ 2];
    c = c - logy;

    const float d = sqrtf( b * b - 4.0f * a * c);

    const float t = ( 2.0f * c) / ( -d - b);

    logx = log10f(C.midPoint.x) + ( t + j) * KNOT_INC_HIGH;

  } else { //if ( logy >= log10f(C.maxPoint.y) ) {

    logx = log10f(C.maxPoint.x);

  }
  
  return pow10f( logx);

}






const SegmentedSplineParams_c9 ODT_48nits =
{
  // coefsLow[10]
  { -1.6989700043f, -1.6989700043f, -1.4779000000f, -1.2291000000f, -0.8648000000f, -0.4480000000f, 0.0051800000f, 0.4511080334f, 0.9113744414f, 0.9113744414f},
  // coefsHigh[10]
  { 0.5154386965f, 0.8470437783f, 1.1358000000f, 1.3802000000f, 1.5197000000f, 1.5985000000f, 1.6467000000f, 1.6746091357f, 1.6878733390f, 1.6878733390f },
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,-6.5f) ),  0.02f},    // minPoint
  {segmented_spline_c5_fwd( 0.18f ),                4.8f},    // midPoint  
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,6.5f) ),   48.0f},    // maxPoint
  0.0f,  // slopeLow
  0.04f  // slopeHigh
};

const SegmentedSplineParams_c9 ODT_1000nits =
{
  // coefsLow[10]
  { -4.9706219331f, -3.0293780669f, -2.1262f, -1.5105f, -1.0578f, -0.4668f, 0.11938f, 0.7088134201f, 1.2911865799f, 1.2911865799f },
  // coefsHigh[10]
  { 0.8089132070f, 1.1910867930f, 1.5683f, 1.9483f, 2.3083f, 2.6384f, 2.8595f, 2.9872608805f, 3.0127391195f, 3.0127391195f },
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,-12.0f) ), 0.0001f},    // minPoint
  {segmented_spline_c5_fwd( 0.18f ),                10.0f},    // midPoint  
  {segmented_spline_c5_fwd( 0.18f*powf(2.,10.0f) ),  1000.0f},    // maxPoint
  3.0f,  // slopeLow
  0.06f  // slopeHigh
};

const SegmentedSplineParams_c9 ODT_2000nits =
{
  // coefsLow[10]
  { -4.9706219331f, -3.0293780669f, -2.1262f, -1.5105f, -1.0578f, -0.4668f, 0.11938f, 0.7088134201f, 1.2911865799f, 1.2911865799f },
  // coefsHigh[10]
  { 0.8019952042f, 1.1980047958f, 1.5943000000f, 1.9973000000f, 2.3783000000f, 2.7684000000f, 3.0515000000f, 3.2746293562f, 3.3274306351f, 3.3274306351f },
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,-12.0f) ), 0.0001f},    // minPoint
  {segmented_spline_c5_fwd( 0.18f ),                10.0f},    // midPoint  
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,11.0f) ),  2000.0f},    // maxPoint
  3.0f,  // slopeLow
  0.12f  // slopeHigh
};

const SegmentedSplineParams_c9 ODT_4000nits =
{
  // coefsLow[10]
  { -4.9706219331f, -3.0293780669f, -2.1262f, -1.5105f, -1.0578f, -0.4668f, 0.11938f, 0.7088134201f, 1.2911865799f, 1.2911865799f },
  // coefsHigh[10]
  { 0.7973186613f, 1.2026813387f, 1.6093000000f, 2.0108000000f, 2.4148000000f, 2.8179000000f, 3.1725000000f, 3.5344995451f, 3.6696204376f, 3.6696204376f },
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,-12.0f) ), 0.0001f},    // minPoint
  {segmented_spline_c5_fwd( 0.18f ),                10.0f},    // midPoint  
  {segmented_spline_c5_fwd( 0.18f*powf(2.0f,12.0f) ),  4000.0f},    // maxPoint
  3.0f,  // slopeLow
  0.3f   // slopeHigh
};















float segmented_spline_c9_fwd
  ( 
    float x,
    SegmentedSplineParams_c9 C = ODT_48nits
  )
{    
  const int N_KNOTS_LOW = 8;
  const int N_KNOTS_HIGH = 8;

  // Check for negatives or zero before taking the log. If negative or zero,
  // set to HALF_MIN.
  float logx = log10f( max(x, HALF_MIN )); 

  float logy;

  if ( logx <= log10f(C.minPoint.x) ) { 

    logy = logx * C.slopeLow + ( log10f(C.minPoint.y) - C.slopeLow * log10f(C.minPoint.x) );

  } else if (( logx > log10f(C.minPoint.x) ) && ( logx < log10f(C.midPoint.x) )) {

    float knot_coord = (N_KNOTS_LOW-1) * (logx-log10f(C.minPoint.x))/(log10f(C.midPoint.x)-log10f(C.minPoint.x));
    int j = knot_coord;
    float t = knot_coord - j;

    array <float, 3> cf = { C.coefsLow[ j], C.coefsLow[ j + 1], C.coefsLow[ j + 2]};
    // NOTE: If the running a version of CTL < 1.5, you may get an 
    // exception thrown error, usually accompanied by "Array index out of range" 
    // If you receive this error, it is recommended that you update to CTL v1.5, 
    // which contains a number of important bug fixes. Otherwise, you may try 
    // uncommenting the below, which is longer, but equivalent to, the above 
    // line of code.
    //
    // float cf[ 3];
    // if ( j <= 0) {
    //     cf[ 0] = C.coefsLow[0];  cf[ 1] = C.coefsLow[1];  cf[ 2] = C.coefsLow[2];
    // } else if ( j == 1) {
    //     cf[ 0] = C.coefsLow[1];  cf[ 1] = C.coefsLow[2];  cf[ 2] = C.coefsLow[3];
    // } else if ( j == 2) {
    //     cf[ 0] = C.coefsLow[2];  cf[ 1] = C.coefsLow[3];  cf[ 2] = C.coefsLow[4];
    // } else if ( j == 3) {
    //     cf[ 0] = C.coefsLow[3];  cf[ 1] = C.coefsLow[4];  cf[ 2] = C.coefsLow[5];
    // } else if ( j == 4) {
    //     cf[ 0] = C.coefsLow[4];  cf[ 1] = C.coefsLow[5];  cf[ 2] = C.coefsLow[6];
    // } else if ( j == 5) {
    //     cf[ 0] = C.coefsLow[5];  cf[ 1] = C.coefsLow[6];  cf[ 2] = C.coefsLow[7];
    // } else if ( j == 6) {
    //     cf[ 0] = C.coefsLow[6];  cf[ 1] = C.coefsLow[7];  cf[ 2] = C.coefsLow[8];
    // } 
    
    array <float, 3> monomials = { t * t, t, 1. };
    logy = dot_f3_f3( monomials, mult_f3_f33( cf, M));

  } else if (( logx >= log10f(C.midPoint.x) ) && ( logx < log10f(C.maxPoint.x) )) {

    float knot_coord = (N_KNOTS_HIGH-1) * (logx-log10f(C.midPoint.x))/(log10f(C.maxPoint.x)-log10f(C.midPoint.x));
    int j = knot_coord;
    float t = knot_coord - j;

    array <float, 3> cf = { C.coefsHigh[ j], C.coefsHigh[ j + 1], C.coefsHigh[ j + 2]};
    // NOTE: If the running a version of CTL < 1.5, you may get an 
    // exception thrown error, usually accompanied by "Array index out of range" 
    // If you receive this error, it is recommended that you update to CTL v1.5, 
    // which contains a number of important bug fixes. Otherwise, you may try 
    // uncommenting the below, which is longer, but equivalent to, the above 
    // line of code.
    //
    // float cf[ 3];
    // if ( j <= 0) {
    //     cf[ 0] = C.coefsHigh[0];  cf[ 1] = C.coefsHigh[1];  cf[ 2] = C.coefsHigh[2];
    // } else if ( j == 1) {
    //     cf[ 0] = C.coefsHigh[1];  cf[ 1] = C.coefsHigh[2];  cf[ 2] = C.coefsHigh[3];
    // } else if ( j == 2) {
    //     cf[ 0] = C.coefsHigh[2];  cf[ 1] = C.coefsHigh[3];  cf[ 2] = C.coefsHigh[4];
    // } else if ( j == 3) {
    //     cf[ 0] = C.coefsHigh[3];  cf[ 1] = C.coefsHigh[4];  cf[ 2] = C.coefsHigh[5];
    // } else if ( j == 4) {
    //     cf[ 0] = C.coefsHigh[4];  cf[ 1] = C.coefsHigh[5];  cf[ 2] = C.coefsHigh[6];
    // } else if ( j == 5) {
    //     cf[ 0] = C.coefsHigh[5];  cf[ 1] = C.coefsHigh[6];  cf[ 2] = C.coefsHigh[7];
    // } else if ( j == 6) {
    //     cf[ 0] = C.coefsHigh[6];  cf[ 1] = C.coefsHigh[7];  cf[ 2] = C.coefsHigh[8];
    // } 

    array <float, 3> monomials = { t * t, t, 1. };
    logy = dot_f3_f3( monomials, mult_f3_f33( cf, M));

  } else { //if ( logIn >= log10f(C.maxPoint.x) ) { 

    logy = logx * C.slopeHigh + ( log10f(C.maxPoint.y) - C.slopeHigh * log10f(C.maxPoint.x) );

  }

  return pow10f(logy);
  
}


float segmented_spline_c9_rev
  ( 
    float y,
    SegmentedSplineParams_c9 C = ODT_48nits
  )
{  
  const int N_KNOTS_LOW = 8;
  const int N_KNOTS_HIGH = 8;

  const float KNOT_INC_LOW = (log10f(C.midPoint.x) - log10f(C.minPoint.x)) / (N_KNOTS_LOW - 1.);
  const float KNOT_INC_HIGH = (log10f(C.maxPoint.x) - log10f(C.midPoint.x)) / (N_KNOTS_HIGH - 1.);
  
  // KNOT_Y is luminance of the spline at each knot
  float KNOT_Y_LOW[ N_KNOTS_LOW];
  for (int i = 0; i < N_KNOTS_LOW; i = i+1) {
    KNOT_Y_LOW[ i] = ( C.coefsLow[i] + C.coefsLow[i+1]) / 2.;
  };

  float KNOT_Y_HIGH[ N_KNOTS_HIGH];
  for (int i = 0; i < N_KNOTS_HIGH; i = i+1) {
    KNOT_Y_HIGH[ i] = ( C.coefsHigh[i] + C.coefsHigh[i+1]) / 2.;
  };

  float logy = log10f( max( y, 1e-10));

  float logx;
  if (logy <= log10f(C.minPoint.y)) {

    logx = log10f(C.minPoint.x);

  } else if ( (logy > log10f(C.minPoint.y)) && (logy <= log10f(C.midPoint.y)) ) {

    unsigned int j;
    array <float, 3> cf;
    if ( logy > KNOT_Y_LOW[ 0] && logy <= KNOT_Y_LOW[ 1]) {
        cf[ 0] = C.coefsLow[0];  cf[ 1] = C.coefsLow[1];  cf[ 2] = C.coefsLow[2];  j = 0;
    } else if ( logy > KNOT_Y_LOW[ 1] && logy <= KNOT_Y_LOW[ 2]) {
        cf[ 0] = C.coefsLow[1];  cf[ 1] = C.coefsLow[2];  cf[ 2] = C.coefsLow[3];  j = 1;
    } else if ( logy > KNOT_Y_LOW[ 2] && logy <= KNOT_Y_LOW[ 3]) {
        cf[ 0] = C.coefsLow[2];  cf[ 1] = C.coefsLow[3];  cf[ 2] = C.coefsLow[4];  j = 2;
    } else if ( logy > KNOT_Y_LOW[ 3] && logy <= KNOT_Y_LOW[ 4]) {
        cf[ 0] = C.coefsLow[3];  cf[ 1] = C.coefsLow[4];  cf[ 2] = C.coefsLow[5];  j = 3;
    } else if ( logy > KNOT_Y_LOW[ 4] && logy <= KNOT_Y_LOW[ 5]) {
        cf[ 0] = C.coefsLow[4];  cf[ 1] = C.coefsLow[5];  cf[ 2] = C.coefsLow[6];  j = 4;
    } else if ( logy > KNOT_Y_LOW[ 5] && logy <= KNOT_Y_LOW[ 6]) {
        cf[ 0] = C.coefsLow[5];  cf[ 1] = C.coefsLow[6];  cf[ 2] = C.coefsLow[7];  j = 5;
    } else if ( logy > KNOT_Y_LOW[ 6] && logy <= KNOT_Y_LOW[ 7]) {
        cf[ 0] = C.coefsLow[6];  cf[ 1] = C.coefsLow[7];  cf[ 2] = C.coefsLow[8];  j = 6;
    }
    
    const array <float, 3> tmp = mult_f3_f33( cf, M);

    float a = tmp[ 0];
    float b = tmp[ 1];
    float c = tmp[ 2];
    c = c - logy;

    const float d = sqrtf( b * b - 4.0f * a * c);

    const float t = ( 2.0f * c) / ( -d - b);

    logx = log10f(C.minPoint.x) + ( t + j) * KNOT_INC_LOW;

  } else if ( (logy > log10f(C.midPoint.y)) && (logy < log10f(C.maxPoint.y)) ) {

    unsigned int j;
    array <float, 3> cf;
    if ( logy > KNOT_Y_HIGH[ 0] && logy <= KNOT_Y_HIGH[ 1]) {
        cf[ 0] = C.coefsHigh[0];  cf[ 1] = C.coefsHigh[1];  cf[ 2] = C.coefsHigh[2];  j = 0;
    } else if ( logy > KNOT_Y_HIGH[ 1] && logy <= KNOT_Y_HIGH[ 2]) {
        cf[ 0] = C.coefsHigh[1];  cf[ 1] = C.coefsHigh[2];  cf[ 2] = C.coefsHigh[3];  j = 1;
    } else if ( logy > KNOT_Y_HIGH[ 2] && logy <= KNOT_Y_HIGH[ 3]) {
        cf[ 0] = C.coefsHigh[2];  cf[ 1] = C.coefsHigh[3];  cf[ 2] = C.coefsHigh[4];  j = 2;
    } else if ( logy > KNOT_Y_HIGH[ 3] && logy <= KNOT_Y_HIGH[ 4]) {
        cf[ 0] = C.coefsHigh[3];  cf[ 1] = C.coefsHigh[4];  cf[ 2] = C.coefsHigh[5];  j = 3;
    } else if ( logy > KNOT_Y_HIGH[ 4] && logy <= KNOT_Y_HIGH[ 5]) {
        cf[ 0] = C.coefsHigh[4];  cf[ 1] = C.coefsHigh[5];  cf[ 2] = C.coefsHigh[6];  j = 4;
    } else if ( logy > KNOT_Y_HIGH[ 5] && logy <= KNOT_Y_HIGH[ 6]) {
        cf[ 0] = C.coefsHigh[5];  cf[ 1] = C.coefsHigh[6];  cf[ 2] = C.coefsHigh[7];  j = 5;
    } else if ( logy > KNOT_Y_HIGH[ 6] && logy <= KNOT_Y_HIGH[ 7]) {
        cf[ 0] = C.coefsHigh[6];  cf[ 1] = C.coefsHigh[7];  cf[ 2] = C.coefsHigh[8];  j = 6;
    }
    
    const array <float, 3> tmp = mult_f3_f33( cf, M);

    float a = tmp[ 0];
    float b = tmp[ 1];
    float c = tmp[ 2];
    c = c - logy;

    const float d = sqrtf( b * b - 4.0f * a * c);

    const float t = ( 2.0f * c) / ( -d - b);

    logx = log10f(C.midPoint.x) + ( t + j) * KNOT_INC_HIGH;

  } else { //if ( logy >= log10f(C.maxPoint.y) ) {

    logx = log10f(C.maxPoint.x);

  }
  
  return pow10f( logx);
}
