
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.Utilities.a1.0.3</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - Utilities</ACESuserName>

//
// Generic functions that may be useful for writing C++ programs
//

#pragma once
#include <cmath>
#include <array>
#include <cstdio>
using namespace std;

#include "ACESlib.CustomDefs.h"



double min( double a, double b)
{
  if (a < b)
    return a;
  else
    return b;
}

double max( double a, double b)
{
  if (a > b)
    return a;
  else
    return b;
}

double min_f3( const array <double, 3> &a)
{
  return min( a[0], min( a[1], a[2]));
}

double max_f3( const array <double, 3> &a)
{
  return max( a[0], max( a[1], a[2]));
}

double clip( double v)
{
  return min(v, 1.0);
}

array <double, 3> clip_f3( const array <double, 3> &in)
{
  array <double, 3> out;
  out[0] = clip( in[0]);
  out[1] = clip( in[1]);
  out[2] = clip( in[2]);

  return out;
}

double clamp( double in, double clampMin, double clampMax)
{
  // Note: Numeric constants can be used in place of a min or max value (i.e. 
  // use HALF_NEG_INF in place of clampMin or HALF_POS_INF in place of clampMax)
  
  return max( clampMin, min(in, clampMax));
}

array <double, 3> clamp_f3(const array <double, 3> &in, double clampMin, double clampMax)
{
  // Note: Numeric constants can be used in place of a min or max value (i.e. 
  // use HALF_NEG_INF in place of clampMin or HALF_POS_INF in place of clampMax)

  array <double, 3> out;
  out[0] = clamp( in[0], clampMin, clampMax);
  out[1] = clamp( in[1], clampMin, clampMax);
  out[2] = clamp( in[2], clampMin, clampMax);
      
  return out;
}

array <double, 3> add_f_f3( double a, const array <double, 3> &b)
{
  array <double, 3> out;
  out[0] = a + b[0];
  out[1] = a + b[1];
  out[2] = a + b[2];
  return out;
}

array <double, 3> powf_f3( const array <double, 3> &a, double b)
{
  array <double, 3> out;
  out[0] = powf(a[0], b);
  out[1] = powf(a[1], b);
  out[2] = powf(a[2], b);
  return out;
}

array <double, 3> pow10f_f3( const array <double, 3> &a)
{
  array <double, 3> out;
  out[0] = pow10f(a[0]);
  out[1] = pow10f(a[1]);
  out[2] = pow10f(a[2]);
  return out;
}

array <double, 3> log10f_f3( const array <double, 3> &a)
{
  array <double, 3> out;
  out[0] = log10f(a[0]);
  out[1] = log10f(a[1]);
  out[2] = log10f(a[2]);
  return out;
}

/* //Use existing roundf instead
double round(double x)
{
  int x1;
 
  if (x < 0.0)
    x1 = (int)(x - 0.5);
  else
    x1 = (int)(x + 0.5);
 
  return (double)x1;
}
*/

int sign( double x)
{
    // Signum function:
    //  sign(X) returns 1 if the element is greater than zero, 0 if it equals zero 
    //  and -1 if it is less than zero

    int y;
    if (x < 0.0) { 
        y = -1;
    } else if (x > 0.0) {
        y = 1;
    } else {
        y = 0;
    }

    return y;
}

void print_f3( const array <double, 3> &m)
{
   printf("%f,\t%f,\t%f\n", m[0], m[1], m[2]);
}

void print_f33( const array <array <double, 3>, 3> &m)
{
   printf( "{ {%f,\t%f,\t%f},\n", m[0][0], m[0][1], m[0][2]);
   printf( "  {%f,\t%f,\t%f},\n", m[1][0], m[1][1], m[1][2]);
   printf( "  {%f,\t%f,\t%f} };\n", m[2][0], m[2][1], m[2][2]);
}

void print_f44( const array <array <double, 4>, 4> &m)
{
    printf("{ {%f,\t%f,\t%f,\t%f},\n", m[0][0], m[0][1], m[0][2], m[0][3]);
    printf("  {%f,\t%f,\t%f,\t%f},\n", m[1][0], m[1][1], m[1][2], m[1][3]);
    printf("  {%f,\t%f,\t%f,\t%f},\n", m[2][0], m[2][1], m[2][2], m[2][3]);
    printf("  {%f,\t%f,\t%f,\t%f} };\n", m[3][0], m[3][1], m[3][2], m[3][3]);
}
