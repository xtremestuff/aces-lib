
// <ACEStransformID>urn:ampas:aces:transformId:v1.5:ACESlib.Transform_Common.a1.0.3</ACEStransformID>
// <ACESuserName>ACES 1.0 Lib - Transform Common</ACESuserName>

//
// Contains functions and constants shared by multiple forward and inverse 
// transforms 
//

#pragma once
#include <cmath>
#include <array>
using namespace std;

#include "ACESlib.CustomDefs.h"
#include "ACESlib.Utilities_Color.h"
#include "ACESlib.Utilities.h"



const array <array <double, 4>, 4> AP0_2_XYZ_MAT = { RGBtoXYZ(AP0, 1.0) };
const array <array <double, 4>, 4> XYZ_2_AP0_MAT = { XYZtoRGB(AP0, 1.0) };

const array <array <double, 4>, 4> AP1_2_XYZ_MAT = { RGBtoXYZ(AP1, 1.0) };
const array <array <double, 4>, 4> XYZ_2_AP1_MAT = { XYZtoRGB(AP1, 1.0) };

const array <array <double, 4>, 4> AP0_2_AP1_MAT = mult_f44_f44( AP0_2_XYZ_MAT, XYZ_2_AP1_MAT);
const array <array <double, 4>, 4> AP1_2_AP0_MAT = mult_f44_f44( AP1_2_XYZ_MAT, XYZ_2_AP0_MAT);

const array <double, 3> AP1_RGB2Y = { AP1_2_XYZ_MAT[0][1],
                             AP1_2_XYZ_MAT[1][1], 
                             AP1_2_XYZ_MAT[2][1] };



const double TINY = (double)1e-10;



double rgb_2_saturation(array <double, 3> rgb)
{
  return ( max( max_f3(rgb), TINY) - max( min_f3(rgb), TINY)) / max( max_f3(rgb), (double)1e-2);
}