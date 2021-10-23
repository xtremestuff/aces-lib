
// Custom Definitions for ACES Libraries
//

#pragma once
#include <cmath>
#include <array>
using namespace std;



#ifndef HALF_MIN
	#define HALF_MIN	5.96046448e-08f	// Smallest positive half
#endif
#ifndef HALF_MAX
	#define HALF_MAX	65504.0f	// Largest positive half
#endif
#ifndef HALF_POS_INF
	#define HALF_POS_INF	0x7c00 // Half Positive Infinity
#endif
#ifndef HALF_NEG_INF
	#define HALF_NEG_INF	0xfc00 // Half Negative Infinity
#endif
#ifndef FLT_NAN
	#define FLT_NAN 0x7fffffff
#endif

struct Chromaticities
{
	array <float, 2>	red;		// CIE xy coordinates of red primary
	array <float, 2>	green;		// CIE xy coordinates of green primary
	array <float, 2>	blue;		// CIE xy coordinates of blue primary
	array <float, 2>	white;		// CIE xy coordinates of white point
};


float pow10f(float x) {
	return powf(10.0f, x);
}

float dot_f3_f3(const array <float, 3> &x, const array <float, 3> &y) {
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

array <float, 3> mult_f_f3(float f, const array <float, 3> &x) {

	array <float, 3> r;

	r[0] = f * x[0];
	r[1] = f * x[1];
	r[2] = f * x[2];

	return r;
}

array <array <float, 3>, 3> transpose_f33(array <array <float, 3>, 3> &a) {

	array <array <float, 3>, 3> r;

	for (int i = 0; i < 3; ++i) {

		for (int j = 0; j < 3; ++j) {
			r[i][j] = a[j][i];
		}
	}

	return r;
}

array <array <float, 4>, 4> transpose_f44(const array <array <float, 4>, 4> &a) {

	array <array <float, 4>, 4> r;

	for (int i = 0; i < 4; ++i) {

		for (int j = 0; j < 4; ++j) {
			r[i][j] = a[j][i];
		}
	}

	return r;
}

array <array <float, 3>, 3> mult_f_f33(float f, const array <array <float, 3>, 3> &a) {

	array <array <float, 3>, 3> r;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			r[i][j] = f * a[i][j];
		}
	}

	return r;
}

array <array <float, 4>, 4> mult_f_f44(float f, const array <array <float, 4>, 4> &a) {

	array <array <float, 4>, 4> r;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			r[i][j] = f * a[i][j];
		}
	}

	return r;
}


array <float, 3> mult_f3_f33(const array <float, 3>& x, const array <array <float, 3>, 3>& a) {
	array <float, 3> r;

	for (int i = 0; i < 3; ++i) {

		r[i] = 0.0f;

		for (int j = 0; j < 3; ++j) {
			r[i] = r[i] + x[j] * a[j][i];
		}
	}
	return r;
}


array <float, 3> mult_f3_f44(const array <float, 3> &x, const array <array <float, 4>, 4> &a) {

	array <float, 3> r;

	for (int i = 0; i < 3; ++i) {

		r[i] = 0.0f;

		for (int j = 0; j < 3; ++j) {
			r[i] = r[i] + x[j] * a[j][i];
		}

		r[i] = r[i] + a[3][i];
	}

	float s = 1.0f / (x[0] * a[0][3] + x[1] * a[1][3] + x[2] * a[2][3] + a[3][3]);

	for (int k = 0; k < 3; ++k) {
		r[k] = r[k] * s;
	}

	return r;
}


array <array <float, 3>, 3> mult_f33_f33(const array <array <float, 3>, 3> &a, const array <array <float, 3>, 3> &b) {

	array <array <float, 3>, 3> r;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {

			r[i][j] = 0.0f;

			for (int k = 0; k < 3; ++k) {
				r[i][j] = r[i][j] + a[i][k] * b[k][j];
			}
		}
	}

	return r;
}


array <array <float, 4>, 4> mult_f44_f44(const array <array <float, 4>, 4> &a, const array <array <float, 4>, 4> &b) {

	array <array <float, 4>, 4> r;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {

			r[i][j] = 0.0f;

			for (int k = 0; k < 4; ++k) {
				r[i][j] = r[i][j] + a[i][k] * b[k][j];
			}
		}
	}

	return r;
}


array <array <float, 3>, 3> add_f33_f33(const array <array <float, 3>, 3> &a, const array <array <float, 3>, 3> &b) {

	array <array <float, 3>, 3> r;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			r[i][j] = a[i][j] + b[i][j];
		}
	}

	return r;
}


array <array <float, 3>, 3> invert_f33(const array <array <float, 3>, 3>& a) {

	float determinant = 0.0f;
	array <array <float, 3>, 3> res = { { {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}} };

	for (int i = 0; i < 3; i++)
		determinant += (a[0][i] * (a[1][(i + 1) % 3] * a[2][(i + 2) % 3] - a[1][(i + 2) % 3] * a[2][(i + 1) % 3]));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			res[i][j] = ((a[(j + 1) % 3][(i + 1) % 3] * a[(j + 2) % 3][(i + 2) % 3]) - (a[(j + 1) % 3][(i + 2) % 3] * a[(j + 2) % 3][(i + 1) % 3])) / determinant;
	}

	return res;

}

array <array <float, 4>, 4> invert_f44(const array <array <float, 4>, 4>& a) {

	float determinant = 0.0f;
	array <array <float, 4>, 4> res = { { {1.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 1.0f} } };
	int len = 3;

	if (a[0][3] != 0 || a[1][3] != 0 || a[2][3] != 0 || a[3][3] != 1)
		len = 4;

	for (int i = 0; i < len; i++)
		determinant += (a[0][i] * (a[1][(i + 1) % len] * a[2][(i + 2) % len] - a[1][(i + 2) % len] * a[2][(i + 1) % len]));

	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++)
			res[i][j] = ((a[(j + 1) % len][(i + 1) % len] * a[(j + 2) % len][(i + 2) % len]) - (a[(j + 1) % len][(i + 2) % len] * a[(j + 2) % len][(i + 1) % len])) / determinant;
	}

	return res;

}


array <array <float, 4>, 4> RGBtoXYZ(const Chromaticities &chroma, float Y)
{
	//
	// X and Z values of RGB value (1, 1, 1), or "white"
	//

	float X = chroma.white[0] * Y / chroma.white[1];
	float Z = (1.0f - chroma.white[0] - chroma.white[1]) * Y / chroma.white[1];

	//
	// Scale factors for matrix rows
	//

	float d = chroma.red[0] * (chroma.blue[1] - chroma.green[1]) +
		chroma.blue[0] * (chroma.green[1] - chroma.red[1]) +
		chroma.green[0] * (chroma.red[1] - chroma.blue[1]);

	float Sr = (X * (chroma.blue[1] - chroma.green[1]) -
		chroma.green[0] * (Y * (chroma.blue[1] - 1.0f) +
			chroma.blue[1] * (X + Z)) +
		chroma.blue[0] * (Y * (chroma.green[1] - 1.0f) +
			chroma.green[1] * (X + Z))) / d;

	float Sg = (X * (chroma.red[1] - chroma.blue[1]) +
		chroma.red[0] * (Y * (chroma.blue[1] - 1.0f) +
			chroma.blue[1] * (X + Z)) -
		chroma.blue[0] * (Y * (chroma.red[1] - 1.0f) +
			chroma.red[1] * (X + Z))) / d;

	float Sb = (X * (chroma.green[1] - chroma.red[1]) -
		chroma.red[0] * (Y * (chroma.green[1] - 1.0f) +
			chroma.green[1] * (X + Z)) +
		chroma.green[0] * (Y * (chroma.red[1] - 1.0f) +
			chroma.red[1] * (X + Z))) / d;

	//
	// Assemble the matrix
	//

	array <array <float, 4>, 4> M;

	M[0][0] = Sr * chroma.red[0];
	M[0][1] = Sr * chroma.red[1];
	M[0][2] = Sr * (1.0f - chroma.red[0] - chroma.red[1]);
	M[0][3] = 0;

	M[1][0] = Sg * chroma.green[0];
	M[1][1] = Sg * chroma.green[1];
	M[1][2] = Sg * (1.0f - chroma.green[0] - chroma.green[1]);
	M[1][3] = 0;

	M[2][0] = Sb * chroma.blue[0];
	M[2][1] = Sb * chroma.blue[1];
	M[2][2] = Sb * (1.0f - chroma.blue[0] - chroma.blue[1]);
	M[2][3] = 0;

	M[3][0] = 0;
	M[3][1] = 0;
	M[3][2] = 0;
	M[3][3] = 1;

	return M;
}

array <array <float, 4>, 4> XYZtoRGB(const Chromaticities &c, float Y) {

	array <array <float, 4>, 4> M = RGBtoXYZ(c, Y);
	M = invert_f44(M);
	return M;
}


template <size_t int1D_N>
float interpolate1D(const array <array <float, int1D_N>, 2> &table, float p)
{
	int size = table.size();

	if (size < 1)
		return 0;

	if (p < table[0][0])
		return table[0][1];

	if (p >= table[size - 1][0])
		return table[size - 1][1];

	int i = 0;
	int j = size;

	while (i < j - 1)
	{
		int k = (i + j) / 2;

		if (table[k][0] == p)
			return table[k][1];
		else if (table[k][0] < p)
			i = k;
		else
			j = k;
	}

	float t = (p - table[i][0]) / (table[i + 1][0] - table[i][0]);
	float s = 1 - t;

	return s * table[i][1] + t * table[i + 1][1];
}
