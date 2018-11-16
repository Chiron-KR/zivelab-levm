// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <complex> 
#include <iostream>

#include "lib\levmar.h"

using namespace std;

#define twopi   6.2831853071795864769252867665590057683943387987502 // 2 * pi

#pragma region Define Rosenbrock

#define ROSD 105.0
double rosenbrock_p[2] = { 1.0, 1.0 }; // best parameters
void rosenbrock(double *p, double *x, int m, int n, void *data)
{
	register int i;
	for (i = 0; i < n; ++i)
		x[i] = ((1.0 - p[0])*(1.0 - p[0]) + ROSD * (p[1] - p[0] * p[0])*(p[1] - p[0] * p[0]));
}
void rosenbrockprime(double *p, double *jac, int m, int n, void *data)
{
	register int i, j;

	for (i = j = 0; i < n; ++i) {
		jac[j++] = (-2 + 2 * p[0] - 4 * ROSD*(p[1] - p[0] * p[0])*p[0]);
		jac[j++] = (2 * ROSD*(p[1] - p[0] * p[0]));
	}
}

#pragma endregion 

#pragma region Define Norris

// https://www.itl.nist.gov/div898/strd/lls/data/Norris.shtml

double norris_p[2] = { -0.262323073774029, 1.00211681802045 }; // best parameters
double norris_y[36] = {
	0.1,	338.8,	118.1,	888,	9.2,	228.1,	668.5,	998.5,	449.1,	778.9,	
	559.2,	0.3,	0.1,	778.1,	668.8,	339.3,	448.9,	10.8,	557.7,	228.3,	
	998,	888.8,	119.6,	0.3,	0.6,	557.6,	339.3,	888,	998.5,	778.9,	
	10.2,	117.6,	228.9,	668.4,	449.2,	0.2
};
double norris_x[36] = {
	0.2,	337.4,	118.2,	884.6,	10.1,	226.5,	666.3,	996.3,	448.6,	777,
	558.2,	0.4,	0.6,	775.5,	666.9,	338,	447.5,	11.6,	556,	228.1,
	995.8,	887.6,	120.2,	0.3,	0.3,	556.8,	339.1,	887.2,	999,	779,
	11.1,	118.3,	229.2,	669.1,	448.9,	0.5
};
void norris(double *p, double *y, int m, int n, void *data)
{
	register int i;
	for (i = 0; i < n; ++i)
		y[i] = p[0] + p[1] * norris_x[i];
}
void norrisprime(double *p, double *jac, int m, int n, void *data)
{
	register int i, j;

	for (i = j = 0; i < n; ++i) {
		jac[j++] = 1.0;
		jac[j++] = norris_x[i];
	}
}

#pragma endregion

#pragma region Define Lanczos1

// https://www.itl.nist.gov/div898/strd/nls/data/lanczos1.shtml
//

double lanczos1_p[6] = { 0.0951, 1.0, 0.8607, 3, 1.5576, 5 }; // best parameters
double lanczos1_y[24] = {
	2.513400000000,	2.044333373291,	1.668404436564,	1.366418021208,	1.123232487372,
	0.926889718004,	0.767933856373,	0.638877552311,	0.533783531740,	0.447936361735,
	0.377584788435,	0.319739319933,	0.272013077375,	0.232496552903,	0.199658954607,
	0.172270412691,	0.149340566017,	0.130070020692,	0.113811932464,	0.100041558756,
	0.088332090845,	0.078335440194,	0.069766937434,	0.062393125367
};
double lanczos1_x[24] = {
	0.000000000000,	0.050000000000,	0.100000000000,	0.150000000000,	0.200000000000,
	0.250000000000,	0.300000000000,	0.350000000000,	0.400000000000,	0.450000000000,
	0.500000000000,	0.550000000000,	0.600000000000,	0.650000000000,	0.700000000000,
	0.750000000000,	0.800000000000,	0.850000000000,	0.900000000000,	0.950000000000,
	1.000000000000,	1.050000000000,	1.100000000000,	1.150000000000
};
void lanczos1(double *p, double *y, int m, int n, void *data)
{
	register int i;

	double *xx = (double*)data;
	for (i = 0; i < n; ++i, xx++)
	{
		double x = *xx;
		y[i] = p[0] * exp(-p[1] * x) + p[2] * exp(-p[3] * x) + p[4] * exp(-p[5] * x);
	}
}

#pragma endregion

#pragma region Define Thurber

// https://www.itl.nist.gov/div898/strd/nls/data/thurber.shtml
//

double thurber_p[7] = {
    1288.139680000000, 1491.079253500000, 583.238368770000, 75.416644291000, 0.966295028640,
    0.397972857970, 0.049727297349 }; // best parameters
double thurber_y[37] = {
	 80.574,	084.248,	087.264,	087.195,	089.076,	
	 089.608,	089.868,	090.101,	092.405,	095.854,
	 100.696,	101.060,	401.672,	390.724,	567.534,	
	 635.316,	733.054,	759.087,	894.206,	990.785,
	1090.109,	1080.914,	1122.643,	1178.351,	1260.531,
	1273.514,	1288.339,	1327.543,	1353.863,	1414.509,
	1425.208,	1421.384,	1442.962,	1464.350,	1468.705,
	1447.894,	1457.628
};
double thurber_x[37] = {
	-3.067,	-2.981,	-2.921,	-2.912,	-2.84,
	-2.797,	-2.702,	-2.699,	-2.633,	-2.481,
	-2.363,	-2.322,	-1.501,	-1.460,	-1.274,
	-1.212,	-1.100,	-1.046,	-0.915,	-0.714,
	-0.566,	-0.545,	-0.400,	-0.309,	-0.109,
	-0.103,	0.01,	0.119,	0.377,	0.79,
	0.963,	1.006,	1.115,	1.572,	1.841,
	2.047,	2.2
};
void thurber(double *p, double *y, int m, int n, void *data)
{
	register int i;

	double *xx = (double*)data;
	for (i = 0; i < n; ++i, xx++)
	{
		double x = *xx;
		y[i] = (p[0] + p[1] * x + p[2] * x*x + p[3] * x*x*x) / (1 + p[4] * x + p[5] * x*x + p[6] * x*x*x);
	}
}

#pragma endregion

#pragma region Define Rat43

// https://www.itl.nist.gov/div898/strd/nls/data/ratkowsky3.shtml

double rat43_p[4] = {
    699.641512700000,   5.277125302500, 0.759629383290, 1.279248385900
}; // best parameters
double rat43_y[15] = {
     16.08, 33.83,  65.80,  97.20,  191.55, 326.20, 386.87, 520.53, 590.03, 651.92,
    724.93, 699.56, 689.96, 637.56, 717.41
};
double rat43_x[15] = {
    1.00,   2.00,   3.00,   4.00,   5.00,   6.00,   7.00,   8.00,   9.00,   10.00,
    11.00,  12.00,  13.00,  14.00,  15.00
};
void rat43(double *p, double *y, int m, int n, void *data)
{
    register int i;

    double *xx = (double*)data;
    for (i = 0; i < n; ++i, xx++)
    {
        double x = *xx;
        y[i] = p[0] / pow(1.0 + exp(p[1] - p[2] * x), 1.0 / p[3]);
    }
}

#pragma endregion

#pragma region Define Randle - Rs-Rct|Cdl

// Model Expression = Rs-Rct|Cdl
// Z = Rs + 1/(1/Rct + Cdl*s)
// Parameters = { Rs, rct, Cdl }
// dZ/dRs = 1
// dZ/dRct = 1 / (Rct ^ 2 * (1 / Rct + Cdl * s) ^ 2)
// dZ/dCdl = -s / (1 / Rct + Cdl * s) ^ 2

double randle_p[3] = {
    100,    100,    1e-6
}; // best parameters
double randle_z1[181] = {
    100.000300,	100.000300,	100.000400,	100.000500,	100.000600,
    100.000800,	100.001000,	100.001300,	100.001600,	100.002000,
    100.002500,	100.003200,	100.004000,	100.005100,	100.006400,
    100.008000,	100.010100,	100.012700,	100.016000,	100.020100,
    100.025300,	100.031900,	100.040100,	100.050500,	100.063600,
    100.080000,	100.100700,	100.126800,	100.159600,	100.200800,
    100.252700,	100.317900,	100.399900,	100.502900,	100.632200,
    100.794600,	100.998300,	101.253600,	101.573100,	101.972400,
    102.470500,	103.090300,	103.859600,	104.810900,	105.982100,
    107.416100,	109.160400,	111.265100,	113.780000,	116.750300,
    120.210800,	124.178600,	128.645700,	133.572700,	138.885300,
    144.475700,	150.209500,	155.937900,	161.512300,	166.800100,
    171.695700,	176.127400,	180.058200,	183.482200,	186.418000,
    188.901400,	190.978100,	192.698200,	194.111500,	195.265300,
    196.202100,	196.959500,	197.569600,	198.059800,	198.452700,
    198.767000,	199.018100,	199.218500,	199.378200,	199.505500,
    199.606800,	199.687400,	199.751500,	199.802500,	199.843100,
    199.875300,	199.900900,	199.921300,	199.937500,	199.950300,
    199.960500,	199.968700,	199.975100,	199.980200,	199.984300,
    199.987500,	199.990100,	199.992100,	199.993700,	199.995000,
    199.996100,	199.996900,	199.997500,	199.998000,	199.998400,
    199.998800,	199.999000,	199.999200,	199.999400,	199.999500,
    199.999600,	199.999700,	199.999800,	199.999800,	199.999800,
    199.999900,	199.999900,	199.999900,	199.999900,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000,	200.000000,	200.000000,	200.000000,	200.000000,
    200.000000
}; // real part of the observed data
double randle_z2[181] = {
    -0.159155,	-0.178574,	-0.200363,	-0.224811,	-0.252242,
    -0.283020,	-0.317553,	-0.356299,	-0.399773,	-0.448551,
    -0.503279,	-0.564685,	-0.633582,	-0.710883,	-0.797613,
    -0.894922,	-1.004099,	-1.126588,	-1.264011,	-1.418185,
    -1.591146,	-1.785179,	-2.002838,	-2.246988,	-2.520832,
    -2.827954,	-3.172360,	-3.558518,	-3.991412,	-4.476589,
    -5.020205,	-5.629080,	-6.310737,	-7.073441,	-7.926211,
    -8.878820,	-9.941744,	-11.126060,	-12.443250,	-13.904920,
    -15.522310,	-17.305620,	-19.263090,	-21.399680,	-23.715420,
    -26.203280,	-28.846640,	-31.616560,	-34.468970,	-37.342450,
    -40.157260,	-42.816550,	-45.210570,	-47.224410,	-48.748990,
    -49.693890,	-49.999560,	-49.646170,	-48.656620,	-47.093060,
    -45.047720,	-42.630520,	-39.956320,	-37.134150,	-34.259740,
    -31.411480,	-28.649490,	-26.016630,	-23.540880,	-21.238040,
    -19.114560,	-17.169990,	-15.399070,	-13.793380,	-12.342620,
    -11.035490,	-9.860397,	-8.805871,	-7.860877,	-7.014985,
    -6.258478,	-5.582389,	-4.978510,	-4.439371,	-3.958201,
    -3.528889,	-3.145933,	-2.804387,	-2.499817,	-2.228251,
    -1.986134,	-1.770287,	-1.577872,	-1.406352,	-1.253463,
    -1.117186,	-0.995719,	-0.887454,	-0.790957,	-0.704950,
    -0.628294,	-0.559972,	-0.499079,	-0.444807,	-0.396436,
    -0.353325,	-0.314902,	-0.280657,	-0.250136,	-0.222935,
    -0.198691,	-0.177084,	-0.157826,	-0.140663,	-0.125366,
    -0.111732,	-0.099582,	-0.088752,	-0.079101,	-0.070498,
    -0.062832,	-0.055999,	-0.049909,	-0.044482,	-0.039644,
    -0.035333,	-0.031491,	-0.028066,	-0.025014,	-0.022294,
    -0.019869,	-0.017708,	-0.015783,	-0.014066,	-0.012537,
    -0.011173,	-0.009958,	-0.008875,	-0.007910,	-0.007050,
    -0.006283,	-0.005600,	-0.004991,	-0.004448,	-0.003964,
    -0.003533,	-0.003149,	-0.002807,	-0.002501,	-0.002229,
    -0.001987,	-0.001771,	-0.001578,	-0.001407,	-0.001254,
    -0.001117,	-0.000996,	-0.000888,	-0.000791,	-0.000705,
    -0.000628,	-0.000560,	-0.000499,	-0.000445,	-0.000396,
    -0.000353,	-0.000315,	-0.000281,	-0.000250,	-0.000223,
    -0.000199,	-0.000177,	-0.000158,	-0.000141,	-0.000125,
    -0.000112,	-9.958178E-5,	-8.875235E-5,	-7.910062E-5,	-7.049850E-5,
    -6.283185E-5
}; // imaginary part of the observed data
double randle_f[181] = {
    1000000.000000,	891250.900000,	794328.200000,	707945.800000,	630957.300000,
    562341.300000,	501187.200000,	446683.600000,	398107.200000,	354813.400000,
    316227.800000,	281838.300000,	251188.600000,	223872.100000,	199526.200000,
    177827.900000,	158489.300000,	141253.800000,	125892.500000,	112201.800000,
    100000.000000,	89125.090000,	79432.820000,	70794.580000,	63095.730000,
    56234.130000,	50118.720000,	44668.360000,	39810.720000,	35481.340000,
    31622.780000,	28183.830000,	25118.860000,	22387.210000,	19952.620000,
    17782.790000,	15848.930000,	14125.380000,	12589.250000,	11220.180000,
    10000.000000,	8912.509000,	7943.282000,	7079.458000,	6309.573000,
    5623.413000,	5011.872000,	4466.836000,	3981.072000,	3548.134000,
    3162.278000,	2818.383000,	2511.886000,	2238.721000,	1995.262000,
    1778.279000,	1584.893000,	1412.538000,	1258.925000,	1122.018000,
    1000.000000,	891.250900,	794.328200,	707.945800,	630.957300,
    562.341300,	501.187200,	446.683600,	398.107200,	354.813400,
    316.227800,	281.838300,	251.188600,	223.872100,	199.526200,
    177.827900,	158.489300,	141.253800,	125.892500,	112.201800,
    100.000000,	89.125090,	79.432820,	70.794580,	63.095730,
    56.234130,	50.118720,	44.668360,	39.810720,	35.481340,
    31.622780,	28.183830,	25.118860,	22.387210,	19.952620,
    17.782790,	15.848930,	14.125380,	12.589250,	11.220180,
    10.000000,	8.912509,	7.943282,	7.079458,	6.309573,
    5.623413,	5.011872,	4.466836,	3.981072,	3.548134,
    3.162278,	2.818383,	2.511886,	2.238721,	1.995262,
    1.778279,	1.584893,	1.412538,	1.258925,	1.122018,
    1.000000,	0.891251,	0.794328,	0.707946,	0.630957,
    0.562341,	0.501187,	0.446684,	0.398107,	0.354813,
    0.316228,	0.281838,	0.251189,	0.223872,	0.199526,
    0.177828,	0.158489,	0.141254,	0.125892,	0.112202,
    0.100000,	0.089125,	0.079433,	0.070795,	0.063096,
    0.056234,	0.050119,	0.044668,	0.039811,	0.035481,
    0.031623,	0.028184,	0.025119,	0.022387,	0.019953,
    0.017783,	0.015849,	0.014125,	0.012589,	0.011220,
    0.010000,	0.008913,	0.007943,	0.007079,	0.006310,
    0.005623,	0.005012,	0.004467,	0.003981,	0.003548,
    0.003162,	0.002818,	0.002512,	0.002239,	0.001995,
    0.001778,	0.001585,	0.001413,	0.001259,	0.001122,
    0.001000
}; // frequency data
void randle(double *p, double *y, int m, int n, void *data)
{
    // PREMISS: the 1st half of y is real part and the other is imaginary part of the observed data
    //          y = { z1[0], z1[1], ..., z1[n/2 - 1], z2[0], z2[1], ..., z2[n/2 - 1] }
    
    register int i;
    int Ns = (int)(n / 2);

    double Rs = p[0];
    double Rct = p[1];
    double Cdl = p[2];

    double *ff = (double*)data; 
    for (i = 0; i < Ns; ++i, ff++)
    {
        double f = *ff; 
        complex <double> s(0.0, twopi * f); // s = j * w, where angular frequency w = 2 * pi * f
        complex <double> z = Rs + Rct / (1.0 + Rct * Cdl * s); //is better than Rs + 1.0 / (1.0 / Rct + Cdl * s);

        y[i] = real(z);
        y[Ns + i] = imag(z);
    }
}
void randleprime(double *p, double *jac, int m, int n, void *data)
{
    register int i;
    int Ns = (int)(n / 2);

    double Rs = p[0];
    double Rct = p[1];
    double Cdl = p[2];    

    double *ff = (double*)data;
    for (i = 0; i < Ns; ++i, ff++)
    {
        double f = *ff;
        complex <double> s (0.0, twopi * f); // s = j * w, where angular frequency w = 2 * pi * f
        complex <double> jac0 = 1.0;
        complex <double> jac1 = 1.0 / ((1.0 + Rct * Cdl * s) * (1.0 + Rct * Cdl * s));
        complex <double> jac2 = - Rct * Rct * s / ((1.0 + Rct * Cdl * s) * (1.0 + Rct * Cdl * s));

        jac[m*i] = real(jac0);
        jac[m*i + 1] = real(jac1);
        jac[m*i + 2] = real(jac2);
        jac[m*(Ns + i)] = imag(jac0);
        jac[m*(Ns + i) + 1] = imag(jac1);
        jac[m*(Ns + i) + 2] = imag(jac2);
    }
}

#pragma endregion

int main()
{
	int i, ret;

	int m; // parameter vector dimension
	int n; // measurement vector dimension

	double p[40]; // initial parameter estimates
    double lb[40]; // lower bound of parameters
    double x[40]; // measurement vector
    
	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used 

	double info[LM_INFO_SZ]; // information regarding the minimization

    int maxiteration = 1000; // max iteration

    //
	// 1. Rosenbrock
    //

	m = 2; n = 2;
	p[0] = -1.0; p[1] = 1.0;
	x[0] = 0.0; x[1] = 0.0;
	ret = dlevmar_der(rosenbrock, rosenbrockprime, p, x, m, n, maxiteration, opts, info, NULL, NULL, NULL);

	printf("Results for Rosenberg\n");
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (i = 0; i < m; ++i)
		printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", rosenbrock_p[i]);
	printf("\n\nMinimization info:\n");
	for (i = 0; i < LM_INFO_SZ; ++i)
		printf("%g ", info[i]);

	printf("\n");
	printf("\n");
	printf("-----------------------------------------\n");
	printf("\n");

    //
	// 2. Norris - Linear Regression
    //

	m = 2; n = 36;
	p[0] = 1.0; p[1] = 1.0;
	
	ret = dlevmar_der(norris, norrisprime, p, norris_y, m, n, maxiteration, opts, info, NULL, NULL, NULL);

	printf("Results for Norris\n");
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (i = 0; i < m; ++i)
		printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", norris_p[i]);
	printf("\n\nMinimization info:\n");
	for (i = 0; i < LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	
	printf("\n");
	printf("\n");
	printf("-----------------------------------------\n");
	printf("\n");

	//
    // 3. Lanczos1 - Non-linear Regression
    //

	m = 6; n = 24;
	p[0] = 1.2; p[1] = 0.3; p[2] = 5.6; p[3] = 5.5; p[4] = 6.5; p[5] = 7.6;

	ret = dlevmar_dif(lanczos1, p, lanczos1_y, m, n, maxiteration, opts, info, NULL, NULL, lanczos1_x);

	printf("Results for Lanczos1\n");
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (i = 0; i < m; ++i)
		printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", lanczos1_p[i]);
	printf("\n\nMinimization info:\n");
	for (i = 0; i < LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	
	printf("\n");
	printf("\n");
	printf("-----------------------------------------\n");
	printf("\n");

    //
	// 4. Thurber - Non-linear Regression
    //

	m = 7; n = 37;
	p[0] = 1000; p[1] = 1000; p[2] = 400; p[3] = 40; p[4] = 0.7; p[5] = 0.3; p[6] = 0.03;

	ret = dlevmar_dif(thurber, p, thurber_y, m, n, maxiteration, opts, info, NULL, NULL, thurber_x);

	printf("Results for Thurber\n");
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (i = 0; i < m; ++i)
		printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", thurber_p[i]);
	printf("\n\nMinimization info:\n");
	for (i = 0; i < LM_INFO_SZ; ++i)
		printf("%g ", info[i]);

    printf("\n");
    printf("\n");
    printf("-----------------------------------------\n");
    printf("\n");

    //
    // 5. Rat43 - Non-linear Regression
    //

    m = 4; n = 15;
    p[0] = 100; p[1] = 10; p[2] = 1; p[3] = 1;
        
    ret = dlevmar_dif(rat43, p, rat43_y, m, n, maxiteration, opts, info, NULL, NULL, rat43_x);

    printf("Results for Rat43\n");
    printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
    for (i = 0; i < m; ++i)
        printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", rat43_p[i]);
    printf("\n\nMinimization info:\n");
    for (i = 0; i < LM_INFO_SZ; ++i)
        printf("%g ", info[i]);

    printf("\n");
    printf("\n");
    printf("-----------------------------------------\n");
    printf("\n");

    //
    // 6. Randle - Complex non-linear Regression
    //

    m = 3; n = 181;
    p[0] = 1; p[1] = 1; p[2] = 1e-3; // initial values of parameters    

    double ydata[362];
    memcpy(ydata, randle_z1, n * sizeof(double));
    memcpy(&ydata[n], randle_z2, n * sizeof(double));

    // Cs may be going to be negative. so set lower bound!!!
    lb[0] = 1e-20; lb[1] = 1e-20; lb[2] = 1e-20; // set lower bound

    //
    // 6.1 analytic Jacobian
    //
    ret = dlevmar_bc_der(randle, randleprime, p, ydata, m, 2 * n, lb, NULL, NULL, maxiteration, opts, info, NULL, NULL, randle_f);

    printf("Results for Randle - with analytic Jacobian\n");
    printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
    for (i = 0; i < m; ++i)
        printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", randle_p[i]);
    printf("\n\nMinimization info:\n");
    for (i = 0; i < LM_INFO_SZ; ++i)
        printf("%g ", info[i]);

    printf("\n");
    printf("\n");
    printf("-----------------------------------------\n");
    printf("\n");

    //
    // 6.2 finite difference approximated Jacobian
    //
    ret = dlevmar_bc_dif(randle, p, ydata, m, 2 * n, lb, NULL, NULL, maxiteration, opts, info, NULL, NULL, randle_f);

    printf("Results for Randle - with finite difference approximated Jacobian \n");
    printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
    for (i = 0; i < m; ++i)
        printf("%12.7g ", p[i]);
    printf("\nExpected: ");
    for (i = 0; i < m; ++i)
        printf("%12.7g ", randle_p[i]);
    printf("\n\nMinimization info:\n");
    for (i = 0; i < LM_INFO_SZ; ++i)
        printf("%g ", info[i]);

    printf("\n");
    printf("\n");
    printf("-----------------------------------------\n");
    printf("\n");

	return 0;
}
