// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include "lib\levmar.h"

//
// 1. Rosenbrock
//
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

//
// 2. Norris
// https://www.itl.nist.gov/div898/strd/lls/data/Norris.shtml
//
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

//
// 3. Lanczos1
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

//
// 4. Thurber
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


int main()
{
	int i, ret;

	int m; // parameter vector dimension
	int n; // measurement vector dimension

	double p[16]; // initial parameter estimates
	double x[16]; // measurement vector

	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used 

	double info[LM_INFO_SZ]; // information regarding the minimization

	// 1. Rosenbrock
	m = 2; n = 2;
	p[0] = -1.0; p[1] = 1.0;
	x[0] = 0.0; x[1] = 0.0;
	ret = dlevmar_der(rosenbrock, rosenbrockprime, p, x, m, n, 10000, opts, info, NULL, NULL, NULL);

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

	// 2. Norris - Linear Regression
	m = 2; n = 36;
	p[0] = 1.0; p[1] = 1.0;
	
	ret = dlevmar_der(norris, norrisprime, p, norris_y, m, n, 10000, opts, info, NULL, NULL, NULL);

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

	// 3. Lanczos1 - Non-linear Regression
	m = 6; n = 24;
	p[0] = 1.2; p[1] = 0.3; p[2] = 5.6; p[3] = 5.5; p[4] = 6.5; p[5] = 7.6;

	ret = dlevmar_dif(lanczos1, p, lanczos1_y, m, n, 10000, opts, info, NULL, NULL, lanczos1_x);

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

	// 4. Thurber - Non-linear Regression
	m = 7; n = 37;
	p[0] = 1000; p[1] = 1000; p[2] = 400; p[3] = 40; p[4] = 0.7; p[5] = 0.3; p[6] = 0.03;

	ret = dlevmar_dif(thurber, p, thurber_y, m, n, 10000, opts, info, NULL, NULL, thurber_x);

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

	return 0;
}
