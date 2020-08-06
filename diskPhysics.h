/*
 *  diskPhysics.h
 *  
 *
 *  Created by Erik on 4/17/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
//#include <gsl/gsl_sf_ellint.h>
//#include <gsl/gsl_integration.h>

#include "vcl/vectorclass.h"
#include "vcl/vectormath_trig.h"
#include "vcl/vectormath_exp.h"

const double gravConst=6.67259e-8;
const double h=6.626e-27; //erg * seconds
const double hbar=1.054e-27; //erg* seconds
const double c=2.99792458e+10; //cm/sec
const double kboltzmann=1.380658e-16; //erg/K
const double sigma=5.67e-5; //Stefan Boltzmann constant in cgs (erg/cm^2/sec/K^4)
const double pi=4*atan(1);
const double AU=1.496e13;
const double amu=1.670539040e-24;//g
const double mSun=1.998e33;//g

//TODO: Move stellar properties to the input file
const double Rstar=6.61e+10; //sun
const double Tstar=2800; //K (brown dwarf)

//coming up are a bunch of long lists of numbers. The first sets are the transition temperatures and frequencies
//of the various CO isotopologues. The last two are the frequencies and opacities of dust in a "typical" 
//protostellar disk. 

static double transitionTemps_12CO[40]={ 5.53, 16.60, 33.19, 55.32, 82.97, 116.16, 154.87, 199.11,
	248.88, 304.16, 364.97, 431.29, 503.13, 580.49, 663.35, 751.72, 845.59, 944.97, 1049.84, 1160.20, 
	1276.05, 1397.38, 1524.19, 1656.47, 1794.23, 1937.44, 2086.12, 2240.24, 2399.82, 2564.83, 2735.28, 
	2911.15, 3092.45, 3279.15, 3471.27, 3668.78, 3871.69, 4079.98, 4293.64, 4512.67};
static double transitionFreqs_12CO[40]={1.152712018e+11, 2.30538e+11, 3.457959899e+11, 4.610407682e+11,
	5.762679305e+11, 6.914730763e+11, 8.06651806e+11, 9.217997e+11, 1.036912393e+12, 1.151985452e+12, 
	1.267014486e+12, 1.381995105e+12, 1.496922909e+12, 1.611793518e+12, 1.7266025057e+12, 1.841345506e+12, 
	1.956018139e+12, 2.070615993e+12, 2.18513468e+12, 2.299569842e+12, 2.413917113e+12, 2.52817206e+12, 
	2.6423303459e+12, 2.756387584e+12, 2.870339407e+12, 2.984181455e+12, 3.097909361e+12, 3.2115187506e+12, 
	3.3250052827e+12, 3.438364611e+12, 3.551592361e+12, 3.66468418e+12, 3.777635728e+12, 3.890442717e+12,
	4.0031007876e+12, 4.115605585e+12, 4.2279527744e+12, 4.340138112e+12, 4.4521571221e+12, 4.5640056399e+12};
static double einsteinA_12CO[40]={7.203e-08,6.910e-07,2.497e-06,6.126e-06,1.221e-05,2.137e-05,3.422e-05,
	5.134e-05,7.330e-05,1.006e-04,1.339e-04,1.735e-04,2.200e-04,2.739e-04,3.354e-04,4.050e-04,4.829e-04,
	5.695e-04,6.650e-04,7.695e-04,8.833e-04,1.006e-03,1.139e-03,1.281e-03,1.432e-03,1.592e-03,1.761e-03,
	1.940e-03,2.126e-03,2.321e-03,2.524e-03,2.735e-03,2.952e-03,3.175e-03,3.404e-03,3.638e-03,3.878e-03,
	4.120e-03,4.365e-03,4.613e-03};
	
static double transitionTemps_13CO[40]={ 5.29, 15.87, 31.73, 52.89, 79.33, 111.05, 148.06, 190.36, 237.93, 
	290.79, 348.92, 412.33, 481.02, 554.97, 634.20, 718.69, 808.44, 903.45, 1003.72, 1109.24, 1220.00, 1336.01,
	1457.27, 1583.75, 1715.47, 1852.42, 1994.58, 2141.96, 2294.55, 2452.35, 2615.35, 2783.54, 2956.91, 3135.47, 
	3319.20, 3508.10, 3702.16, 3901.37, 4105.73, 4315.23 };
static double transitionFreqs_13CO[40]={ 1.102013542798e+11, 2.203986841281e+11, 3.305879652218e+11, 4.407651734547e+11, 
	5.509262850456e+11, 6.610672766472e+11, 7.711841254539e+11, 8.812728093107e+11, 9.913293068214e+11, 
	1.1013495974571e+12, 1.2113296616644e+12, 1.3212654809740e+12, 1.4311530381090e+12, 1.5409883170934e+12, 
	1.6507673033603e+12, 1.7604859838606e+12, 1.8701403471712e+12, 1.9797263836035e+12, 2.0892400853116e+12, 
	2.1986774464010e+12, 2.3080344630367e+12, 2.4173071335521e+12, 2.5264914585567e+12, 2.6355834410451e+12, 
	2.7445790865050e+12, 2.8534744030259e+12, 2.9622654014074e+12, 3.0709480952674e+12, 3.1795185011510e+12, 
	3.2879726386382e+12, 3.3963065304531e+12, 3.5045162025716e+12, 3.6125976843302e+12, 3.7205470085342e+12, 
	3.8283602115665e+12, 3.9360333334954e+12, 4.0435624181834e+12, 4.1509435133956e+12, 4.2581726709079e+12, 4.3652459466157e+12 };
static double einsteinA_13CO[40]={6.294E-08,6.038E-07,2.181E-06,5.353E-06,1.067E-05,1.868E-05,2.991E-05,
	4.487E-05,6.406E-05,8.797E-05,1.170E-04,1.517E-04,1.924E-04,2.395E-04,2.934E-04,3.543E-04,4.225E-04,
	4.984E-04,5.821E-04,6.737E-04,7.735E-04,8.815E-04,9.978E-04,1.122E-03,1.255E-03,1.396E-03,1.545E-03,
	1.702E-03,1.866E-03,2.038E-03,2.217E-03,2.403E-03,2.594E-03,2.792E-03,2.995E-03,3.203E-03,3.414E-03,
	3.630E-03,3.848E-03,4.068E-03};
				
static double transitionTemps_C18O[40]={5.27, 15.81, 31.61, 52.68, 79.02, 110.63, 147.50, 189.63, 237.03, 289.68, 347.60, 
	410.77, 479.19, 552.86, 631.79, 715.95, 805.36, 900.02, 999.90, 1105.02, 1215.37, 1330.94, 1451.73, 1577.74, 1708.96, 
	1845.38, 1987.01, 2133.83, 2285.85, 2443.05, 2605.43, 2772.98, 2945.70, 3123.58, 3306.62, 3494.81, 3688.13, 3886.60, 4090.18, 4298.89};
static double transitionFreqs_C18O[40]={1.097821734e11, 2.195603541e11, 3.293305525e11, 4.390887658e11, 5.488310055e11, 
	6.585532782e11, 7.682515933e11, 8.779219553e11, 9.875603822e11, 1.0971628753e12, 1.2067254487e12, 1.3162441143e12, 
	1.4257148854e12, 1.5351337765e12, 1.6444968033e12, 1.7537999810e12, 1.8630393590e12, 1.9722108700e12, 2.0813106385e12, 
	2.1903346365e12, 2.2992788998e12, 2.4081394565e12, 2.5169123372e12, 2.6255935745e12, 2.7341792037e12, 2.8426652625e12, 
	2.9510477910e12, 3.0593228324e12, 3.1674864325e12, 3.2755346398e12, 3.3834635063e12, 3.4912690868e12, 3.5989474393e12, 
	3.7064946254e12, 3.8139067097e12, 3.9211797609e12, 4.0283098507e12, 4.1352930551e12, 4.2421254535e12, 4.3488031296e12};
static double einsteinA_C18O[40]={6.266e-08,6.011e-07,2.172e-06,5.330e-06,1.062e-05,1.860e-05,2.978e-05,
	4.468e-05,6.380e-05,8.762e-05,1.166e-04,1.512e-04,1.918e-04,2.388e-04,2.925e-04,3.533e-04,4.214e-04,
	4.972e-04,5.808e-04,6.725e-04,7.723e-04,8.803e-04,9.969e-04,1.122e-03,1.255e-03,1.396e-03,1.546e-03,
	1.704e-03,1.869e-03,2.042e-03,2.222e-03,2.409e-03,2.603e-03,2.802e-03,3.008e-03,3.219e-03,3.433e-03,
	3.652e-03,3.875e-03,4.100e-03};
				
static double freqvals[206]={3000000000, 1e+10, 1.5e+10, 1.875e+10, 3e+10, 4.285714286e+10, 9.375e+10, 1.714285714e+11,
	1.973684211e+11, 2.272727273e+11, 2.608695652e+11, 3.009027081e+11, 3.460207612e+11, 3.978779841e+11, 4.573170732e+11, 
	5.263157895e+11, 6.048387097e+11, 6.960556845e+11, 8.571428571e+11, 1.034482759e+12, 1.214574899e+12, 1.401869159e+12, 
	1.538461538e+12, 1.694915254e+12, 2e+12, 2.255639098e+12, 2.479338843e+12, 2.727272727e+12, 3e+12, 3.300330033e+12, 
	3.627569528e+12, 3.989361702e+12, 4.39238653e+12, 4.830917874e+12, 5.309734513e+12, 5.847953216e+12, 6.423982869e+12, 
	7.075471698e+12, 7.772020725e+12, 8.522727273e+12, 9.404388715e+12, 1.034482759e+13, 1.140684411e+13, 1.255230126e+13, 
	1.265822785e+13, 1.282051282e+13, 1.293103448e+13, 1.298701299e+13, 1.304347826e+13, 1.315789474e+13, 1.333333333e+13, 
	1.382488479e+13, 1.401869159e+13, 1.41509434e+13, 1.421800948e+13, 1.428571429e+13, 1.442307692e+13, 1.463414634e+13, 
	1.52284264e+13, 1.595744681e+13, 1.621621622e+13, 1.639344262e+13, 1.648351648e+13, 1.657458564e+13, 1.675977654e+13, 
	1.704545455e+13, 1.764705882e+13, 1.796407186e+13, 1.818181818e+13, 1.829268293e+13, 1.840490798e+13, 1.863354037e+13, 
	1.886792453e+13, 1.898734177e+13, 1.910828025e+13, 1.923076923e+13, 1.948051948e+13, 1.973684211e+13, 1.986754967e+13, 
	2e+13, 2.027027027e+13, 2.054794521e+13, 2.083333333e+13, 2.097902098e+13, 2.112676056e+13, 2.127659574e+13, 2.158273381e+13, 
	2.173913043e+13, 2.189781022e+13, 2.205882353e+13, 2.222222222e+13, 2.23880597e+13, 2.255639098e+13, 2.272727273e+13, 
	2.290076336e+13, 2.325581395e+13, 2.34375e+13, 2.362204724e+13, 2.380952381e+13, 2.4e+13, 2.419354839e+13, 2.459016393e+13, 
	2.479338843e+13, 2.5e+13, 2.521008403e+13, 2.542372881e+13, 2.564102564e+13, 2.586206897e+13, 2.608695652e+13, 2.631578947e+13, 
	2.654867257e+13, 2.678571429e+13, 2.702702703e+13, 2.777777778e+13, 2.830188679e+13, 2.884615385e+13, 2.97029703e+13, 
	3.06122449e+13, 3.157894737e+13, 3.260869565e+13, 3.370786517e+13, 3.448275862e+13, 3.488372093e+13, 3.529411765e+13, 
	3.571428571e+13, 3.658536585e+13, 3.797468354e+13, 3.846153846e+13, 3.896103896e+13, 3.947368421e+13, 4e+13, 4.109589041e+13, 
	4.166666667e+13, 4.225352113e+13, 4.285714286e+13, 4.347826087e+13, 4.411764706e+13, 4.545454545e+13, 4.6875e+13, 
	4.761904762e+13, 4.838709677e+13, 4.918032787e+13, 5.084745763e+13, 5.263157895e+13, 5.454545455e+13, 5.660377358e+13, 
	5.769230769e+13, 5.882352941e+13, 6.060606061e+13, 6.25e+13, 6.666666667e+13, 7.317073171e+13, 8.108108108e+13, 8.450704225e+13, 
	8.695652174e+13, 8.823529412e+13, 8.955223881e+13, 9.090909091e+13, 9.230769231e+13, 9.375e+13, 9.677419355e+13, 1.034482759e+14, 
	1.153846154e+14, 1.234567901e+14, 1.363636364e+14, 1.492537313e+14, 1.639344262e+14, 1.818181818e+14, 1.986754967e+14, 
	2.189781022e+14, 2.4e+14, 2.654867257e+14, 2.912621359e+14, 3.201707577e+14, 3.333333333e+14, 3.521126761e+14, 3.875968992e+14,
	4.285714286e+14, 4.6875e+14, 5.154639175e+14, 5.454545455e+14, 5.671077505e+14, 6.237006237e+14, 6.818181818e+14, 
	7.556675063e+14, 8.219178082e+14, 9.146341463e+14, 1.003344482e+15, 1.10701107e+15, 1.214574899e+15, 1.339285714e+15, 
	1.376146789e+15, 1.470588235e+15, 1.621621622e+15, 1.775147929e+15, 1.960784314e+15, 2.158273381e+15, 2.362204724e+15, 
	2.608695652e+15, 2.857142857e+15, 3.154574132e+15, 3.289473684e+15, 3.468208092e+15, 3.816793893e+15, 4.195804196e+15, 4.615384615e+15};
	
	static double opacvals[206]={4.794e-06, 1.56e-05, 2.937e-05, 3.967e-05, 7.057e-05, 0.0001189, 0.0003892, 0.001064, 0.001363, 
	0.001778, 0.002357, 0.003118, 0.00449, 0.006935, 0.01123, 0.0186, 0.03325, 0.06494, 0.1219, 0.2009, 0.2775, 0.3822, 
	0.4548, 0.5367, 0.6828, 0.8127, 0.9177, 1.024, 1.131, 1.24, 1.344, 1.448, 1.548, 1.644, 1.738, 1.829, 1.905, 1.992, 
	2.101, 2.198, 2.33, 2.426, 2.434, 2.402, 2.402, 2.402, 2.402, 2.402, 2.402, 2.407, 2.415, 2.435, 2.443, 2.448, 2.451, 
	2.453, 2.464, 2.48, 2.524, 2.574, 2.584, 2.591, 2.594, 2.597, 2.601, 2.603, 2.604, 2.603, 2.602, 2.592, 2.582, 2.56, 
	2.537, 2.526, 2.514, 2.502, 2.478, 2.453, 2.439, 2.426, 2.415, 2.405, 2.394, 2.388, 2.382, 2.376, 2.364, 2.358, 2.351, 
	2.345, 2.338, 2.328, 2.318, 2.308, 2.297, 2.273, 2.26, 2.247, 2.263, 2.278, 2.292, 2.321, 2.334, 2.347, 2.371, 2.394, 
	2.417, 2.439, 2.462, 2.483, 2.496, 2.51, 2.523, 2.553, 2.57, 2.578, 2.574, 2.542, 2.459, 2.275, 1.945, 1.89, 1.873, 
	1.855, 1.836, 1.802, 1.802, 1.8, 1.798, 1.797, 1.81, 1.856, 1.878, 1.902, 1.925, 1.948, 1.969, 1.994, 2.012, 2.021, 
	2.03, 2.041, 2.06, 1.943, 1.819, 1.693, 1.685, 1.676, 1.668, 1.648, 1.642, 1.704, 1.829, 1.935, 2.1, 2.208, 2.302, 
	2.281, 2.258, 2.233, 2.176, 2.153, 1.947, 1.848, 1.821, 1.77, 1.819, 1.875, 1.936, 1.946, 1.89, 1.936, 1.967, 2.027, 
	2.056, 2.095, 2.093, 2.123, 2.19, 2.334, 2.415, 2.468, 2.603, 2.744, 2.958, 3.092, 3.214, 3.255, 3.275, 3.278, 3.276, 
	3.274, 3.267, 3.248, 3.218, 3.175, 3.131, 3.102, 3.073, 3.034, 2.982, 2.939, 2.867, 2.699, 2.501, 2.292};

//blackbody distribution
double blackBody(const double temp, const double freq);

Vec4d blackBodyAVX(const double temp, const Vec4d freq);

double tempFromSurfaceBrightness(const double intensity, const double freq);

double doppler(double frequency, double velocity);

Vec4d dopplerAVX(Vec4d frequency, double velocity);

//planck mean opacity
double PMO(const double temp, const std::vector<double>& frequencies, const std::vector<double>& opacities);

//surface temperature based on star properties
double Tsurf(double r, const std::vector<double>& frequencies,const std::vector<double>& opacities);

double temperature(double r, double z, double freq, 
const std::vector<double>& frequencies,const std::vector<double>& opacities);

//Given the mass and other density parameters, gives the surface mass density coefficient
double sigma0FromMass(double mass, double P, double r0, double rmin, double rmax);

double gaussian(double x, double mu, double sigma);

double gaussianNotNorm(double x, double mu, double sigma);

double dustOpac(const double frequency);

Vec4d dustOpacAVX(const Vec4d frequency);

struct expCutoffDensity {
        double Sigma0;
        double rc;
        double r0;
        double h0;
        double P;//density index
        double S;//scale height index
        double p1, p2, p3, p4; //ring positions (au)
        double d1, d2, d3, d4; //ring depths (0 to 1)
        double w1, w2, w3, w4; //ring widths (au)

        expCutoffDensity(): Sigma0(0), rc(0), r0(0), h0(0), P(0), S(0), p1(0), p2(0), p3(0), p4(0), d1(0), d2(0), d3(0), d4(0), w1(0), w2(0), w3(0), w4(0) { }
        expCutoffDensity(const expCutoffDensity& other): Sigma0(other.Sigma0), rc(other.rc), r0(other.r0), h0(other.h0),
                P(other.P), S(other.S), p1(other.p1), p2(other.p2), p3(other.p3), p4(other.p4), d1(other.d1), d2(other.d2), 
				d4(other.d4), d3(other.d3), w1(other.w1), w2(other.w2), w3(other.w3), w4(other.w4) {  }
        expCutoffDensity(double Sigma0, double rc, double r0, double h0, double P, double S, double p1, double p2, double p3, double p4, double d1, double d2, double d3, double d4, double w1, double w2, double w3, double w4):
                Sigma0(Sigma0), rc(rc), r0(r0), h0(h0), P(P), S(S), p1(p1), p2(p2), p3(p3), p4(p4), d1(d1), d2(d2), d3(d3), d4(d4), w1(w1), w2(w2), w3(w3), w4(w4) {    }

        expCutoffDensity& operator= (const expCutoffDensity& other){
                Sigma0=other.Sigma0; rc=other.rc; r0=other.r0; h0=other.h0; P=other.P; S=other.S;
                p1=other.p1; p2=other.p2; d1=other.d1; d2=other.d2; w1=other.w1; w2=other.w2;
                return *this;
        }

        double surfaceMassDensity(const double r) const{
                return (Sigma0*pow(r/rc, -P)) * exp(-(pow(r/rc, 2-P)));
        }

        double scaleHeight(const double r) const{
                return (h0*pow(r/r0,S));
        }
		
        /*double scaleHeight(const double r) const{
			double T100=14.2243;
        	double prefactor=sqrt(kboltzmann*T100/gravConst/1.12/mSun/3.819239518e-24)/pow(100*AU,-0.25);
			return prefactor*pow(r,1.25);
        }*/

        double operator()(double r, double theta, double phi) const{
                double r_cyl=r*sin(theta);
                double z=r*cos(theta);
                double h=scaleHeight(r_cyl);
                double ring1=(1-d1*(gaussianNotNorm(r_cyl,p1,w1)));
                double ring2=(1-d2*(gaussianNotNorm(r_cyl,p2,w2)));
                double ring3=(1-d3*(gaussianNotNorm(r_cyl,p3,w3)));
				double ring4=(1-d4*(gaussianNotNorm(r_cyl,p4,w4)));

                return((surfaceMassDensity(r_cyl)/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h))*ring1*ring2*ring3*ring4);
        }
};


struct expCutoffDensityNoRings {
	double Sigma0;
	double rc;
	double r0;
	double h0;
	double P;//density index
	double S;//scale height index
	
	expCutoffDensityNoRings(): Sigma0(0), rc(0), r0(0), h0(0), P(0), S(0) {	}
	expCutoffDensityNoRings(const expCutoffDensity& other): Sigma0(other.Sigma0), rc(other.rc), r0(other.r0), h0(other.h0), P(other.P), S(other.S) {	}
	expCutoffDensityNoRings(double Sigma0, double rc, double r0, double h0, double P, double S): Sigma0(Sigma0), rc(rc), r0(r0), h0(h0), P(P), S(S) {	}
		
	expCutoffDensityNoRings& operator= (const expCutoffDensityNoRings& other){
		Sigma0=other.Sigma0; rc=other.rc; r0=other.r0; h0=other.h0; P=other.P; S=other.S; 
		return *this;
	}
	
	double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/rc, -P)) * exp(-(pow(r/rc, 2-P)));
	}
	
	double scaleHeight(const double r) const{
		return (h0*pow(r/r0,S));
	}
	
	double operator()(double r, double theta, double phi) const{
		double r_cyl=r*sin(theta);
		double z=r*cos(theta);
		double h=scaleHeight(r_cyl);
		
		return((surfaceMassDensity(r_cyl)/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h)));
	}
};

//Modified density for testing nonaxisymmetric disks
struct densityPert {
	double Sigma0;
	double r0;
	double h0;
	double P;//density index
	double S;//scale height index
		
	densityPert(): Sigma0(0), r0(0), h0(0), P(0), S(0) {	}
	densityPert(double Sigma0, double r0, double h0, double P, double S): Sigma0(Sigma0), r0(r0), h0(h0), P(P), S(S) {	}
	
	/*double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/r0, -P));
	}*/
	double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/r0, -P));//*exp(-pow(r/r0, 2-P)));
	}
	
	double scaleHeight(const double r) const{
		return (h0*pow(r/r0,S));
	}
	
	/*double operator()(double r, double theta, double phi) const{
		//double mass=3.978e+30;
		double r_cyl=r*sin(theta);
		double z=r*cos(theta);
		double h=scaleHeight(r_cyl);
		//double perturbation=0;
		double SMD=surfaceMassDensity(r_cyl)/(sqrt(2*pi)*h) * exp(-z*z/(2*h*h));
		//perturbation=mass*gaussian(phi,pi,1)*gaussian(r,50*AU,8*AU)/r/r;
		return(SMD);
		//return 5e-19;
	}*/
	
	double operator()(double r, double theta, double phi) const{
		double mass=3.978e+30;
		double r_cyl=r*sin(theta);
		//double z=r*cos(theta);
		//double h=scaleHeight(r_cyl);
		double perturbation=0;
		//double SMD=surfaceMassDensity(r_cyl)/(sqrt(2*pi)*h) * exp(-z*z/(2*h*h));
		perturbation=mass*gaussian(phi,pi/2,0.5)*gaussian(r,150*AU,10*AU)/r_cyl/r_cyl;
		//perturbation=(1.0-(sqrt(2*pi)*10*AU*gaussian(r,150*AU,10*AU)));
		/*double rmax=84*AU; double pitch=0.138;	double stddev=4*AU;
		perturbation=exp(-(r_cyl-rmax*exp(pitch*(phi)))*(r_cyl-rmax*exp(pitch*(phi)))/(stddev*stddev));	//first spiral
		//the second spiral requires some finagling
		if(phi>=pi){
			double phi2=phi-pi;
			perturbation+=exp(-(r_cyl-rmax*exp(pitch*(phi2)))*(r_cyl-rmax*exp(pitch*(phi2)))/(stddev*stddev));
		} else if(phi<pi){
			double phi2=phi+pi;
			perturbation+=exp(-(r_cyl-rmax*exp(pitch*(phi2)))*(r_cyl-rmax*exp(pitch*(phi2)))/(stddev*stddev));
		}
		perturbation*=SMD;
		perturbation*= 0.08236*r_cyl/AU - 2.396;*/
		return(perturbation);
	}
};

//used for velocity calculation. See Bertin&Lodato eq 4 for details of the integral
//double integrand(double rprime, void* p);

enum isotopologue_type { iso_12CO, iso_13CO, iso_C17O, iso_C18O };

//Opacity of C12O16. See Isella et al 2007
//This is a functor so that I can do as much of the calculation ahead of time as possible
//The calculation of the partition function is based on a Taylor expansion. For details,
//see Mangum and Shirley 2015.
//Later I'll add capability to choose which isotopologue you want. For now C12O16 is good.
struct COopacFast {
	double COFraction, dipoleMoment, COmass,B;
	double densityRatio; //ratio of chosen isotopologue to 12CO
	int lowerEnergyLevel;
	isotopologue_type isot;
	double* transitionTemps;
	double* restFrequencies;
	double* einsteinAs;
		
	COopacFast(int lowerEnergyLevel, isotopologue_type isot) : lowerEnergyLevel(lowerEnergyLevel){
		setTransitionTemps(isot);
		COFraction = 1e-5; //ratio of 12CO to H2
		COFraction *= 14; //convert from density ratio to number ratio
		dipoleMoment = 1.1e-19; //dipole moment of CO
		B=5.763596e+10; //Hz. The rigid rotor rotation constant for CO
	}
	
	COopacFast(){
		setTransitionTemps(iso_12CO);
		lowerEnergyLevel = 0;
		COFraction = 1e-5; //density ratio of CO to H2
		COFraction *= 14; //convert from density ratio to number ratio
		dipoleMoment = 1.12e-19; //dipole moment of CO
		B=5.763596e+10; //Hz. The rigid rotor rotation constant for CO
	}
	
	double operator()(double temperature, double frequency, double vTurb, bool freezeout) const{
		double freezeOutFactor=1;
		if(freezeout){
			if(temperature < 30)
				freezeOutFactor=temperature/20.0-0.5;
			if(temperature < 10)
				return 0;
		}
		double nu0=restFrequencies[lowerEnergyLevel];
		double partitionFunction = (kboltzmann*temperature)/(h*B) + 1.0/3.0 + (h*B)/(15.0*kboltzmann*temperature)
			+4.0/315.0 *(h*B)*(h*B)/(kboltzmann*temperature)/(kboltzmann*temperature);
		double boltzFactor= multiplicity(lowerEnergyLevel)*exp(-energy(lowerEnergyLevel)/(kboltzmann*temperature))/partitionFunction;
		double deltaV=sqrt(2*kboltzmann*temperature/COmass + vTurb*vTurb);
		double deltaNu= c/nu0 * (frequency - nu0);
		//std::cout.precision(12);
		//std::cout << nu0 << "\t" << frequency << std::endl;
		double profile= c/(nu0*sqrt(pi)*deltaV)*exp(-deltaNu*deltaNu/deltaV/deltaV);
		//std::cout << "nu0: " << nu0 << "\tfreq: " << frequency << "\tdeltaNu: "<< deltaNu << "\tdeltaV: " << deltaV <<std::endl;
		double einsteinA=einsteinAs[lowerEnergyLevel];
		double opacity=c*c/(8*pi*nu0*nu0);
		//std::cout << "c*c/(8*pi*nu0*nu0): " <<  opacity << std::endl;
		opacity*=multiplicity(lowerEnergyLevel+1)/multiplicity(lowerEnergyLevel);
		//std::cout << "mult factor: " << multiplicity(lowerEnergyLevel+1)/multiplicity(lowerEnergyLevel) << std::endl;
		opacity*=einsteinA;
		//std::cout << "Einstein A: " << einsteinA << std::endl;		
		//double profile=2.60076e-06;
		//std::cout << profile << std::endl;
		opacity*=profile;
		//std::cout << "profile: " << profile << std::endl;
		opacity*=boltzFactor;
		opacity*=(1-exp(-h*nu0/(kboltzmann*temperature)));
		//std::cout << "exp factor: " << (1-exp(-h*nu0/(kboltzmann*temperature))) << std::endl;
		//std::cout << "opac? " << c*c/(8*pi*nu0*nu0) * multiplicity(lowerEnergyLevel+1)/multiplicity(lowerEnergyLevel) * einsteinA * profile * (1-exp(-h*nu0/(kboltzmann*temperature))) << std::endl;
		//std::cout << "nu0: " << nu0 << "\tTemp: " << temperature << std::endl;
		//std::cout << "Boltz: " << boltzFactor << "\tprof: " << profile << "\tA21: " << einsteinA << "\t" << (1-exp(-h*nu0/(kboltzmann*temperature))) << std::endl;
		
		//the integrator works with the mass density, we need to convert to number density here
		opacity*=(COFraction*densityRatio/COmass);
		
		return opacity*freezeOutFactor;
	}
	
	//AVX version
	Vec4d operator()(double temperature, Vec4d frequency, double vTurb, bool freezeout) const{
		double freezeOutFactor=1;
		if(freezeout){
			if(temperature < 30)
				freezeOutFactor=temperature/20.0-0.5;
			if(temperature < 10)
				return 0;
		}
		double nu0=restFrequencies[lowerEnergyLevel];
		double partitionFunction = (kboltzmann*temperature)/(h*B) + 1.0/3.0 + (h*B)/(15.0*kboltzmann*temperature)
			+4.0/315.0 *(h*B)*(h*B)/(kboltzmann*temperature)/(kboltzmann*temperature);
		double boltzFactor= multiplicity(lowerEnergyLevel)*exp(-energy(lowerEnergyLevel)/(kboltzmann*temperature))/partitionFunction;
		double deltaV=sqrt(2*kboltzmann*temperature/COmass + vTurb*vTurb);
		
		Vec4d deltaNu= (frequency - nu0) * c/nu0;
		Vec4d profile= c/(nu0*sqrt(pi)*deltaV)*exp(-deltaNu*deltaNu/deltaV/deltaV);
		double einsteinA=einsteinAs[lowerEnergyLevel];
		Vec4d opacity=c*c/(8*pi*nu0*nu0);
		opacity*=multiplicity(lowerEnergyLevel+1)/multiplicity(lowerEnergyLevel);
		opacity*=einsteinA;
		opacity*=profile;
		opacity*=boltzFactor;
		opacity*=(1-exp(-h*nu0/(kboltzmann*temperature)));
		
		//the integrator works with the mass density, we need to convert to number density here
		opacity*=(COFraction*densityRatio/COmass);
		
		return opacity*freezeOutFactor;
	}
	
	double multiplicity(int level) const {
		return 2*level+1;
	}
	
	double energy(int level) const{
		double T1=transitionTemps[0];
		return 0.5*kboltzmann*level*(level+1)*T1;
	}
		
	//Transition temperatures, energies, and frequencies come from the LAMDA database
	//http://home.strw.leidenuniv.nl/~moldata/CO.html
	void setTransitionTemps(isotopologue_type isot){
		switch (isot){
			case iso_12CO:
			densityRatio=1;
			transitionTemps = transitionTemps_12CO;
			restFrequencies=transitionFreqs_12CO;
			einsteinAs=einsteinA_12CO;
			COmass=28*amu;
			break;
			case iso_13CO:
			densityRatio=(1.0/70.0);
			transitionTemps = transitionTemps_13CO;
			restFrequencies = transitionFreqs_13CO;
			einsteinAs=einsteinA_13CO;
			COmass=29*amu;
			break;
			case iso_C17O:
			//I'll add the numbers when I have to.
			std::cout << "This CO isotopologue isn't in yet. Add the numbers" << std::endl;
			exit(1);
			break;
			case iso_C18O:
			densityRatio=(1.0/500.0);
			transitionTemps = transitionTemps_C18O;
			restFrequencies = transitionFreqs_C18O;
			einsteinAs=einsteinA_C18O;
			COmass=30*amu;
			break;	
		}
	}
};

//H2(16)O
//This molecule is a lot more complicated than CO.
//So this only covers the 3_13 - 2_20 transition. There are enough transitions that you can't just specify
//the lower level and be done. For now the data is hardcoded to be stuff I need. We probably won't ever do 
//other water lines, so it doesn't need to be done like the CO.
//The partition function comes from Herzberg (1945) eq V28
struct H2Oopac {
	double H2Omass,A,B,C, H2OFraction,nu0,energy;
	
	/*H2Oopac(int test){
		H2Omass= 18*amu;
		A=8.358403e11;
		B=4.353517e11;
		C=2.781387e11;
		H2OFraction=1e-4; //from Hogerheijde et al 2011 (Science)
		H2OFraction *= 9; //convert from density ratio to number ratio
		nu0=1.83310087e11; //hardcoded for now.
		energy=142.278460 * c*h; // energy of the upper state
	}*/
	
	H2Oopac(int test){
		H2Omass= 20*amu;
		A=8.358403e11;
		B=4.353517e11;
		C=2.781387e11;
		H2OFraction=1e-4; //from Hogerheijde et al 2011 (Science)
		H2OFraction/=500.0;
		H2OFraction *= 10; //convert from density ratio to number ratio
		nu0=2.03407e11; //hardcoded for now.
		energy=142.278460 * c*h; // energy of the upper state
	}
	
	double operator()(double temperature, double frequency, double vTurb, bool freezeout) const{
		double freezeOutFactor=1;
		if(freezeout){
			if(temperature < 150)
				freezeOutFactor=1e-3;
		}
		double partitionFunction = exp(h*sqrt(B*C)/4/kboltzmann/temperature);
		partitionFunction *= sqrt(pi)/2;
		partitionFunction *= sqrt(kboltzmann*kboltzmann*kboltzmann*temperature*temperature*temperature/h/h/h/A/B/C);
		partitionFunction *= (1.0 + (1.0 - sqrt(B*C)/A)*h*sqrt(B*C)/kboltzmann/temperature/12.0);
		double deltaV=sqrt(2*kboltzmann*temperature/H2Omass + vTurb*vTurb);
		double deltaNu= c/nu0 * (frequency - nu0);
		double profile= c/(nu0*sqrt(pi)*deltaV)*exp(-deltaNu*deltaNu/deltaV/deltaV);
		double einsteinA=3.5755e-6;
		double opacity=c*c/(8*pi*nu0*nu0);
		opacity*=7*exp(-energy/kboltzmann/temperature)/partitionFunction; //7 from multiplicity
		opacity*=einsteinA;
		//double profile=2.60076e-06;
		opacity*=profile;
		opacity*=(1-exp(-h*nu0/(kboltzmann*temperature)));
		
		//the integrator works with the mass density, we need to convert to number density here
		opacity*=(H2OFraction/H2Omass);
		
		return opacity*freezeOutFactor;
	}
};

//basic power law. Temperature along the midplane is (roughly) a power law.
struct powerLaw {
	double T0;
	double r0;
	double q; //index
	double offset;
	
	powerLaw(): T0(0), r0(0), q(0), offset(0) {	}
	powerLaw(double T0, double r0, double q, double offset):
	 T0(T0), r0(r0), q(q), offset(offset) {	}
	
	double operator()(const double r) const{
		return(T0*pow((r/r0),-q)+offset);
	}
};

//A simple parametric temperature from Rosenfeld et al (2013)
/*struct paramTemp{
	double midptemp;	//midplane temp at reference pt
	double surftemp;	//surface temp at reference pt
	double mprefpt, surfrefpt;	//reference points
	
	paramTemp():midptemp(0),surftemp(0), mprefpt(0),surfrefpt(0){	}
	paramTemp(double midptemp, double surftemp, double mprefpt, double surfrefpt) :
		midptemp(midptemp),surftemp(surftemp),mprefpt(mprefpt),surfrefpt(surfrefpt)	{	}
	
	//I'm not sure where these values came from, so I'm leaving them hardwired
	double zq(double r) const{
		return 63*AU*pow(r/(200*AU),1.3)*exp(-(r/(800*AU))*(r/(800*AU)));
	}
	double delta(double r)const{
		return 0.0034*(r-200*AU)/AU + 2.5;
	}
	double Tmidplane(double r)const{
		return midptemp*pow(r/mprefpt,-0.3);
	}
	double Tsurf(double r, double z)const{
		return surftemp*sqrt(surfrefpt)*pow(r*r+z*z,-0.25);
	}
	double operator ()(double rSphere, double theta, double phi) const{
		//note that here I use the spherical radius
		double r=rSphere*cos(theta);
		double z=rSphere*sin(theta);
		if(z<zq(r)){
			return Tsurf(r,z)+(Tmidplane(r)-Tsurf(r,z))*pow(cos(pi*z/(2*zq(r))),2*delta(r));
		}else{
			return Tsurf(r,z);
		}
	}
};*/

//very simple temperature model
//vertically isothermal, radially goes as (r)^-1/2
struct paramTemp{
	double midptemp;	//midplane temp at reference pt
	double mprefpt;
	double dummy1, dummy2;	//dummies so that the constructor matches the more advanced case
	
	paramTemp():midptemp(0), mprefpt(0), dummy1(0), dummy2(0){	}
	paramTemp(double midptemp, double mprefpt, double dummy1, double dummy2) :
		midptemp(midptemp),mprefpt(mprefpt),dummy1(dummy1),dummy2(dummy2)	{	}

	double operator ()(double rSphere, double theta, double phi) const{
		double r=rSphere*cos(theta-(pi/2));
		return midptemp*pow(r/mprefpt,-0.5);
	}
};

//Opacity of the dust is complicated and is best calculated elsewhere.
//It can be read in from a file.
struct dustOpacity{
	std::vector<std::pair<double,double> > freqsAndOpacs;
	
	dustOpacity(){	}
	
	dustOpacity(std::string filename){	
		std::ifstream infile(filename.c_str());
		double temp1, temp2;
		if(!infile){
			std::cout << "Could not open file " << filename << std::endl;
			exit(1);
		}
		while(!infile.eof()){
			infile >> temp1 >> temp2;
			temp1/=1e4; //micron to cm
			temp1=3e10/temp1;
			freqsAndOpacs.push_back(std::pair<double,double>(temp1,temp2));
		}
		std::sort(freqsAndOpacs.begin(),freqsAndOpacs.end(),
					[&](const std::pair<double,double>& first, const std::pair<double,double>& second)->bool{
					return(first.first<second.first);
				});
		infile.close();
	}
	
	double operator () (double inputfreq) const{
		if(freqsAndOpacs.size()==0){
			std::cout << "Can't interpolate on empty opacity list" << std::endl;
			exit(1);
		}
		int index=0;
		//std::cout << "input: " << inputfreq << std::endl;
		while(inputfreq > freqsAndOpacs[index].first){
			index++;
			if(index>=freqsAndOpacs.size())
				break;
		}
		double x1,x2,y1,y2;
		x1=freqsAndOpacs[index-1].first;	x2=freqsAndOpacs[index].first;
		y1=freqsAndOpacs[index-1].second;	y2=freqsAndOpacs[index].second;
		if(index==0){
			x1=freqsAndOpacs[0].first;	x2=freqsAndOpacs[1].first;
			y1=freqsAndOpacs[0].second;	y2=freqsAndOpacs[1].second;
		}
		if(index==freqsAndOpacs.size()-1){
			int lastindex=freqsAndOpacs.size()-1;
			x1=freqsAndOpacs[lastindex-1].first;	x2=freqsAndOpacs[lastindex].first;
			y1=freqsAndOpacs[lastindex-1].second;	y2=freqsAndOpacs[lastindex].second;
		}
		double outputOpac=y1+((y2-y1)/(x2-x1))*(inputfreq-x1);
		//std::cout << "dust opacity: " << outputOpac << " (for frequency " << inputfreq << ")\n";
		return outputOpac;
	}
	
//this one really doesn't benifit from being vectorized since it's literally just a 
//table lookup, but I need an AVX version for compatibility with the other stuff.
	Vec4d operator () (Vec4d freqs) const{
		if(freqsAndOpacs.size()==0){
			std::cout << "Can't interpolate on empty opacity list" << std::endl;
			exit(1);
		}
		
		Vec4d outputs(0);
		for(int i=0;i<4;i++){
			double frequency=freqs[i];
			int index=0;
			//std::cout << "input: " << inputfreq << std::endl;
			while(frequency > freqsAndOpacs[index].first){
				index++;
				if(index>=freqsAndOpacs.size())
					break;
			}
			double x1,x2,y1,y2;
			x1=freqsAndOpacs[index-1].first;	x2=freqsAndOpacs[index].first;
			y1=freqsAndOpacs[index-1].second;	y2=freqsAndOpacs[index].second;
			if(index==0){
				x1=freqsAndOpacs[0].first;	x2=freqsAndOpacs[1].first;
				y1=freqsAndOpacs[0].second;	y2=freqsAndOpacs[1].second;
			}
			if(index==freqsAndOpacs.size()-1){
				int lastindex=freqsAndOpacs.size()-1;
				x1=freqsAndOpacs[lastindex-1].first;	x2=freqsAndOpacs[lastindex].first;
				y1=freqsAndOpacs[lastindex-1].second;	y2=freqsAndOpacs[lastindex].second;
			}
			double outputOpac=y1+((y2-y1)/(x2-x1))*(frequency-x1);
			outputs.insert(i,outputOpac);
		}
		return outputs;
	}
};