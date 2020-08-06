/*
 *  diskPhysics.cpp
 *
 *  Created by Erik on 4/17/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "diskPhysics.h"

//blackbody distribution
double blackBody(const double temp, const double freq){
	double exponential= exp(h*freq/(kboltzmann*temp));
	exponential= 1.0/(exponential -1);
	exponential *= 2*h*freq*freq*freq/(c*c);
	return exponential;
}

Vec4d blackBodyAVX(const double temp, const Vec4d freq){
	Vec4d arg= freq*h/(kboltzmann*temp);
	Vec4d exponential= exp(arg);
	exponential-=1.0;
	exponential=pow(exponential,-1.0);
	exponential *= 2*h*freq*freq*freq/(c*c);
	return exponential;
}

//inverse of the blackbody
double tempFromSurfaceBrightness(const double intensity, const double freq){
	double intens = 2*h*freq*freq*freq/c/c/intensity;
	return h*freq/kboltzmann/log(intens+1);
}

//Doppler shift
double doppler(double frequency, double velocity){
	return frequency/(1+velocity/c);
}

Vec4d dopplerAVX(Vec4d frequency, double velocity){
	return frequency/(1+velocity/c);
}

//Planck Mean Opacity, approximated by using a table of dust opacities
double PMO(const double temp, const std::vector<double>& frequencies, const std::vector<double>& opacities){
	double range=frequencies[0];
	double intensity=range*blackBody(temp,frequencies[0])*opacities[0];
	for(int i=1;i<frequencies.size();i++){
		//std::cout << range << "\t" << range*blackBody(temp,frequencies[i]) << std::endl;
		range=frequencies[i]-frequencies[i-1];
		intensity+=range*blackBody(temp,frequencies[i])*opacities[i];
	}
	return intensity*pi/sigma/temp/temp/temp/temp;
}

//surface temperature based on star properties
double Tsurf(double r, const std::vector<double>& frequencies,const std::vector<double>& opacities){
	double tol=1e-6;
	double Tstart=300; //A reasonable-ish temperature we start the iterator at.
	double planckMean=PMO(Tstart,frequencies,opacities);
	double first=pow(planckMean,-0.25)*Tstar*sqrt(Rstar/2/r);
	double delta=1;
	double next;
	while(delta>tol){
		double newpmo=PMO(first,frequencies,opacities);
		next=pow(newpmo,-0.25)*Tstar*sqrt(Rstar/2/r);
		delta=std::abs(next-first);
		first=next;
		//std::cout << next << "\t" << delta << std::endl;
	}
	return next;
}

/*double sigma0FromMass(double mass, double P, double r0, double rmin, double rmax){
	return ((2-P)*mass/(2*pi*pow(r0,P)*(pow(rmax,2-P)-pow(rmin,2-P))) );
}*/

//power law with cutoff version of sigma 0. //rmin and rmax aren't used here, but I'm keeping them for compatibility
double sigma0FromMass(double mass, double P, double r0, double rmin, double rmax){
	return mass*(2-P)/(2*pi*r0*r0);
}

double gaussian(double x, double mu, double sigma){
	return (sqrt(1/(2*pi))/sigma * exp(-(x-mu)*(x-mu)/(2*sigma*sigma)));
}

double gaussianNotNorm(double x, double mu, double sigma){
	return (exp(-(x-mu)*(x-mu)/(2*sigma*sigma)));
}

double dustOpac(const double frequency){
	int index=0;
	while(frequency > freqvals[index]){
		//std::cout << index << "\t" << inputfreq << "\t" << freqvals[index] << std::endl;
		index++;
		if(index>=205)
			break;
	}
	double x1,x2,y1,y2;
	x1=freqvals[index-1];	x2=freqvals[index];
	y1=opacvals[index-1];	y2=opacvals[index];
	if(index==0){
		x1=freqvals[0];	x2=freqvals[1];
		y1=opacvals[0];	y2=opacvals[1];
	}
	if(index==205){
		x1=freqvals[204];	x2=freqvals[205];
		y1=opacvals[204];	y2=opacvals[205];
	}
	double outputOpac=y1+((y2-y1)/(x2-x1))*(frequency-x1);
	//multiply by 100 to correct for the dust to gas ratio
	//TODO: Make the dust/gas ratio an adjustable parameter
	outputOpac*=100;
	return outputOpac;
}

//this one really doesn't benifit from being vectorized since it's literally just a 
//table lookup, but I need an AVX version for compatibility with the other stuff.
Vec4d dustOpacAVX(const Vec4d freqs){
	Vec4d outputs(0);
	for(int i=0;i<4;i++){
		double frequency=freqs[i];
		int index=0;
		while(frequency > freqvals[index]){
			//std::cout << index << "\t" << inputfreq << "\t" << freqvals[index] << std::endl;
			index++;
			if(index>=205)
				break;
		}
		double x1,x2,y1,y2;
		x1=freqvals[index-1];	x2=freqvals[index];
		y1=opacvals[index-1];	y2=opacvals[index];
		if(index==0){
			x1=freqvals[0];	x2=freqvals[1];
			y1=opacvals[0];	y2=opacvals[1];
		}
		if(index==205){
			x1=freqvals[204];	x2=freqvals[205];
			y1=opacvals[204];	y2=opacvals[205];
		}
		double outputOpac=y1+((y2-y1)/(x2-x1))*(frequency-x1);
		//multiply by 100 to correct for the dust to gas ratio
		//TODO: Make the dust/gas ratio an adjustable parameter
		outputOpac*=100;
		outputs.insert(i,outputOpac);
	}

	return outputs;
}