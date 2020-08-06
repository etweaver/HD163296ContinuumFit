//modelEval.cpp
//given a set of parameters, make and evaluate a model disk

#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <utility>
#include <vector>
#include "vcl/vectorclass.h"
#include "vcl/vectormath_trig.h"

#include "diskPhysics.h"
#include "geometry.h"
#include "grid.h"
#include "image.h"
#include "ParameterSet.h"
#include "worker.h"

struct diskParams{
	double inc, PA, logSig0, rc, P, h0dust, Sdust, h0gas, Sgas, p1,p2,p3, d1, d2, d3, w1, w2, w3, mStar, deltaF, dx, dy;
	diskParams()=default;
	diskParams(double inc, double PA, double logSig0,double rc,double P,double h0dust,double Sdust,
	double h0gas,double Sgas, double p1, double p2, double p3, double d1,
	double d2,double d3,double w1,double w2,double w3,double mStar,double deltaF,double dx,double dy):
	inc(inc), PA(PA), logSig0(logSig0), rc(rc), P(P), h0dust(h0dust), Sdust(Sdust), 
	h0gas(h0gas), Sgas(Sgas), p1(p1), p2(p2), p3(p3), d1(d1), d2(d2), 
	d3(d3), w1(w1), w2(w2), w3(w3), mStar(mStar), deltaF(deltaF), dx(dx), dy(dy)	{	}
};

struct newDensity {
	double Sigma0;
	double rc;
	double h0;
	double P;//density index
	double S;//scale height index
	double p1, p2, p3, p4; //ring positions (au)
	double d1, d2, d3, d4; //ring depths (0 to 1)
	double w1, w2, w3, w4; //ring widths (au)

	newDensity(): Sigma0(0), rc(0), h0(0), P(0), S(0), p1(0), p2(0), p3(0), p4(0), d1(0), d2(0), d3(0), d4(0), w1(0), w2(0), w3(0), w4(0) { }
	newDensity(const newDensity& other): Sigma0(other.Sigma0), rc(other.rc), h0(other.h0),
	P(other.P), S(other.S), p1(other.p1), p2(other.p2), p3(other.p3), p4(other.p4), d1(other.d1), d2(other.d2), 
	d3(other.d3), d4(other.d4), w1(other.w1), w2(other.w2), w3(other.w3), w4(other.w4) {  }
	newDensity(double Sigma0, double rc, double h0, double P, double S, double p1, double p2, double p3, double p4,
	double d1, double d2, double d3, double d4, double w1, double w2, double w3, double w4):
	Sigma0(Sigma0), rc(rc), h0(h0), P(P), S(S), p1(p1), p2(p2), p3(p3), p4(p4), d1(d1), d2(d2), d3(d3), 
	d4(d4), w1(w1), w2(w2), w3(w3), w4(w4) {	}

	newDensity& operator= (const newDensity& other){
		Sigma0=other.Sigma0; rc=other.rc; h0=other.h0; P=other.P; S=other.S;
		p1=other.p1; p2=other.p2; p3=other.p3; p4=other.p4; d1=other.d1; d2=other.d2; 
		d3=other.d3; d4=other.d4; w1=other.w1; w2=other.w2; w3=other.w3; w4=other.w4;
		return *this;
	}

	double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/rc, -P)) * exp(-(pow(r/rc, 2-P)));
	}

	double scaleHeight(const double r) const{
		return (h0*pow(r/rc,S));
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

		return(((surfaceMassDensity(r_cyl))/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h))*ring1*ring2*ring3*ring4);
	}
};

struct uvTable{
	unsigned int nfreqs;
	double** u;
	double** v;
	double** real;
	double** imag;
	double** weight;
	std::vector<std::unique_ptr<double[]> > u_store;
	std::vector<std::unique_ptr<double[]> > v_store;
	std::vector<std::unique_ptr<double[]> > real_store;
	std::vector<std::unique_ptr<double[]> > imag_store;
	std::vector<std::unique_ptr<double[]> > weight_store;
	unsigned int size; //size per channel

	uvTable()=default;

	uvTable(uvTable&& other)=default;

	double* makeAligned(size_t vectorSize, double* unalignedAddress){
		const int alignBy = vectorSize*sizeof(double);
		double* alignedAddress = (double*)(((size_t)unalignedAddress + alignBy - 1) & (-alignBy));
		return alignedAddress;
	}

	void readFile(std::string inFileName, int nchans, unsigned int nLines){
		const size_t vectorSize=4;
		const size_t padding=vectorSize-1;
		nfreqs=nchans;
		u=new double*[nchans];
		v=new double*[nchans];
		real=new double*[nchans];
		imag=new double*[nchans];
		weight=new double*[nchans];
		u_store.resize(nchans);
		v_store.resize(nchans);
		real_store.resize(nchans);
		imag_store.resize(nchans);
		weight_store.resize(nchans);
		size=nLines;
		
		for(int f=0;f<nchans; f++){
			u_store[f].reset(new double[nLines+2*padding]);
			v_store[f].reset(new double[nLines+2*padding]);
			real_store[f].reset(new double[nLines+2*padding]);
			imag_store[f].reset(new double[nLines+2*padding]);
			weight_store[f].reset(new double[nLines+2*padding]);

			u[f]=makeAligned(vectorSize,u_store[f].get());
			v[f]=makeAligned(vectorSize,v_store[f].get());
			real[f]=makeAligned(vectorSize,real_store[f].get());
			imag[f]=makeAligned(vectorSize,imag_store[f].get());
			weight[f]=makeAligned(vectorSize,weight_store[f].get());

			std::ifstream infile(inFileName);
			std::string line;
			double tempu, tempv, tempReal, tempImag, tempWeight;
			for(int i=0; i<nLines; i++){
				std::getline(infile, line);
				std::stringstream stream;
				stream << line;
				stream >> tempu >> tempv >> tempReal >> tempImag >> tempWeight;
				u[f][i]=tempu;
				v[f][i]=tempv;
				real[f][i]=tempReal;
				imag[f][i]=tempImag;
				weight[f][i]=tempWeight;
				//std::cout << tempu << "\t" << tempv << "\t" << tempReal << "\t" << tempImag << "\t" << tempWeight << std::endl;
			}
			infile.close();
			if(nLines%vectorSize){
				for(int i=0;i<vectorSize-nLines%vectorSize;i++){
					u[f][i+nLines]=v[f][i+nLines]=real[f][i+nLines]=imag[f][i+nLines]=weight[f][i+nLines]=0;
				}
			}
		}
		//don't forget to adjust the size per channel based on the additional padding
		//I'm assuming that all channels have the same number of points
		if(nLines%vectorSize)
			size+=vectorSize-nLines%vectorSize;
		
	}
};

//given a final set of parameters, print the file to a fits file and make the uv tables
void printModel(const grid<newDensity>& g1, image& im, const double PA, const std::vector<uvPoint>& data, ThreadPool& pool){
	im.propagate(g1,-PA,grid<newDensity>::continuum, pool);
	
	double totalFlux1=0;
	for(int i=0; i<im.hpix; i++){
		for(int j=0; j<im.vpix; j++){
			totalFlux1+=im.data[0][i][j];
		}
	}
	//im.data[0][1400][850]=0.05; //indeces are data[f][y][x]
	std::cout << "total flux: " << totalFlux1 << " jy" << std::endl;
	im.printToFits("model.fits");
	fourierImage FFTs=FFTDifferent(im);
	//FFTs.printToFits("fft.fits");
	double result=0;
	
	//std::ofstream UVoutfile("uvtableSim.txt");
	std::ofstream UVoutfile1("uvtableData.txt");
	std::ofstream UVoutfile2("uvtableResid.txt");
	//UVoutfile1.precision(12);
	//UVoutfile2.precision(12);
	//double inc=(90-22.004)*(pi/180); double PA=(110.3-90)*pi/180;
	for(int i=0;i<data.size();i++){
		std::pair<double,double> interpPair = FFTs.interp(0,data[i].u, data[i].v);
		UVoutfile1 << /*data[i].u << "\t" << data[i].v << "\t" <<*/ interpPair.first << "\t" << interpPair.second << std::endl;
		double diffReal=data[i].real-interpPair.first;
		double diffImag=data[i].imaginary-interpPair.second;
		UVoutfile2 << diffReal << "\t" << diffImag << std::endl;
		//temporary: need to deproject the uv data.
		//double rotU=data[i].u*cos(PA)-data[i].v*sin(PA);
		//double rotV=data[i].u*sin(PA)+data[i].v*cos(PA);
		//rotV*=cos(inc);
		//UVoutfile << sqrt(rotU*rotU + rotV*rotV) << "\t" << sqrt(data[i].u*data[i].u + data[i].v*data[i].v) << "\t" << interpPair.first << "\t" << interpPair.second << std::endl;
		//UVoutfile << data[i].u << "\t" << data[i].v << "\t" << sqrt(data[i].u*data[i].u+data[i].v*data[i].v) << "\t" << interpPair.first << "\t" << interpPair.second << std::endl;
		result+=(diffReal*diffReal + diffImag*diffImag)*data[i].weight;
	}
	UVoutfile1.close();
	UVoutfile2.close();
}

double chiSquaredAVX(image& model, const uvTable& data, unsigned int index, const double dx, const double dy){
	fourierImage FFT=FFTDifferent(model);
	FFT.offset(dx,dy);
	double result=0;
	int totalNumber=data.size;

	//std::cout << "starting chi^2 " << index << "\t" << dx << "\t" << dy << std::endl;

	for(int i=0;i<data.size;i+=4){
		Vec4d uPointVec(0);
		Vec4d vPointVec(0);
		Vec4d realVec(0);
		Vec4d imagVec(0);
		Vec4d weightVec(0);
		uPointVec.load_a(data.u[index]+i);
		vPointVec.load_a(data.v[index]+i);
		realVec.load_a(data.real[index]+i);
		imagVec.load_a(data.imag[index]+i);
		weightVec.load_a(data.weight[index]+i);
			
		std::pair<Vec4d,Vec4d> interpPair = FFT.interpAVX(0,uPointVec, vPointVec);
		Vec4d diffReal=realVec-interpPair.first;
		Vec4d diffImag=imagVec-interpPair.second;

		Vec4d resultVec=(diffReal*diffReal + diffImag*diffImag)*weightVec;

		for(int i=0;i<resultVec.size();i++){
			//std::cout << uPointVec[i] << "\t" << vPointVec[i] << "\t" << realVec[i] << "\t" << imagVec[i] << "\t" << weightVec[i] << "\t\t" << diffReal[i] <<std::endl;
			result+=resultVec[i];
		}
	}

	return result;
}

class modelWorker : public Worker{
public:
	modelWorker(std::string serverAddress, uvTable&& data):
		Worker(serverAddress), data(std::move(data)),
		dens(-1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0),	//-1 is just a placeholder for things we'll put in once we get the parameters
		g(0.01*AU,400*AU,60*(pi/180),120*(pi/180),0, 6.28318530717, -1, false, true, "diskmodels/HD163296/amr_grid.inp",
			"diskmodels/HD163296/dust_temperature_phi0.ascii", "diskmodels/HD163296/dust_temperature_phi0.ascii",
			"diskmodels/HD163296/dustopac.txt", dens, false, 0),
		im(1500, 1500, 2000*AU, 2000*AU , origin, origin, 2.30538e+11, 4,{0,0,0,0}, 101*3.0857e+18)
	{
		//the grid has separate dust and gas density structures, but it didn't used to, and the constructors don't yet reflect this
		//now the two structures are nearly the same, so we only need to change a few things
		g.dustDens=g.dens;
		g.dustDens.h0=-1; g.dustDens.S=-1;
		im.RA = 246.5987125; im.DEC = -24.720525; //these also need to go in the constructor
		vect pos(-1,0,-1);
        	im.position=pos;
	}
			
protected:
	virtual double compute(const std::vector<double>& params) override {
		ThreadPool pool(1); //this will be taken out once things are working

		//first we need to extract the full set of disk parameters from the input
		unsigned int chainIndex, fIndex;
		double inc, PA, logSig0, rc, P, h0dust, Sdust, h0gas, Sgas, p1,p2,p3, d1, d2, d3, w1, w2, w3, mStar, deltaF, dx, dy;
		//some are hardcoded for now:
		//inc=-0.768; PA=2.307;
        //	p1=48.23*AU; p2=85.37*AU; p3=98.89*AU;
		chainIndex=(unsigned int)params[0]; fIndex=(unsigned int)params[1];
		inc=-params[2]; PA=params[3];
		logSig0=params[4]; rc=params[5]*AU; P=params[6]; h0dust=params[7]*AU; Sdust=params[8];
		p1=params[9]*AU; p2=params[10]*AU; p3=params[11]*AU; d1=params[12]; d2=params[13]; d3=params[14]; 
		w1=params[15]*AU; w2=params[16]*AU; w3=params[17]*AU;
		
		//now we need to reset the disk parameters with these values
		dens.Sigma0=pow(10,logSig0); dens.rc=rc; dens.P=P; dens.h0=h0dust; dens.S=Sdust; dens.p1=p1; dens.p2=p2; dens.p3=p3;
		dens.w1=w1; dens.w2=w2; dens.w3=w3; dens.d1=d1; dens.d2=d2; dens.d3=d3;
		g.dustDens=dens; //dust structure
		dens.h0=h0gas; dens.S=Sgas;
		g.dens=dens;	//gas structure
		g.starMass=2.43*mSun;
		vect pos(10000*AU*cos(inc),0,10000*AU*sin(inc));
		im.position=pos;
		//std::cout << "disk set up" << std::endl;

		//now we need to set up the frequencies
		unsigned int nFreqs=1;
		double deltaV[1]={0};
		//double deltaV[6]={6.8,6.4,6,5.6,5.2,0};
		/*std::vector<double> globalFrequencies;
		for(int i=0; i< nFreqs;i++){
                	double frequency=im.centfreq+(deltaV[i]*1e5/c*im.centfreq);
                	globalFrequencies.push_back(frequency);
        	}
		std::vector<double> frequencies;
		double startfreq=globalFrequencies[fIndex];
		double freqStep=2e5/4;
		//double freqStep=0;
		frequencies.push_back(startfreq-(freqStep*1.5));
		frequencies.push_back(startfreq-(freqStep*0.5));
		frequencies.push_back(startfreq+(freqStep*0.5));
		frequencies.push_back(startfreq+(freqStep*1.5));
		std::cout.precision(10);
		for(int i=0;i<frequencies.size();i++){ 
			frequencies[i]+=deltaF;
			//std::cout << frequencies[i] << std::endl;
		}*/
		im.frequencies={2.314936e11};
		im.freqbins=1;
		//std::cout << "frequencies set up" << std::endl;

		im.propagate(g,-PA,grid<newDensity>::continuum, pool);
		//std::cout << "propagation done" << std::endl;
		//im.printToFits("junk.fits");

		//now we need to flatten the image along frequency
		/*std::vector<double> finalfreq; finalfreq.push_back(im.centfreq);
		image finalIm(im.vpix, im.hpix, im.width, im.height , origin, origin, im.centfreq, 1, finalfreq, im.distance);
		for(int i=0; i<im.hpix; i++){
			for(int j=0; j<im.vpix; j++){
				finalIm.data[0][j][i]=0;
				for(int f=0; f<4; f++){
					finalIm.data[0][j][i]+=im.data[f][j][i]/4;
				}
			}
		}*/
		//im.printToFits("junk.fits");
		//std::cout << "frequency average done" << std::endl;
		chi2=chiSquaredAVX(im, data,fIndex,0,0);
		std::cout << "chi squared: " << chi2 << std::endl;
		return chi2;
	}

private:
	newDensity dens;
	grid<newDensity> g;
	image im, finalImg;
	uvTable data;
	double chi2;

};

int main(int argc, char* argv[]){
	uvTable data;
	data.readFile("diskmodels/HD163296_cont_final_spav_noflags_nohighfreq.txt.cor",1,1549897);
	std::cout << "uvtable read in" << std::endl;
	
	std::string serverAddress="palladium-ii.cnweaver.staticcling.org:14000";
	modelWorker testWorker(serverAddress, std::move(data));
	testWorker.run();

	/*image residual=im;
	image realFrame=fitsExtract("DoAr25frame.fits",frequencies,dist);
	for(int i=0;i<realFrame.hpix;i++){
	for(int j=0; j<realFrame.vpix;j++){
	residual.data[0][j][i]=realFrame.data[0][j][i]-im.data[0][j][i]*8178.6086/4;
	}
	}
	residual.printToFits("res.fits");
	*/
	//std::cout << inc << "\t" << PA << "\t" << logSig0 << "\t" << rc << "\t" << P << "\t" << h0 << "\t" << S << "\t" << p1 << "\t" 
	//	<< p2 << "\t" << p3 << "\t" << d1 << "\t" << d2 << "\t" << d3 << "\t" << w1 << "\t" << w2 << "\t" << w3 << std::endl;
	//std::cout << p3 << "\t" << d3 << "\t" << w3 << std::endl;
	
	return 0;
}
