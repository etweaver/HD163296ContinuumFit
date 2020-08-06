/*
 *  image.h
 *  
 *  Created by Erik on 3/19/15.
 */

//the image that will be displayed
//each pixel's value will be the result of a ray that goes out, normal
//to the surface, which propagates through the disk
//The disk is always centered at the origin.

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <chrono>

#include "ThreadPool.h"
#include "fitsio.h"
#include "fftw3.h"
#include "marray.h"

#include "vcl/vectorclass.h"
#include "vcl/vectormath_trig.h"
#include "aligned_alloc.h"

//print cfitsio error message and exit
void fitserror(int status){
	if	(status){
		fits_report_error(stderr,	status);	/*	print	error	report	*/
		exit(status);	/*	terminate the program, returning error status	*/
	}
	return;
}

class image{
public:
	unsigned int vpix, hpix; //vertical and horizontal pixel numbers
	double width, height; //centimeters
	vect position, target; //position and target of image normal
	double centfreq;
	unsigned int freqbins;
	std::vector<double> frequencies;
	double rayTraceTime;
	marray<double,3> data;
	double RA, DEC;
	double distance; //the physical distance the object is away, for determining the pixel size. Not necessarily the distance away the camera is.
	
	//the usual constructors/destructors
	image():vpix(0), hpix(0), width(0), height(0), position(), target(), centfreq(0),freqbins(0), distance(0), frequencies(), data() {}
	image(unsigned int vpix, unsigned int hpix,double w, double h , vect p, vect t, double centfreq, unsigned int freqbins, std::vector<double> frequencies, double distance):
		vpix(vpix),hpix(hpix),width(w),height(h), position(p), target(t), centfreq(centfreq), freqbins(freqbins), distance(distance), frequencies(frequencies), data({freqbins,hpix,vpix}){}
	image(const image& other):
		vpix(other.vpix),hpix(other.hpix), width(other.width),height(other.height),position(other.position),target(other.target),
		centfreq(other.centfreq), freqbins(other.freqbins), distance(other.distance), frequencies(other.frequencies), data(other.data){}
	//~image(){ delete[] data; }
	image& operator=(const image& other){
		if(this!=&other){
			vpix=other.vpix;
			hpix=other.hpix;
			width=other.width;
			height=other.height;
			position=other.position;
			target=other.target;
			centfreq=other.centfreq;
			freqbins=other.freqbins;
			distance=other.distance;
			RA=other.RA;
			DEC=other.DEC;
			frequencies=other.frequencies;
			data=other.data;
		}
		return(*this);
	}
	
	//actually determine the values of each pixel of an image, for a given disk
	//for now, we assume parallel rays, so the center of each pixel casts a
	//ray through the disk
	template<class density>
	void propagate(const grid<density>& g, const double rotAngle, typename grid<density>::prop_type type, ThreadPool& pool){
		//Galario setup stuff
		//for Marco's FFT system, we need the center of the image to be at pixel hpix/2,vpix/2
		//I'm not quite convinced this is necessary, given that the ray tracing will probably 
		//always take more time than the FFT, but I need to make sure this works.
		//So, if the image has an even number of pixels in either dimension, it needs to be padded 
		//with an extra row/column. This will then be stripped off at the end.
		//If the size dimensions are odd, then we just need to strip off the edges at the end.
		//And suitably warn the user to pad the image.
		unsigned int hpixSpecial, vpixSpecial;
		hpixSpecial = (hpix%2) ? hpix : hpix+1; //if even, add another column
		vpixSpecial = (vpix%2) ? vpix : vpix+1; 
		
		unsigned int hpixFinal, vpixFinal;
		hpixFinal = (hpix%2) ? hpix-1 : hpix; //if odd, cut the last column, for the final size
		vpixFinal = (vpix%2) ? vpix-1 : vpix; 
		
		//the width and height need to be adjusted also
		double widthSpecial, heightSpecial;
		widthSpecial = (hpix%2) ? width : width*(hpixSpecial/hpix); //if we added a column, extend the width to match that.
		heightSpecial = (vpix%2) ? height : height*(vpixSpecial/vpix);
		
		double widthFinal, heightFinal;
		widthFinal = (hpix%2) ? width : width*(hpixSpecial/hpix); //if we added a column, extend the width to match that.
		heightFinal = (vpix%2) ? height : height*(vpixSpecial/vpix);
		
		/*std::vector<double> frequencies;
		double freq=centfreq;
		freq-=freqrange/2;
		for(int i=0;i<freqbins;i++){
			freq+=(freqrange/freqbins);
			frequencies.push_back(freq);
			//std::cout << freq << std::endl;
		}
		assert(frequencies.size() > 0);*/
		
		//now we need to repackage the frequencies into groups of four
		unsigned int nfreqsAVX=frequencies.size()/4;
		unsigned int extra=frequencies.size()%4;
		std::vector<Vec4d,aligned_allocator<Vec4d> > freqsAVX (aligned_allocator<Vec4d>(5));
		//std::cout << frequencies.size() << "\t" << nfreqsAVX << "\t" << extra << "\t" << frequencies[29] << std::endl;
		for(int f=0;f<nfreqsAVX;f++){
			Vec4d batch(frequencies[f*4],frequencies[f*4+1],frequencies[f*4+2],frequencies[f*4+3]);
			//Vec4d batch(0);
			freqsAVX.push_back(batch);
		}
		Vec4d stragglers(-1.0);//the last couple left over, padded with an obvious error value.
		for(int f=nfreqsAVX*4;f<frequencies.size();f++){
			stragglers.insert(f-nfreqsAVX*4,frequencies[f]);
		}
		freqsAVX.push_back(stragglers);
		
		std::chrono::high_resolution_clock::time_point start, finish;
		start = std::chrono::high_resolution_clock::now();
		
		//float* longdataOld;
		marray<double,3> longdata({frequencies.size(),hpixFinal,vpixFinal});
		//std::cout << hpixSpecial << "\t" << vpixSpecial << "\t" << hpixFinal << "\t" << vpixFinal << std::endl;
		//std::cout << "Allocating " << frequencies.size()*hpixFinal*vpixFinal << " pixels for final image" << std::endl;
		//TODO: Change back to double when done with CASA
		//longdataOld=new float[hpix*vpix*frequencies.size()];
		
		//this will be done on a pixel by pixel basis, for now
		//but first we calculate the normal vector for the pixels, since it will be the same
		//for all of them.
		vect normal= target-position;
		if(normal.mag()==0){
			std::cout << "Error in image propagation: Normal vector is 0." << std::endl;
			std::cout << "position and target vectors must not be the same." << std:: endl;
			exit(1);
		}
		//g.propagateRayAVX(line(target,normal),frequencies,position,type); //just to print diagnostic info.
		
		const double radDist= position.r(); //radial distance of the image center from disk center
		const double theta= position.theta(); //angle between image center and disk verticle (z axis)
		const double phi= position.phi();	//angle between image center and y axis
		//now, this is flux coming out of the pixel. We want the intensity viewed by the observer
		//To get this, we need to know how far away the observer is, and the area of each pixel.
		double pixelArea=(heightSpecial*widthSpecial)/(vpixSpecial*hpixSpecial); //area per pixel
		
		//currently, the disk position is always the origin. At some point in the future, this may change,
		//and we'll need to change the distance to (position-g.pos).mag(). For now, this works.
		double realDist=position.mag();
		
		double cgsToJansky=1e23;
		double solidAngle = pixelArea/(distance*distance);
		//double cgsToJansky=1;
		//double solidAngle = 1;
		//std::cout << "solid angle: " << solidAngle << std::endl;
		int numpix=0;//temporary. Number of pixels that intersect the disk
		
		std::vector<std::future<void>> results;
		std::mutex storeMut;
		std::vector<double> totalFluxes(frequencies.size(),0.0); //total flux at each frequency, in units of janskys
		
		for(int i=0;i<hpixFinal;i++){
			results.emplace_back(pool.enqueue([&,i](){
				std::vector<std::vector<double>> values;
				values.reserve(vpixFinal);
				for(int j=0;j<vpixFinal;j++){
					//we begin with the array of pixels in the xy plane, and calculate the
					//"displacement vector" from the center of the image to the pixel
					//The ray for each pixel used to come from the center. I'm changing it to be the lower left corner.
					//This is so that the image center will actually be the center, not offset by half a pixel.
					//We need this so that the FFT can actually use the symmetry of the transforms properly later.
					double xpos=(i+0.5)*(widthSpecial/hpixSpecial)-widthSpecial/2.0;	
					double ypos=(j+0.5)*(heightSpecial/vpixSpecial)-heightSpecial/2.0;
					vect displacement(xpos,ypos,0);
					//this is the euler angle stuff. Don't touch it
					displacement.rotateZ(rotAngle);
					displacement.rotateY(theta);
					displacement.rotateZ(phi);
					vect radial(0,0,radDist);
					radial.rotateY(theta);
					radial.rotateZ(phi);
					displacement+=radial;
					
					vect normal=-radial;
					line l(displacement, normal);
					std::vector<double> pix=g.propagateRayAVX(l,freqsAVX,position,type);
					//std::vector<double> pix(frequencies.size(),0);
					//std::cout << "Not configured for images right now. Fix return type of grid::propagate" << std::endl;
					values.push_back(pix);//list of values at each freqency for a row
					bool all0=true;
					for(auto p : pix){
						if(p!=0.0)
							all0=false;
					}
					if(!all0)
						numpix++;
				}
				std::unique_lock<std::mutex> lock(storeMut);
				for(int j=0;j<vpixFinal;j++){
					for(int k=0;k<frequencies.size();k++){
						longdata[k][j][i]=values[j][k]*solidAngle*cgsToJansky;
					}
				}
				for(int k=0;k<frequencies.size();k++){
					double rowval=0;
					for(int j=0;j<vpixFinal;j++){
						rowval+=values[j][k]*solidAngle*cgsToJansky;
					}
					totalFluxes[k]+=rowval;
				}
			}));
		}
		for(auto& r : results) //wait for all jobs to finish
			r.get();
		
		//std::cout << "starting averaging" << std::endl;
		/*int centfreqnum=freqbins/2;
		double binwidth=freqrange/freqbins;
		//optional spectral averaging.
		//TODO: Make this an adjustable parameter
		int numToStack=10; //Average by stacking every n images
		unsigned int numFrequenciesAveraged=ceil((double) freqbins/numToStack);
		marray<double,3> longdataAveraged({numFrequenciesAveraged,hpixFinal,vpixFinal});
		//std::cout << frequencies.size() << "\t" << numFrequenciesAveraged << std::endl;
		//std::cout << "longdata size: " << hpix*vpix*frequencies.size() << "\tlongdataAveraged size: " << hpix*vpix*numFrequenciesAveraged << std::endl;
		for(int i=0; i<hpix; i++){
			for(int j=0; j<vpix; j++){
		    	for(int kk=0; kk<numFrequenciesAveraged; kk++){
		        	double value=0;
		            for(int k=0; k<numToStack && kk*numToStack+k<frequencies.size(); k++)
		                value+=longdata[kk*numToStack+k][i][j];
		            value/=numToStack;
					longdataAveraged[kk][i][j]=value;
				}
			}
		}
		longdata=longdataAveraged;
		//std::cout << "made longdataAveraged" << std::endl;
		freqbins=numFrequenciesAveraged;
		centfreqnum=numFrequenciesAveraged/2;
		binwidth=freqrange/numFrequenciesAveraged;*/
		
		finish = std::chrono::high_resolution_clock::now();
		rayTraceTime=std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
		
		//clean up
		//std::cout << numpix << " pixels intersect the disk" << std::endl; //temporary
		//if we stripped a row/column, we need to adjust the width and height of the image.
		width=widthFinal;
		height=heightFinal;
		hpix=hpixFinal;
		vpix=vpixFinal;
		//std::cout << "Warning: Width and height have been adjusted for Galario compatibility." << std::endl;
		//delete[] longdataAveraged;
		//std::cout << "all done " << std::endl;
		
		data=longdata;
		//return longdata;
	}
	
	void printToFits(std::string outFileName){
		//RA=239.1666917;
		//DEC=-22.027885833;
		std::cout << "printing file " << outFileName << std::endl;
		//Write fits file
		fitsfile* fptr;
		int status = 0;	//ye olde C style error code
		long naxis = 4;
		long nelements;
		int centfreqnum=freqbins/2 + 1;
		double freqrange=frequencies.back()-frequencies.front();
		double binwidth=freqrange/frequencies.size();
		//int naxis4=freqbins;
		/*int centfreqnum=freqbins/2+1;
		double binwidth=freqrange/freqbins;
		int naxis4=freqbins;*/
		int naxis4=(int)data.extent(0);
		long fpixel[4] = {1,1,1,1};
		//long naxes[4] = {hpix, vpix, (int)frequencies.size(),1};
		long naxes[4] = {vpix, hpix, naxis4,1};
		
		std::string fitsFileName(outFileName);
		fitsFileName.erase(outFileName.rfind('.'));
		fitsFileName.append(".fits");
		fitsFileName.insert(fitsFileName.begin(),'!'); //this is a safety feature that will make cfitsio overwrite an existing file of the same name
		if (fits_create_file(&fptr, fitsFileName.c_str(), &status))  
			fitserror( status );
		//There is an exclamation point because that tells cfitsio to overwrite the existing file
		//TODO: Set back to DOUBLE_IMG when doen with CASA
		if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) // create new FITS file 
			fitserror( status );           // call printerror if error occurs
		nelements=naxes[0]*naxes[1]*naxes[2];
					
		char bunit[10] = "JY/PIXEL";
		char btype[10] = "INTENSITY";
		fits_write_key(fptr , TSTRING , "BUNIT"    , &bunit   , "units of flux"   , &status);
		fits_write_key(fptr , TSTRING , "BTYPE"    , &btype   , "units of flux"   , &status);
		
		double temp=0;
		char RADESYS[10] = "FK5";
		fits_write_key(fptr , TSTRING , "RADESYS", &RADESYS, "", &status);
		double LONPOLE=180;
		fits_write_key(fptr , TDOUBLE , "LONPOLE", &LONPOLE, "", &status);
		double LATPOLE=-2.195626763250E+01;
		fits_write_key(fptr , TDOUBLE , "LATPOLE", &LATPOLE, "", &status);
		double EQUINOX = 2000.0;
		fits_write_key(fptr , TDOUBLE , "EQUINOX"    , &EQUINOX   , ""   , &status);
		char ctype1[10] = "RA---SIN";
		fits_write_key(fptr , TSTRING , "CTYPE1"    , &ctype1   , "axis 1 coordinate type"   , &status);
		//std::cout << RA << "\t" << DEC << std::endl;
		fits_write_key(fptr , TDOUBLE , "CRVAL1", &RA, "value at reference pixel", &status);
		double degperpixh = -(atan(width/(distance))/hpix) * (180/pi);
		fits_write_key(fptr , TDOUBLE , "CDELT1", &degperpixh, "increment per pixel (deg)", &status);
		int refpix1=hpix/2 + 1;
		fits_write_key(fptr , TINT , "CRPIX1", &refpix1, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA1", &temp, "rotation angle", &status);
		char ctype2[10] = "DEC--SIN";
		fits_write_key(fptr , TSTRING , "CTYPE2"    , &ctype2   , "axis 2 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL2", &DEC, "value at reference pixel", &status);
		double degperpixv = (atan(height/(distance))/vpix)* (180/pi);
		fits_write_key(fptr , TDOUBLE , "CDELT2", &degperpixv, "increment per pixel (deg)", &status);
		int refpix2=vpix/2 + 1;
		fits_write_key(fptr , TINT , "CRPIX2", &refpix2, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA2", &temp, "rotation angle", &status);
		char ctype3[10] = "FREQ-LSR";
		fits_write_key(fptr , TSTRING , "CTYPE3"    , &ctype3   , "axis 3 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL3", &centfreq, "value at reference pixel", &status);
		if(binwidth==0)
			binwidth=1e9;//set it to 1GHz if we just took a single frame.
		fits_write_key(fptr , TDOUBLE , "CDELT3", &binwidth, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX3", &centfreqnum, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA3", &temp, "rotation angle", &status);
				
		char cunit[10] = "HZ";
		char ctype[10] = "STOKES";
		fits_write_key(fptr , TSTRING , "CUNIT3"    , &cunit   , "units"   , &status);
		fits_write_key(fptr , TSTRING , "CTYPE4"    , &ctype   , "axis 3 coordinate type"   , &status);
		double temp1=1; double temp2=0; int temp3=1;
		fits_write_key(fptr , TDOUBLE , "CRVAL4", &temp1, "value at reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CDELT4", &temp1, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX4", &temp3, "references pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA4", &temp2, "rotation angle", &status);
		
		int velref=257;
		char specsys[10] = "LSRK";
		fits_write_key(fptr , TDOUBLE , "RESTFREQ", &centfreq, "Rest Frequency (Hz)", &status);
		fits_write_key(fptr , TSTRING , "SPECSYS", &specsys, "Spectral reference frame", &status);
		fits_write_key(fptr , TINT , "VELREF", &velref, "1 LSR, 2 HEL, 3 OBS, +256 Radio", &status);
		
		fits_write_key(fptr , TDOUBLE , "RA", &RA, "Right Ascention, J2000", &status);
		fits_write_key(fptr , TDOUBLE , "DEC", &DEC, "Declination, J2000", &status);
		
		//temporary: gotta flip the actual data around the vertical
		//indeces are data[f][y][x]
		marray<double,3> tempData=data;
		for(int f=0;f<data.extent(0);f++){
			for(int i=1;i<vpix;i++){
				for(int j=0;j<hpix;j++){
					tempData[f][i][j]=data[f][vpix-i][j];
				}
			}
		}
		
		//TODO: set back to TDOUBLE when done with CASA
		if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, &tempData.front(), &status))  
			fitserror( status );
		if (fits_close_file(fptr, &status))  
			fitserror( status );
	}
	
	void readInputFile(std::string filename){
		std::cout << "Reading file " << filename << std::endl;
		std::string line;
		std::ifstream infile(filename.c_str());
		int linenum=0;
		std::string dummy;
		while(!infile.eof()){
			std::getline(infile, line);
			if(line[0]=='#')
				continue;
			std::stringstream stream;
			stream << line;
			switch (linenum) {

				case 13:
					stream >> dummy >> dummy >> height;
					if(height <= 0.0){
						std::cout << "Error reading input file" << std::endl << "height <= 0 AU" << std::endl;
						exit(1);
					}
					height*=AU;
					break;
				case 14:
					stream >> dummy >> dummy >> width;
					if(width <= 0.0){
						std::cout << "Error reading input file" << std::endl << "width <= 0 AU" << std::endl;
						exit(1);
					}
					width*=AU;
					break;
				case 15:
					stream >> dummy >> dummy >> vpix;
					if(vpix < 1){
						std::cout << "Error reading input file" << std::endl << "number of vertical pixels < 1" << std::endl;
						exit(1);
					}
					break;
				case 16:
					stream >> dummy >> dummy >> hpix;
					if(hpix < 1){
						std::cout << "Error reading input file" << std::endl << "number of horizontal pixels < 1" << std::endl;
						exit(1);
					}
					break;
				case 17:
					stream >> dummy >> dummy >> centfreq;
					if(centfreq <= 0.0){
						std::cout << "Error reading input file" << std::endl << "bad central frequency" << std::endl;
						exit(1);
					}
					break;
				case 18:
					stream >> dummy >> dummy >> freqbins;
					if(freqbins < 1){
						std::cout << "Error reading input file" << std::endl << "number of frequency bins < 1" << std::endl;
						exit(1);
					}
					break;
					
				default:
					break;
			}
			linenum++;
		}
		//delete[] data;
		data=marray<double,3>({freqbins,hpix,vpix});
		std::cout << "image size " << hpix << "x" << vpix << std::endl;
		std::cout << "image physical size " << height << "x" <<  width << std::endl;
		std::cout << "image position " << position << std::endl;
		std::cout << "image target " << target << std::endl;
		std::cout << "central frequency " << centfreq << std::endl;
		std::cout << "frequency bins " << freqbins << std::endl;
	}
	
	//special version of propagate() for the case of finding/characterizing the matching temperature point
	//and emission region. Instead of grid::propagate returning a vector of intensities at each frequency,
	//this version works at only a single frequency, and returns a vector of width infos, of the form:
	//[5%pt,33%pt,matchpt,67%pt,95%pt]
	//no fits file written, just a text file for plotting with gnuplot
	/*void propagateSpecial(grid& g, double rotAngle, double frequency, std::string outFileName){
		
		std::chrono::high_resolution_clock::time_point start, finish;
		start = std::chrono::high_resolution_clock::now();
		
		grid::widthInfo* longdata;
		longdata=new grid::widthInfo[hpix*vpix*5]; //propagate returns a vector of 5 width infos.
		const size_t nThreadsMax=std::thread::hardware_concurrency();
		const size_t nThreads=(nThreadsMax>15 ? 15 : nThreadsMax);
		//const size_t nThreads=1;
		std::cout << "Using " << nThreads << " threads" << std::endl;
		
		//this will be done on a pixel by pixel basis, for now
		//but first we calculate the normal vector for the pixels, since it will be the same
		//for all of them.
		vect normal= target-position;
		if(normal.mag()==0){
			std::cout << "Error in image propagation: Normal vector is 0." << std::endl;
			std::cout << "position and target vectors must not be the same." << std:: endl;
			exit(1);
		}
		const double radDist= position.r(); //radial distance of the image center from disk center
		const double theta= position.theta(); //angle between image center and disk verticle (z axis)
		const double phi= position.phi();	//angle between image center and y axis
		
		ThreadPool pool(nThreads);
		std::vector<std::future<void>> results;
		std::mutex storeMut;
		
		for(int i=0;i<hpix;i++){
			results.emplace_back(pool.enqueue([&,i](){
				std::vector<std::vector<grid::widthInfo> > values;
				values.reserve(vpix);
				for(int j=0;j<vpix;j++){
					//we begin with the array of pixels in the xy plane, and calculate the
					//"displacement vector" from the center of the image to the pixel
					double xpos=(i+0.5)*(width/hpix)-width/2.0;	//x position of the center of pixel i,j
					double ypos=(j+0.5)*(height/vpix)-height/2.0;
					vect displacement(xpos,ypos,0);
					//Don't touch the Euler angle stuff here
					displacement.rotateZ(rotAngle);
					displacement.rotateY(theta);
					displacement.rotateZ(phi);
					vect radial(0,0,radDist);
					radial.rotateY(theta);
					radial.rotateZ(phi);
					displacement+=radial;
					
					vect normal=-radial;
					line l(displacement, normal);
					std::vector<double> frequencies;
					frequencies.push_back(frequency);
					std::vector<grid::widthInfo> pix=g.propagateRaySpecial(l,frequency,position);
					values.push_back(pix);
				}
				std::unique_lock<std::mutex> lock(storeMut);
				for(int j=0;j<vpix;j++){
 					//(*this)[i][j]=0;//just set this to 0 
					for(int k=0;k<5;k++){
						longdata[k*hpix*vpix + j*vpix + i]=values[j][k];
					}
				}
			}));
		}
		for(auto& r : results) //wait for all jobs to finish
			r.get();
		
		//This prints out the most recent image to something you can plot
		std::ofstream outfile(outFileName.c_str());
		outfile.precision(12);
		for(int i=0;i<hpix;i++){
			for(int j=0;j<vpix;j++){
				double xpos=(i+0.5)*(width/hpix)-width/2.0;
				double ypos=(j+0.5)*(height/vpix)-height/2.0;
				double width=(longdata[0*hpix*vpix + j*vpix + i].pos - longdata[4*hpix*vpix + j*vpix + i].pos).mag();
				vect centroid=(longdata[0*hpix*vpix + j*vpix + i].pos + longdata[4*hpix*vpix + j*vpix + i].pos)/2/AU;
				double intens=longdata[4*hpix*vpix + j*vpix + i].intensity;
				if(intens > 4e-13){
					outfile << xpos/AU << "\t" << ypos/AU << "\t" << width/AU << "\t" << centroid.x << "\t" << centroid.y << "\t" << centroid.z << "\t" << intens << "\n";
				}
			}
		}
		outfile.close();
		
		
		//clean up
		delete[] longdata;
	}*/

	//to reduce fourier noise, we need to be able to pad the outsides of the image with zeroes, which the then
	//crop out again after the FFTs.
	//padnum is an integer number of pixels to pad each side with
	image pad(unsigned int padNum){
		unsigned int newhpix=hpix+2*padNum;
		unsigned int newvpix=vpix+2*padNum;
		double widthScaleFactor=(double)newhpix/(double)hpix;
		double heightScaleFactor=(double)newvpix/(double)vpix;
		double newWidth=width*widthScaleFactor;
		double newHeight=height*heightScaleFactor;
		image result(newvpix, newhpix, width, height, position, target, centfreq, freqbins, frequencies, distance);
		
		std::cout << newhpix << "\t" << newvpix << std::endl;
		//the padded image is now set up, but now we need to copy over the original image data,
		//so that it is properly in the middle.
		//the image.data constructor fills in the entire area with zeroes, so all we need to do
		//is copy over the real data
		for(int f=0; f<freqbins; f++){
			for(int i=0;i<newhpix;i++){
				for(int j=0;j<newvpix;j++){
					int iprime=i-padNum;
					int jprime=j-padNum;
					//std::cout << i << "\t" << iprime << "\t" << j << "\t" << jprime << "\t";
					if(iprime>=0 && iprime<newhpix-padNum*2){
						if(jprime>=0 && jprime<newvpix-padNum*2){
							//std::cout << "check";
							result.data[f][j][i]=data[f][jprime][iprime];
						}
					}
					//std::cout << std::endl;
				}
				//std::cout << std::endl;
			}
		}
		return result;
	}
	
};

image fitsExtract(std::string filename, std::vector<double> frequencies){
	//first read out the image properties from the header
	int status = 0;	//ye olde C style error code
	int nkeys;
	fitsfile* fptr;
	char card[FLEN_CARD];
	
	fits_open_file(&fptr, filename.c_str(), READONLY, &status);
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	unsigned int xpix, ypix, nfreq;

	//read in x,y,and freq axis lengths.
	//technically, there's polarization too, but I have never needed it ever.
	fits_read_key(fptr,TINT, "naxis1", &xpix,NULL,&status);
	fits_read_key(fptr,TINT, "naxis2", &ypix,NULL,&status);
	fits_read_key(fptr,TINT, "naxis3", &nfreq,NULL,&status);
	if(status)
		fits_report_error(stderr, status);
	//std::cout << xpix << "\t" << ypix << std::endl;
	
	double* data=new double[ypix*xpix*nfreq];
	long start[4]={1,1,1,1};
	fits_read_pix(fptr, TDOUBLE, start, ypix*xpix*nfreq, 0, data, NULL, &status);
	if(status)
		fits_report_error(stderr, status);
	
	//now everything is set up for the image except for the width, height, position, and frequencies
	int xrefpix, yrefpix, centfreqnum;
	fits_read_key(fptr,TINT, "crpix1", &xrefpix, NULL,&status);
	fits_read_key(fptr,TINT, "crpix2", &yrefpix, NULL,&status);
	fits_read_key(fptr,TINT, "crpix3", &centfreqnum, NULL,&status);
	double xrefval, yrefval, centfreq;
	fits_read_key(fptr,TDOUBLE, "crval1", &xrefval, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "crval2", &yrefval, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "crval3", &centfreq, NULL,&status);
	//std::cout << xrefval << "\t" << yrefval << std::endl;
	double xinc, yinc, binwidth;
	fits_read_key(fptr,TDOUBLE, "cdelt1", &xinc, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "cdelt2", &yinc, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "cdelt3", &binwidth, NULL,&status);
	//double width=xinc*xpix; //note that these are in units of degrees
	//double height=yinc*ypix;
	
	//the standard fits files measure position relative to some pixel, often the center
	//It's easier for me to shift that down to the corner
	double cornerRA=xrefval-((xrefpix)*xinc);
	double cornerdec=yrefval-((yrefpix)*yinc);

	//std::cout << xrefval << "\t" << xrefpix << "\t" << xinc << std::endl;
	//std::cout << yrefval << "\t" << yrefpix << "\t" << yinc << "\t" << yrefval-(yrefpix)*yinc << std::endl;
	//std::cout << cornerRA << "\t" << cornerdec << std::endl;

	image img(ypix, xpix, 0, 0 , origin, origin, centfreq, nfreq, frequencies, 0);
	//cfitsio gives you an old c array, so we get to do pointer arithmatic
	std::copy(data,data+(ypix*xpix*nfreq),img.data.begin());
	img.RA=cornerRA; img.DEC=cornerdec;

	fits_close_file(fptr, &status);
	if(status)
		fits_report_error(stderr, status);
	
	return img;
}

//there are two versions because this one includes a distanc as an argument
//I can put the distance as a keyword in fits files I make, but there isn't much point
//since real data will never come with it. So, if you want to properly initialize the 
//disk, you need to read it in, and specify it's distance.
image fitsExtract(std::string filename, std::vector<double> frequencies, double distance){
	//first read out the image properties from the header
	int status = 0;	//ye olde C style error code
	int nkeys;
	fitsfile* fptr;
	char card[FLEN_CARD];
	
	fits_open_file(&fptr, filename.c_str(), READONLY, &status);
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	unsigned int xpix, ypix, nfreq;

	//read in x,y,and freq axis lengths.
	//technically, there's polarization too, but I have never needed it ever.
	fits_read_key(fptr,TINT, "naxis1", &xpix,NULL,&status);
	fits_read_key(fptr,TINT, "naxis2", &ypix,NULL,&status);
	fits_read_key(fptr,TINT, "naxis3", &nfreq,NULL,&status);
	if(status)
		fits_report_error(stderr, status);
	//std::cout << xpix << "\t" << ypix << std::endl;
	
	double* data=new double[ypix*xpix*nfreq];
	long start[4]={1,1,1,1};
	fits_read_pix(fptr, TDOUBLE, start, ypix*xpix*nfreq, 0, data, NULL, &status);
	if(status)
		fits_report_error(stderr, status);
	
	//now everything is set up for the image except for the width, height, position, and frequencies
	int xrefpix, yrefpix, centfreqnum;
	fits_read_key(fptr,TINT, "crpix1", &xrefpix, NULL,&status);
	fits_read_key(fptr,TINT, "crpix2", &yrefpix, NULL,&status);
	fits_read_key(fptr,TINT, "crpix3", &centfreqnum, NULL,&status);
	double xrefval, yrefval, centfreq;
	fits_read_key(fptr,TDOUBLE, "crval1", &xrefval, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "crval2", &yrefval, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "crval3", &centfreq, NULL,&status);
	//std::cout << xrefval << "\t" << yrefval << std::endl;
	double xinc, yinc, binwidth;
	fits_read_key(fptr,TDOUBLE, "cdelt1", &xinc, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "cdelt2", &yinc, NULL,&status);
	fits_read_key(fptr,TDOUBLE, "cdelt3", &binwidth, NULL,&status);
	double width=abs(tan(xinc*pi/180)*distance*xpix); 
	double height=abs(tan(yinc*pi/180)*distance*ypix);
	
	//the standard fits files measure position relative to some pixel, often the center
	//It's easier for me to shift that down to the corner
	double cornerRA=xrefval-((xrefpix)*xinc);
	double cornerdec=yrefval-((yrefpix)*yinc);

	//std::cout << xrefval << "\t" << xrefpix << "\t" << xinc << std::endl;
	//std::cout << yrefval << "\t" << yrefpix << "\t" << yinc << "\t" << yrefval-(yrefpix)*yinc << std::endl;
	//std::cout << cornerRA << "\t" << cornerdec << std::endl;
	
	image img(ypix, xpix, width, height, origin, origin, centfreq, nfreq, frequencies, distance);
	//cfitsio gives you an old c array, so we get to do pointer arithmatic
	std::copy(data,data+(ypix*xpix*nfreq),img.data.begin());
	img.RA=cornerRA; img.DEC=cornerdec; img.distance=distance;

	fits_close_file(fptr, &status);
	if(status)
		fits_report_error(stderr, status);
	
	return img;
}

//the FFTs are similar format to the image class, but are different enough (no sky coordinates, 
//centered differently, different units, etc) that they merit their own class.
//The fourierImage class contains the height, width, and pixel information, and both the real and imaginary
//parts of the FFT.
//It's worth mentioning that the real and imaginary parts are highly redundant in this form.
//Because the starting image is real, the FFT components will be symmetric, and so could be stored in half the
//space. Ideally, in the future this could be compressed down to the form that FFTW releases them in, but
//since I have to write all this in order to verify that everything is correct, we'll work with this easier 
//version for now.
//The loss in efficiency shouldn't be large. I think there's only going to be one FFT extent at any given time,
//and the time required to move memory and suchlike is insignificant compared to the ray tracing.
class fourierImage{
private:
	double xpixsize, ypixsize;
public:
	unsigned int vpix, hpix; //vertical and horizontal pixel numbers
	double width, height; //uv coordinates - units of wavelength
	double centfreq;
	unsigned int freqbins;
	std::vector<double> frequencies;
	double distance;
	
	marray<double,3> realPart;
	marray<double,3> imaginaryPart;
	
	fourierImage():vpix(0), hpix(0), width(0), height(0), centfreq(0), freqbins(0), distance(0), frequencies(), realPart(), imaginaryPart() {}
	
	fourierImage(unsigned int vpix, unsigned int hpix,double w, double h, double centfreq, unsigned int freqbins, std::vector<double> frequencies, double distance):
		vpix(vpix),hpix(hpix),width(w),height(h), centfreq(centfreq), freqbins(freqbins), distance(distance), frequencies(frequencies), realPart({freqbins,hpix,vpix}), imaginaryPart({freqbins,hpix,vpix}){
			xpixsize=(atan(width/distance)/hpix);
			ypixsize=(atan(height/distance)/vpix);
		}
	
	fourierImage(const fourierImage& other):
		vpix(other.vpix),hpix(other.hpix), width(other.width),height(other.height),centfreq(other.centfreq), freqbins(other.freqbins), distance(other.distance), 
		xpixsize(other.xpixsize), ypixsize(other.ypixsize), frequencies(other.frequencies), realPart(other.realPart), imaginaryPart(other.imaginaryPart){}
		
	fourierImage& operator=(const fourierImage& other){
		width=other.width; height=other.height;
		centfreq=other.centfreq; freqbins=other.freqbins;
		distance=other.distance;
		xpixsize=other.xpixsize;
		ypixsize=other.ypixsize;
		realPart=other.realPart; imaginaryPart=other.imaginaryPart;
		return(*this);
	}

	
	//write out both the real and imaginary parts to fits files
	void printToFits(std::string outFileName){
		fitsfile* fptr;
		int status = 0;	//ye olde C style error code
		long naxis = 4;
		long nelements;
		double freqrange=frequencies.back()-frequencies.front();
		int centfreqnum=frequencies.size()/2 + 1;
		double binwidth=freqrange/frequencies.size();
		int naxis4=freqbins;
		long fpixel[4] = {1,1,1,1};
		long naxes[4] = {vpix, hpix, naxis4,1};
		
		std::string fitsFileNameReal(outFileName);
		std::string fitsFileNameIm(outFileName);
		fitsFileNameReal.erase(outFileName.rfind('.'));
		fitsFileNameIm.erase(outFileName.rfind('.'));
		fitsFileNameReal.append("Real.fits");
		fitsFileNameIm.append("Im.fits");
		std::cout << "printing file " << fitsFileNameReal << std::endl;
		fitsFileNameReal.insert(fitsFileNameReal.begin(),'!'); //this is a safety feature that will make cfitsio overwrite an existing file of the same name
		fitsFileNameIm.insert(fitsFileNameIm.begin(),'!');
		if (fits_create_file(&fptr, fitsFileNameReal.c_str(), &status))  
			fitserror( status );
		
		if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) // create new FITS file 
			fitserror( status );           // call printerror if error occurs
		nelements=naxes[0]*naxes[1]*naxes[2];
		
		char bunit[10] = "Lambda";
		char btype[10] = "uv dist";
		fits_write_key(fptr , TSTRING , "BUNIT"    , &bunit   , "units of wavelength"   , &status);
		fits_write_key(fptr , TSTRING , "BTYPE"    , &btype   , "units of wavelength"   , &status);
		
		char ctype1[10] = "u";
		fits_write_key(fptr , TSTRING , "CTYPE1"    , &ctype1   , "axis 1 coordinate type"   , &status);
		double centerPos=0;
		double temp=0;
		fits_write_key(fptr , TDOUBLE , "CRVAL1", &centerPos, "value at reference pixel", &status);
		double distperpixh = width/hpix;
		fits_write_key(fptr , TDOUBLE , "CDELT1", &distperpixh, "increment per pixel", &status);
		int refpix1=hpix/2 + 1;
		fits_write_key(fptr , TINT , "CRPIX1", &refpix1, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA1", &temp, "rotation angle", &status);
		char ctype2[10] = "v";
		fits_write_key(fptr , TSTRING , "CTYPE2"    , &ctype2   , "axis 2 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL2", &centerPos, "value at reference pixel", &status);
		double distperpixv = height/vpix;
		fits_write_key(fptr , TDOUBLE , "CDELT2", &distperpixv, "increment per pixel", &status);
		int refpix2=vpix/2 + 1;
		fits_write_key(fptr , TINT , "CRPIX2", &refpix2, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA2", &temp, "rotation angle", &status);
		char ctype3[10] = "FREQ";
		fits_write_key(fptr , TSTRING , "CTYPE3"    , &ctype3   , "axis 3 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL3", &centfreq, "value at reference pixel", &status);
		if(binwidth==0)
			binwidth=1e9;//set it to 1GHz if we just took a single frame.
		fits_write_key(fptr , TDOUBLE , "CDELT3", &binwidth, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX3", &centfreqnum, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA3", &temp, "rotation angle", &status);
		char cunit[10] = "HZ";
		char ctype[10] = "STOKES";
		fits_write_key(fptr , TSTRING , "CUNIT3"    , &cunit   , "units"   , &status);
		fits_write_key(fptr , TSTRING , "CTYPE4"    , &ctype   , "axis 3 coordinate type"   , &status);
		double temp1=1; double temp2=0; int temp3=1;
		fits_write_key(fptr , TDOUBLE , "CRVAL4", &temp1, "value at reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CDELT4", &temp1, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX4", &temp3, "references pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA4", &temp2, "rotation angle", &status);
		
		if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, &realPart.front(), &status))  
			fitserror( status );
		if (fits_close_file(fptr, &status))  
			fitserror( status );
		
		//now write the imaginary part
		std::cout << "printing file " << fitsFileNameIm << std::endl;
		if (fits_create_file(&fptr, fitsFileNameIm.c_str(), &status))  
			fitserror( status );
		if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) // create new FITS file 
			fitserror( status );           // call printerror if error occurs
		//none of the values for the axes changed, so we can reuse them all
		fits_write_key(fptr , TSTRING , "BUNIT"    , &bunit   , "units of wavelength"   , &status);
		fits_write_key(fptr , TSTRING , "BTYPE"    , &btype   , "units of wavelength"   , &status);
		fits_write_key(fptr , TSTRING , "CTYPE1"    , &ctype1   , "axis 1 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL1", &centerPos, "value at reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CDELT1", &distperpixh, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX1", &refpix1, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA1", &temp, "rotation angle", &status);
		fits_write_key(fptr , TSTRING , "CTYPE2"    , &ctype2   , "axis 2 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL2", &centerPos, "value at reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CDELT2", &distperpixv, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX2", &refpix2, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA2", &temp, "rotation angle", &status);
		fits_write_key(fptr , TSTRING , "CTYPE3"    , &ctype3   , "axis 3 coordinate type"   , &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL3", &centfreq, "value at reference pixel", &status);
		if(binwidth==0)
			binwidth=1e9;//set it to 1GHz if we just took a single frame.
		fits_write_key(fptr , TDOUBLE , "CDELT3", &binwidth, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX3", &centfreqnum, "reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA3", &temp, "rotation angle", &status);
		fits_write_key(fptr , TDOUBLE , "CRVAL4", &temp1, "value at reference pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CDELT4", &temp1, "increment per pixel", &status);
		fits_write_key(fptr , TINT , "CRPIX4", &temp3, "references pixel", &status);
		fits_write_key(fptr , TDOUBLE , "CROTA4", &temp2, "rotation angle", &status);
		
		if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, &imaginaryPart.front(), &status))  
			fitserror( status );
		if (fits_close_file(fptr, &status))  
			fitserror( status );
		
	}
	
	//interpolate in uv space. Given u and v in units of meters (the default for things
	//like casa), return the real and imaginary value of the FFT. freqIndex tells which 
	//channel to go to.
	//bilinear interpolation isn't hard, but I keep having to rederive the whole thing, so, 
	//for my own future reference, here is the schema I'll try to stick to for the rest of this:
	//		
	//		|
	//		|
	//	  y2|---*-------*-------*-----		R1=(x2-x)/(x2-x1)Q11 + (x-x1)/(x2-x1)Q21
	//		|	|Q12	|R2		|Q22		R2=(x2-x)/(x2-x1)Q12 + (x-x1)/(x2-x1)Q22
	//		|	|		|		|
	//	   y|---|-------*-------|-----		P=(y2-y)/(y2-y1)R1 + (y-y1)/(y2-y1)R2
	//		|	|		|P		|
	//		|	|		|		|
	//	  y1|---*-------*-------*-----
	//		|	|Q11	|R1		|Q21
	//		|__________________________
	//			x1		x		x2	
	//		
	//		
	std::pair<double, double> interp(unsigned int freqIndex, double u, double v){
		if(freqIndex >= frequencies.size()){
			std::cout << "Error: Bad channel number in interp." << std::endl;
			exit(1);
		}
		double frequency=frequencies[freqIndex];
		
		//double kilolambda=c/frequency*10; //we also convert to mks here, hence the *10, not *1000
		double pixSizex=(atan(width/distance)/hpix);
		double pixSizey=(atan(height/distance)/vpix);
		
		//output from FFTW ranges from -0.5 to 0.5, in units of freq/npix
		//convert u and v into x and y coordinates
		double x=u*pixSizex*hpix;
		double y=v*pixSizey*vpix;
		
		x+=(int)hpix/2;
		y+=(int)vpix/2;
		double xindex=floor(x); double yindex=floor(y);
		//now we need to check the bounds
		int xnextIndex, ynextIndex;
		xnextIndex=xindex+1; ynextIndex=yindex+1;
		if(xindex<0){
			xindex=0;xnextIndex=1;
		}
		if(xindex>=hpix-1){	//if we're at the edge, back up one and we'll extrapolate
			xindex=hpix-2; xnextIndex=hpix-1;
		}
		if(yindex<0){
			yindex=0;ynextIndex=1;
		}
		if(yindex>=vpix-1){	//if we're at the edge, back up one and we'll extrapolate
			yindex=hpix-2; ynextIndex=hpix-1;
		}
		
		double q11, q12, q21, q22;
		q11=realPart[freqIndex][yindex][xindex];
		q12=realPart[freqIndex][ynextIndex][xindex];
		q21=realPart[freqIndex][yindex][xnextIndex];
		q22=realPart[freqIndex][ynextIndex][xnextIndex];
		//std::cout << xindex << "\t" << yindex << "\t" << q11 << std::endl;
		//std::cout << xindex << "\t" << ynextIndex << "\t" << q12 << std::endl;
		//std::cout << xnextIndex << "\t" << yindex << "\t" << q21 << std::endl;
		//std::cout << xnextIndex << "\t" << ynextIndex << "\t" << q22 << std::endl;
		
		double r1,r2,p;
		r1=(xnextIndex-x)*q11 + (x-xindex)*q21;
		r2=(xnextIndex-x)*q12 + (x-xindex)*q22;
		double interpReal=(ynextIndex-y)/(ynextIndex-yindex)*r1 + (y-yindex)/(ynextIndex-yindex)*r2;
		
		//now the imaginary part
		q11=imaginaryPart[freqIndex][yindex][xindex];
		q12=imaginaryPart[freqIndex][ynextIndex][xindex];
		q21=imaginaryPart[freqIndex][yindex][xnextIndex];
		q22=imaginaryPart[freqIndex][ynextIndex][xnextIndex];
		
		r1=(xnextIndex-x)*q11 + (x-xindex)*q21;
		r2=(xnextIndex-x)*q12 + (x-xindex)*q22;
		double interpIm=(ynextIndex-y)/(ynextIndex-yindex)*r1 + (y-yindex)/(ynextIndex-yindex)*r2;
		
		//std::cout << x << "\t" << y << "\t" << interpReal  << "\t" << interpIm << "\t "<< sqrt(u*u+v*v) << std::endl;
		
		std::pair<double,double> result=std::make_pair(interpReal,interpIm);
		return result;
	}
	
	std::pair<Vec4d, Vec4d> interpAVX(unsigned int freqIndex, Vec4d u, Vec4d v){
		if(freqIndex >= freqbins){
			std::cout << "Error: Bad channel number in interp." << std::endl;
			exit(1);
		}
		double frequency=frequencies[freqIndex];
				
		//output from FFTW ranges from -0.5 to 0.5, in units of freq/npix
		//convert u and v into x and y coordinates
		Vec4d x=u*xpixsize*hpix;
		Vec4d y=v*ypixsize*vpix;
		
		x+=hpix/2;
		y+=vpix/2;
		Vec4i xindex=truncate_to_int32(x); Vec4i yindex=truncate_to_int32(y);
		//now we need to check the bounds
		Vec4i xnextIndex, ynextIndex;
		xnextIndex=xindex+1; ynextIndex=yindex+1;
		
		for(int i=0;i<xindex.size();i++){
			if(xindex[i]<0){
				xindex.insert(i,0);xnextIndex.insert(i,1);
			}
			if(xindex[i]>=hpix-1){	//if we're at the edge, back up one and we'll extrapolate
				xindex.insert(i,hpix-2); xnextIndex.insert(i,hpix-1);
			}
			if(yindex[i]<0){
				yindex.insert(i,0); ynextIndex.insert(i,1);
			}
			if(yindex[i]>=vpix-1){	//if we're at the edge, back up one and we'll extrapolate
				yindex.insert(i,hpix-2); ynextIndex.insert(i,hpix-1);
			}
		}
		Vec4d q11(realPart[freqIndex][yindex[0]][xindex[0]], realPart[freqIndex][yindex[1]][xindex[1]], realPart[freqIndex][yindex[2]][xindex[2]], realPart[freqIndex][yindex[3]][xindex[3]]);
		Vec4d q12(realPart[freqIndex][ynextIndex[0]][xindex[0]], realPart[freqIndex][ynextIndex[1]][xindex[1]], realPart[freqIndex][ynextIndex[2]][xindex[2]], realPart[freqIndex][ynextIndex[3]][xindex[3]]);
		Vec4d q21(realPart[freqIndex][yindex[0]][xnextIndex[0]], realPart[freqIndex][yindex[1]][xnextIndex[1]], realPart[freqIndex][yindex[2]][xnextIndex[2]], realPart[freqIndex][yindex[3]][xnextIndex[3]]);
		Vec4d q22(realPart[freqIndex][ynextIndex[0]][xnextIndex[0]], realPart[freqIndex][ynextIndex[1]][xnextIndex[1]], realPart[freqIndex][ynextIndex[2]][xnextIndex[2]], realPart[freqIndex][ynextIndex[3]][xnextIndex[3]]);
		
		Vec4d r1,r2, interpReal, interpIm;
		Vec4d xindexDouble=to_double(xindex); Vec4d xnextIndexDouble=to_double(xnextIndex);
		Vec4d yindexDouble=to_double(yindex); Vec4d ynextIndexDouble=to_double(ynextIndex);
		r1=(xnextIndexDouble-x)*q11 + (x-xindexDouble)*q21;
		r2=(xnextIndexDouble-x)*q12 + (x-xindexDouble)*q22;
		interpReal=(ynextIndexDouble-y)/(ynextIndexDouble-yindexDouble)*r1 + (y-yindexDouble)/(ynextIndexDouble-yindexDouble)*r2;
		
		//now the imaginary part
		q11=(imaginaryPart[freqIndex][yindex[0]][xindex[0]], imaginaryPart[freqIndex][yindex[1]][xindex[1]], imaginaryPart[freqIndex][yindex[2]][xindex[2]], imaginaryPart[freqIndex][yindex[3]][xindex[3]]);
		q12=(imaginaryPart[freqIndex][ynextIndex[0]][xindex[0]], imaginaryPart[freqIndex][ynextIndex[1]][xindex[1]], imaginaryPart[freqIndex][ynextIndex[2]][xindex[2]], imaginaryPart[freqIndex][ynextIndex[3]][xindex[3]]);
		q21=(imaginaryPart[freqIndex][yindex[0]][xnextIndex[0]], imaginaryPart[freqIndex][yindex[1]][xnextIndex[1]], imaginaryPart[freqIndex][yindex[2]][xnextIndex[2]], imaginaryPart[freqIndex][yindex[3]][xnextIndex[3]]);
		q22=(imaginaryPart[freqIndex][ynextIndex[0]][xnextIndex[0]], imaginaryPart[freqIndex][ynextIndex[1]][xnextIndex[1]], imaginaryPart[freqIndex][ynextIndex[2]][xnextIndex[2]], imaginaryPart[freqIndex][ynextIndex[3]][xnextIndex[3]]);
		
		r1=(xnextIndexDouble-x)*q11 + (x-xindexDouble)*q21;
		r2=(xnextIndexDouble-x)*q12 + (x-xindexDouble)*q22;
		interpIm=(ynextIndexDouble-y)/(ynextIndexDouble-yindexDouble)*r1 + (y-yindexDouble)/(ynextIndexDouble-yindexDouble)*r2;
		
		//std::cout << x << "\t" << y << "\t" << interpReal  << "\t" << interpIm << "\t "<< sqrt(u*u+v*v) << std::endl;
		
		std::pair<Vec4d,Vec4d> result=std::make_pair(interpReal,interpIm);
		return result;
	}
	
	//offset the image in the sky by applying a phase offset to the fourier transform
	//we could technically do this more easily during the interpolation stage above,
	//but for testing it's best to have it separate.
	//the units of the offset supplied are in pixels
	void offset(double xOffset, double yOffset){
		//a translation in real space is easy to apply in fourier space as a phase offset.
		//this is honestly easier than calculating how to repoint the camera or adding and stripping
		//rows of pixels to the final image, as we just need to multiply all the points by a phase factor.
		for(int k=0;k<freqbins;k++){
			double frequency=frequencies[k];
			unsigned int start=k*hpix*vpix;
			//for each point, we need to calculate u and v, and then the phase factor.
			//the FFT is defined on a -0.5 to 0.5 grid (I guess that's in wavenumber space)
			for(size_t i=0; i<vpix; i++){
				for(size_t j=0;j<hpix;j++){
					//first get u and v position
					double u=double(j)/hpix - 0.5;
					double v=double(i)/vpix - 0.5;
					double real=realPart[k][j][i];
					double im=imaginaryPart[k][j][i];
					double phaseOffset=-2*pi*(xOffset*u+yOffset*v); //the real and imaginary parts will rotated by this amount
					double offsetReal = real*cos(phaseOffset) - im*sin(phaseOffset);
					double offsetIm = real*sin(phaseOffset) + im*cos(phaseOffset);
					realPart[k][j][i]=offsetReal;
					imaginaryPart[k][j][i]=offsetIm;
				}
				//std::cout << std::endl;
			}
		}
	}
	
};

//the image class represents the sky, and so is entirely real. So to do a backtransform for something that is complex,
//the input data needs to be either a fourier class, or a pair of images.
fourierImage backTransform(const fourierImage& im){
	unsigned int hpix=im.hpix;
	unsigned int vpix=im.vpix;
	
	fourierImage results(vpix,hpix,im.width,im.height,im.centfreq,im.freqbins, im.frequencies, im.distance);
	//std::cout << im.centfreq << "\t" << im.freqrange << "\t" << im.freqbins << "\t" << hpix << "\t" << vpix << std::endl;
	
	fftw_complex *out1;
	fftw_plan p1;
	//note about out1: It should probably be an std::unique_ptr, like the following things, but because the fftw_complex is itself
	//a c style array of length 2, things get... ugly. There may be a way to do this properly, but it's going to be a mess.
	unsigned int out1size=vpix * hpix;
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
	
	const unsigned int dataSize=hpix*vpix*results.frequencies.size();
	std::unique_ptr<double[]> FFTdataREAL (new double[dataSize]);
	std::unique_ptr<double[]> FFTdataCMPLX (new double[dataSize]);
	fftw_complex* inputdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
	//std::cout << "out1size: " << out1size << std::endl;
	unsigned int length = hpix*vpix; //number of pixels in each image
	for(int k=0; k<im.frequencies.size();k++){
		unsigned int start=k*length;
		int toggle = 1;
		
		for(int i=0; i<hpix;i++,toggle*=-1){
			for(int j=0; j<vpix;j++,toggle*=-1){
				inputdata[j+i*vpix][0]= im.realPart[k][i][j]*toggle;
				inputdata[j+i*vpix][1]= im.imaginaryPart[k][i][j]*toggle;
			}
		}
		p1=fftw_plan_dft_2d(hpix, vpix, &inputdata[0], out1,-1, FFTW_ESTIMATE);
		fftw_execute(p1);

		//now we need to extract the data from out1 and put it into the real and complex arrays
		toggle = 1;
		for(size_t i=0; i<vpix; i++,toggle*=-1){
			for(size_t j=0;j<hpix;j++,toggle*=-1){
				FFTdataREAL[start+j+i*vpix]=out1[j+i*vpix][0]*toggle/hpix/vpix;
				FFTdataCMPLX[start+j+i*vpix]=out1[j+i*vpix][1]*toggle/hpix/vpix;
			}
			//std::cout << std::endl;
		}
		//now to put the FFT data into the final images
		std::copy(FFTdataREAL.get(),FFTdataREAL.get()+length,results.realPart.begin()+start);
		std::copy(FFTdataCMPLX.get(),FFTdataCMPLX.get()+length,results.imaginaryPart.begin()+start);
		
	}
	fftw_free(out1);
	fftw_free(inputdata);
	return results;
}

//other signiture
fourierImage backTransform(const image& realInput, const image& imaginaryInput){
	//make sure the images are compatible
	if(realInput.hpix != imaginaryInput.hpix ||
		realInput.vpix != imaginaryInput.vpix ||
		realInput.freqbins != imaginaryInput.freqbins){
			std::cout << "ERROR: cannot do fourier transform given two images of different size" << std::endl;
			exit(1);
		}
	image RealPart=realInput;	
	image ImaginaryPart=imaginaryInput;	//initialize both of these to the starting image to get most of the parameters right
	//we'll overwrite the data, and need to change the axes/units later

	unsigned int hpix=realInput.hpix;
	unsigned int vpix=realInput.vpix;
	
	fourierImage results(vpix,hpix,realInput.width,realInput.height,realInput.centfreq,realInput.freqbins, realInput.frequencies, realInput.distance);
	//std::cout << im.centfreq << "\t" << im.freqrange << "\t" << im.freqbins << "\t" << hpix << "\t" << vpix << std::endl;
	
	fftw_complex *out1;
	fftw_plan p1;
	//note about out1: It should probably be an std::unique_ptr, like the following things, but because the fftw_complex is itself
	//a c style array of length 2, things get... ugly. There may be a way to do this properly, but it's going to be a mess.
	unsigned int out1size=vpix * hpix;
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
	
	const unsigned int dataSize=hpix*vpix*results.frequencies.size();
	std::unique_ptr<double[]> FFTdataREAL (new double[dataSize]);
	std::unique_ptr<double[]> FFTdataCMPLX (new double[dataSize]);
	fftw_complex* inputdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
	//std::cout << "out1size: " << out1size << std::endl;
	unsigned int length = hpix*vpix; //number of pixels in each image
	for(int k=0; k<results.frequencies.size();k++){
		unsigned int start=k*length;
		int toggle = 1;
		
		for(int i=0; i<hpix;i++,toggle*=-1){
			for(int j=0; j<vpix;j++,toggle*=-1){
				inputdata[j+i*vpix][0]= realInput.data[k][i][j]*toggle;
				inputdata[j+i*vpix][1]= imaginaryInput.data[k][i][j]*toggle;
			}
		}
		p1=fftw_plan_dft_2d(hpix, vpix, &inputdata[0], out1,-1, FFTW_ESTIMATE);
		fftw_execute(p1);

		//now we need to extract the data from out1 and put it into the real and complex arrays
		toggle = 1;
		for(size_t i=0; i<vpix; i++,toggle*=-1){
			for(size_t j=0;j<hpix;j++,toggle*=-1){
				FFTdataREAL[start+j+i*vpix]=out1[j+i*vpix][0]*toggle;
				FFTdataCMPLX[start+j+i*vpix]=out1[j+i*vpix][1]*toggle;
			}
			//std::cout << std::endl;
		}
		//now to put the FFT data into the final images
		std::copy(FFTdataREAL.get(),FFTdataREAL.get()+length,results.realPart.begin()+start);
		std::copy(FFTdataCMPLX.get(),FFTdataCMPLX.get()+length,results.imaginaryPart.begin()+start);
		
	}
	fftw_free(out1);
	fftw_free(inputdata);
	return results;
}

fourierImage FFTDifferent(const image& im){
	image RealPart=im;	
	image ImaginaryPart=im;	//initialize both of these to the starting image to get most of the parameters right
	//we'll overwrite the data, and need to change the axes/units later

	unsigned int hpix=im.hpix;
	unsigned int vpix=im.vpix;
	fourierImage results(vpix,hpix,im.width,im.height,im.centfreq,im.freqbins, im.frequencies, im.distance);
	
	for(int f=0;f<results.frequencies.size();f++){
		unsigned int start=f*hpix*vpix;
		//note about out1: It should probably be an std::unique_ptr, like the following things, but because the fftw_complex is itself
		//a c style array of length 2, things get... ugly. There may be a way to do this properly, but it's going to be a mess.
		fftw_complex *out1;
		fftw_plan p1;
		unsigned int out1size=vpix * hpix;
		out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
	
		fftw_complex* inputdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out1size);
		unsigned int length = hpix*vpix; //number of pixels in each image
		int toggle = 1;
		
		for(int i=0; i<vpix;i++,toggle*=-1){
			for(int j=0; j<hpix;j++,toggle*=-1){
				inputdata[j+i*hpix][0]= im.data[f][i][j]*toggle;
				inputdata[j+i*hpix][1]= 0;
			}
		}
		p1=fftw_plan_dft_2d(hpix, vpix, &inputdata[0], out1,-1, FFTW_ESTIMATE);
		fftw_execute(p1);

		//now we need to extract the data from out1 and put it into the real and complex arrays
		toggle = 1;
		for(size_t i=0; i<hpix; i++,toggle*=-1){
			for(size_t j=0;j<vpix;j++,toggle*=-1){
				results.realPart[f][j][i]=out1[i+j*hpix][0]*toggle;
				results.imaginaryPart[f][j][i]=out1[i+j*hpix][1]*toggle;
			}
		}
		//now to put the FFT data into the final images
		fftw_free(out1);
		fftw_free(inputdata);
	}
	return results;
}

fourierImage FFTmultiThread(const image& im, ThreadPool& pool){
	const unsigned int hpix=im.hpix;
	const unsigned int vpix=im.vpix;
	const unsigned int length = hpix*vpix; //number of pixels in each image
	fourierImage results(vpix,hpix,im.width,im.height,im.centfreq,im.freqbins, im.frequencies, im.distance);
	std::vector<std::future<void> > slices;

	std::vector<fftw_complex*> inputs(results.frequencies.size()), outputs(results.frequencies.size());
	std::vector<fftw_plan> plans(results.frequencies.size());
	for(int f=0;f<results.frequencies.size();f++){
		inputs[f] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
		outputs[f] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
		plans[f]=fftw_plan_dft_2d(hpix, vpix, inputs[f], outputs[f], FFTW_FORWARD, FFTW_MEASURE);
	}

	for(int f=0;f<results.frequencies.size();f++){
		slices.emplace_back(pool.enqueue([&im,&results,&inputs,&outputs,&plans,f,length,vpix,hpix](){
			auto input=inputs[f];
			auto output=outputs[f];
			auto plan=plans[f];

			int toggle = 1;
			for(int i=0; i<vpix;i++,toggle*=-1){
				for(int j=0; j<hpix;j++,toggle*=-1){
					input[j+i*hpix][0]= im.data[f][i][j]*toggle;
					input[j+i*hpix][1]= 0;
				}
			}
			fftw_execute(plan);

			//now we need to extract the data from out1 and put it into the real and complex arrays
			toggle = 1;
			for(size_t i=0; i<hpix; i++,toggle*=-1){
				for(size_t j=0;j<vpix;j++,toggle*=-1){
					results.realPart[f][j][i]=output[i+j*hpix][0]*toggle;
					results.imaginaryPart[f][j][i]=output[i+j*hpix][1]*toggle;
				}
			}
		}));
	}
	
	for(auto& slice : slices) //wait for all jobs to finish
		slice.get();

	for(int f=0;f<results.frequencies.size();f++){
		fftw_free(inputs[f]);
		fftw_free(outputs[f]);
		fftw_destroy_plan(plans[f]);
	}
	return results;
}

//a point in uv space, including real, imaginary components, and a weight
//u and v are taken to be in units of *meters*.
//everything else is in cgs, but this needs to match casa, which gives this in meters.
//real and imaginary are in units of janskys.
struct uvPoint{
	double u,v,real,imaginary,weight;	
};

//beam smearing is done via convlution with a (usually gaussian) beam.
//the convolution is done by just pointwise multiplying the fourier transforms
//together, so it makes sense to have a separate image type for the beam.
//Since the beam is generally a fairly simple 2d gaussian, we can calculate it
//and its fourier transform when it's initialized.
struct beam{
	unsigned int vpix, hpix; //vertical and horizontal pixel numbers
	double width, height; //centimeters
	marray<double,2> data;	//actual beam, in the image plane
	marray<double,2> realPart;	//FFT
	marray<double,2> imaginaryPart;	//FFT
	double distance;
	
	//the usual constructors. These leave the actual arrays empty
	beam():vpix(0),hpix(0),width(0),height(0),distance(0) {	}
	beam(unsigned int vpix, unsigned int hpix, double width, double height, double distance):vpix(vpix),hpix(hpix),width(width),height(height),distance(distance){	}
	beam(const beam& other):vpix(other.vpix),hpix(other.hpix),width(other.width),height(other.height),distance(other.distance) {	}
	beam& operator=(const beam& other){
		hpix=other.hpix;
		vpix=other.vpix;
		width=other.width;
		height=other.height;
		distance=other.distance;
		data=other.data;
		realPart=other.realPart; imaginaryPart=other.imaginaryPart;
		return(*this);
	}
	
	//special constructor:
	//This makes the beam associated with a given image. It takes an image, the major and minor beam axes, and the position angle
	//Note that the axes are specified in arcseconds. This is because that's what comes out of CASA and the actual observations.
	beam(const image& im, double bmaj, double bmin, double PA):vpix(im.vpix),hpix(im.hpix),width(im.width),height(im.height),distance(im.distance), 
		data({hpix,vpix}),realPart({hpix,vpix}),imaginaryPart({hpix,vpix}){
		//first we need to convert the beam dimensions from arcsec to pixels
		double degperpixh = (atan(width/(distance))/hpix) * (180/pi);
		double degperpixv = (atan(height/(distance))/vpix)* (180/pi);
		double majorAxis=bmaj/3600/degperpixh;
		double minorAxis=bmin/3600/degperpixh;
		
		//now we need to construct the image plane beam and its FFTs
		for(int i=0;i<hpix;i++){
			for(int j=0;j<vpix;j++){
				double xpos=(double)i; //plus a half pixel to offset
				double ypos=(double)j;

				//rotate
				double newXpos=xpos*cos(PA)+ypos*sin(PA);
				double newYpos=xpos*sin(PA)-ypos*cos(PA);
				double newXcen=(double)hpix/2*cos(PA)+(double)vpix/2*sin(PA);
				double newYcen=(double)hpix/2*sin(PA)-(double)vpix/2*cos(PA);
				data[j][i]=gaussian(newXpos,newXcen, majorAxis)*gaussian(newYpos,newYcen, minorAxis);
				double prefactor=1/(2*pi)/majorAxis/minorAxis;
				
				//we can do the FFTs here too
				//because we're using a gaussian, the formulas are all analytic
				//see Isella 2013 for the full derivation
				//xpos and ypos now correspond to u and v in the fourier plane
				double xinc=1.0/hpix; double yinc=1.0/vpix;
				double xpos2=xinc*(i+0.5)-0.5; double ypos2=yinc*(j+0.5)-0.5;
				realPart[j][i]=2*prefactor*pi*majorAxis*minorAxis*exp(-2*pi*pi*( (xpos2*cos(PA) + ypos2*sin(PA))*(xpos2*cos(PA) + ypos2*sin(PA))*majorAxis*majorAxis + 
					(-xpos2*sin(PA) + ypos2*cos(PA))*(-xpos2*sin(PA) + ypos2*cos(PA))*minorAxis*minorAxis ));
				imaginaryPart[j][i]=0; //for a function centered at 0,0, the imaginary part is 0 because of that phase term
			}
		}
		
	}
	
};

struct radialImage{
	struct point{
		double radius;
		double value; //the value can be either surface brightness or brightness temperature
		point(double radius, double value):radius(radius),value(value) {}
	};
	std::vector<point> data;
	double rmin, rmax;
	unsigned int npoints;
	
	radialImage():rmin(0),rmax(0),npoints(0) {}
	radialImage(double rmin, double rmax, unsigned int npoints):rmin(rmin),rmax(rmax),npoints(npoints),data(){}
	radialImage(const radialImage& other):rmin(other.rmin),rmax(other.rmax),npoints(other.npoints),data(other.data){}
	
	template<class density>
	void propagate(const grid<density>& g, typename grid<density>::prop_type type){
		data.clear();
		std::vector<double> frequencies; //for now, this will just be moonochromatic
		frequencies.push_back(230.358e9);
		double rpos;
		double range=rmax-rmin;
		double stepsize=range/npoints;
		//initialize position
		rpos=rmin+stepsize/2;
		vect down(0,0,-1);
		for(int i=0;i<npoints;i++,rpos+=stepsize){
			vect pos(rpos,0,1e18); //keep y at 0, keep z way above the disk
			line path(pos,down);
			std::vector<double> pix=g.propagateRay(path,frequencies,pos,type);
			double btemp=tempFromSurfaceBrightness(pix[0], 230.358e9);
			double midptemp=g.temperature(rpos,pi/2,0);
			data.push_back(point(rpos,btemp/midptemp));
		}
	}
	
	
};
