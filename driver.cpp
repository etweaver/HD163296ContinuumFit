#include "driver.h"
#include "ParameterSet.h"
#include <random>
#include <sstream>
#include <string>

const double pi=4*atan(1);

struct moveProposal{
        std::vector<double> coordinates;
        double logAcceptanceFactor;

        moveProposal()=default;
        moveProposal(std::size_t n):coordinates(n){}
        moveProposal(moveProposal&& mp):
        coordinates(std::move(mp.coordinates)),logAcceptanceFactor(mp.logAcceptanceFactor){}
};

struct ensembleMember{
	std::vector<double> coords;
	double currentValue;
	std::vector<double> proposedCoords;
	std::vector<double> partialValues;
	double logAcceptanceFactor;
	
	ensembleMember()=default;
	ensembleMember(unsigned int nFrequencies):partialValues(nFrequencies,0) {	}

	void proposeJump(moveProposal mp){
		proposedCoords=std::move(mp.coordinates);
		logAcceptanceFactor=mp.logAcceptanceFactor;
	}

	template<typename RNG>
	bool decideJump(std::uniform_real_distribution<double>& acceptDist, RNG& rng, ParameterSet& params){
		double proposedValue=0;
		for(auto v : partialValues)
			proposedValue+=v;
		double logRatio=(currentValue-proposedValue)+logAcceptanceFactor;
		bool accept=(logRatio>=0) || (logRatio>log(acceptDist(rng)));
		if(!params.inBounds(proposedCoords)){
			accept=false;
			std::cout << "step rejected: out of bounds" << std::endl;
		}
		if(accept){
			std::cout << "move accepted: old value: " << currentValue << ", new value: " << proposedValue << ", logAcceptanceFactor: " << logAcceptanceFactor << std::endl;
			std::copy(proposedCoords.begin(),proposedCoords.end(),coords.begin());
			currentValue=proposedValue;
		}
		return accept;
	}
};

struct stretchMove{

        //square root distribution
        struct sqrtDist{
                double range; //extends from 1/r to r
                double norm;
                mutable std::uniform_real_distribution<double> uniform; //mutable is needed to override a const later because the
                                                                                                                                //uniform distribution has no const call operator. Ugly but necessary.

                sqrtDist(double range): range(range),norm(1/((sqrt(range)-sqrt(1/range)))),uniform(norm*sqrt(1/range),norm*sqrt(range)){
                        if(range<=1)
                                throw std::domain_error("square_root_distribution requires range>1");
                }

                template<typename RNG>
                double operator()(RNG& rng) const{
                        double v=uniform(rng)/(norm);
                        return(v*v);
                }
        };

        sqrtDist jumpDist;
        stretchMove():jumpDist(2){}

        template<typename RNG>
        moveProposal operator()(const std::vector<double>& coords, const std::vector<ensembleMember>& ensemble, RNG& rng) const {
                assert(!coords.empty());
                assert(ensemble.size() > 1);
                moveProposal mp(coords.size());

                //choose a random member of the ensemble, but not the one we are already using
                std::uniform_int_distribution<int> idxDist(0,ensemble.size()-1);
                int idx=idx=idxDist(rng);;
                unsigned int maxTrials=100; //if we need to keep trying
                //std::cout << coords.size() << "\t" << ensemble.size() << "\t" << idx << std::endl;
                while(std::equal(coords.begin(),coords.end(),ensemble[idx].coords.begin())){
                        idx=idxDist(rng);
                        if(!--maxTrials)
                                throw std::runtime_error("StretchMove failed too many times to find a distinct ensemble member. "
                                                         "Ensmeble may have become degenerate.");
                }

                //jump distance
                double z=jumpDist(rng);

                //direction
                for(std::size_t i=0; i<coords.size(); i++)
                        mp.coordinates[i]=ensemble[idx].coords[i]+z*(coords[i]-ensemble[idx].coords[i]);
                //and the penalty associated with going there
                mp.logAcceptanceFactor=(coords.size()-1)*log(z);

                return mp;
        }
};

template<typename RNG>
double randInRange(double min, double max, RNG& rng){
	std::uniform_real_distribution<double> dist(min,max);
	return dist(rng);
}

template<typename RNG, class Jumper>
class modelingServer : public WorkServerImpl{
	public:
		modelingServer()=default;
		modelingServer(unsigned int nSamples, unsigned int nFrequencies, unsigned int ensembleSize,
			       	const ParameterSet& params, std::size_t rngSeed, bool readFromFile):
			nSamples(nSamples),nFrequencies(nFrequencies),ensembleSize(ensembleSize),counter(0),
			ensembleCounter(0),accepted(0),rejected(0),params(params),rng(rngSeed),firstGeneration(true){
				if(readFromFile){
					std::cout << "You haven't updated this for CO yet" << std::endl;
					exit(1);
					std::string lineIn;
					firstGeneration=false;
					std::ifstream statefile("finalState.txt");
					int numRead=0;
					while(!statefile.eof() && numRead < ensembleSize){
						ensembleMember member(nFrequencies);
						std::getline(statefile, lineIn);
						if(lineIn[0]=='#')
							continue;
						if(lineIn[0]=='\n')
							continue;
						std::vector<double> ensembleMember;
						double inc, PA, logSig, rc, P, S, h0, r1, r2, r3, d1, d2, d3, w1, w2, w3, like;
						std::stringstream stream;
						stream << lineIn;
						stream >> inc>> PA>> logSig >> rc >> P>> h0>> S>> r1 >> r2 >> r3>> d1>> d2>> d3>> w1>> w2>> w3>> like;

						member.coords = {inc, PA, logSig, rc, P, h0, S, r1, r2, r3, d1, d2, d3, w1, w2, w3, like};
						member.currentValue=like;
						ensemble.push_back(member);
						numRead++;
						//for(auto i : ensembleMember)		
							//std::cout << i << "\t";
						//std::cout << std::endl;
					}
					statefile.close();
					std::cout << "Ensemble loaded from file" << std::endl;
				}else{
					for(int i=0;i<ensembleSize;i++){
						ensembleMember member(nFrequencies);
						member.coords.push_back(randInRange(0.77,0.78,rng));//inc
						member.coords.push_back(randInRange(2.31,2.33,rng));//PA
						member.coords.push_back(randInRange(0.34,0.345,rng));//logSigma0
						member.coords.push_back(randInRange(79.91,71.92,rng));//rc
						member.coords.push_back(randInRange(-0.1,-0.11,rng));//P
						member.coords.push_back(randInRange(0.45,0.55,rng));//h0
        				//member.coords.push_back(randInRange(9,10,rng));//h0gas
						member.coords.push_back(randInRange(0.7,1.1,rng));//S
        				//member.coords.push_back(randInRange(0.7,1.2,rng));//Sgas*/
        				member.coords.push_back(randInRange(45.6,45.7,rng));//r1
        				member.coords.push_back(randInRange(90.4,91.2,rng));//r2
        				member.coords.push_back(randInRange(98,99,rng));//r3
        				member.coords.push_back(randInRange(0.99,0.995,rng));//d1
        				member.coords.push_back(randInRange(0.95,1,rng));//d2
        				member.coords.push_back(randInRange(-1.4,1.3,rng));//d3
        				member.coords.push_back(randInRange(11.1,11.2,rng));//w1
        				member.coords.push_back(randInRange(3.42,3.5,rng));//w2
        				member.coords.push_back(randInRange(1,2,rng));//w3
        				//member.coords.push_back(randInRange(1.2,1.28,rng));//mStar
        				//member.coords.push_back(randInRange(-50000,5000,rng));//deltaF*/
        				//member.coords.push_back(randInRange(-3,3,rng));//dx
        				//member.coords.push_back(randInRange(-3,3,rng));//dy
        				ensemble.push_back(member);
					}
				}
				/*std::cout << "starting ensemble: " << std::endl;
				for(int i=0;i<ensembleSize;i++){
					for(int j=0;j<17;j++){
						std::cout << ensemble[i].coords[j] << "\t";	
					}
					std::cout << std::endl;
				}*/
		}
		
	protected:
		virtual void processResult(WorkItem item, double value) override{
			//std::cout << "processing result" << std::endl;
			std::cout << item << std::endl;
			unsigned int chainIndex=(unsigned int)item.parameters[0];
			unsigned int freqIndex=(unsigned int)item.parameters[1];
			//std::cout << "got indeces: " << chainIndex << "\t" << freqIndex << std::endl;
			if(chainIndex > ensembleSize)
				throw std::runtime_error("Error: Bad chain index");
			if(freqIndex > nFrequencies)
				throw std::runtime_error("Error: Bad frequency index");
			//to construct the ensemble from each result, we need to go to
			//the correct spot for the member piece, and add its value to the 
			//value list. The coordinates are already there from when the proposed
			//ensemble was generated.
			ensemble[chainIndex].partialValues[freqIndex]=value;
			ensembleCounter++;
			
			//std::cout << "value placed" << std::endl;
			//if we're done with the set, we need to process it
			if(ensembleCounter.load()==ensembleSize*nFrequencies){
				//std::cout << "processing finished ensemble" << std::endl;
				if(firstGeneration){
					for(int i=0;i<ensembleSize;i++){
						ensemble[i].currentValue=0;
						for(auto v : ensemble[i].partialValues)
							ensemble[i].currentValue+=v;
					}
					firstGeneration=false;
				}
				else{
					unsigned int acceptedThisRun=0;
					unsigned int rejectedThisRun=0;
					std::uniform_real_distribution<double> acceptDist(0,1);
					for(int i=0;i<ensembleSize;i++){
						if(ensemble[i].decideJump(acceptDist,rng,params)){
							acceptedThisRun++;
						}else{
							rejectedThisRun++;
						}
					}
					accepted+=acceptedThisRun;
					rejected+=rejectedThisRun;
					counter++;
					std::ofstream outfile;
					outfile.open("samples.txt", std::ofstream::app);
					for(int i=0;i<ensembleSize;i++){
						for(int j=0;j<ensemble[i].coords.size();j++){
							outfile << ensemble[i].coords[j] << "\t";
						}
						outfile << ensemble[i].currentValue << std::endl;
					}
					outfile.close();
				}
				ensembleCounter=0;
			}

		}

		virtual bool done() const override{
			return counter.load()>=nSamples;
		}
			
		//I think that this should generate ensembleSize*nFrequencies work units
		//starting from an ensemble generated by generateNextStep(). Each needs 
		//its ensemble number and frequency number and then to be sent out.
		virtual std::vector<WorkItem> generateNextWorkBlock() override{
			std::cout << "generating work block.";
			if(firstGeneration)
				std::cout << "\tFirst Generation.";
			std::cout << std::endl;
			std::vector<WorkItem> workBlock;
			for(int i=0;i<ensembleSize;i++){
				auto proposal=jumper(ensemble[i].coords,ensemble,rng);
				ensemble[i].proposeJump(std::move(proposal));
				for(int j=0; j<nFrequencies; j++){
					WorkItem w({(double)i,(double)j});
					for(int k=0; k<ensemble[i].coords.size(); k++)
						w.parameters.push_back(firstGeneration ? ensemble[i].coords[k] : ensemble[i].proposedCoords[k]);
					workBlock.push_back(w);
				}
			}
			return workBlock;
			//return {WorkItem{{0,5,0.129439, 137.253, 0.044272, 1.44035, 0.0695371, 12.7746, 0.989993, 0.976257, 0.4996, 0.8097, 5.30054, 8.2923, 12.8918, 1.21569, -391105, -0.280648, 14.6624}}};
		}

	private:
		unsigned int nSamples;
		unsigned int nFrequencies;
		unsigned int ensembleSize;
		std::atomic<unsigned int> counter;		//from 0 to nSamples
		std::atomic<unsigned int> ensembleCounter;	//from 0 to ensembleSize*nFrequencies. Just for tracking when each ensemble finishes
		std::vector<ensembleMember> ensemble;
		unsigned int accepted;				//steps that were accepted
		unsigned int rejected;
		ParameterSet params;
		RNG rng;
		Jumper jumper;
		bool firstGeneration;
		//list of all states of the ensemble over time
	//	std::vector<std::vector<std::vector<double>>> ensembleHistory;

};

int main(int argc, char* argv[]){
	std::mt19937 rng(137);

	auto randInRange=[&rng](const double min, const double max){
		double range = max-min;
		double randNum=rng()/(double)rng.max();
		return min+randNum*range;
	};

	ParameterSet params;
	params.addParameter("inclination");
	params.setParameterLowerLimit("inclination",0); params.setParameterUpperLimit("inclination",pi/2);
	params.addParameter("PA");
	params.setParameterLowerLimit("PA",0); params.setParameterUpperLimit("PA",2*pi);
	params.addParameter("logSigma0");
	params.setParameterLowerLimit("logSigma0",log10(1e-1)); params.setParameterUpperLimit("logSigma0",log10(1e2));
	params.addParameter("rc"); //units of au
	params.setParameterLowerLimit("rc",50); params.setParameterUpperLimit("rc",250);
	params.addParameter("P"); //unitless
	params.setParameterLowerLimit("P",-1); params.setParameterUpperLimit("P",1.999);
	params.addParameter("h0dust"); //au
	params.setParameterLowerLimit("h0dust",0.01); params.setParameterUpperLimit("h0dust",15);
	//params.addParameter("h0gas"); //au
	//params.setParameterLowerLimit("h0gas",0.01); params.setParameterUpperLimit("h0gas",15);
	params.addParameter("Sdust"); //unitless
	params.setParameterLowerLimit("Sdust",-1); params.setParameterUpperLimit("Sdust",1.25);
	//params.addParameter("Sgas"); //unitless
	//params.setParameterLowerLimit("Sgas",-1); params.setParameterUpperLimit("Sgas",1.25);
	params.addParameter("r1"); //au
	params.setParameterLowerLimit("r1",30); params.setParameterUpperLimit("r1",65);
	params.addParameter("r2"); //au
	params.setParameterLowerLimit("r2",70); params.setParameterUpperLimit("r2",100);
	params.addParameter("r3"); //au
	params.setParameterLowerLimit("r3",120); params.setParameterUpperLimit("r3",180);
	params.addParameter("d1"); //unitless
	params.setParameterLowerLimit("d1",0); params.setParameterUpperLimit("d1",1);
	params.addParameter("d2"); //unitless
	params.setParameterLowerLimit("d2",0); params.setParameterUpperLimit("d2",1);
	params.addParameter("d3"); //unitless
	params.setParameterLowerLimit("d3",0); params.setParameterUpperLimit("d3",1);
	params.addParameter("w1"); //au
	params.setParameterLowerLimit("w1",0); params.setParameterUpperLimit("w1",25);
	params.addParameter("w2"); //au
	params.setParameterLowerLimit("w2",0); params.setParameterUpperLimit("w2",25);
	params.addParameter("w3"); //au
	params.setParameterLowerLimit("w3",0); params.setParameterUpperLimit("w3",50);
	/*params.addParameter("mStar"); //mSun
	params.setParameterLowerLimit("mStar",0.1); params.setParameterUpperLimit("mStar",4);
	params.addParameter("deltaF"); //Hz
	params.setParameterLowerLimit("deltaF",-1e6); params.setParameterUpperLimit("deltaF",1e6);
	params.addParameter("dx"); //pixels
	params.setParameterLowerLimit("dx",-25); params.setParameterUpperLimit("dx",25);
	params.addParameter("dy"); //pixels
	params.setParameterLowerLimit("dy",-25); params.setParameterUpperLimit("dy",25);*/
	
	modelingServer<std::mt19937,stretchMove> testServer(500,1,100,params,137,false);
	testServer.setWatchdogInterval(std::chrono::seconds(1));
	testServer.setMaxWorkerSilenceTime(std::chrono::seconds(2));
	testServer.run();

	return 0;
}
