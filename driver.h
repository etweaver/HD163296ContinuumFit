#ifndef DRIVER_H
#define DRIVER_H

#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <thread>

#include <libcuckoo/cuckoohash_map.hh>

#include <capnp/ez-rpc.h>

#include "build/distribute.capnp.h"

///A record of a remote worker thread
struct Worker{
	Worker():id(0){}
	Worker(uint64_t id):id(id){
		updateSeenTime();
	}
	
	uint64_t id;
	std::chrono::steady_clock::time_point lastSeen;
	
	void updateSeenTime(){
		lastSeen=std::chrono::steady_clock::now();
	}
};

bool operator==(const Worker& w1, const Worker& w2){
	return w1.id==w2.id;
}

struct WorkerHasher{
	using result_type=std::size_t;
	using argument_type=Worker;
	
	result_type operator()(const Worker& w) const{
		return hasher(w.id);
	}
	
	std::hash<uint64_t> hasher;
};

///A record for a computation to be done, which is described by some number of 
///floating point parameters. 
struct WorkItem{
	WorkItem(std::vector<double> params):
	parameters(std::move(params)),assignedWorkerID(0)
	{}
	
	void start(uint64_t workerID){
		assignedWorkerID=workerID;
		startTime=std::chrono::steady_clock::now();
	}
	
	std::vector<double> parameters;
	uint64_t assignedWorkerID;
	std::chrono::steady_clock::time_point startTime;
};

std::ostream& operator<<(std::ostream& is, const WorkItem& item){
	is << '[';
	for(double p : item.parameters)
		is << p << ' ';
	return is << ']';
}

bool operator==(const WorkItem& i1, const WorkItem& i2){
	return i1.parameters==i2.parameters;
}

struct WorkItemHasher{
	using result_type=std::size_t;
	using argument_type=WorkItem;
	
	result_type operator()(const WorkItem& item) const{
		result_type result=0;
		for(double p : item.parameters)
			result^=hasher(p);
		return result;
	}
	
	std::hash<double> hasher;
};

class WorkServerImpl: public WorkServer::Server {
public:
	WorkServerImpl():
	listeningPort(14000),
	doneFlag(kj::newPromiseAndFulfiller<void>()),
	watchdogInterval(std::chrono::seconds(5)),
	maxWorkerSilenceTime(std::chrono::seconds(5)){
		token=0;
		std::ifstream rfile("/dev/urandom");
		while(token==0){
			rfile.read((char*)&token,sizeof(token));
			if(!rfile)
				throw std::runtime_error("Unable to read random data from /dev/urandom");
		}
	}
	
	~WorkServerImpl(){
		if(watchdog.joinable())
			watchdog.join();
	}
	
	kj::Promise<void> registerWorker(RegisterWorkerContext context) override{
		auto workerID = context.getParams().getId().getHash();
		workers.insert(workerID,Worker(workerID));
		std::cout << "Worker " << workerID << " registered" << std::endl;
		context.getResults().getResp().setType(WorkServer::RegistrationResponse::Type::ACKNOWLEDGED);
		return kj::READY_NOW;
	}
	kj::Promise<void> heartbeat(HeartbeatContext context) override{
		auto workerID = context.getParams().getId().getHash();
		//std::cout << "Worker " << workerID << " sent heartbeat" << std::endl;
		workers.update_fn(workerID,[](auto& worker){ worker.updateSeenTime(); });
		if(done()){
			context.getResults().getResp().setType(WorkServer::HeartbeatResponse::Type::SHUT_DOWN);
			std::cout << "Directed worker " << workerID << " to shut down" << std::endl;
			workers.erase(workerID);
		}
		else
			context.getResults().getResp().setType(WorkServer::HeartbeatResponse::Type::ACKNOWLEDGED);
		return kj::READY_NOW;
	}
	kj::Promise<void> requestWork(RequestWorkContext context) override{
		auto workerID = context.getParams().getId().getHash();
		workers.insert(workerID,Worker(workerID));
		//std::cout << "Worker " << workerID << " requested more work" << std::endl;
		if(done()){ //we're done; instruct the worker to exit
			context.getResults().getResp().setType(WorkServer::WorkResponse::Type::SHUT_DOWN);
			std::cout << "Directed worker " << workerID << " to shut down" << std::endl;
			workers.erase(workerID);
		}
		else{
			workers.update_fn(workerID,[](auto& worker){ worker.updateSeenTime(); });
			std::lock_guard<std::mutex> lock(readyWorkMutex);
			if(readyWork.empty()){
				std::cout << "Had to ask worker " << workerID << " to wait" << std::endl;
				context.getResults().getResp().setType(WorkServer::WorkResponse::Type::WAIT);
			}
			else{
				auto item=readyWork.front();
				readyWork.pop();
				item.start(workerID);
				outstandingWork.insert(item,item);
				//std::cout << "Assigning work item " << item << " to worker " << workerID << std::endl;
				context.getResults().getResp().setType(WorkServer::WorkResponse::Type::WORK_ITEM);
				auto respParams=context.getResults().getResp().initParameters(item.parameters.size());
				for(std::size_t i=0; i!=item.parameters.size(); i++)
					respParams.set(i,item.parameters[i]);
			}
		}
		return kj::READY_NOW;
	}
	kj::Promise<void> workResult(WorkResultContext context) override{
		auto workerID = context.getParams().getId().getHash();
		auto workResult = context.getParams().getResult();
		std::vector<double> resultParams;
		resultParams.reserve(workResult.getParameters().size());
		for(double p : workResult.getParameters())
			resultParams.push_back(p);
		WorkItem resultItem(std::move(resultParams));
		WorkItem inProgress(resultItem);
		try{
			inProgress=outstandingWork.find(resultItem);
		}catch(...){}
		if(inProgress.assignedWorkerID!=workerID){
			std::cout << "Work item " << resultItem << " is assigned to worker " 
			<< inProgress.assignedWorkerID << ", not " << workerID << std::endl;
		}
		else if(workResult.isResult()){
			//std::cout << "Worker " << workerID << " completed work item "
			//	<< resultItem << std::endl;
			outstandingWork.erase(inProgress);
			processResult(inProgress,workResult.getResult());
			//TODO: getting more work generated should happen asynchronously
			if(outstandingWork.empty() && !done()){
				std::lock_guard<std::mutex> lock(readyWorkMutex);
				if(readyWork.empty()){
					auto newWork=generateNextWorkBlock();
					for(WorkItem& item : newWork)
						readyWork.emplace(std::move(item));
				}
			}
		}
		else if(workResult.isErrMsg()){
			outstandingWork.erase(inProgress);
			std::cout << "Worker " << workerID << " failed to compute work item "
			<< resultItem << ": " << workResult.getErrMsg().cStr() << std::endl;
			std::lock_guard<std::mutex> lock(readyWorkMutex);
			readyWork.push(inProgress);
		}
		if(!done())
			context.getResults().getResp().setType(WorkServer::WorkResultResponse::Type::ACKNOWLEDGED);
		else{
			context.getResults().getResp().setType(WorkServer::WorkResultResponse::Type::SHUT_DOWN);
		}
		return kj::READY_NOW;
	}
	kj::Promise<void> workerShutdown(WorkerShutdownContext context) override{
		auto workerID = context.getParams().getId().getHash();
		std::cout << "Worker " << workerID << " sent shutdown notification" << std::endl;
		workers.erase(workerID);
		return kj::READY_NOW;
	}
	kj::Promise<void> serverShutdown(ServerShutdownContext context) override{
		auto sentToken = context.getParams().getToken();
		if(sentToken==token) //ignore unauthorized shutdown requests
			doneFlag.fulfiller->fulfill();
		return kj::READY_NOW;
	}
	
	///Run the server to completion of all work
	void run(){
		auto newWork=generateNextWorkBlock();
		for(WorkItem& item : newWork)
			readyWork.emplace(std::move(item));
		startWatchdog();
		
		kj::NullDisposer disp;
		kj::Own<WorkServerImpl> wp(this,disp);
		capnp::EzRpcServer server(std::move(wp), kj::StringPtr("0.0.0.0:"+std::to_string(listeningPort)));
		doneFlag.promise.wait(server.getWaitScope());
	}
	
	///Set the port on which the server should listen for workers
	void setPort(unsigned int port){
		listeningPort=port;
	}
	///Get the port on which the server will listen for workers
	unsigned int getPort() const{ return listeningPort; }
	
	///Set how often the watchdog thread should wake up.
	///\pre may not be called after run() has been called
	void setWatchdogInterval(std::chrono::steady_clock::duration interval){
		watchdogInterval=interval;
	}
	///Get how often the watchdog thread will wake up.
	std::chrono::steady_clock::duration getWtachdogInterval() const{
		return watchdogInterval; 
	}
	
	///Set the maximum time a worker may be out of contact before it is assumed 
	///missing and its work reassigned.
	///\pre may not be called after run() has been called
	void setMaxWorkerSilenceTime(std::chrono::steady_clock::duration time){
		maxWorkerSilenceTime=time;
	}
	///Get the maximum time a worker may be out of contact before it is assumed 
	///missing and its work reassigned.
	std::chrono::steady_clock::duration getMaxWorkerSilenceTime() const{
		return maxWorkerSilenceTime; 
	}
	
private:
	///The port on which the server will listen for workers
	unsigned int listeningPort;
	libcuckoo::cuckoohash_map<uint64_t,Worker> workers;
	///Protects the readWork queue
	std::mutex readyWorkMutex;
	///Should only be manipulated while readyWorkMutex is held
	std::queue<WorkItem> readyWork;
	libcuckoo::cuckoohash_map<WorkItem,WorkItem,WorkItemHasher> outstandingWork;
	std::thread watchdog;
	kj::PromiseFulfillerPair<void> doneFlag;
	uint64_t token;
	///The period the watchdog thread should sleep between checks
	std::chrono::steady_clock::duration watchdogInterval;
	///The maximum time a worker can go without checking in before being 
	///considered missing and having its work reassigned
	std::chrono::steady_clock::duration maxWorkerSilenceTime;
	
	///Check for work items which have been outstanding for a long time, and
	///put them back in the ready pool to be assigned to more responsive workers
	void reassignMIAWork(){
		std::chrono::steady_clock::time_point now=std::chrono::steady_clock::now();
		std::unordered_set<WorkItem,WorkItemHasher> missingItems;
		for(const auto& item : outstandingWork.lock_table()){
			//just collect the ones we want to look at more closely, in order to minimize the time we keep the table locked
			//std::cout << "Work item " << item.second << " has been outstanding for " 
			//<< std::chrono::duration_cast<std::chrono::duration<double>>(now-item.second.startTime).count() 
			//<< " seconds" << std::endl; 
			if(now-item.second.startTime > maxWorkerSilenceTime)
				missingItems.insert(item.second);
		}
		for(WorkItem item : missingItems){
			Worker worker;
			bool reassign=false;
			try{
				worker=workers.find(item.assignedWorkerID);
			}catch(...){
				reassign=true; //if we can find the worker supposedly doing the work, we need to reassign
			}
			if(worker.id){
				if(now-worker.lastSeen > maxWorkerSilenceTime)
					reassign=true;
			}
			if(reassign){
				item.assignedWorkerID=0;
				std::cout << "Putting work item " << item << " back in the ready queue " << std::endl;
				std::lock_guard<std::mutex> lock(readyWorkMutex);
				readyWork.push(item);
				outstandingWork.erase(item);
			}
		}
	}
	
	///Clean up records of workers which have become unresponsive. 
	void reapMIAWorkers(){
		std::chrono::steady_clock::time_point now=std::chrono::steady_clock::now();
		std::unordered_set<Worker,WorkerHasher> missingWorkers;
		for(const auto& worker : workers.lock_table()){
			//just collect the ones we want to look at more closely, in order to minimize the time we keep the table locked
			if(now-worker.second.lastSeen > maxWorkerSilenceTime)
				missingWorkers.insert(worker.second);
		}
		for(Worker worker : missingWorkers){
			std::cout << "Reaping missing worker " << worker.id << std::endl;
			workers.erase(worker.id);
		}
	}
	
	///Run the watchdog thread which will handle missing work items/workers and
	///shut down the server when all work is done. 
	void startWatchdog(){
		watchdog=std::thread([this](){
			while(!done()){
				std::this_thread::sleep_for(watchdogInterval);
				std::cout << "Watchdog waking up" << std::endl;
				reassignMIAWork();
				reapMIAWorkers();
			}
			while(!workers.empty()){
				std::this_thread::sleep_for(watchdogInterval);
				reapMIAWorkers();
			}
			
			std::cout << "Shutting down" << std::endl;
			//this is the incredibly stupid set of hoops we have to jump 
			//through to notify the main thread to exit the 'event loop'
			capnp::EzRpcClient client("localhost:"+std::to_string(listeningPort));
			auto shutdown = client.getMain<WorkServer>().serverShutdownRequest();
			shutdown.setToken(token);
			try{
				shutdown.send().wait(client.getWaitScope());
			}catch(...){/*don't care*/}
		});
	}
	
protected:
	///Derived classes should override this function. 
	///The implementation of this function must be thread-safe.
	///\return Whether all necessary work items have been completed
	virtual bool done() const{
		return true;
	}
	
	///Process the result of a completed work item.
	///Derived classes should override this function.
	///\param item the original work item
	///\param value the value which was computed for the work item
	virtual void processResult(WorkItem item, double value){}
	
	///Generate more work items to be computed. This function will only be 
	///called if done() returns false.
	///Derived classes should override this function.
	///\return a collextion of new work items to be sent out to the workers. 
	virtual std::vector<WorkItem> generateNextWorkBlock(){
		return {};
	}
};

#endif //DRIVER_H
