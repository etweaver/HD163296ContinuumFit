#ifndef WORKER_H
#define WORKER_H

#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include <capnp/ez-rpc.h>

#include "build/distribute.capnp.h"

class Worker{
public:
	Worker(std::string serverAddress):
	id(generateWorkerID()),
	serverAddress(serverAddress),
	client(serverAddress),
	workServer(client.getMain<WorkServer>()),
	waitScope(client.getWaitScope()),
	doneFlag(false),
	heartbeatInterval(std::chrono::seconds(1)),
	maxRuntime(std::chrono::seconds(0)),
	startTime(std::chrono::steady_clock::now())
	{
		std::cout << "Worker ID is " << id << std::endl;
		registerWorker();
	}
	
	~Worker(){
		doneFlag=true;
		auto shutDown = workServer.workerShutdownRequest();
		shutDown.getId().setHash(id);
		try{
			auto shutDownResult = shutDown.send().wait(waitScope);
		}catch(...){
			//too late to care about errors now
		}
		if(heartbeat.joinable())
			heartbeat.join();
	}
	
	uint64_t getID() const{ return id; }
	
	///Set the maximum time for which this worker should run
	///\param time the maximum time for the worker ot run, or zero to indicate 
	///            no limit
	void setMaxRunTime(std::chrono::steady_clock::duration time){
		maxRuntime=time;
	}
	///Set how often to send heartbeat messages to the server
	void setHeartbeatInterval(std::chrono::steady_clock::duration interval){
		heartbeatInterval=interval;
	}
	
	void registerWorker(){
		auto registration = workServer.registerWorkerRequest();
		registration.getId().setHash(id);
		auto registerResult = registration.send().wait(waitScope);
	}
	
	void run(){
		startHeartbeat();
		std::vector<double> parameters;
		while(true){
			auto getWork=workServer.requestWorkRequest();
			getWork.getId().setHash(id);
			auto workResponse = getWork.send().wait(waitScope);
			switch(workResponse.getResp().getType()){
				case WorkServer::WorkResponse::Type::WORK_ITEM:
				{
					std::size_t nParam=workResponse.getResp().getParameters().size();
					try{
						parameters.resize(nParam);
						for(std::size_t i=0; i<nParam; i++)
							parameters[i]=workResponse.getResp().getParameters()[i];
						double calcResult=compute(parameters);
						
						auto result = workServer.workResultRequest();
						result.getId().setHash(id);
						auto resultParams=result.getResult().initParameters(nParam);
						for(std::size_t i=0; i<nParam; i++)
							resultParams.set(i,workResponse.getResp().getParameters()[i]);
						result.getResult().setResult(calcResult);
						auto resultResult = result.send().wait(waitScope);
					}catch(std::exception& ex){
						auto result = workServer.workResultRequest();
						result.getId().setHash(id);
						auto resultParams=result.getResult().initParameters(nParam);
						for(std::size_t i=0; i<nParam; i++)
							resultParams.set(i,workResponse.getResp().getParameters()[i]);
						result.getResult().setErrMsg(ex.what());
						auto resultResult = result.send().wait(waitScope);
					}
					break;
				}
				case WorkServer::WorkResponse::Type::WAIT:
				{
					std::this_thread::sleep_for(std::chrono::seconds(1));
					break;
				}
				case WorkServer::WorkResponse::Type::SHUT_DOWN:
					std::cout << "Got shut down command from server" << std::endl;
					doneFlag=true;
					return;
			}
			if(maxRuntime>std::chrono::steady_clock::duration(0)){
				if(std::chrono::steady_clock::now()-startTime > maxRuntime){
					doneFlag=true;
					return;
				}
			}
		}
	}
	
protected:
	virtual double compute(const std::vector<double>& parameters){
		return 0;
	}
	
private:
	uint64_t id;
	std::string serverAddress;
	capnp::EzRpcClient client;
	WorkServer::Client workServer;
	kj::WaitScope& waitScope;
	std::atomic<bool> doneFlag;
	std::thread heartbeat;
	std::chrono::steady_clock::duration heartbeatInterval;
	std::chrono::steady_clock::duration maxRuntime;
	std::chrono::steady_clock::time_point startTime;
	
	void startHeartbeat(){
		heartbeat=std::thread([this](){
			while(!doneFlag){
				std::this_thread::sleep_for(heartbeatInterval);
				if(doneFlag)
					return;
				//Create and destroy the client on every iteration to avoid keeping 
				//a connection open and doubling the number of fds the server must 
				//handle. This assumedly has some cost of its own, but as long as 
				//it's small compared to heartbeatInterval we win out. 
				capnp::EzRpcClient client(serverAddress);
				auto workServer=client.getMain<WorkServer>();
				auto& waitScope=client.getWaitScope();
				auto heartbeat = workServer.heartbeatRequest();
				heartbeat.getId().setHash(id);
				try{
					auto result = heartbeat.send().wait(waitScope);
					if(result.getResp().getType()==WorkServer::HeartbeatResponse::Type::SHUT_DOWN)
						return;
				}catch(...){
					//TODO: possibly take action if server does not respond?
				}
			}
		});
	}
	
	static uint64_t generateWorkerID(){
		uint64_t hash=0;
		std::ifstream rfile("/dev/urandom");
		while(hash==0){
			rfile.read((char*)&hash,sizeof(hash));
			if(!rfile)
				throw std::runtime_error("Unable to read random data from /dev/urandom");
		}
		return hash;
	}
};

#endif //WORKER_H
