#random file ID, $(capnp id)
@0xecade5795dee9c66;

interface WorkServer {
	struct WorkerID{
		hash @0: UInt64;
	}
	struct RegistrationResponse{
		enum Type{
			acknowledged @0;
			shutDown @1;
		}
		type @0: Type;
	}

	registerWorker @0 (id:WorkerID) -> (resp:RegistrationResponse);

	struct WorkResponse{
		enum Type{
			workItem @0;
			wait @1;
			shutDown @2;
		}
		type @0: Type;
        parameters @1: List(Float64);
	}

	requestWork @1 (id:WorkerID) -> (resp:WorkResponse);

	struct HeartbeatResponse{
		enum Type{
			acknowledged @0;
			shutDown @1;
		}
		type @0: Type;
	}

	heartbeat @2 (id:WorkerID) -> (resp:HeartbeatResponse);

	struct WorkResult{
        parameters @0: List(Float64);
		union{
			result @1: Float64;
			errMsg @2: Text;
		}
	}
	struct WorkResultResponse{
		enum Type{
			acknowledged @0;
			shutDown @1;
		}
		type @0: Type;
	}

	workResult @3 (id:WorkerID, result:WorkResult) -> (resp:WorkResultResponse);

	workerShutdown @4 (id:WorkerID) -> (resp:Void);
	
	serverShutdown @5 (token:UInt64) -> (resp:Void); 
}
