
#CC=clang
#CXX=clang++
AR=ar
LD=clang++
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib
PREFIX=/usr/local/
CXXFLAGS+= -g -fPIC -O3 -mavx -std=c++17
LDFLAGS+= -lcfitsio -lfftw3 -lcapnp-rpc -lcapnp -lkj-async -lkj -pthread


.PHONY: all clean
all : modelEval driver

clean : 
	rm -rf build/*.o
	rm -rf build/distribute.capnp.h
	rm -rf build/distribute.capnp.c++

modelEval : build/modelEval.o build/geometry.o build/diskPhysics.o build/ParameterSet.o build/distribute.capnp.o
	$(CXX) -o modelEval $^ $(LDFLAGS)

sampler : build/sampler.o build/geometry.o build/diskPhysics.o build/ParameterSet.o
	$(CXX) -o sampler $^ $(LDFLAGS)

driver: build/driver.o build/ParameterSet.o  build/distribute.capnp.o
	$(CXX) -o driver $^ $(LDFLAGS)

build/modelEval.o : modelEval.cpp geometry.h grid.h image.h worker.h build/distribute.capnp.h
	$(CXX) $(CXXFLAGS) modelEval.cpp -c -o build/modelEval.o
build/sampler.o : sampler.cpp geometry.h grid.h image.h mcmc.h
	$(CXX) $(CXXFLAGS) sampler.cpp -c -o build/sampler.o
build/geometry.o : geometry.cpp diskPhysics.h geometry.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) geometry.cpp -c -o build/geometry.o
build/image.o : image.cpp image.h diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) image.cpp -c -o build/image.o
build/diskPhysics.o : diskPhysics.cpp diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) diskPhysics.cpp -c -o build/diskPhysics.o
build/ParameterSet.o : ParameterSet.cpp ParameterSet.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) ParameterSet.cpp -c -o build/ParameterSet.o
build/distribute.capnp.h : distribute.capnp Makefile
	capnp compile -oc++ distribute.capnp
	mv distribute.capnp.h build/

build/distribute.capnp.c++ : build/distribute.capnp.h Makefile
	mv distribute.capnp.c++ build/

build/distribute.capnp.o : build/distribute.capnp.c++ Makefile
	$(CXX) $(CXXFLAGS) -c build/distribute.capnp.c++ -o build/distribute.capnp.o

build/driver.o : driver.cpp driver.h build/distribute.capnp.h ParameterSet.h
	$(CXX) $(CXXFLAGS) -c driver.cpp -o build/driver.o

install :
	cp disk.so $(PREFIX)/include/
