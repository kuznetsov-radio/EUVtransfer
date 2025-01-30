EUVtransfer		:	EUVmain.o ExtMath.o IDLinterface.o getparms.o
				g++ $^ -shared -fopenmp -o EUVtransfer.so
EUVmain.o		:	EUVmain.cpp EUVmain.h ExtMath.h IDLinterface.h
				g++ -c -std=c++11 -O3 -fPIC -fopenmp -D LINUX EUVmain.cpp
ExtMath.o		:	ExtMath.cpp EUVmain.h ExtMath.h IDLinterface.h
				g++ -c -std=c++11 -O3 -fPIC -fopenmp -D LINUX ExtMath.cpp
IDLinterface.o		:	IDLinterface.cpp EUVmain.h ExtMath.h IDLinterface.h
				g++ -c -std=c++11 -O3 -fPIC -fopenmp -D LINUX IDLinterface.cpp
getparms.o		:	getparms.cpp
				g++ -c -std=c++11 -O3 -fPIC -fopenmp -D LINUX getparms.cpp