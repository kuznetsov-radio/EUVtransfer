# -D LINUX in sources means "not Windows" (OpenMP vs MSVC PPL). It is set on
# Darwin and Linux only; Windows builds use EUVtransfer.vcxproj, not this makefile.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
BINDIR = ../binaries

CXX ?= g++
COMMON_COMPILE = -c -std=c++11 -O3 -fPIC

DEPS = EUVmain.h ExtMath.h IDLinterface.h

ifeq ($(UNAME_S),Darwin)
LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /opt/homebrew/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /usr/local/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
$(error Homebrew libomp not found (need omp.h). Install with: brew install libomp)
endif
ifeq ($(UNAME_M),arm64)
ONAME = $(BINDIR)/EUVtransfer_arm64.so
else
ONAME = $(BINDIR)/EUVtransfer_x86_64.so
endif
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -Xpreprocessor -fopenmp -I$(LIBOMP)/include
LINK_FLAGS = -shared
LDLIBS = -L$(LIBOMP)/lib -lomp -Wl,-rpath,$(LIBOMP)/lib
else ifeq ($(UNAME_S),Linux)
ONAME = $(BINDIR)/EUVtransfer.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LINK_FLAGS = -shared
LDLIBS = -fopenmp
else ifeq ($(OS),Windows_NT)
$(error Windows builds use EUVtransfer.vcxproj; this makefile supports Darwin and Linux only)
else
ONAME = $(BINDIR)/EUVtransfer.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LINK_FLAGS = -shared
LDLIBS = -fopenmp
endif

.DEFAULT_GOAL := $(ONAME)

OBJS = EUVmain.o ExtMath.o IDLinterface.o getparms.o

$(ONAME): $(OBJS) | $(BINDIR)
	$(CXX) $(LINK_FLAGS) -o $(ONAME) $(OBJS) $(LDLIBS)

EUVmain.o: EUVmain.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) EUVmain.cpp
ExtMath.o: ExtMath.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) ExtMath.cpp
IDLinterface.o: IDLinterface.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) IDLinterface.cpp
getparms.o: getparms.cpp
	$(CXX) $(CXXFLAGS) getparms.cpp

$(BINDIR):
	mkdir -p $@

.PHONY: install clean EUVtransfer

EUVtransfer: $(ONAME)

install: $(ONAME)

clean:
	rm -f *.o $(ONAME)
