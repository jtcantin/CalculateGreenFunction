MKLROOT = /global/software/intel/composerxe-2011.0.084/mkl

CC = icc 

#use openmp, use c++11 standard, -DNDEBUG: to disable bound check in eigen to speed up
CFLAGS = -O3 -openmp -Wall -DNDEBUG


#include eigen, mkl, boost
CINCLUDE = -I/home/pxiang/eigen -I/home/pxiang/boost_1_55_0 -I/home/pxiang/gtest-1.7.0/include  -I$(MKLROOT)/include 


#link with mkl and google test (static library)
FLAGSLIB = /home/pxiang/gtest-1.7.0/libgtest.a -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm



SOURCES = $(wildcard *.cpp) $(wildcard */*.cpp) 


green: $(SOURCES) Makefile
	$(CC) $(CFLAGS) $(CINCLUDE)   $(SOURCES)  $(FLAGSLIB) -o $@
# Tab before $(CC)

clean:
	rm -f green
# Tab before "rm"
