#### Makefile to build the C++ wrapper example.
#### You can use this file as a skeleton to create the Makefile build your
#### software using the SLMSeg C++ wrapper.

R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 			$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 			$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)
MY_INCLUDE := 		/stornext/snfs130/NGIRD/scratch/phuang/repo/wgs_cnv/slmseg
MY_DIR :=			/stornext/snfs130/NGIRD/scratch/phuang/repo/wgs_cnv/connector

## if you need to set an rpath to R itself, also uncomment
RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CXX := 			$(shell $(R_HOME)/bin/R CMD config CXX)
CPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS := 		$(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LDLIBS := 		$(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

## source files
SRC :=			slmseg_connector.cpp

## program name
BIN :=			slmcc

all:
	$(CXX) -I$(MY_INCLUDE) -I$(MY_DIR) $(CPPFLAGS) $(CXXFLAGS) -Wall -O2 -g $(LDLIBS) -L$(MY_INCLUDE) -L$(MY_DIR) -lslm $(SRC) -fpic -shared -o libslmcc.so

clean:
	rm -rf *.o libslmcc.so

