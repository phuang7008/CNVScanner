IDIR=/opt/rh/devtoolset-8/root/usr/include
CODEDIR=/stornext/snfs130/NGIRD/scratch/phuang/repo/scandium/src/sequenceStats
HTSDIR=/stornext/snfs130/NGIRD/scratch/phuang/htslib/v1.10.2/htslib
MYSQLDIR=/usr/include/mysql
MYSQLLIB=/usr/lib64
MYLIB=/stornext/snfs130/NGIRD/scratch/phuang/repo/common_libs
CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
#CC=/hgsc_software/gcc/gcc-7.3.0/bin/gcc
#CFLAGS= -g -Wall -O2 -fopenmp
CFLAGS= -g -Wall -Wextra -O0 -fopenmp
#CFLAGS= -g3 -gdwarf-2 -fvar-tracking -fvar-tracking-assignments -Wall -O0 -fopenmp

# define any directories containing header files other than /usr/bin
INCLUDES= -I$(IDIR) -I$(HTSDIR) -I$(MYSQLDIR) -I$(MYLIB)

# define library paths in addition to /usr/lib
LFLAGS = -L$(HTSDIR) -L$(MYSQLLIB) -L$(MYLIB)

# define any libraries to be link into executables:
#LIBS =  $(HTSDIR)/libhts.a -lm -lz -lpthread -lmysqlclient -lrt -lssl -lcrypto -ldl
LIBS =  $(HTSDIR)/libhts.a $(MYLIB)/libcov_tracking.a -lm -lz -lpthread -lmysqlclient_r -lrt -lssl -lcrypto -ldl -lcurl -llzma -lbz2

# define dependent header files
_DEPS = analysis.h excluded_regions.h stats.h terms.h user_inputs.h utility.h data_structure.h utils.h
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C ssource files
SRCS = main.c analysis.c excluded_regions.c stats.c user_inputs.c utils.c

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = wgs_cnv

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:	$(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
#$(CODEDIR)/%.o: %.c $(DEPS)
#.c.o:
.c.o: $(DEPS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
