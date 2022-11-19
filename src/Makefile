IDIR=/opt/rh/devtoolset-8/root/usr/include
CODEDIR=/stornext/snfs130/NGIRD/scratch/phuang/repo/wgs_cnv/src
HTSDIR=/stornext/snfs130/NGIRD/scratch/phuang/htslib/
MYLIB=/stornext/snfs130/NGIRD/scratch/phuang/repo/common_libs
#CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
CC=gcc

# no compiler optimization
CFLAGS= -g -Wall -Wextra -fopenmp

# with compiler optimization
#CFLAGS= -g -Wall -Wextra -O0 -fopenmp
#CFLAGS= -g3 -gdwarf-2 -fvar-tracking -fvar-tracking-assignments -Wall -O0 -fopenmp

# define any directories containing header files other than /usr/bin
INCLUDES = -I$(IDIR) -I$(HTSDIR) -I$(MYLIB)

# define library paths in addition to /usr/lib
LFLAGS = -L$(HTSDIR) -L$(MYLIB)

# define any libraries to be link into executables:
LIBS =  $(HTSDIR)/libhts.a $(MYLIB)/libcov_tracking.a -lm -lz -lpthread -lrt -lssl -lcrypto -ldl -lcurl -llzma -lbz2 -lmcheck

# define dependent header files
_DEPS = $(wildcard *.h)
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C source files
SRCS = $(wildcard *.c)

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
#.c.o:		this one will give me the warning: 
#			ignoring prerequisites on suffix rule definition stackoverflow
#.c.o: $(DEPS)
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
