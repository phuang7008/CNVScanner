noinst_LIBRARIES = libcov_tracking.a
libcov_tracking_a_SOURCES=coverage_tracking.c coverage_tracking.h data_structure.h utility.c utility.h

libcov_tracking_a_CFLAGS = -Wall -g -fopenmp
libcov_tracking_a_CPPFLAGS=-I$(includedir) -I$(abs_top_builddir)/htslib
libcov_tracking_a_LIBADD = $(abs_top_builddir)/htslib/lib/libhts.a

