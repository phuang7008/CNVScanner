bin_PROGRAMS=cnvscanner
cnvscanner_SOURCES=analysis.c analysis.h breakpoints.c breakpoints.h calculate_stdev.c calculate_stdev.h excluded_regions.c excluded_regions.h fileProcessing.c fileProcessing.h improperly_paired_reads.c improperly_paired_reads.h main.c reports.c reports.h stats.c stats.h std_cnv_calling.c std_cnv_calling.h terms.h user_inputs.c user_inputs.h utils.c utils.h

cnvscanner_LDFLAGS = -lm -lz -lpthread -llzma -lcrypto -lbz2 -lcurl

cnvscanner_CFLAGS = -Wall -g -fopenmp
cnvscanner_CPPFLAGS=-I$(includedir) -I$(abs_top_builddir)/htslib -I$(abs_top_builddir)/lib
cnvscanner_LDADD = $(abs_top_builddir)/htslib/lib/libhts.a $(abs_top_builddir)/lib/libcov_tracking.a
