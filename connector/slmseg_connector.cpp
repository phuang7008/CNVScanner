/*****************************************************************************
This is the connector between C++ and C
*/

#include <math.h>

#include "slmseg_connector.h"
#include "slmseg.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void slmseg_call(unsigned int chr, char* file_in, char* out_file_name, Segment_Array* segment_array, double omega, double eta, double stepeta, double fw) {
    unsigned int tmp_chr = 0;
    if (strcmp(segment_array->chrom_id, "X") == 0 || strcmp(segment_array->chrom_id, "chrX") == 0) {
        tmp_chr = 23;
    } else if (strcmp(segment_array->chrom_id, "Y") == 0 || strcmp(segment_array->chrom_id, "chrY") == 0) {
        tmp_chr = 24;
    } else {
        tmp_chr = atoi(segment_array->chrom_id);
    }

    if (tmp_chr != chr) {
        fprintf(stderr, "Something went wrong as chr %d and chrom_id %s not the same!\n", chr, segment_array->chrom_id);
        exit(EXIT_FAILURE);
    }

    static SLMSeg segdata(omega, eta, stepeta, fw);

    fprintf(stderr, "Loading file...\n");
    segdata.load_signal_file(file_in, chr);

    fprintf(stderr, "Performing analysis...\n");
    segdata.SLM();

    fprintf(stderr, "Data output for chromosome %d:\n", chr);

    FILE * out_fhd;
    out_fhd = fopen(out_file_name, "w");

    std::vector<double> data_seg = segdata.data_seg();
    std::vector<unsigned int> data_pos = segdata.pos_data();

    segment_array->size = 0;
    segment_array->capacity = Array_SIZE;
    segment_array->segments = (Segment_Details*) calloc(Array_SIZE, sizeof(Segment_Details));

    double prev_data = 0.0;
    unsigned int prev_pos = 0;
    bool beginning_flag = true;
    uint32_t counter=0;

    for (unsigned int i = 0; i < data_seg.size(); i++)
    {
	    //fprintf(out_fhd, "%d\t%" PRIu32 "\t%.2f\n",chr, data_pos.at(i), data_seg.at(i));
	    if (beginning_flag || (int(data_seg.at(i) * 10000) != int(prev_data * 10000)) ) {
	        if (!beginning_flag) {
		        fprintf(out_fhd, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.2f\n", segment_array->chrom_id, prev_pos, data_pos.at(i), prev_data);
		        //fprintf(out_fhd, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.2f\n",chr, prev_pos, data_pos.at(i), prev_data);
                segment_array->segments[counter].start = prev_pos;
                segment_array->segments[counter].end   = data_pos.at(i);
                segment_array->segments[counter].log2R_mean    = prev_data;
                segment_array->segments[counter].ave_coverage  = exp2(prev_data);
                counter++;
                segment_array->size++;

                if (segment_array->capacity == segment_array->size + 5) {
                    // dynamic expand the array size
                    //
                    segment_array->capacity = segment_array->capacity + 10000;
                    segment_array->segments = (Segment_Details*) realloc(segment_array->segments, \
                            segment_array->capacity * sizeof(Segment_Details));
                    if (segment_array->segments == NULL) {
                        fprintf(stderr, "ERROR: Dynamic Memory allocation for segment_array->segments failed!\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
	        prev_data = data_seg.at(i);
	        prev_pos  = data_pos.at(i);

	        if (beginning_flag) beginning_flag = false;
	    }
    }

    // output the last line
    fprintf(out_fhd, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.2f\n", segment_array->chrom_id, prev_pos, data_pos.at(data_seg.size()-1), prev_data);
    segment_array->segments[counter].start = prev_pos;
    segment_array->segments[counter].end   = data_pos.at(data_seg.size()-1);
    segment_array->segments[counter].log2R_mean  = prev_data;
    segment_array->segments[counter].ave_coverage  = exp2(prev_data);
    counter++;
    segment_array->size++;

    // remove all elements in the vector data of SLMSeg
    //
    segdata.erase_all_vector_data();

    fclose(out_fhd);
}

#ifdef __cplusplus
}
#endif
