/*****************************************************************************
This is the connector between C++ and C
*/

#ifndef SLMSEG_CONNECTOR_H
#define SLMSEG_CONNECTOR_H

#include <string.h>
#include <inttypes.h>   // for PRIu32 and PRIu64 

#include "storage.h"

#ifdef __cplusplus
extern "C" {
#endif


// Use type 'unsigned int' for chromosome ID, for X will be 23, while Y 24
//
void slmseg_call(unsigned int chr, char* file_in, char* out_file_name, Segment_Array* segments, double omega, double eta, double stepeta, double fw);

#ifdef __cplusplus
}
#endif

#endif
