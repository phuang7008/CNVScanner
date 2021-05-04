/*
 * =====================================================================================
 *
 *       Filename:  fileProcessing.h
 *
 *    Description:  to process text file such as mappability and gc content files
 *
 *        Version:  1.0
 *        Created:  04/21/2021 02:42:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston TX
 *
 * =====================================================================================
 */

#ifndef FILE_PROCESSING_H
#define FILE_PROCESSING_H

#include <stdlib.h>
#include <stdio.h>

#include "terms.h"
#include "utility.h"

void processFile(char* chrom_id, char* file_name, Binned_Data_Wrapper* binned_data_wrapper);

#endif
