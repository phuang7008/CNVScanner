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

#include "reports.h"
#include "terms.h"
#include "utility.h"

uint32_t processFile(char* chrom_id, char* file_name, Binned_Data_Wrapper *binned_data_wraper);

void khashInsertion(khash_t(khIntStr) *khash_in, uint32_t key, char* value);

void outputHashTable(khash_t(khIntStr) * khash_in, int type, User_Input *user_inputs, char *chrom_id);

void forDebug();

#endif
