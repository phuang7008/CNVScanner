/*
 * =====================================================================================
 *
 *       Filename:  storage.h
 *
 *    Description:  header file for some data structure for storages
 *
 *        Version:  1.0
 *        Created:  02/06/2024 03:46:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, pemhuang@gmail.com
 *        Company:  Burnaby, BC, Canada
 *
 * =====================================================================================
 */

#ifndef STORAGE_H
#define STORAGE_H

#define Array_SIZE 50000

typedef struct {
    uint32_t start;
    uint32_t end;
    double log2R_mean;
    double ave_coverage;
} Segment_Details;

typedef struct {
    char *chrom_id;
    uint32_t size;
    uint32_t capacity;
    Segment_Details *segments;
} Segment_Array;

#endif
