

#ifndef DSF_H
#define DSF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsf_defines.h"


typedef struct complex_nr
{
    INPRECISION im;
    INPRECISION re; 

} complex_nr;


typedef struct dsf
{
    INPRECISION sr_inv;
    INPRECISION frequency;
    INPRECISION distance;
    float sr;

    complex_nr *phasor_a;
    complex_nr *increment_a;

    complex_nr *phasor_b;
    complex_nr *increment_b;

} dsf;



void dsf_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size);

complex_nr *new_complex_nr();

dsf *dsf_new();

void dsf_free(dsf *x);

void set_increment_to_freq(complex_nr* increment, INPRECISION freq, INPRECISION sr_inv);

void set_phasor_to_argument(complex_nr *phasor, INPRECISION argument);

void normalize_phasor(complex_nr *phasor);

void multiply_complex(complex_nr *a, complex_nr *b, complex_nr *result);

void divide_phasors(complex_nr *numerator, complex_nr *denominator, complex_nr *result);

void divide_complex_by_length(complex_nr *numerator, 
        complex_nr *denominator, float len, complex_nr *result);

void divide_complex(complex_nr *numerator, complex_nr *denominator, complex_nr *result);

void geometric_series(complex_nr *a, complex_nr *b, complex_nr *result);

void dsf_set_frequency(dsf *x, float frequency);

void dsf_set_distance(dsf *x, float distance);

#endif
