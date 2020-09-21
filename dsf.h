

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
    INPRECISION frequency;
    INPRECISION sr_inv;
    float sr;

    complex_nr *phasor;
    complex_nr *increment;

} dsf;



void dsf_run(dsf *x, OUTPRECISION *out, int vector_size);

complex_nr *new_complex_nr();

dsf *dsf_new();

void dsf_free(dsf *x);

void set_increment_to_freq(dsf *x);

void set_phasor_to_argument(complex_nr *phasor, INPRECISION argument);

void normalize_phasor(complex_nr *phasor);

void multiply_complex(complex_nr *a, complex_nr *b, complex_nr *result);

void dsf_set_frequency(dsf *x, float frequency);

#endif
