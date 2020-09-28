

#ifndef DSF_H
#define DSF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsf_defines.h"
#include "m_pd.h"


typedef struct complex_nr
{
    double im;
    double re; 

} complex_nr;


typedef struct dsf
{
    INPRECISION sr_inv;
    INPRECISION frequency;
    INPRECISION distance;

    INPRECISION weight;
    int usr_num_of_sines;
    int num_of_sines;
    int norm_counter;
    float sr;

    complex_nr *phasor_a;
    complex_nr *increment_a;

    complex_nr *phasor_b;
    complex_nr *increment_b;

} dsf;


int min(int a, int b);

void dsf_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size);

complex_nr *new_complex_nr();

dsf *dsf_new();

void dsf_free(dsf *x);

void set_increment_to_freq(complex_nr* increment, INPRECISION freq, INPRECISION sr_inv);

void set_phasor_to_argument(complex_nr *phasor, double argument);

void normalize_phasor(complex_nr *phasor);

void adjust_phasor(complex_nr *phasor);

void multiply_complex(complex_nr *a, complex_nr *b, complex_nr *result);

void power_complex(complex_nr *x, int power, complex_nr *result);

void power_complex_naiv(complex_nr *x, int power, complex_nr *result); 

void divide_complex(complex_nr *numerator, complex_nr *denominator, complex_nr *result);

void one_minus_complex(complex_nr *x, complex_nr *result);

void scale_complex(complex_nr *x, double scalar, complex_nr *result);

void calculate_series(dsf *x, complex_nr *result, double *norm_factor);

void geometric_series_denominator(dsf *x, complex_nr *result);

void geometric_series_numerator(dsf *x, complex_nr *result);

void geometric_series_factor(dsf *x, complex_nr *result);

void geometric_series(dsf *x, complex_nr *result);

INPRECISION norm_factor(INPRECISION len, int num_of_sines);

void dsf_set_fundamental(dsf *x, float frequency);

void dsf_set_distance(dsf *x, float distance);

void dsf_set_weight(dsf *x, float weight);

void dsf_set_num_of_sines(dsf *x, int num_of_sines);

void dsf_set_frequency(dsf *x, float frequency);

void dsf_set_ratio(dsf *x, float ratio);


#endif
