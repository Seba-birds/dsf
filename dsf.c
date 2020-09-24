
#include "dsf.h"



complex_nr *new_complex_nr()
{
    complex_nr *new_nr = (complex_nr *)malloc(sizeof(complex_nr));
    new_nr->im = 0.03141;
    new_nr->re = 0.999506;

    return new_nr;
}


dsf *dsf_new()
{
    dsf *x = (dsf *)malloc(sizeof(dsf));

    x->phasor_a = new_complex_nr();
    x->increment_a = new_complex_nr();
    x->phasor_b = new_complex_nr();
    x->increment_b = new_complex_nr();

    return x; 
}


void dsf_free(dsf *x)
{
    free(x->increment_a);
    free(x->phasor_a);
    free(x->increment_b);
    free(x->phasor_b);

    free(x);
} 


void set_increment_to_freq(complex_nr* increment, INPRECISION freq, INPRECISION sr_inv)
{
    INPRECISION argument = (freq * sr_inv) * PI; 
    set_phasor_to_argument(increment, argument);
}


void set_phasor_to_argument(complex_nr *phasor, INPRECISION argument)
{
    phasor->im = sin(argument);
    phasor->re = cos(argument);
} 


void multiply_complex(complex_nr *a, complex_nr *b, complex_nr *result)
{
    INPRECISION re = (a->re * b->re) - (a->im * b->im);
    INPRECISION im = (a->re * b->im) + (a->im * b->re);
    result->re = re;
    result->im = im; 
}


void divide_phasors(complex_nr *numerator, complex_nr *denominator, complex_nr *result)
{ 
    // because length of phasors is 1,
    // denominator becomes 1 after expansion with its conjugate
    complex_nr conjugate;
    conjugate.im = -denominator->im;
    conjugate.re = denominator->re; 
    multiply_complex(numerator, &conjugate, result); 
}


void divide_complex_by_length(complex_nr *numerator, 
        complex_nr *denominator, float len, complex_nr *result)
{ 
    // if we know the length of denominator
    // we can just use that to divide instead
    // of another complex multiplication
    complex_nr conjugate;
    conjugate.im = -denominator->im;
    conjugate.re = denominator->re;

    complex_nr expanded_num;
    multiply_complex(numerator, &conjugate, &expanded_num);

    float denom = len * len;

    result->im = expanded_num.im / denom;
    result->re = expanded_num.re / denom;
}


void divide_complex(complex_nr *numerator, complex_nr *denominator, complex_nr *result)
{ 
    complex_nr conjugate;
    conjugate.im = -denominator->im;
    conjugate.re = denominator->re;

    complex_nr expanded_num;
    multiply_complex(numerator, &conjugate, &expanded_num);

    complex_nr expanded_denom;
    multiply_complex(denominator, &conjugate, &expanded_denom); 

    result->im = expanded_num.im / expanded_denom.re;
    result->re = expanded_num.re / expanded_denom.re; 
}


void normalize_phasor(complex_nr *phasor)
{
    int too_big = phasor->im > 1.0;
    int too_small = phasor->im < -1.0;
    double in_range = !(too_big || too_small) * (double)phasor->im;
    
    double x = too_big - too_small + in_range; 


    double argument = asin(x);
    phasor->re = (INPRECISION) cos(argument); 
    phasor->im = (INPRECISION) x;
}


void dsf_set_frequency(dsf *x, float frequency)
{
    x->frequency = frequency;
    set_increment_to_freq(x->increment_a, x->frequency, x->sr_inv);
}

void dsf_set_distance(dsf *x, float distance)
{ 
    x->distance = distance;
    set_increment_to_freq(x->increment_b, x->distance, x->sr_inv);
}


int phasor_close_to_one(complex_nr* x)
{
    return x->re > 0.995; 
}


void geometric_series(complex_nr *a, complex_nr *b, complex_nr *result)
{
    complex_nr factor;

    if(b->re > 0.9995)
    {
        // following rule of l'hopital for "0/0":
        factor.im = 0;
        factor.re = 1;
    }
    else
    {
        complex_nr power;
        multiply_complex(b, b, &power); 

        // assuming w = 0.5, w^2 = 0.25 
        complex_nr numerator;
        numerator.im = -power.im * 0.25;
        numerator.re = 1 - power.re * 0.25;

        complex_nr denominator;
        denominator.im = -b->im;
        denominator.re = 1 - b->re; 

        divide_complex_by_length(&numerator, &denominator, 0.25, &factor); 
    }

    multiply_complex(a, &factor, result); 
}


INPRECISION norm_factor(INPRECISION len, int num_of_sines)
{
    INPRECISION norm_factor = 0;
    if(len > 0.9995 && len < 1.0005)
    {
        norm_factor = 1;
    }
    else
    { 
        norm_factor = (1 - pow(len, num_of_sines)) / (1 - len);
    }

    return norm_factor; 
}


void dsf_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size)
{
    for(int i = 0; i < vector_size; i++) { 

        multiply_complex(x->phasor_a, x->increment_a, x->phasor_a);
        multiply_complex(x->phasor_b, x->increment_b, x->phasor_b);

        complex_nr result; 
        geometric_series(x->phasor_a, x->phasor_b, &result); 

        // 0.6 to compensate for 1.5 as max signal strength
        out1[i] = result.re * norm_factor(0.5, 2);
        out2[i] = x->phasor_b->re;
    }

}

