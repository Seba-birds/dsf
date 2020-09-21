
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
    x->phasor = new_complex_nr();
    x->increment = new_complex_nr();
    return x; 
}


void dsf_free(dsf *x)
{
    free(x->increment);
    free(x->phasor);
    free(x);
} 


void set_increment_to_freq(dsf *x)
{
    INPRECISION argument = (x->frequency * x->sr_inv) * PI; 
    set_phasor_to_argument(x->increment, argument);
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
    set_increment_to_freq(x);
}

void dsf_run(dsf *x, OUTPRECISION *out, int vector_size)
{
    for(int i = 0; i < vector_size; i++) {
        multiply_complex(x->phasor, x->increment, x->phasor);

        out[i] = x->phasor->re;
    }

}

