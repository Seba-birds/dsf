
#include "dsf.h"



complex_nr *new_complex_nr()
{
    complex_nr *new_nr = (complex_nr *)malloc(sizeof(complex_nr));
    new_nr->im = 0.03141;
    new_nr->re = 0.999506;

    adjust_phasor(new_nr);

    return new_nr;
}


dsf *dsf_new()
{
    dsf *x = (dsf *)malloc(sizeof(dsf));

    x->phasor_a = new_complex_nr();
    x->increment_a = new_complex_nr();
    x->phasor_b = new_complex_nr();
    x->increment_b = new_complex_nr();


    x->weight = 0.5;
    x->num_of_sines = 2;
    x->norm_counter = 0;

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
    double argument = (freq * sr_inv) * PI; 
    set_phasor_to_argument(increment, argument);
}


void set_phasor_to_argument(complex_nr *phasor, double argument)
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


void adjust_phasor(complex_nr *phasor)
{
    double length = sqrt(phasor->im * phasor->im + phasor->re * phasor->re);
    double adjustment = 1.0 / length;
    phasor->im *= adjustment;
    phasor->re *= adjustment;
}


void dsf_set_frequency(dsf *x, float frequency)
{
    x->frequency = frequency;
    set_increment_to_freq(x->increment_a, x->frequency, x->sr_inv);

    adjust_phasor(x->phasor_a);
    adjust_phasor(x->increment_a);
}

void dsf_set_distance(dsf *x, float distance)
{ 
    x->distance = distance;
    set_increment_to_freq(x->increment_b, x->distance, x->sr_inv);

    adjust_phasor(x->phasor_b);
    adjust_phasor(x->increment_b);
}


void dsf_set_weight(dsf *x, float weight)
{
    x->weight = weight; 
}


void dsf_set_num_of_sines(dsf *x, int num_of_sines)
{
    if(num_of_sines > 0)
    {
        x->num_of_sines = num_of_sines;
    } 
}




int phasor_close_to_one(complex_nr* x)
{
    return x->re > 0.995; 
}


void power_complex(complex_nr *x, int power, complex_nr *result)
{
    result->im = x->im;
    result->re = x->re;
    for(int i = 1; i < power; i++)
    {
        multiply_complex(x, result, result);
    } 
}


//void geometric_series(complex_nr *a, complex_nr *b, complex_nr *result)
void geometric_series(dsf *x, complex_nr *result)
{
    complex_nr factor;

    if(x->phasor_b->re == 1)
    {
        // following rule of l'hopital for "0/0":
        factor.im = 0;
        factor.re = 1;
    }
    else
    {
        complex_nr power;
        //multiply_complex(x->phasor_b, x->phasor_b, &power); 
        power_complex(x->phasor_b, x->num_of_sines, &power);

        // assuming w = 0.5, w^2 = 0.25 
        double scale = pow(x->weight, x->num_of_sines);
        //double scale = pow(0.5, 2);
        complex_nr numerator;
        numerator.im = -power.im * scale;
        numerator.re = 1 - power.re * scale;

        complex_nr denominator;
        denominator.im = -x->phasor_b->im;
        denominator.re = 1 - x->phasor_b->re; 

        divide_complex_by_length(&numerator, &denominator, 1.0 - x->weight, &factor); 
    }

    multiply_complex(x->phasor_a, &factor, result); 
}


INPRECISION norm_factor(INPRECISION len, int num_of_sines)
{
    INPRECISION norm_factor = 0;
    if(len == 1.0)
    {
        norm_factor = 1.0;
    }
    else
    { 
        norm_factor = (1.0 - len) / (1.0 - pow(len, num_of_sines));
    }

    return norm_factor; 
}


void dsf_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size)
{
    for(int i = 0; i < vector_size; i++) { 

        multiply_complex(x->phasor_a, x->increment_a, x->phasor_a);
        multiply_complex(x->phasor_b, x->increment_b, x->phasor_b);

        complex_nr result; 
        if(x->num_of_sines > 1)
        geometric_series(x, &result); 


        out1[i] = result.re * norm_factor(x->weight, x->num_of_sines);
        out2[i] = x->phasor_b->re;
    }
    
}

