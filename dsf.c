/**
  \file dsf.c
  \author Sebastian Zimmermann
  \date September 2020
  \brief dsf.c contains the main functionality of the dsf synthesis.

  \see James A. Moorer, "The Synthesis of Complex Audio Spectra by Means of Discrete Summation Formulas", 1976
  \see https://ccrma.stanford.edu/files/papers/stanm5.pdf
  \see Tim Stilson,  Julius Smith, "Alias-Free Digital Synthesis of Classic Analog Waveforms", 1996
  \see https://ccrma.stanford.edu/~stilti/papers/blit.pdf
  */

#include "dsf.h"

/**
  \brief get the smaller value of two integers
  \param a first value
  \param b second value
  \return smaller of a and b
  */
int min(int a, int b)
{
    int x = a < b ? a : b;
    return x; 
}


/**
  \brief initialize a complex_nr 

  Allocates necessary memory and sets the 
  complex number to a degree so that it could
  be used as an increment to generate an oscillation
  at 440 Hz. This enables fast sanity checking.

  \return Pointer to initialized complex_nr
  */
complex_nr *new_complex_nr()
{
    complex_nr *new_nr = (complex_nr *)malloc(sizeof(complex_nr));
    new_nr->im = 0.03141;
    new_nr->re = 0.999506;

    adjust_phasor(new_nr);

    return new_nr;
}


/**
  \brief initialize dsf memory

  The dsf struct holds all variables necessary for
  calculating the dsf.

  \return Pointer to allocated and initialized memory.
  \see dsf
  */ 
dsf *dsf_new()
{
    dsf *x = (dsf *)malloc(sizeof(dsf));

    x->phasor_a = new_complex_nr();
    x->increment_a = new_complex_nr();
    x->phasor_b = new_complex_nr();
    x->increment_b = new_complex_nr();


    x->weight = 0.5;
    x->usr_num_of_sines = 2;
    x->num_of_sines = 2;
    x->norm_counter = 0;
    x->frequency = 1;
    x->distance = 1;

    dsf_set_fundamental(x, 220);
    dsf_set_distance(x, 220);

    return x; 
}

/**
  \brief free dsf struct memory
  */
void dsf_free(dsf *x)
{
    free(x->increment_a);
    free(x->phasor_a);
    free(x->increment_b);
    free(x->phasor_b);

    free(x);
} 


/**
  \brief set angle of increment to frequency

  Rotating phasors are being rotated by multiplications with
  a complex number called increment. set_increment_to_freq()
  sets the angle of
  an increment with respect to the sample frequency
  so that those multiplications produce
  the correct oscillation frequency of the rotating phasor.

  \param increment The complex number that is to be adjusted
  for the given frequency
  \param freq The frequency to which the increment is to be set
  \param sr_inv The inverted sample frequency 

  */
void set_increment_to_freq(complex_nr* increment, INPRECISION freq, INPRECISION sr_inv)
{
    // multiply by 2.0 to get fraction of Nyquist frequency
    double argument = (freq * sr_inv * 2.0) * PI; 
    set_phasor_to_argument(increment, argument);
}

/**
  \brief set angle of a phasor

  \param phasor The phasor whose angle is to be set
  \param argument The angle to which to set the phasor to.
  The range \f$[0 \dots 2\pi]\f$ describes a full rotation.
  */
void set_phasor_to_argument(complex_nr *phasor, double argument)
{
    phasor->im = sin(argument);
    phasor->re = cos(argument);
} 

/**
  \brief multiplying two complex numbers

  \param a Factor a
  \param b Factor b
  \param result Memory to which the result is being written
  */
void multiply_complex(complex_nr *a, complex_nr *b, complex_nr *result)
{
    INPRECISION re = (a->re * b->re) - (a->im * b->im);
    INPRECISION im = (a->re * b->im) + (a->im * b->re);
    result->re = re;
    result->im = im; 
}

/**
  \brief divide complex numbers

  \param numerator Numerator of division
  \param denominator Denominator of division
  \param result Pointer to memory where result of division is written to
  */
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


/**
  \brief adjust phasor degeneration caused by digital noise

  Due to finite machine precision, over time phasors might
  diverge from the unit circle. In certain intervals
  adjust_phasor(complex_nr *phasor) is being called to 
  re-adjust the length of a phasor to 1.

  \param phasor Phasor to be adjusted.
  */
void adjust_phasor(complex_nr *phasor)
{
    double length = sqrt(phasor->im * phasor->im + phasor->re * phasor->re);
    double adjustment = 1.0 / length;
    phasor->im *= adjustment;
    phasor->re *= adjustment;
}

/**
  \brief set frequency of dsf

  Sets fundamental and distances relative to a new frequency.
  This keeps the overtone structure and corresponds with
  setting the oscillating frequency of a conventional oscillator.

  \param x dsf variables
  \param frequency Frequency to set the dsf to
  */ 
void dsf_set_frequency(dsf *x, float frequency)
{
    if(frequency != 0.0)
    { 
        float ratio = frequency / x->frequency;
        float new_distance = x->distance * ratio;
        dsf_set_fundamental(x, frequency);
        dsf_set_distance(x, new_distance); 
    }
}

/**
  \brief set ratio of overtones

  All overtones are spaced with equal distance.
  This function derives that distance as a ratio
  from the fundamental frequency.

  \param x dsf variables
  \param ratio Ratio of fundamental frequency to
  overtone spacing
  \see dsf_set_fundamental(dsf *x, float frequency)
  \see dsf_set_distance(dsf *x, float distance)
  */ 
void dsf_set_ratio(dsf *x, float ratio)
{
    float new_distance = ratio * x->frequency; 
    dsf_set_distance(x, new_distance); 
}


/**
  \brief set fundamental frequency

  This sets the fundamental frequency without changing
  the distance between the overtones.
  
  \param x dsf variables
  \param frequency New fundamental frequency
  \see dsf_set_ratio(dsf *x, float ratio)
  \see dsf_set_frequency(dsf *x, float frequency)
  */
void dsf_set_fundamental(dsf *x, float frequency)
{
    x->frequency = frequency; 
    // new frequency: adjust number of generated sines
    // to avoid aliasing by increased frequency
    dsf_set_num_of_sines(x, x->usr_num_of_sines); 

    set_increment_to_freq(x->increment_a, x->frequency, x->sr_inv);

    adjust_phasor(x->phasor_a);
    adjust_phasor(x->increment_a);
}


/**
  \brief set distance between overtones

  This sets the distance between overtones
  without changing the fundamental frequency.  
  
  \param x dsf variables
  \param distance New distance between overtones
  \see dsf_set_ratio(dsf *x, float ratio)
  \see dsf_set_frequency(dsf *x, float frequency)
  */
void dsf_set_distance(dsf *x, float distance)
{ 
    x->distance = distance; 
    // new distance: adjust number of generated sines
    // to avoid aliasing by increased gaps between sines
    dsf_set_num_of_sines(x, x->usr_num_of_sines); 

    set_increment_to_freq(x->increment_b, x->distance, x->sr_inv);

    adjust_phasor(x->phasor_b);
    adjust_phasor(x->increment_b);
}


/**
  \brief set diminishing factor of overtones

  The overtone series is tailing of by a weight
  factor \f$ w \f$, where the nth overtone is scaled by
  \f$ w^n \f$.

  \param x dsf variables
  \param weight New value for \f$ w \f$
  */ 
void dsf_set_weight(dsf *x, float weight)
{
    x->weight = weight; 
}

/**
  \brief set total number of partials to generate

  This function takes user input x and 
  calculates the maximum number of 
  partials y that can be produced without aliasing
  and sets the number of actually produced partials
  to the lower of x and y.

  \param x dsf variables
  \param num_of_sines Number of partials to produce
  */ 
void dsf_set_num_of_sines(dsf *x, int num_of_sines)
{
    if(num_of_sines > 0)
    {
        x->usr_num_of_sines = num_of_sines;
        if(x->distance)
        {
            int max_num_of_sines = (int) (((x->sr / 2.0) - x->frequency) / x->distance) + 1;
            x->num_of_sines = min(max_num_of_sines, num_of_sines);
        }
        else
        {
            // x->distance == 0: only one frequency
            x->num_of_sines = 1; 
        }

        char msg[90];
        sprintf(msg, "sr: %f, usr: %d, actual: %d\n", x->sr, x->usr_num_of_sines, x->num_of_sines);
        //post(msg);
    } 
} 


void power_complex(complex_nr *x, int power, complex_nr *result)
{
    // stop recursion:
    if(power <= 1)
    {
        result->im = x->im;
        result->re = x->re; 
    }
    else
    {
        complex_nr intermediate;
        power_complex(x, power/2, &intermediate);

        multiply_complex(&intermediate, &intermediate, result);

        if(power % 2)
        {
            // uneven exponent: multiply result one more time with x
            multiply_complex(x, result, result); 
        } 
    } 
}



void power_complex_naiv(complex_nr *x, int power, complex_nr *result)
{
    result->im = x->im;
    result->re = x->re;
    for(int i = 1; i < power; i++)
    {
        multiply_complex(x, result, result);
    } 
}



void calculate_series(dsf *x, complex_nr *result, double *norm_factor)
{
    INPRECISION powered_weight;

    complex_nr running_sum;
    complex_nr powered_b;

    // w^0 * b^0
    powered_b.im = 0.0;
    powered_b.re = 1.0; 

    running_sum.im = 0.0;
    running_sum.re = 1.0;

    powered_weight = 1.0;
    *norm_factor = 1.0;

    // increasing powers
    for(int i = 1; i < x->num_of_sines; i++)
    {
        // b^i
        multiply_complex(x->phasor_b, &powered_b, &powered_b); 
        // w^i
        powered_weight *= x->weight;

        // track norm factor
        *norm_factor += powered_weight;

        // add w^i * b^i to sum
        running_sum.im += powered_weight * powered_b.im;
        running_sum.re += powered_weight * powered_b.re; 
    } 

    multiply_complex(&running_sum, x->phasor_a, result); 
    // track length for phasor a:
    *norm_factor += 1.0; 
}



/*
 *

 SUM(0..k) a * (wb)^k = a * (1 - (wb)^k) / (1 - wb)

 (wb)^k = w^k * b^k

 *
 * 
 */


void one_minus_complex(complex_nr *x, complex_nr *result)
{
    result->im = 0.0 - x->im;
    result->re = 1.0 - x->re; 
}


void scale_complex(complex_nr *x, double scalar, complex_nr *result)
{ 
    result->im = scalar * x->im;
    result->re = scalar * x->re; 
}


void geometric_series_denominator(dsf *x, complex_nr *result)
{
    // calculating (1 - w * b)
    // where w is x->weight
    // and b is current position of distance phasor x->phasor_b

    double w = x->weight;

    complex_nr wb;
    scale_complex(x->phasor_b, w, &wb); 

    one_minus_complex(&wb, result); 
}


void geometric_series_numerator(dsf *x, complex_nr *result)
{ 
    // calculating (1 - (w * b)^n) = (1 - (w^n * b^n))
    // where w is x->weight,
    // b is current position of distance phasor x->phasor_b,
    // and n is number of sines (or number of sines + 1?)
    int n = x->num_of_sines;
    double wn = pow(x->weight, n);

    complex_nr bn;
    power_complex(x->phasor_b, n, &bn);

    complex_nr wnbn;
    scale_complex(&bn, wn, &wnbn);

    one_minus_complex(&wnbn, result); 
}


void geometric_series_factor(dsf *x, complex_nr *result)
{
    // calculating (1 - (wb)^k) / (1 - wb)
    complex_nr numerator;
    geometric_series_numerator(x, &numerator);

    complex_nr denominator;
    geometric_series_denominator(x, &denominator); 

    divide_complex(&numerator, &denominator, result);
}


void geometric_series(dsf *x, complex_nr *result)
{
    complex_nr factor;

    if(x->phasor_b->re == 1 && x->weight == 1)
    {
        // following rule of l'hopital for "0/0":
        factor.im = 0;
        factor.re = 1;
    }
    else
    {
        geometric_series_factor(x, &factor); 
    }

    multiply_complex(x->phasor_a, &factor, result); 
}


INPRECISION norm_factor(INPRECISION len, int num_of_sines)
{
    INPRECISION norm_factor = 0;
    if(len == 1.0)
    {
        // l'hopital "0/0":
        norm_factor = 1.0 / num_of_sines;
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
        geometric_series(x, &result); 

        out1[i] = result.re * norm_factor(x->weight, x->num_of_sines); 
        out2[i] = x->phasor_b->re;
    }
    
}

