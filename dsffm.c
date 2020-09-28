
#include "dsffm.h"
#include "dsflib.h"
#include "dsf.h"

void dsffm_set_fundamental(dsf *x, float frequency)
{
    set_increment_to_freq(x->increment_a, frequency, x->sr_inv); 
}

void dsffm_run(dsf *x, 
        INPRECISION *in1, INPRECISION *in2, INPRECISION *in3,
        OUTPRECISION *out1, OUTPRECISION *out2, 
        int vector_size)
{
    for(int i = 0; i < vector_size; i++) { 

        // process all inputs:
        float freq = mtof(in1[i]); 
        dsf_set_frequency(x, freq); 
        dsf_set_ratio(x, in2[i] + 1.0); 
        float weight = clip(0.0, 1.5, x->weight + in3[i]); 
        dsf_set_weight(x, weight);


        multiply_complex(x->phasor_a, x->increment_a, x->phasor_a); 
        multiply_complex(x->phasor_b, x->increment_b, x->phasor_b); 

        complex_nr result; 
        geometric_series(x, &result); 

        out1[i] = result.re * norm_factor(x->weight, x->num_of_sines); 
        out2[i] = x->phasor_a->re;
    } 
}
