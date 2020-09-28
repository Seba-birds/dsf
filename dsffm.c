
#include "dsffm.h"
#include "dsflib.h"

void dsffm_set_fundamental(dsf *x, float frequency)
{
    set_increment_to_freq(x->increment_a, frequency, x->sr_inv); 
}

void dsffm_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size)
{
    for(int i = 0; i < vector_size; i++) { 

        multiply_complex(x->phasor_a, x->increment_a, x->phasor_a); 

        out1[i] = x->phasor_a->im;
        out2[i] = x->phasor_a->re;
    } 
}
