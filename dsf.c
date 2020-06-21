
#include "dsf.h"



dsf *dsf_new()
{
    dsf *x = (dsf *)malloc(sizeof(dsf));
    x->phasor = 0;
    x->increment = 0.01;
    return x; 
}

void dsf_free(dsf *x)
{
    free(x);
}

void dsf_set_increment(dsf *x, float inc)
{
    x->increment = inc; 
}

void dsf_run(dsf *x, OUTPRECISION *out, int vector_size)
{
    for(int i = 0; i < vector_size; i++) {
        x->phasor += x->increment;
        x->phasor -= (int)x->phasor;

        out[i] = x->phasor;
    }

}

