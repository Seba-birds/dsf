

#ifndef DSF_H
#define DSF_H

#include <stdio.h>
#include <stdlib.h>
#include "dsf_defines.h"


typedef struct dsf
{
    OUTPRECISION phasor;
    OUTPRECISION increment;
} dsf;


//void dsf_run(dsf *x, OUTPRECISION *out, int vector_size);

dsf *dsf_new();

//void dsf_free(dsf *x);

#endif
