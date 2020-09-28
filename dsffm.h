
#ifndef DSFFM_H
#define DSFFM_H

#include "dsflib.h"
#include "m_pd.h"



void dsffm_run(dsf *x, OUTPRECISION *out1, OUTPRECISION *out2, int vector_size); 

void dsffm_set_fundamental(dsf *x, float frequency);


#endif
