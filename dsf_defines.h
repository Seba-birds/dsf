/**
 * @file stp_defines.h
 * @author Thomas Resch <br>
 * Audiocommunication Group, Technical University Berlin <br>
 * University of Applied Sciences Nordwestschweiz (FHNW), Music-Academy, Research and Development <br>
 * Preprocessor Instructions for the stp library
 */

#ifndef DSF_DEFINES_H
#define DSF_DEFINES_H

/** Data Type Macro for flexible input vector single or double floating point precision */
#define DSF_INPUTVECTOR_USE_FLOAT
#define DSF_OUTPUTVECTOR_USE_FLOAT

#ifdef DSF_INPUTVECTOR_USE_FLOAT
typedef float INPRECISION;
#else
typedef double INPRECISION;
#endif

/** Data Type Macro for flexible output vector single or double floating point precision */
#ifdef DSF_OUTPUTVECTOR_USE_FLOAT
typedef float OUTPRECISION;
#else
typedef double OUTPRECISION;
#endif

#endif /* DSF_DEFINES_H */
