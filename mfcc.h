/* mfcc.h */

#ifndef HEADER_MFCC
#define HEADER_MFCC
#endif

#include "utl.h"

#define NCHANNELS (20)
void DCT_IIE_20(DTYPE *pIn,DTYPE *pOut);
void fbank(DTYPE *pEnergy,DTYPE *pFBank);
