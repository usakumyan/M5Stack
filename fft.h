/* fft.h */
#ifndef HEADER_FFT
#include "utl.h"
#define HEADER_FFT
extern const DTYPE pWindow[256];
#endif

/* FFT for real sequence                          */
/* input                                          */
/*   pRe[128] : Real part of input sequence       */
/*   pIm[128] : Imaginary part of input sequence  */
/* output                                         */
/*   pRe[65] : Real part of output sequence       */
/*   pIm[65] : Imaginary part of output sequence  */
int FFTExec128(DTYPE *pRe,DTYPE *pIm);
/* FFT for real sequence                          */
/* input                                          */
/*   pRe[256] : Real part of input sequence       */
/*   pIm[256] : don't care                        */
/* output                                         */
/*   pRe[129] : Real part of output sequence      */
/*   pIm[129] : Imaginary part of output sequence */
void RealFFTExec256(DTYPE *pRe,DTYPE *pIm);

