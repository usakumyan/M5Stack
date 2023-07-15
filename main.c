#include <stdio.h>
#include <stdlib.h>

#include "utl.h"
#include "fft.h"
#include "mfcc.h"

#define FFTSIZE (256)
#define NCHANNELS (20)
int main(int argc,char **argv[])
{
  DTYPE pRe[FFTSIZE];
  DTYPE pIm[FFTSIZE];
  DTYPE pFBank[NCHANNELS];
  DTYPE pMFCC[NCHANNELS];
  

  for (IDX i=0;i<256;i++) {
    pRe[i]=pRe[i]*pWindow[i]; // pWindow[] is define in fft.h
  }
  RealFFTExec256(pRe,pIm);
  for (IDX i=0;i<FFTSIZE/2;i++) {
    pRe[i]=pRe[i]*pRe[i]+pIm[i]*pIm[i];
  }
  // calculate Filterbank output (FBank)
  fbank(pRe,pFBank);
  // calculate MFCCs
  DCT_IIE_20(pFBank,pMFCC);
  
  return 0;
}
