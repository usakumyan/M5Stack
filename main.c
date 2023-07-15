#include <stdio.h>
#include <stdlib.h>

typedef float DTYPE;
typedef short IDX;
#include "fft.h"

#define FFTSIZE (256)
int main(int argc,char **argv[])
{
  DTYPE pRe[FFTSIZE];
  DTYPE pIm[FFTSIZE];

  for (IDX i=0;i<256;i++) {
    pRe[i]=pRe[i]*pWindow[i]; // pWindow[] is define in fft.h
  }
  RealFFTExec256(pRe,pIm);
  
  return 0;
}
