//#define M5StackCore2
#ifdef M5StackCore2
#include <M5Core2.h>
#include <drive/i2s.h>
#include "SD.h"
#else
#include <stdio.h>
#include <stdlib.h>
#endif

typedef float DTYPE;
typedef short IDX;
#include "fft.h"

/* FFT for real sequence                          */
/* input                                          */
/*   pRe[128] : Real part of input sequence       */
/*   pIm[128] : Imaginary part of input sequence  */
/* output                                         */
/*   pRe[65] : Real part of output sequence       */
/*   pIm[65] : Imaginary part of output sequence  */
int FFTExec128(DTYPE *pRe,DTYPE *pIm)
{
  const IDX nRealFFTSize=256;
  const IDX nFFTSize=128;
  IDX       j, lmx, li;
  DTYPE     *_pX=NULL;
  DTYPE     *_pY=NULL;
  DTYPE     *pSin=NULL;
  DTYPE     *pCos=NULL;
  IDX       lf, lix;
  IDX       mv2, mm1;
  DTYPE     t1, t2;

  lf  = nRealFFTSize/nFFTSize;
  lmx = nFFTSize;
  for(;;) {
    lix =  lmx;
    lmx /= 2;
    if (lmx <= 1) break;
    pSin = _sintbl;
    pCos = _sintbl + nRealFFTSize/4;
    for (j=0; j<lmx; j++) {
      _pX = &pRe[j];
      _pY = &pIm[j];
      for (li=lix; li<=nFFTSize ; li+=lix) {
	t1           =  *(_pX) - *(_pX + lmx);
	t2           =  *(_pY) - *(_pY + lmx);
	*(_pX)       += *(_pX + lmx);
	*(_pY)       += *(_pY + lmx);
	*(_pX + lmx) =  *pCos * t1 + *pSin * t2;
	*(_pY + lmx) =  *pCos * t2 - *pSin * t1;
	_pX          += lix;
	_pY          += lix;
      }
      pSin += lf;
      pCos += lf;
    }
    lf += lf;
  }
  _pX = pRe;
  _pY = pIm;
  for (li=nFFTSize/2; li--; _pX+=2,_pY+=2) {
    t1         =  *(_pX) - *(_pX + 1);
    t2         =  *(_pY) - *(_pY + 1);
    *(_pX)     += *(_pX + 1);
    *(_pY)     += *(_pY + 1);
    *(_pX + 1) = t1;
    *(_pY + 1) = t2;
  }
  j   = 0;
  _pX = pRe;
  _pY = pIm;
  mv2 = nFFTSize / 2;
  mm1 = nFFTSize - 1;
  for (lmx=0; lmx<mm1; lmx++) {
    if ((li=lmx-j)<0) {
      t1          = *(_pX);
      t2          = *(_pY);
      *(_pX)      = *(_pX + li);
      *(_pY)      = *(_pY + li);
      *(_pX + li) = t1;
      *(_pY + li) = t2;
    }
    li = mv2;
    while (li<=j) {
      j  -= li;
      li /= 2;
    }
    j   += li;
    _pX =  pRe + j;
    _pY =  pIm + j;
  }

  return(0);
}
/* FFT for real sequence                          */
/* input                                          */
/*   pRe[256] : Real part of input sequence       */
/*   pIm[256] : don't care                        */
/* output                                         */
/*   pRe[129] : Real part of output sequence      */
/*   pIm[129] : Imaginary part of output sequence */
void RealFFTExec256(DTYPE *pRe,DTYPE *pIm)
{
  const IDX nRealFFTSize=256;
  const IDX nFFTSize=128;
  DTYPE _r[128];
  DTYPE _i[128];
  DTYPE dTempRe1=0;
  DTYPE dTempIm1=0;
  DTYPE dTempRe2=0;
  DTYPE dTempIm2=0;
  DTYPE dTempRe3=0;
  DTYPE dTempIm3=0;
  DTYPE *pSin=NULL;
  DTYPE *pCos=NULL;
  
  pSin = _sintbl;
  pCos = _sintbl + nRealFFTSize/4;

  for (IDX i=0;i<nFFTSize;i++) {
    _r[i]=pRe[2*i];
    _i[i]=pRe[2*i+1];
  }
  FFTExec128(_r,_i);
  
  for (IDX k=1;k<nFFTSize/2;k++) {
    dTempRe1=0.5+0.5*pSin[k];
    dTempIm1=    0.5*pCos[k];
    dTempRe2=_r[k]-_r[nFFTSize-k];
    dTempIm2=_i[k]+_i[nFFTSize-k];
    dTempRe3=dTempRe1*dTempRe2-dTempIm1*dTempIm2;
    dTempIm3=dTempRe1*dTempIm2+dTempRe2*dTempIm1;
    pRe[k         ]=_r[k         ]-dTempRe3;
    pIm[k         ]=_i[k         ]-dTempIm3;
    pRe[nFFTSize-k]=_r[nFFTSize-k]+dTempRe3;
    pIm[nFFTSize-k]=_i[nFFTSize-k]-dTempIm3;
  }
  pRe[0         ]= _r[0]+_i[0];
  pIm[0         ]= 0;
  pRe[nFFTSize/2]= _r[nFFTSize/2];
  pIm[nFFTSize/2]=-_i[nFFTSize/2];  
  pRe[nFFTSize  ]= _r[0]-_i[0];
  pIm[nFFTSize  ]= 0;  
}
