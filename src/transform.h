/*
  MDCT transform
*/
#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

int Transform_Real(FFTWFContext_t *fftctx,
                   float *mdctDelayBuffer,
                   float *timeSignal,
                   int chIncrement,
                   float *realOut,
                   int windowSequence);
#endif
