/*
  CFFTN header file
*/

#ifndef __cfftn_h
#define __cfftn_h

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _FFTW3
#include <fftw3.h>
typedef struct {
        fftwf_plan plan4, plan8, plan64, plan512;
} FFTWFContext_t;
void init_plans(FFTWFContext_t *ctx);
void destroy_plans(FFTWFContext_t *ctx);
#else
#define init_plans(c);
#define destroy_plans(c);
typedef int FFTWFContext_t;
#endif


int CFFTN(FFTWFContext_t *ctx, float *afftData,int len, int isign);

#endif
