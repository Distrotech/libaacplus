/*
  IIR resampler tool box
*/

#include <memory.h>
#include <assert.h>
#include <stdio.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define MAX_COEFF 32


struct IIR_PARAM{
  const float *coeffIIRa;
  const float *coeffIIRb;
  int noOffCoeffs;
  int transitionFactor;
  int delay;
};

static const float set1_a[] = {
  0.004959855f,   0.025814206f,   0.080964205f,   0.182462303f,
  0.322109621f,   0.462709927f,   0.552160404f,   0.552160404f,
  0.462709927f,   0.322109621f,   0.182462303f,   0.080964205f,
  0.025814206f,   0.004959855f
};

static const float set1_b[] = {
  0.0f,          -1.038537170f,   2.627279635f,  -1.609574122f,
  2.205922661f,  -0.751928739f,   0.787128253f,  -0.105573173f,
  0.131638380f,   0.003884641f,   0.010544805f,   0.001232040f,
  0.000320798f,   0.000023031f
};
  
  
static struct IIR_PARAM const set1 = {
  set1_a,
  set1_b,
  14,
  218,
  6
};



/*
  Reset downsampler instance and clear delay lines     
  
  returns status
*/
int 
InitIIR21_Resampler(IIR21_RESAMPLER *ReSampler)
     
{
  

   
  ReSampler->iirFilter.ptr   =   0;
 
   
  ReSampler->iirFilter.coeffIIRa = set1.coeffIIRa;
  ReSampler->iirFilter.coeffIIRb = set1.coeffIIRb;
  ReSampler->iirFilter.noOffCoeffs = set1.noOffCoeffs;
  ReSampler->delay=set1.delay;

  assert(ReSampler->iirFilter.noOffCoeffs <= BUFFER_SIZE);

   
  ReSampler->ratio = 2;

   
  ReSampler->pending = 1;

  

  return 1;
}



/* 
   NOTE: enabling NEWIIR would save us some wMOPS, but result in 1 LSB diffs, it is worth a CR
*/

#ifdef NEWIIR
static float
AdvanceMAFilter( IIR_FILTER *iirFilter
                 )
{
  float y;
  int j;
  int ptr = iirFilter->ptr;
  int i = ptr + (BUFFER_SIZE-1);

  

   /*  --> ptr isn't needed */  /* counting previous operations */

   
  y = (iirFilter->coeffIIRa[0] * iirFilter->ring_buf_2[i & (BUFFER_SIZE-1)]);

   /* iirFilter->noOffCoeffs
                  iirFilter->coeffIIRa[]
                  iirFilter->ring_buf_2[]
               */
  
  for (j=1; j<iirFilter->noOffCoeffs; j++) {
    i--;

    
    y += (iirFilter->coeffIIRa[j] * iirFilter->ring_buf_2[i & (BUFFER_SIZE-1)]);
  }

  

  return y;
}


static void
AdvanceARFilter( IIR_FILTER *iirFilter,
                 float input
                 )

{
  int j;
  float y;
  int ptr = iirFilter->ptr;
  int i = ptr + (BUFFER_SIZE-1);

  

     /* counting previous operations */

    
  y = input + (iirFilter->coeffIIRb[1] * (-iirFilter->ring_buf_2[i & (BUFFER_SIZE-1)]));

   /* iirFilter->noOffCoeffs
                  iirFilter->coeffIIRb[]
                  iirFilter->ring_buf_2[i]
                  iirFilter->ring_buf_2[ptr]
               */

  
  for (j=2; j<iirFilter->noOffCoeffs; j++) {
    i--;

     
    y += (iirFilter->coeffIIRb[j] * (-iirFilter->ring_buf_2[i & (BUFFER_SIZE-1)]));
  }

  
  iirFilter->ring_buf_2[ptr] = y;

  /* pointer update */
  iirFilter->ptr = (ptr+1) & (BUFFER_SIZE-1);

  

}
#else  //NEWIIR



/*
  faster simple folding operation    

  returns filtered value
*/


static float
AdvanceIIRFilter(IIR_FILTER *iirFilter,
                 float input
                 )

{
  float y = 0.0f;
  int j = 0;
  int i;

  

   /* counting previous operations */

   
  iirFilter->ring_buf_1[iirFilter->ptr] = input;

   /* pointer for iirFilter->ring_buf_1,
                              iirFilter->ring_buf_2,
                              iirFilter->coeffIIRa,
                              iirFilter->coeffIIRb
               */
   
  for (i = iirFilter->ptr; i > iirFilter->ptr - iirFilter->noOffCoeffs; i--, j++) {
     
    y += iirFilter->coeffIIRa[j] * iirFilter->ring_buf_1[i & (BUFFER_SIZE - 1)] - iirFilter->coeffIIRb[j] * iirFilter->ring_buf_2[i & (BUFFER_SIZE - 1)];
  }

  
  iirFilter->ring_buf_2[(iirFilter->ptr) & (BUFFER_SIZE - 1)] = y;

  iirFilter->ptr = ++iirFilter->ptr & (BUFFER_SIZE - 1);

  

  return y;
}
#endif  //NEWIIR






/*
  Downsample numInSamples of type short

  returns  success of operation
*/

int 
IIR21_Downsample(IIR21_RESAMPLER *DownSampler,
                 float *inSamples,
                 int numInSamples,
                 int inStride,
                 float *outSamples,
                 int *numOutSamples,
                 int outStride
               )
{
  int i;
  *numOutSamples=0;

  

   /* counting previous operations */

   /* pointer for inSamples[],
                              outSamples[]
               */
   
  for(i=0;i<numInSamples;i++){
    float iirOut;

#ifdef NEWIIR
    
    AdvanceARFilter(&(DownSampler->iirFilter), inSamples[i*inStride]);
#else
    
    iirOut = AdvanceIIRFilter(&(DownSampler->iirFilter), inSamples[i*inStride]);
#endif

    
    DownSampler->pending++;

     
    if(DownSampler->pending == DownSampler->ratio){

#ifdef NEWIIR
      
      outSamples[(*numOutSamples)*outStride] = AdvanceMAFilter(&(DownSampler->iirFilter));;
#else
      
      outSamples[(*numOutSamples)*outStride] = iirOut;
#endif

      (*numOutSamples)++;

      
      DownSampler->pending=0;
    }
  }

  

  return 1;
}


int 
IIR21_Upsample(IIR21_RESAMPLER *UpSampler,
               float *inSamples,
               int numInSamples,
               int inStride,
               float *outSamples,
               int *numOutSamples,
               int outStride
               )
{
  int i,k;
  int idxOut=0;

  

   /* counting previous operations */

   /* pointer for inSamples[],
                              outSamples[]
               */
   
  for(i=0;i<numInSamples;i++){

      
    outSamples[idxOut] = AdvanceIIRFilter(&(UpSampler->iirFilter), inSamples[i*inStride] * UpSampler->ratio);

    idxOut += outStride;

    
    for (k=0; k<UpSampler->ratio-1; k++) {

       
      outSamples[idxOut] = AdvanceIIRFilter(&(UpSampler->iirFilter), 0.0f);

      idxOut += outStride;
    }      
  }

   
  *numOutSamples=numInSamples*UpSampler->ratio;

  

  return 1;
}
