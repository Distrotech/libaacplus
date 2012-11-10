/*
  MDCT transform
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define LS_TRANS ((FRAME_LEN_LONG-FRAME_LEN_SHORT)/2) /* 448 */



static void preModulationDCT(float *x,int m,const float *sineWindow)
{
  int i;
  float wre, wim, re1, re2, im1, im2;

  

   /* pointers for x[2*i],
                               x[m-1-2*i],
                               sineWindow[i*2],
                               sineWindow[m-1-2*i]
               */
  
  for(i=0;i<m/4;i++){

    
    re1 = x[2*i];
    im2 = x[2*i+1];
    re2 = x[m-2-2*i];
    im1 = x[m-1-2*i];

    
    wim = sineWindow[i*2];
    wre = sineWindow[m-1-2*i];

      
    x[2*i]   = im1*wim + re1* wre;

      
    x[2*i+1] = im1*wre - re1* wim;

    
    wim = sineWindow[m-2-2*i];
    wre = sineWindow[2*i+1];

      
    x[m-2-2*i] = im2*wim + re2*wre;

      
    x[m-1-2*i] = im2*wre - re2*wim;
  }

  
}


static void postModulationDCT(float *x,int m, const float *trigData, int step,int trigDataSize)
{
  int i;
  float wre, wim, re1, re2, im1, im2;
  const float *sinPtr = trigData;
  const float *cosPtr = trigData+trigDataSize;

  

   /* counting previous operations */

  
  wim = *sinPtr;
	wre = *cosPtr;

   /* pointers for x[2*i],
                               x[m-1-2*i]
               */
  
  for(i=0;i<m/4;i++){

    
    re1=x[2*i];
    im1=x[2*i+1];
    re2=x[m-2-2*i];
    im2=x[m-1-2*i];

      
    x[2*i]=(re1*wre + im1*wim);

      
    x[m-1-2*i]=(re1*wim - im1*wre);

    sinPtr+=step;
    cosPtr-=step;

    
    wim=*sinPtr;
    wre=*cosPtr;

      
    x[m-2-2*i] = (re2*wim + im2* wre);

      
    x[2*i+1]   = (re2*wre - im2* wim);
  }

  
}


static void mdct(FFTWFContext_t *fftctx,
                 float *dctdata,
                 const float *trigData,
                 const float *sineWindow,
                 int n,
                 int ld_n)
{
  
  

  
  preModulationDCT(dctdata,n,sineWindow);

   
  CFFTN(fftctx, dctdata,n/2,-1);

  assert (LD_FFT_TWIDDLE_TABLE_SIZE >= ld_n-1);

    
  postModulationDCT(dctdata, n, trigData, 1 << (LD_FFT_TWIDDLE_TABLE_SIZE - (ld_n - 1)), FFT_TWIDDLE_TABLE_SIZE);

  
}




/*!
    \brief 
    the mdct delay buffer has a size of 1600,
    so the calculation of LONG,START,STOP must be  spilt in two 
    passes with 1024 samples and a mid shift,
    the SHORT transforms can be completed in the delay buffer,
    and afterwards a shift
*/

static void shiftMdctDelayBuffer(
                                 float *mdctDelayBuffer, /*! start of mdct delay buffer */
                                 float *timeSignal,      /*! pointer to new time signal samples, interleaved */
                                 int chIncrement         /*! number of channels */
                                 )
{
  int i;

  

      
  memmove(mdctDelayBuffer,mdctDelayBuffer+FRAME_LEN_LONG,(BLOCK_SWITCHING_OFFSET-FRAME_LEN_LONG)*sizeof(float));

   /* pointers for mdctDelayBuffer[],
                               timeSignal[]
               */
  
  for(i=0;i<FRAME_LEN_LONG ;i++) {

    
    mdctDelayBuffer[BLOCK_SWITCHING_OFFSET-FRAME_LEN_LONG+i] = timeSignal[i*chIncrement];
  }

  
}

int Transform_Real(FFTWFContext_t *fftctx, float *mdctDelayBuffer,float *timeSignal,int chIncrement,float *realOut,int blockType)
{
  int i,w;
  float timeSignalSample;
  float ws1, ws2;

  float *dctIn;

  

  
  switch(blockType){


    case LONG_WINDOW:
      
      dctIn = realOut;

       /* pointer for mdctDelayBuffer[i],
                                  mdctDelayBuffer[FRAME_LEN_LONG-i-1],
                                  LongWindowKBD[i],
                                  LongWindowKBD[FRAME_LEN_LONG-i-1],
                                  dctIn[FRAME_LEN_LONG/2+i] 
                   */
      
      for(i=0;i<FRAME_LEN_LONG/2;i++){

        
        timeSignalSample = mdctDelayBuffer[i];

        
        ws1 = timeSignalSample * LongWindowKBD[i];

        
        timeSignalSample = mdctDelayBuffer[(FRAME_LEN_LONG-i-1)];

        
        ws2 = timeSignalSample * LongWindowKBD[FRAME_LEN_LONG-i-1];

         
        dctIn[FRAME_LEN_LONG/2+i] = ws1 - ws2;
      }

      
      shiftMdctDelayBuffer(mdctDelayBuffer,timeSignal,chIncrement);

       /* pointer for mdctDelayBuffer[i],
                                  mdctDelayBuffer[FRAME_LEN_LONG-i-1],
                                  LongWindowKBD[i],
                                  LongWindowKBD[FRAME_LEN_LONG-i-1],
                                  dctIn[FRAME_LEN_LONG/2+i] 
                   */
      
      for(i=0;i<FRAME_LEN_LONG/2;i++){
      
        
        timeSignalSample = mdctDelayBuffer[i];

        
        ws1 = timeSignalSample * LongWindowKBD[FRAME_LEN_LONG-i-1];

        
        timeSignalSample = mdctDelayBuffer[(FRAME_LEN_LONG-i-1)];

        
        ws2 = timeSignalSample * LongWindowKBD[i];

          
        dctIn[FRAME_LEN_LONG/2-i-1] = -(ws1 + ws2);

      }

      
      mdct(fftctx, dctIn, fftTwiddleTab, LongWindowSine, FRAME_LEN_LONG, 10);

    break;

    case START_WINDOW:
      
      dctIn = realOut;

     /* pointer for mdctDelayBuffer[i],
                                mdctDelayBuffer[FRAME_LEN_LONG-i-1],
                                LongWindowKBD[i],
                                LongWindowKBD[FRAME_LEN_LONG-i-1],
                                dctIn[FRAME_LEN_LONG/2+i] 
                   */
    
    for(i=0;i<FRAME_LEN_LONG/2;i++){

      
      timeSignalSample = mdctDelayBuffer[i];

      
      ws1 = timeSignalSample * LongWindowKBD[i];

      
      timeSignalSample = mdctDelayBuffer[(FRAME_LEN_LONG-i-1)];

      
      ws2 = timeSignalSample * LongWindowKBD[FRAME_LEN_LONG-i-1];

       
      dctIn[FRAME_LEN_LONG/2+i] = ws1 - ws2;
    }

    
    shiftMdctDelayBuffer(mdctDelayBuffer,timeSignal,chIncrement);
 

     /* pointer for mdctDelayBuffer[i],
                                dctIn[FRAME_LEN_LONG/2+i] 
                 */
    
    for(i=0;i<LS_TRANS;i++)
    {
      
      timeSignalSample = mdctDelayBuffer[i];
      ws1 = timeSignalSample;
      ws2 = 0.0f;

        
      dctIn[FRAME_LEN_LONG/2-i-1] = -(ws1 + ws2);
    }
    
     /* pointer for mdctDelayBuffer[i+LS_TRANS],
                                mdctDelayBuffer[FRAME_LEN_LONG-i-1-LS_TRANS],
                                ShortWindowSine[i],
                                ShortWindowSine[FRAME_LEN_SHORT-i-1],
                                dctIn[FRAME_LEN_LONG/2-i-1-LS_TRANS]
                 */
    
    for(i=0;i<FRAME_LEN_SHORT/2;i++){

      
      timeSignalSample= mdctDelayBuffer[i+LS_TRANS];

      
      ws1 = timeSignalSample * ShortWindowSine[FRAME_LEN_SHORT-i-1];

      
      timeSignalSample= mdctDelayBuffer[(FRAME_LEN_LONG-i-1-LS_TRANS)];

      
      ws2 = timeSignalSample * ShortWindowSine[i];

        
      dctIn[FRAME_LEN_LONG/2-i-1-LS_TRANS] = -(ws1 + ws2);
    }

    
    mdct(fftctx, dctIn, fftTwiddleTab, LongWindowSine, FRAME_LEN_LONG, 10);

    break;

    case STOP_WINDOW:

    
    dctIn = realOut;

     /* pointer for mdctDelayBuffer[FRAME_LEN_LONG-i-1],
                                dctIn[FRAME_LEN_LONG/2+i] 
                   */
    
    for(i=0;i<LS_TRANS;i++){

      
      ws1 = 0.0f;
      timeSignalSample= mdctDelayBuffer[(FRAME_LEN_LONG-i-1)];
      ws2 = timeSignalSample;

       
      dctIn[FRAME_LEN_LONG/2+i] = ws1 - ws2;
    }
    
     /* pointer for mdctDelayBuffer[i+LS_TRANS],
                                mdctDelayBuffer[FRAME_LEN_LONG-i-1-LS_TRANS],
                                ShortWindowSine[i],
                                ShortWindowSine[FRAME_LEN_SHORT-i-1],
                                dctIn[FRAME_LEN_LONG/2-i-1-LS_TRANS]
                 */
    
    for(i=0;i<FRAME_LEN_SHORT/2;i++){

      
      timeSignalSample= mdctDelayBuffer[(i+LS_TRANS)];

      
      ws1 = timeSignalSample * ShortWindowSine[i];

      
      timeSignalSample= mdctDelayBuffer[(FRAME_LEN_LONG-LS_TRANS-i-1)];

      
      ws2 = timeSignalSample * ShortWindowSine[FRAME_LEN_SHORT-i-1];

       
      dctIn[FRAME_LEN_LONG/2+i+LS_TRANS] = ws1 - ws2;
    }

    
    shiftMdctDelayBuffer(mdctDelayBuffer,timeSignal,chIncrement);
 

     /* pointer for mdctDelayBuffer[i],
                                mdctDelayBuffer[FRAME_LEN_LONG-i-1],
                                LongWindowKBD[i],
                                LongWindowKBD[FRAME_LEN_LONG-i-1],
                                dctIn[FRAME_LEN_LONG/2+i] 
                 */
    
    for(i=0;i<FRAME_LEN_LONG/2;i++){

      
      timeSignalSample= mdctDelayBuffer[i];

      
      ws1 = timeSignalSample * LongWindowKBD[FRAME_LEN_LONG-i-1];

      
      timeSignalSample= mdctDelayBuffer[(FRAME_LEN_LONG-i-1)];

      
      ws2 = timeSignalSample * LongWindowKBD[i];

        
      dctIn[FRAME_LEN_LONG/2-i-1] = -(ws1 + ws2);
    }

    
    mdct(fftctx, dctIn, fftTwiddleTab, LongWindowSine, FRAME_LEN_LONG, 10);

    break;

    case SHORT_WINDOW:

    
    for(w=0;w<TRANS_FAC;w++){

      
      dctIn = realOut+w*FRAME_LEN_SHORT;

       /* pointer for mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+i],
                                  mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT-i-1],
                                  mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT+i],
                                  mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT*2-i-1],
                                  ShortWindowSine[i],
                                  ShortWindowSine[FRAME_LEN_SHORT-i-1],
                                  dctIn[FRAME_LEN_SHORT/2+i],
                                  dctIn[FRAME_LEN_SHORT/2-i-1]
                   */
      
      for(i=0;i<FRAME_LEN_SHORT/2;i++){

        
        timeSignalSample= mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+i];

        
        ws1 = timeSignalSample * ShortWindowSine[i];

        
        timeSignalSample= mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT-i-1];

        
        ws2 = timeSignalSample * ShortWindowSine[FRAME_LEN_SHORT-i-1];

         
        dctIn[FRAME_LEN_SHORT/2+i] = ws1 - ws2;
        
        
        timeSignalSample= mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT+i];

        
        ws1 = timeSignalSample * ShortWindowSine[FRAME_LEN_SHORT-i-1];

        
        timeSignalSample= mdctDelayBuffer[TRANSFORM_OFFSET_SHORT+w*FRAME_LEN_SHORT+FRAME_LEN_SHORT*2-i-1];

        
        ws2 = timeSignalSample * ShortWindowSine[i];

          
        dctIn[FRAME_LEN_SHORT/2-i-1] = -(ws1 + ws2);
      }

      
      mdct(fftctx, dctIn, fftTwiddleTab, ShortWindowSine, FRAME_LEN_SHORT, 7);
    }

    
    shiftMdctDelayBuffer(mdctDelayBuffer,timeSignal,chIncrement);
    break;

  }
  

  return 0;
}
