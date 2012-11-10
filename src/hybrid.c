/*
  Hybrid Filter Bank
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static void fourChannelFiltering( FFTWFContext_t *fftctx,
                                  const float *pQmfReal,
                                  const float *pQmfImag,
                                  float **mHybridReal,
                                  float **mHybridImag,
                                  int chOffset)
{
  int i, k, n;
  int midTap = HYBRID_FILTER_DELAY;
  float real, imag;
  float cum[8];

  

  
   
  for(i = 0; i < QMF_TIME_SLOTS; i++) {

    
    cum[5] = cum[4] = 0;

     
    for(k = 0; k < 4; k++) {

      
      cum[5] -= p4_13[4*k] * pQmfReal[i+4*k];

      
      cum[4] += p4_13[4*k] * pQmfImag[i+4*k];
    }
    
    

    
    real = imag = 0;

     
    for(k = 0; k < 3; k++) {

      
      real += p4_13[4*k+3] * pQmfReal[i+4*k+3];
      imag += p4_13[4*k+3] * pQmfImag[i+4*k+3];
    }

      
    cum[6] = (imag + real) * 0.70710678118655f;
    cum[7] = (imag - real) * 0.70710678118655f;

     
    cum[0] = p4_13[midTap] * pQmfReal[i+midTap];
    cum[1] = p4_13[midTap] * pQmfImag[i+midTap];

    
    real = imag = 0;

     
    for(k = 0; k < 3; k++) {

      
      real += p4_13[4*k+1] * pQmfReal[i+4*k+1];
      imag += p4_13[4*k+1] * pQmfImag[i+4*k+1];
    }

      
    cum[2] = (real - imag ) * 0.70710678118655f;
    cum[3] = (real + imag ) * 0.70710678118655f;

    
    CFFTN(fftctx, cum, 4, 1);

     
    for(n = 0; n < 4; n++) {

      
      mHybridReal[i][n + chOffset] = cum[2*n];
      mHybridImag[i][n + chOffset] = cum[2*n+1];
    }
  }
  
}

static void eightChannelFiltering( FFTWFContext_t *fftctx,
                                   const float *pQmfReal,
                                   const float *pQmfImag,
                                   float **mHybridReal,
                                   float **mHybridImag)
{
  int i, n;
  float real, imag;
  int midTap = HYBRID_FILTER_DELAY;
  float cum[16];

  


  
   
  for(i = 0; i < QMF_TIME_SLOTS; i++) {

     
    real = p8_13[4]  * pQmfReal[4+i] +
           p8_13[12] * pQmfReal[12+i];

     
    imag = p8_13[4]  * pQmfImag[4+i] +
           p8_13[12] * pQmfImag[12+i];

      
    cum[4] =  (imag - real) * 0.70710678118655f;

      
    cum[5] = -(imag + real) * 0.70710678118655f;

     
    real = p8_13[3]  * pQmfReal[3+i] +
           p8_13[11] * pQmfReal[11+i];

     
    imag = p8_13[3]  * pQmfImag[3+i] +
           p8_13[11] * pQmfImag[11+i];

      
    cum[6] =   imag * 0.92387953251129f - real * 0.38268343236509f;

      
    cum[7] = -(imag * 0.38268343236509f + real * 0.92387953251129f);

       
    cum[9] = -( p8_13[2]  * pQmfReal[2+i] +
                p8_13[10] * pQmfReal[10+i] );

      
    cum[8] =    p8_13[2]  * pQmfImag[2+i] +
                p8_13[10] * pQmfImag[10+i];

     
    real = p8_13[1]  * pQmfReal[1+i] +
           p8_13[9] * pQmfReal[9+i];

     
    imag = p8_13[1]  * pQmfImag[1+i] +
           p8_13[9] * pQmfImag[9+i];

      
    cum[10] = imag * 0.92387953251129f + real * 0.38268343236509f;

      
    cum[11] = imag * 0.38268343236509f - real * 0.92387953251129f;

     
    real = p8_13[0]  * pQmfReal[i] +
           p8_13[8] * pQmfReal[8+i];

     
    imag = p8_13[0]  * pQmfImag[i] +
           p8_13[8] * pQmfImag[8+i];

      
    cum[12] = (imag + real) * 0.70710678118655f;

      
    cum[13] = (imag - real) * 0.70710678118655f;

    
    real = p8_13[7]  * pQmfReal[7+i];
    imag = p8_13[7]  * pQmfImag[7+i];

      
    cum[14] = imag * 0.38268343236509f + real * 0.92387953251129f;

      
    cum[15] = imag * 0.92387953251129f - real * 0.38268343236509f;

     
    cum[0] = p8_13[midTap]  * pQmfReal[midTap+i];
    cum[1] = p8_13[midTap]  * pQmfImag[midTap+i];

    
    real = p8_13[5]  * pQmfReal[5+i];
    imag = p8_13[5]  * pQmfImag[5+i];

      
    cum[2] = real * 0.92387953251129f - imag * 0.38268343236509f;

      
    cum[3] = real * 0.38268343236509f + imag * 0.92387953251129f;

    
    CFFTN(fftctx, cum, 8, 1);

     
    for(n = 0; n < 8; n++) {

      
      mHybridReal[i][n] = cum[2*n];
      mHybridImag[i][n] = cum[2*n+1];
    }
  }
  
}


/**************************************************************************/
/*!
  \brief   HybridAnalysis

  \return none.

*/
/**************************************************************************/
void
HybridAnalysis ( FFTWFContext_t *fftctx,
                 const float **mQmfReal,
                 const float **mQmfImag,
                 float **mHybridReal,
                 float **mHybridImag,
                 HANDLE_HYBRID hHybrid )   /*!< Handle to HYBRID struct. */

{
  int  n, band;
  HYBRID_RES hybridRes;
  int  chOffset = 0;

  

  

  
  for(band = 0; band < NO_QMF_BANDS_IN_HYBRID; band++) {

     
    hybridRes = (HYBRID_RES)aHybridResolution[band];

     
         
    memcpy(hHybrid->pWorkReal, hHybrid->mQmfBufferReal[band],
           QMF_BUFFER_MOVE * sizeof(float));

         
    memcpy(hHybrid->pWorkImag, hHybrid->mQmfBufferImag[band],
           QMF_BUFFER_MOVE * sizeof(float));

     
    for(n = 0; n < QMF_TIME_SLOTS; n++) {
      
      hHybrid->pWorkReal [QMF_BUFFER_MOVE + n] = mQmfReal [n] [band];
      hHybrid->pWorkImag [QMF_BUFFER_MOVE + n] = mQmfImag [n] [band];
    }

         
    memcpy(hHybrid->mQmfBufferReal[band], hHybrid->pWorkReal + QMF_TIME_SLOTS,
           QMF_BUFFER_MOVE * sizeof(float));

         
    memcpy(hHybrid->mQmfBufferImag[band], hHybrid->pWorkImag + QMF_TIME_SLOTS,
           QMF_BUFFER_MOVE * sizeof(float));

    switch(hybridRes) {
    case HYBRID_4_CPLX:
       

      /* filtering. */
      
      fourChannelFiltering( fftctx,
                            hHybrid->pWorkReal,
                            hHybrid->pWorkImag,
                            mHybridReal,
                            mHybridImag,
                            chOffset );

     break;
    case HYBRID_8_CPLX:
       

      /* filtering. */
      
      eightChannelFiltering( fftctx,
                             hHybrid->pWorkReal,
                             hHybrid->pWorkImag,
                             mHybridReal,
                             mHybridImag);
      break;
    default:
      assert(0);
    }

    
    chOffset += hybridRes;
  }
  
}

/**************************************************************************/
/*!
  \brief  HybridSynthesis

  \return none.

*/
/**************************************************************************/
void
HybridSynthesis ( const float **mHybridReal,
                  const float **mHybridImag,
                  float **mQmfReal,
                  float **mQmfImag,
                  HANDLE_HYBRID hHybrid )      /*!< Handle to HYBRID struct. */
{
  int  k, n, band;
  HYBRID_RES hybridRes;
  int  chOffset = 0;

  
  

  
  for(band = 0; band < NO_QMF_BANDS_IN_HYBRID; band++) {

     
    hybridRes = (HYBRID_RES)aHybridResolution[band];

    
    for(n = 0; n < QMF_TIME_SLOTS; n++) {

      
      mQmfReal [n] [band] = mQmfImag [n] [band] = 0;

       
      for(k = 0; k < hybridRes; k++) {

        
        mQmfReal [n] [band] += mHybridReal [n] [chOffset + k];
        mQmfImag [n] [band] += mHybridImag [n] [chOffset + k];
      }
      
    }

    
    chOffset += hybridRes;
  }
  
}


/**************************************************************************/
/*!
  \brief    CreateHybridFilterBank


  \return   errorCode, noError if successful.

*/
/**************************************************************************/
int
CreateHybridFilterBank ( HANDLE_HYBRID hs,         /*!< Handle to HYBRID struct. */
                         float **pPtr)
{
  int i;
  float *ptr = *pPtr;

  

    /* counting previous operations */

    
  hs->pWorkReal = ptr;ptr+= QMF_TIME_SLOTS + QMF_BUFFER_MOVE;
  hs->pWorkImag = ptr;ptr+= QMF_TIME_SLOTS + QMF_BUFFER_MOVE;

   
  hs->mQmfBufferReal = (float **)ptr;

  
  ptr+= NO_QMF_BANDS_IN_HYBRID * sizeof (float *)/sizeof(float);

   
  hs->mQmfBufferImag = (float **)ptr;

  
  ptr+= NO_QMF_BANDS_IN_HYBRID * sizeof (float *)/sizeof(float);

   /* hs->mQmfBufferReal[]
                  hs->mQmfBufferImag[]
               */
  
  for (i = 0; i < NO_QMF_BANDS_IN_HYBRID; i++) {

     
    hs->mQmfBufferReal[i] = ptr;ptr+= QMF_BUFFER_MOVE;

     
    hs->mQmfBufferImag[i] = ptr;ptr+= QMF_BUFFER_MOVE;

  }

  
  *pPtr=ptr;

  

  return 0;
}


/**************************************************************************/
/*!
  \brief

  \return   none

*/
/**************************************************************************/
void
DeleteHybridFilterBank ( HANDLE_HYBRID *phHybrid ) /*!< Pointer to handle to HYBRID struct. */
{
  
  
  return;
}
