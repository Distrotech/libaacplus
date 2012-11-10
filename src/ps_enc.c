/*
  parametric stereo encoding
*/
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


#define FLOAT_EPSILON                   0.0001f
#define INV_LOG2                        1.4427f
#define NO_QMFSB_IN_SUBQMF              ( 3 )
#define SUBQMF_BINS_ENERGY              ( 8 )
#define SUBQMF_GROUPS_MIX               ( 16 )

#define MAX_MAGNITUDE_CORR              ( 2.0f )

#define PS_MODE_LOW_FREQ_RES_IID_ICC    ( 0x00020000 )

static void downmixToMono( HANDLE_PS_ENC   pms,
                           float            **qmfLeftReal,
                           float            **qmfLeftImag,
                           float            **qmfRightReal,
                           float            **qmfRightImag );

int GetPsMode(int bitRate)
{
  int psMode = 0;

  

   /* counting previous operation */

   
  if(bitRate < 21000)
  {
    
    psMode =
      PS_MODE_LOW_FREQ_RES_IID_ICC;
  }
  else {
   
  if(bitRate < 128000)
  {
    
    psMode = 0;
  }
  }

  

  return psMode;
}

int
CreatePsEnc(SBRRam_t *sbrram,
            HANDLE_PS_ENC h_ps_e,
            int psMode)
{
  int i;
  int err;

  float *ptr1, *ptr2, *ptr3, *ptr4;

  

  
  ptr1 = &sbrram->sbr_envYBuffer[QMF_TIME_SLOTS * QMF_CHANNELS];
  ptr2 = &sbrram->PsBuf2[0];
  ptr3 = &sbrram->PsBuf3[0];
  ptr4 = &sbrram->PsBuf4[0];

  
  if (h_ps_e == NULL)
  {
    
    return 1;
  }

   
  h_ps_e->psMode = psMode;

   
  h_ps_e->bPrevZeroIid = 0;
  h_ps_e->bPrevZeroIcc = 0;

     
  h_ps_e->bHiFreqResIidIcc = ( ( psMode & PS_MODE_LOW_FREQ_RES_IID_ICC ) != 0 )? 0: 1;

    
  h_ps_e->iidIccBins = ( h_ps_e->bHiFreqResIidIcc )? NO_IID_BINS: NO_LOW_RES_IID_BINS;

   
  h_ps_e->aaaICCDataBuffer = (float **)ptr1;

  
  ptr1+= NO_BINS * sizeof (float *)/sizeof(float);

   
  h_ps_e->aaaIIDDataBuffer = (float **)ptr1;

  
  ptr1+= NO_BINS * sizeof (float *)/sizeof(float);

   /* h_ps_e->aaaICCDataBuffer[][]
                  h_ps_e->aaaIIDDataBuffer[][]
                 */
  
  for (i=0 ; i < NO_BINS ; i++) {

    
    h_ps_e->aaaICCDataBuffer[i] = ptr1;

    
    ptr1+= SYSTEMLOOKAHEAD + 1;

    
    h_ps_e->aaaIIDDataBuffer[i] = ptr1;

    
    ptr1+= SYSTEMLOOKAHEAD + 1;
  }

   
  h_ps_e->hHybridLeft = &h_ps_e->hybridLeft;
  h_ps_e->hHybridRight = &h_ps_e->hybridRight;


  
  err = CreateHybridFilterBank ( h_ps_e->hHybridLeft,
                                 &ptr4);

  
  if(err != 0) {

       
      DeletePsEnc (&h_ps_e);

      
      return 1; /* "Failed to create hybrid filterbank." */
  }

  
  err = CreateHybridFilterBank ( h_ps_e->hHybridRight,
                                 &ptr4);

  
  if(err != 0) {
       
      DeletePsEnc (&h_ps_e);

      
      return 1; /* "Failed to create hybrid filterbank." */
  }

   /* h_ps_e->mHybridRealLeft[]
                  h_ps_e->mHybridImagLeft[]
                  h_ps_e->mHybridRealRight[]
                  h_ps_e->mHybridImagRight[]
                */
  
  for (i=0 ; i < NO_SUBSAMPLES ; i++) {
    
    h_ps_e->mHybridRealLeft[i] = ptr3;

    
    ptr3+=NO_HYBRID_BANDS;

    
    h_ps_e->mHybridImagLeft[i] = ptr3;

    
    ptr3+=NO_HYBRID_BANDS;

    
    h_ps_e->mHybridRealRight[i] = ptr1;

    
    ptr1+=NO_HYBRID_BANDS;

    
    h_ps_e->mHybridImagRight[i] = ptr1;

    
    ptr1+=NO_HYBRID_BANDS;
  }

    
  h_ps_e->tempQmfLeftReal = (float **)ptr1;
  ptr1+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

    
  h_ps_e->tempQmfLeftImag= (float **)ptr1;
  ptr1+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

    
  h_ps_e->histQmfLeftReal = (float **)ptr2;
  ptr2+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

    
  h_ps_e->histQmfLeftImag= (float **)ptr2;
  ptr2+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

    
  h_ps_e->histQmfRightReal = (float **)ptr2;
  ptr2+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

    
  h_ps_e->histQmfRightImag= (float **)ptr2;
  ptr2+= HYBRID_FILTER_DELAY * sizeof (float *)/sizeof(float);

   /* h_ps_e->tempQmfLeftReal[]
                  h_ps_e->tempQmfLeftImag[]
                  h_ps_e->histQmfLeftReal[]
                  h_ps_e->histQmfLeftImag[]
                  h_ps_e->histQmfRightReal[]
                  h_ps_e->histQmfRightImag[]
               */
  
  for (i=0 ; i < HYBRID_FILTER_DELAY ; i++) {
    
    h_ps_e->tempQmfLeftReal[i] = ptr1;

    
    ptr1+=NO_QMF_BANDS;

    
    h_ps_e->tempQmfLeftImag[i] = ptr1;

    
    ptr1+=NO_QMF_BANDS;

    
    h_ps_e->histQmfLeftReal[i] = ptr2;

    
    ptr2+=NO_QMF_BANDS;

    
    h_ps_e->histQmfLeftImag[i] = ptr2;

    
    ptr2+=NO_QMF_BANDS;

    
    h_ps_e->histQmfRightReal[i] = ptr2;

    
    ptr2+=NO_QMF_BANDS;

    
    h_ps_e->histQmfRightImag[i] = ptr2;

    
    ptr2+=NO_QMF_BANDS;

  }

   
  if ( ( h_ps_e->histQmfLeftReal     == NULL ) || ( h_ps_e->histQmfLeftImag     == NULL ) ||
       ( h_ps_e->histQmfRightReal    == NULL ) || ( h_ps_e->histQmfRightImag    == NULL ) )
  {
     
    DeletePsEnc( &h_ps_e );

    
    return 1; /* "Memory allocation in function CreatePsEnc() failed." */
  }

   
  if (h_ps_e->psBitBuf.isValid == 0) {
      
    CreateBitBuffer(&h_ps_e->psBitBuf,
                    (unsigned char*)ptr1,
                    255+15);

     /* h_ps_e->aaaICCDataBuffer[][] */
    for (i=0; i<h_ps_e->iidIccBins; i++) {
      h_ps_e->aaaICCDataBuffer[i][0] = -1.0f; 
    }
  }
  

  return 0;
}

void
DeletePsEnc(HANDLE_PS_ENC *h_ps_e)
{
  
  

  return;
}


/***************************************************************************/
/*!

  \brief  Extracts the parameters (IID, ICC) of the current frame.

****************************************************************************/
void
EncodePsFrame(FFTWFContext_t *fftctx,
              HANDLE_PS_ENC pms,
              float **iBufferLeft,
              float **rBufferLeft,
              float **iBufferRight,
              float **rBufferRight)
{
  int     env;
  int     i;
  int     bin;
  int     subband, maxSubband;
  int     startSample, stopSample;
  float **hybrLeftImag;
  float **hybrLeftReal;
  float **hybrRightImag;
  float **hybrRightReal;

  

   
  HybridAnalysis ( fftctx,
                   (const float**) rBufferLeft,
                   (const float**) iBufferLeft,
                   pms->mHybridRealLeft,
                   pms->mHybridImagLeft,
                   pms->hHybridLeft);

   
  HybridAnalysis ( fftctx,
                   (const float**) rBufferRight,
                   (const float**) iBufferRight,
                   pms->mHybridRealRight,
                   pms->mHybridImagRight,
                   pms->hHybridRight);

   /*  pms->iidIccBins
                   pms->bHiFreqResIidIcc
                   pms->aaaIIDDataBuffer[][]
                   pms->aaaICCDataBuffer[][]
                */

  
  for ( bin = 0; bin < pms->iidIccBins; bin++ ) {
    
    pms->aaaIIDDataBuffer[bin][1] = pms->aaaIIDDataBuffer[bin][0];
    pms->aaaICCDataBuffer[bin][1] = pms->aaaICCDataBuffer[bin][0];
  }

  
  for ( env = 0; env < 2; env++ ) {

     
    hybrLeftReal  = pms->mHybridRealLeft;
    hybrLeftImag  = pms->mHybridImagLeft;
    hybrRightReal = pms->mHybridRealRight;
    hybrRightImag = pms->mHybridImagRight;

    
    if ( env == 0  ) {

      
      startSample   = 0;

      
      stopSample    = NO_SUBSAMPLES/2;
    }
    else {
      
      startSample   = NO_SUBSAMPLES/2;

      
      stopSample    = NO_SUBSAMPLES;
    }

     /* hiResBandBorders[]
                    pms->powerLeft[]
                    pms->powerRight[]
                    pms->powerCorrReal[]
                    pms->powerCorrImag[]
                    pms->histQmfLeftReal;
                    pms->histQmfLeftImag;
                    pms->histQmfRightReal;
                    pms->histQmfRightImag;
                 */
    
    for ( bin = 0; bin < NO_BINS; bin++ ) {

       
      if ( bin >= SUBQMF_BINS_ENERGY ) {

        
        if ( env == 0 ) {

          
          hybrLeftReal  = pms->histQmfLeftReal;
          hybrLeftImag  = pms->histQmfLeftImag;
          hybrRightReal = pms->histQmfRightReal;
          hybrRightImag = pms->histQmfRightImag;
        }
        else {

          
          hybrLeftReal  = rBufferLeft -HYBRID_FILTER_DELAY;
          hybrLeftImag  = iBufferLeft -HYBRID_FILTER_DELAY;
          hybrRightReal = rBufferRight-HYBRID_FILTER_DELAY;
          hybrRightImag = iBufferRight-HYBRID_FILTER_DELAY;
        }
      }

        
      maxSubband = ( bin < SUBQMF_BINS_ENERGY )? hiResBandBorders[bin] + 1: hiResBandBorders[bin + 1];

      
      for ( i = startSample; i < stopSample; i++ ) {

          
        if ( bin >= SUBQMF_BINS_ENERGY && i == HYBRID_FILTER_DELAY ) {

          
          hybrLeftReal  = rBufferLeft -HYBRID_FILTER_DELAY;
          hybrLeftImag  = iBufferLeft -HYBRID_FILTER_DELAY;
          hybrRightReal = rBufferRight-HYBRID_FILTER_DELAY;
          hybrRightImag = iBufferRight-HYBRID_FILTER_DELAY;
        }

         /* hybrLeftReal[][]
                        hybrLeftImag[][]
                        hybrRightReal[][]
                        hybrRightImag[][]
                     */
        
        for ( subband = hiResBandBorders[bin]; subband < maxSubband; subband++ ) {

           
          pms->powerLeft    [bin] += hybrLeftReal [i][subband] * hybrLeftReal [i][subband] +
                                     hybrLeftImag [i][subband] * hybrLeftImag [i][subband];

           
          pms->powerRight   [bin] += hybrRightReal[i][subband] * hybrRightReal[i][subband] +
                                     hybrRightImag[i][subband] * hybrRightImag[i][subband];

           
          pms->powerCorrReal[bin] += hybrLeftReal [i][subband] * hybrRightReal[i][subband] +
                                     hybrLeftImag [i][subband] * hybrRightImag[i][subband];

             
          pms->powerCorrImag[bin] += hybrLeftImag [i][subband] * hybrRightReal[i][subband] -
                                     hybrLeftReal [i][subband] * hybrRightImag[i][subband];
        }
      }   /* for (i) */

      
      if ( env == 0 ) {

         
        pms->powerLeft    [bin] += FLOAT_EPSILON;
        pms->powerRight   [bin] += FLOAT_EPSILON;
        pms->powerCorrReal[bin] += FLOAT_EPSILON;
        pms->powerCorrImag[bin] += FLOAT_EPSILON;
      }
    } /* bins loop */

    
    if (env == 0) {
      float  tempLeft;
      float  tempRight;
      float  tempCorrR;
      float  tempCorrI;

       /* pms->powerLeft[]
                      pms->powerRight[]
                      pms->powerCorrReal[]
                      pms->powerCorrImag[]
                      pms->aaaICCDataBuffer[][][]
                   */
      
      for ( bin = 0; bin < pms->iidIccBins; bin++ ) {

        
        if ( pms->bHiFreqResIidIcc ) {
          
          tempLeft  = pms->powerLeft    [bin];
          tempRight = pms->powerRight   [bin];
          tempCorrR = pms->powerCorrReal[bin];
          tempCorrI = pms->powerCorrImag[bin];
        }
        else {
           
          tempLeft  = pms->powerLeft    [2 * bin] + pms->powerLeft    [2 * bin + 1];
          tempRight = pms->powerRight   [2 * bin] + pms->powerRight   [2 * bin + 1];
          tempCorrR = pms->powerCorrReal[2 * bin] + pms->powerCorrReal[2 * bin + 1];
          tempCorrI = pms->powerCorrImag[2 * bin] + pms->powerCorrImag[2 * bin + 1];
        }

         
        if (bin > NO_IPD_BINS) {
           
          tempCorrR = tempCorrR * tempCorrR + tempCorrI * tempCorrI;

             
          pms->aaaICCDataBuffer[bin][0] = ( float )sqrt(tempCorrR / ( tempLeft * tempRight ));
        }
        else
        {
             
          pms->aaaICCDataBuffer[bin][0] = tempCorrR / ( float )sqrt( tempLeft * tempRight );
        }

         
        if ( pms->aaaICCDataBuffer[bin][0] > 1.0f ) {

          
          pms->aaaICCDataBuffer[bin][0] = 0;
        }
        else
        {
             
          pms->aaaICCDataBuffer[bin][0] = ( float )sqrt( ( 1.0f - pms->aaaICCDataBuffer[bin][0] ) * 0.5f );
        }

           
        pms->aaaIIDDataBuffer[bin][0] = INV_LOG2 * ( float )log( sqrt( tempLeft / tempRight ) );
      }

          
      memset( pms->powerLeft    , 0, sizeof( pms->powerLeft     ) );

          
      memset( pms->powerRight   , 0, sizeof( pms->powerRight    ) );

          
      memset( pms->powerCorrReal, 0, sizeof( pms->powerCorrReal ) );

          
      memset( pms->powerCorrImag, 0, sizeof( pms->powerCorrImag ) );
    } /* if (env==0) */
  } /* envelopes loop */

  
  downmixToMono( pms, rBufferLeft, iBufferLeft, rBufferRight, iBufferRight );

  
}

/***************************************************************************/
/*!
  \brief  generates mono downmix

****************************************************************************/
static void
downmixToMono( HANDLE_PS_ENC pms,
               float **qmfLeftReal,
               float **qmfLeftImag,
               float **qmfRightReal,
               float **qmfRightImag )
{
  int     i;
  int     group;
  int     subband;
  int     maxSubband;
  int     bin;
  float   temp;
  float   temp2;
  float   tempLeftReal;
  float   tempLeftImag;
  float   tempRightReal;
  float   tempRightImag;
  float **hybrLeftReal;
  float **hybrLeftImag;
  float **hybrRightReal;
  float **hybrRightImag;

  

  
  for ( i = 0; i < HYBRID_FILTER_DELAY; i++ ) {
        
    memcpy( pms->tempQmfLeftReal [i], qmfLeftReal [NO_SUBSAMPLES-HYBRID_FILTER_DELAY+i], NO_QMF_BANDS * sizeof( **qmfLeftReal  ) );

        
    memcpy( pms->tempQmfLeftImag [i], qmfLeftImag [NO_SUBSAMPLES-HYBRID_FILTER_DELAY+i], NO_QMF_BANDS * sizeof( **qmfLeftImag  ) );
  }

   
  hybrLeftReal  = pms->mHybridRealLeft;
  hybrLeftImag  = pms->mHybridImagLeft;
  hybrRightReal = pms->mHybridRealRight;
  hybrRightImag = pms->mHybridImagRight;

   /* bins2groupMap[group]
                  groupBordersMix[group]
               */
  
  for ( group = 0; group < NO_IPD_GROUPS; group++ ) {

    
    bin = ( ~NEGATE_IPD_MASK ) & bins2groupMap[group];

     
    if ( group >= SUBQMF_GROUPS_MIX ) {

      
      hybrLeftReal  = qmfLeftReal -HYBRID_FILTER_DELAY;
      hybrLeftImag  = qmfLeftImag -HYBRID_FILTER_DELAY;
      hybrRightReal = qmfRightReal-HYBRID_FILTER_DELAY;
      hybrRightImag = qmfRightImag-HYBRID_FILTER_DELAY;
    }

      /* ... ? */ 
    maxSubband = ( group < SUBQMF_GROUPS_MIX )? groupBordersMix[group] + 1: groupBordersMix[group + 1];

    
    for ( i = NO_SUBSAMPLES-1; i >= 0; i-- ) {

        
      if ( group >= SUBQMF_GROUPS_MIX && i == HYBRID_FILTER_DELAY-1 ) {

        
        hybrLeftReal  = pms->histQmfLeftReal;
        hybrLeftImag  = pms->histQmfLeftImag;
        hybrRightReal = pms->histQmfRightReal;
        hybrRightImag = pms->histQmfRightImag;
      }

       /* hybrLeftReal [i][subband]
                      hybrLeftImag [i][subband]
                      hybrRightReal[i][subband]
                      hybrRightImag[i][subband]
                   */
      
      for ( subband = groupBordersMix[group]; subband < maxSubband; subband++ ) {

        
        tempLeftReal  = hybrLeftReal [i][subband];
        tempLeftImag  = hybrLeftImag [i][subband];
        tempRightReal = hybrRightReal[i][subband];
        tempRightImag = hybrRightImag[i][subband];

          
        temp = ( tempLeftReal  * tempLeftReal  + tempLeftImag  * tempLeftImag  +
                 tempRightReal * tempRightReal + tempRightImag * tempRightImag ) * 0.5f + FLOAT_EPSILON;

          
        temp2 = temp + ( tempLeftReal * tempRightReal +
                         tempLeftImag * tempRightImag );

          
        if (temp > temp2 * 2.0f * MAX_MAGNITUDE_CORR * MAX_MAGNITUDE_CORR) {
          
          temp = MAX_MAGNITUDE_CORR;
        }
        else {
            
          temp = ( float )sqrt( 0.5f * temp / temp2 );
        }

         
        tempLeftReal = ( tempLeftReal + tempRightReal ) * temp;
        tempLeftImag = ( tempLeftImag + tempRightImag ) * temp;

         
        if ( group < SUBQMF_GROUPS_MIX ) {
          
          hybrLeftReal[i][subband] = tempLeftReal;
          hybrLeftImag[i][subband] = tempLeftImag;
        }
        else {

          
          qmfLeftReal [i][subband] = tempLeftReal;
          qmfLeftImag [i][subband] = tempLeftImag;
        }
      } /* subbands loop */
    } /* loop (i) */
  } /* groups loop */

  
  for ( i = 0; i < HYBRID_FILTER_DELAY; i++ ) {
        
    memcpy( pms->histQmfLeftReal [i], pms->tempQmfLeftReal [i], NO_QMF_BANDS * sizeof( **qmfLeftReal  ) );

        
    memcpy( pms->histQmfLeftImag [i], pms->tempQmfLeftImag [i], NO_QMF_BANDS * sizeof( **qmfLeftImag  ) );

        
    memcpy( pms->histQmfRightReal[i], qmfRightReal[NO_SUBSAMPLES-HYBRID_FILTER_DELAY+i], NO_QMF_BANDS * sizeof( **qmfRightReal ) );

        
    memcpy( pms->histQmfRightImag[i], qmfRightImag[NO_SUBSAMPLES-HYBRID_FILTER_DELAY+i], NO_QMF_BANDS * sizeof( **qmfRightImag ) );
  }

   
  HybridSynthesis( ( const float** )pms->mHybridRealLeft,
                   ( const float** )pms->mHybridImagLeft,
                   qmfLeftReal,
                   qmfLeftImag,
                   pms->hHybridLeft );

  
}

