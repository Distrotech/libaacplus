/*
  Envelope estimation
*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif



#define QUANT_ERROR_THRES 200


/***************************************************************************/
/*!

  \brief  Quantisation of the panorama value

  \return the quantized pan value

****************************************************************************/
static int
mapPanorama (int nrgVal,
             int ampRes,
             int *quantError
             )
{
  int i;
  int min_val, val;
  int panTable[2][10] = {{0,2,4,6,8,12,16,20,24},
                         { 0, 2, 4, 8, 12 }};
  int maxIndex[2] = {9,5};

  int panIndex;
  int sign;

  

   
  sign = nrgVal > 0 ? 1 : -1;

  
  nrgVal = sign * nrgVal;

  
  min_val = INT_MAX;
  panIndex = 0;

   /* pointer for maxIndex[ampRes],
                              panTable[ampRes][i]
               */
  
  for (i = 0; i < maxIndex[ampRes]; i++) {

     
    val = abs (nrgVal - panTable[ampRes][i]);

     
    if (val < min_val) {

      
      min_val = val;
      panIndex = i;
    }
  }

  
  *quantError=min_val;

     /* counting post-operations */

  

  return panTable[ampRes][maxIndex[ampRes]-1] + sign * panTable[ampRes][panIndex];
}


/***************************************************************************/
/*!

  \brief  Quantisation of the noise floor levels

  \return void

****************************************************************************/
static void
sbrNoiseFloorLevelsQuantisation (int   *iNoiseLevels,
                                 float *NoiseLevels,
                                 int   coupling
                                 )
{
  int i;
  int dummy;

  

   /* pointer for NoiseLevels[i],
                              iNoiseLevels[i]
               */
  
  for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {
    int tmp;

     
    tmp = NoiseLevels[i] > (float)30.0f ? 30: (int) (NoiseLevels[i] + (float)0.5);

    
    if (coupling) {

        
      tmp = tmp < -30 ? -30  : tmp;

      
      tmp = mapPanorama (tmp,1,&dummy);
    }

    
    iNoiseLevels[i] = tmp;

  }

  
}

/***************************************************************************/
/*!

  \brief  Calculation of noise floor for coupling

  \return void

****************************************************************************/
static void
coupleNoiseFloor(float *noise_level_left,
                 float *noise_level_right
                 )
{
  int i;

  

   /* pointer for noise_level_left[i]
                              noise_level_right[i]
               */
  
  for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {
    float temp1, temp2;

      
    temp1 = (float) pow (2.0f, (- noise_level_right[i] + NOISE_FLOOR_OFFSET));
    temp2 = (float) pow (2.0f, (- noise_level_left[i]  + NOISE_FLOOR_OFFSET));

        
    noise_level_left[i] =  (float)(NOISE_FLOOR_OFFSET - log((temp1 + temp2) / 2) * ILOG2);

       
    noise_level_right[i] = (float) (log(temp2 / temp1) * ILOG2);
  }

  
}


/***************************************************************************/
/*!

  \brief  calculates the SBR envelope values

  \return void

****************************************************************************/
static void
calculateSbrEnvelope (float **YBufferLeft,
                      float **YBufferRight,
                      const SBR_FRAME_INFO *frame_info,
                      int *sfb_nrgLeft,
                      int *sfb_nrgRight,
                      HANDLE_SBR_CONFIG_DATA h_con,
                      HANDLE_ENV_CHANNEL h_sbr,
                      SBR_STEREO_MODE stereoMode,
                      int* maxQuantError)

{
  int i, j, k, l, count, m = 0;
  int no_of_bands, start_pos, stop_pos, li, ui;
  FREQ_RES freq_res;

  int ca = 2 - h_sbr->encEnvData.init_sbr_amp_res;
  int quantError;
  int nEnvelopes = frame_info->nEnvelopes;
  int short_env = frame_info->shortEnv - 1;
  int timeStep = h_sbr->sbrExtractEnvelope.time_step;
  int missingHarmonic = 0;

  

     /* counting previous operations */

   
  if (stereoMode == SBR_COUPLING) {

    
    *maxQuantError = 0;
  }

   /* pointers for frame_info->borders[i]
                               frame_info->freqRes[i]
                               h_sbr->encEnvData.addHarmonic[i]
                               sfb_nrgLeft[m]
                               sfb_nrgRight[m]
               */
  
  for (i = 0; i < nEnvelopes; i++) {

    
    start_pos = timeStep * frame_info->borders[i];
    stop_pos = timeStep * frame_info->borders[i + 1];

    
    freq_res = frame_info->freqRes[i];

     
    no_of_bands = h_con->nSfb[freq_res];

     
    if (i == short_env) {
      
      stop_pos = stop_pos - timeStep;
    }


     /* pointers for h_con->freqBandTable[freq_res][j],
                                 h_sbr->encEnvData.addHarmonic[j]
                 */
    
    for (j = 0; j < no_of_bands; j++) {

      float nrgRight = 0;
      float nrgLeft = 0;
      float temp;

       /* counting previous operations */

      
      li = h_con->freqBandTable[freq_res][j];
      ui = h_con->freqBandTable[freq_res][j + 1];

       
      if(freq_res == FREQ_RES_HIGH){

          
        if(j == 0 && ui-li > 1){
          
          li++;
        }
      }
      else{
          
        if(j == 0 && ui-li > 2){
          
          li++;
        }
      }


      
      missingHarmonic = 0;

       
      if(h_sbr->encEnvData.addHarmonicFlag){

         
        if(freq_res == FREQ_RES_HIGH){

          
          if(h_sbr->encEnvData.addHarmonic[j]){

            
            missingHarmonic = 1;
          }
        }
        else{
          int i;
          int startBandHigh = 0;
          int stopBandHigh = 0;

           /* counting previous operations */

           /* pointer for h_con->freqBandTable[FREQ_RES_HIGH][startBandHigh] */
          
          while(h_con->freqBandTable[FREQ_RES_HIGH][startBandHigh] < h_con->freqBandTable[FREQ_RES_LOW][j]) {
             /* while() condition */

            
            startBandHigh++;
          }

           /* pointer for h_con->freqBandTable[FREQ_RES_HIGH][stopBandHigh] */
          
          while(h_con->freqBandTable[FREQ_RES_HIGH][stopBandHigh] < h_con->freqBandTable[FREQ_RES_LOW][j + 1]) {
             /* while() condition */

            
            stopBandHigh++;
          }

           /* pointer for h_sbr->encEnvData.addHarmonic[] */
          
          for(i = startBandHigh; i<stopBandHigh; i++){

            
            if(h_sbr->encEnvData.addHarmonic[i]){

              
              missingHarmonic = 1;
            }
          }
        }

      }


      
      if(missingHarmonic){

        float tmpNrg = 0;

         /* counting previous operation */

        
        count = (stop_pos - start_pos);

         /* pointer for YBufferLeft[][] */
        
        for (l = start_pos; l < stop_pos; l++) {

          
          nrgLeft += YBufferLeft[l/2][li];
        }

        
        for (k = li+1; k < ui; k++){

          
          tmpNrg = 0;

           /* pointer for YBufferLeft[][] */
          
          for (l = start_pos; l < stop_pos; l++) {

            
            tmpNrg += YBufferLeft[l/2][k];
          }

           
          if(tmpNrg > nrgLeft){

            
            nrgLeft = tmpNrg;
          }
        }

         
        if(ui-li > 2){

          
          nrgLeft = nrgLeft*0.398107267f;

        }
        else{

           
          if(ui-li > 1){

            
            nrgLeft = nrgLeft*0.5f;
          }
        }

         
        if (stereoMode == SBR_COUPLING) {

           /* pointer for YBufferRight[][] */
          
          for (l = start_pos; l < stop_pos; l++) {

            
            nrgRight += YBufferRight[l/2][li];
          }


          
          for (k = li+1; k < ui; k++){

            
            tmpNrg = 0;

             /* pointer for YBufferRight[][] */
            
            for (l = start_pos; l < stop_pos; l++) {

              
              tmpNrg += YBufferRight[l/2][k];
            }

             
            if(tmpNrg > nrgRight){

              
              nrgRight = tmpNrg;
            }
          }


           
          if(ui-li > 2){

            
            nrgRight = nrgRight*0.398107267f;

          }
          else{

             
            if(ui-li > 1){

              
              nrgRight = nrgRight*0.5f;
            }
          }

        }


         
        if (stereoMode == SBR_COUPLING) {

          
          temp = nrgLeft;

           
          nrgLeft = (nrgRight + nrgLeft) * (float)0.5f;

           
          nrgRight = (temp + 1) / (nrgRight + 1);
        }


    } /* end missingHarmonic */
    else{

         
        count = (stop_pos - start_pos) * (ui - li);


        
        for (k = li; k < ui; k++) {

           /* pointer for YBufferLeft[][] */
          
          for (l = start_pos; l < stop_pos; l++) {

            
            nrgLeft += YBufferLeft[l/2][k];
          }
        }

         
        if (stereoMode == SBR_COUPLING) {

          
          for (k = li; k < ui; k++) {

             /* pointer for YBufferRight[][] */
            
            for (l = start_pos; l < stop_pos; l++) {

              
              nrgRight += YBufferRight[l/2][k];
            }
          }
        }

        
       if (stereoMode == SBR_COUPLING) {

        
        temp = nrgLeft;

         
        nrgLeft = (nrgRight + nrgLeft) * 0.5f;

         
        nrgRight = (temp + 1) / (nrgRight + 1);
      }

    }


          
      nrgLeft = (float) (log (nrgLeft / (count * 64) + EPS) * ILOG2);

       
      nrgLeft = max(nrgLeft,0);

        
      sfb_nrgLeft[m] = (int) (ca * nrgLeft + 0.5);

       
      if (stereoMode == SBR_COUPLING) {

         
        nrgRight = (float) (log (nrgRight) * ILOG2);

            
        sfb_nrgRight[m] = mapPanorama ((int)((float)ca * nrgRight),h_sbr->encEnvData.init_sbr_amp_res,&quantError);

         
        if(quantError > *maxQuantError) {
          
          *maxQuantError = quantError;
        }
      }
      m++;

    } /* j */


     
    if (h_con->useParametricCoding) {

      
      m-=no_of_bands;

       /* pointers for h_sbr->sbrExtractEnvelope.envelopeCompensation[j]
                                   sfb_nrgLeft[m]
                                   sfb_nrgRight[m]
                   */
      
      for (j = 0; j < no_of_bands; j++) {

          
        if (freq_res==FREQ_RES_HIGH && h_sbr->sbrExtractEnvelope.envelopeCompensation[j]){

             
          sfb_nrgLeft[m] -= (int)(ca * abs(h_sbr->sbrExtractEnvelope.envelopeCompensation[j]));
        }

        
        if(sfb_nrgLeft[m] < 0) {

          
          sfb_nrgLeft[m] = 0;
        }
        m++;
      }
    }
  } /* i*/

  
}



/***************************************************************************/
/*!

  \brief  calculates the noise floor and the envelope values.

****************************************************************************/
void
extractSbrEnvelope (FFTWFContext_t *fftctx,
                    float *timeInPtr,
                    float *pCoreBuffer,
                    unsigned int timeInStride,
                    HANDLE_SBR_CONFIG_DATA h_con,
                    HANDLE_SBR_HEADER_DATA sbrHeaderData,
                    HANDLE_SBR_BITSTREAM_DATA sbrBitstreamData,
                    HANDLE_ENV_CHANNEL h_envChan[MAX_CHANNELS],
                    struct PS_ENC *hPsEnc,
                    HANDLE_SBR_QMF_FILTER_BANK hSynthesisQmfBank,
                    HANDLE_COMMON_DATA hCmonData)
{
  int ch, i, j, c;
  int nEnvelopes[MAX_CHANNELS];
  int transient_info[MAX_CHANNELS][2];

  const SBR_FRAME_INFO *frame_info[MAX_CHANNELS];

  int nChannels = h_con->nChannels;
  int nInChannels = (hPsEnc) ? 2:nChannels;

  SBR_STEREO_MODE stereoMode = h_con->stereoMode;

  FREQ_RES res[MAX_NUM_NOISE_VALUES];
  int v_tuning[6] = { 0, 2, 4, 0, 0, 0 };

  int   sfb_nrg    [MAX_CHANNELS][MAX_NUM_ENVELOPE_VALUES];

  float noiseFloor[MAX_CHANNELS][MAX_NUM_NOISE_VALUES];
  int   noise_level[MAX_CHANNELS][MAX_NUM_NOISE_VALUES];

#if (MAX_CHANNELS>1)
  int   sfb_nrg_coupling[MAX_CHANNELS][MAX_NUM_ENVELOPE_VALUES];
  int   noise_level_coupling[MAX_CHANNELS][MAX_NUM_NOISE_VALUES];
  int   maxQuantError;
#endif

  

     

   /* pointer for res[i] */
  
  for(i=0; i<MAX_NUM_NOISE_VALUES; i++) {

    
    res[i] = FREQ_RES_HIGH;
  }

      
  memset(noiseFloor,0,sizeof(noiseFloor));

      
  memset(transient_info,0,sizeof(transient_info));


   /* pointers for h_envChan[ch]->sbrExtractEnvelope,
                               h_envChan[ch]->sbrExtractEnvelope,
                               h_envChan[ch]->sbrQmf
               */
  
  for(ch = 0; ch < nInChannels;ch++){

    
    if (hPsEnc) {
        
      sbrAnalysisFiltering  (timeInPtr ? timeInPtr+ch:NULL,
                             MAX_CHANNELS,
                             h_envChan[ch]->sbrExtractEnvelope.rBuffer,
                             h_envChan[ch]->sbrExtractEnvelope.iBuffer,
                             &h_envChan[ch]->sbrQmf);
    }
    else
    {
        
      sbrAnalysisFiltering  (timeInPtr ? timeInPtr+ch:NULL,
                             timeInStride,
                             h_envChan[ch]->sbrExtractEnvelope.rBuffer,
                             h_envChan[ch]->sbrExtractEnvelope.iBuffer,
                             &h_envChan[ch]->sbrQmf);
    }
  } /* ch */

#ifndef MONO_ONLY

   
  if (hPsEnc && hSynthesisQmfBank) {

     
    EncodePsFrame(fftctx, hPsEnc,
                  h_envChan[0]->sbrExtractEnvelope.iBuffer,
                  h_envChan[0]->sbrExtractEnvelope.rBuffer,
                  h_envChan[1]->sbrExtractEnvelope.iBuffer,
                  h_envChan[1]->sbrExtractEnvelope.rBuffer);

    
    SynthesisQmfFiltering (h_envChan[0]->sbrExtractEnvelope.rBuffer,
                           h_envChan[0]->sbrExtractEnvelope.iBuffer,
                           pCoreBuffer,
                           (HANDLE_SBR_QMF_FILTER_BANK)hSynthesisQmfBank);

     /* pointer for pCoreBuffer[] */
      
    for (i = 0; i < 1024; i++) {

       
      pCoreBuffer[i] = -pCoreBuffer[i];
    }
  }

#endif /* #ifndef MONO_ONLY */

   /* pointers for h_envChan[ch]->sbrExtractEnvelope,
                               h_envChan[ch]->sbrQmf,
                               h_envChan[ch]->TonCorr,
                               h_envChan[ch]->SbrEnvFrame,
                               h_envChan[ch]->sbrTransientDetector,
                               transient_info[ch],
                               h_con->freqBandTable[][],
                               h_con->nSfb[]
               */
  
  for(ch = 0; ch < nChannels;ch++) {

     
    getEnergyFromCplxQmfData (h_envChan[ch]->sbrExtractEnvelope.YBuffer + h_envChan[ch]->sbrExtractEnvelope.YBufferWriteOffset,
                              h_envChan[ch]->sbrExtractEnvelope.rBuffer,
                              h_envChan[ch]->sbrExtractEnvelope.iBuffer);




     
    CalculateTonalityQuotas(&h_envChan[ch]->TonCorr,
                            h_envChan[ch]->sbrExtractEnvelope.rBuffer,
                            h_envChan[ch]->sbrExtractEnvelope.iBuffer,
                            h_con->freqBandTable[HI][h_con->nSfb[HI]]);

    
    transientDetect (h_envChan[ch]->sbrExtractEnvelope.YBuffer,
                     &h_envChan[ch]->sbrTransientDetector,
                     transient_info[ch],
                     h_envChan[ch]->sbrExtractEnvelope.time_step
                     );


     
    frameSplitter(h_envChan[ch]->sbrExtractEnvelope.YBuffer,
                  &h_envChan[ch]->sbrTransientDetector,
                  h_con->freqBandTable[1],
                  h_con->nSfb[1],
                  h_envChan[ch]->sbrExtractEnvelope.time_step,
                  h_envChan[ch]->sbrExtractEnvelope.no_cols,
                  transient_info[ch]);
  } /* ch */


#if (MAX_CHANNELS>1)
   
  if (stereoMode == SBR_COUPLING) {

     
    if (transient_info[0][1] && transient_info[1][1]) {

        
      transient_info[0][0] =
        min (transient_info[1][0], transient_info[0][0]);

      
      transient_info[1][0] = transient_info[0][0];
    }
    else {

       
      if (transient_info[0][1] && !transient_info[1][1]) {

        
        transient_info[1][0] = transient_info[0][0];
      }
      else {

         
        if (!transient_info[0][1] && transient_info[1][1]) {

          
          transient_info[0][0] = transient_info[1][0];
        }
        else {

            
          transient_info[0][0] =
            max (transient_info[1][0], transient_info[0][0]);

          
          transient_info[1][0] = transient_info[0][0];
        }
      }
    }
  }
#endif


    
  frame_info[0] = frameInfoGenerator (&h_envChan[0]->SbrEnvFrame,
                                      h_envChan[0]->sbrExtractEnvelope.pre_transient_info,
                                      transient_info[0],
                                      v_tuning);

   
  h_envChan[0]->encEnvData.hSbrBSGrid = &h_envChan[0]->SbrEnvFrame.SbrGrid;


#if (MAX_CHANNELS>1)

  
  switch (stereoMode) {
  case SBR_LEFT_RIGHT:
  case SBR_SWITCH_LRC:

      
    frame_info[1] = frameInfoGenerator (&h_envChan[1]->SbrEnvFrame,
                                        h_envChan[1]->sbrExtractEnvelope.pre_transient_info,
                                        transient_info[1],
                                        v_tuning);

     
    h_envChan[1]->encEnvData.hSbrBSGrid = &h_envChan[1]->SbrEnvFrame.SbrGrid;

      
    if (frame_info[0]->nEnvelopes != frame_info[1]->nEnvelopes) {

      
      stereoMode = SBR_LEFT_RIGHT;

    } else {

       /* pointer for frame_info[0]->borders[i],
                                  frame_info[1]->borders[i]
                   */
       
      for (i = 0; i < frame_info[0]->nEnvelopes + 1; i++) {

         
        if (frame_info[0]->borders[i] != frame_info[1]->borders[i]) {

          
          stereoMode = SBR_LEFT_RIGHT;
          break;
        }
      }

       /* pointer for frame_info[0]->freqRes[i],
                                  frame_info[1]->freqRes[i]
                   */
       
      for (i = 0; i < frame_info[0]->nEnvelopes; i++) {

         
        if (frame_info[0]->freqRes[i] != frame_info[1]->freqRes[i]) {

          
          stereoMode = SBR_LEFT_RIGHT;
          break;
        }
      }

        
      if (frame_info[0]->shortEnv != frame_info[1]->shortEnv) {

        
        stereoMode = SBR_LEFT_RIGHT;
      }
    }
    break;
  case SBR_COUPLING:

    
    frame_info[1] = frame_info[0];

     
    h_envChan[1]->encEnvData.hSbrBSGrid = &h_envChan[0]->SbrEnvFrame.SbrGrid;
    break;
  case SBR_MONO:
    /* nothing to do */
    break;
  default:
    assert (0);
  }

#endif /* (MAX_CHANNELS>1) */


   /* pointers for h_envChan[ch]->sbrExtractEnvelope,
                                h_envChan[ch]->encEnvData,
                                h_envChan[ch]->encEnvData.hSbrBSGrid->frameClass,
                                h_envChan[ch]->sbrQmf,
                                h_envChan[ch]->sbrCodeEnvelope,
                                h_envChan[ch]->sbrCodeNoiseFloor,
                                h_envChan[ch]->sbrExtractEnvelope,
                                h_envChan[ch]->TonCorr,
                                h_con->freqBandTable[][],
                                h_con->nSfb[],
                                frame_info[ch]->nEnvelopes,
                                frame_info[ch]->freqRes[i],
                                nEnvelopes[ch],
                                noiseFloor[ch],
                                transient_info[ch]
                */
  
  for(ch = 0; ch < nChannels;ch++){

    
    h_envChan[ch]->sbrExtractEnvelope.pre_transient_info[0] = transient_info[ch][0];
    h_envChan[ch]->sbrExtractEnvelope.pre_transient_info[1] = transient_info[ch][1];
    h_envChan[ch]->encEnvData.noOfEnvelopes = nEnvelopes[ch] = frame_info[ch]->nEnvelopes;


    
    for (i = 0; i < nEnvelopes[ch]; i++){

         
      h_envChan[ch]->encEnvData.noScfBands[i] =
      (frame_info[ch]->freqRes[i] == FREQ_RES_HIGH ? h_con->nSfb[FREQ_RES_HIGH] : h_con->nSfb[FREQ_RES_LOW]);
    }


     
   if( ( h_envChan[ch]->encEnvData.hSbrBSGrid->frameClass == FIXFIX ) &&
        ( nEnvelopes[ch] == 1 ) ) {

       
      if( h_envChan[ch]->encEnvData.init_sbr_amp_res != SBR_AMP_RES_1_5 ) {

        
        InitSbrHuffmanTables (&h_envChan[ch]->encEnvData,
                              &h_envChan[ch]->sbrCodeEnvelope,
                              &h_envChan[ch]->sbrCodeNoiseFloor,
                              SBR_AMP_RES_1_5);
      }
    }
    else {

        
      if(sbrHeaderData->sbr_amp_res != h_envChan[ch]->encEnvData.init_sbr_amp_res ) {

        
        InitSbrHuffmanTables (&h_envChan[ch]->encEnvData,
                              &h_envChan[ch]->sbrCodeEnvelope,
                              &h_envChan[ch]->sbrCodeNoiseFloor,
                              sbrHeaderData->sbr_amp_res);
      }
    }


     
    TonCorrParamExtr(&h_envChan[ch]->TonCorr,
                     h_envChan[ch]->encEnvData.sbr_invf_mode_vec,
                     noiseFloor[ch],
                     &h_envChan[ch]->encEnvData.addHarmonicFlag,
                     h_envChan[ch]->encEnvData.addHarmonic,
                     h_envChan[ch]->sbrExtractEnvelope.envelopeCompensation,
                     frame_info[ch],
                     transient_info[ch],
                     h_con->freqBandTable[HI],
                     h_con->nSfb[HI],
                     h_envChan[ch]->encEnvData.sbr_xpos_mode);

    
    h_envChan[ch]->encEnvData.sbr_invf_mode = h_envChan[ch]->encEnvData.sbr_invf_mode_vec[0];
    h_envChan[ch]->encEnvData.noOfnoisebands = h_envChan[ch]->TonCorr.sbrNoiseFloorEstimate.noNoiseBands;




  } /* ch */


  
  switch (stereoMode) {

  case SBR_MONO:
     
    calculateSbrEnvelope (h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL, frame_info[0], sfb_nrg[0],
                          NULL,h_con,h_envChan[0], SBR_MONO,NULL);
    break;

#if (MAX_CHANNELS>1)

  case SBR_LEFT_RIGHT:

     
    calculateSbrEnvelope (h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL, frame_info[0], sfb_nrg[0],
                          NULL, h_con, h_envChan[0], SBR_MONO,NULL);

     
    calculateSbrEnvelope (h_envChan[1]->sbrExtractEnvelope.YBuffer, NULL, frame_info[1], sfb_nrg[1],
                          NULL, h_con,h_envChan[1], SBR_MONO,NULL);
    break;

  case SBR_COUPLING:
     
    calculateSbrEnvelope (h_envChan[0]->sbrExtractEnvelope.YBuffer, h_envChan[1]->sbrExtractEnvelope.YBuffer, frame_info[0],
                          sfb_nrg[0], sfb_nrg[1], h_con,h_envChan[0], SBR_COUPLING,&maxQuantError);

    break;

  case SBR_SWITCH_LRC:

     
    calculateSbrEnvelope (h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL, frame_info[0], sfb_nrg[0],
                          NULL, h_con,h_envChan[0], SBR_MONO,NULL);

     
    calculateSbrEnvelope (h_envChan[1]->sbrExtractEnvelope.YBuffer, NULL, frame_info[1], sfb_nrg[1],
                          NULL, h_con,h_envChan[1], SBR_MONO,NULL);

     
    calculateSbrEnvelope (h_envChan[0]->sbrExtractEnvelope.YBuffer, h_envChan[1]->sbrExtractEnvelope.YBuffer, frame_info[0],
                          sfb_nrg_coupling[0], sfb_nrg_coupling[1], h_con,h_envChan[0], SBR_COUPLING,&maxQuantError);
    break;

#endif /* (MAX_CHANNELS>1) */

  default:
    assert(0);
  }



  
  switch (stereoMode) {

  case SBR_MONO:

    
    sbrNoiseFloorLevelsQuantisation (noise_level[0], noiseFloor[0], 0);

       
    codeEnvelope (noise_level[0], res,
                  &h_envChan[0]->sbrCodeNoiseFloor,
                  h_envChan[0]->encEnvData.domain_vec_noise, 0,
                  (frame_info[0]->nEnvelopes > 1 ? 2 : 1), 0,
                 sbrBitstreamData->HeaderActive);

    break;

#if (MAX_CHANNELS>1)

  case SBR_LEFT_RIGHT:

    
    sbrNoiseFloorLevelsQuantisation (noise_level[0],noiseFloor[0], 0);

       
    codeEnvelope (noise_level[0], res,
                  &h_envChan[0]->sbrCodeNoiseFloor,
                  h_envChan[0]->encEnvData.domain_vec_noise, 0,
                  (frame_info[0]->nEnvelopes > 1 ? 2 : 1), 0,
                  sbrBitstreamData->HeaderActive);


    
    sbrNoiseFloorLevelsQuantisation (noise_level[1],noiseFloor[1], 0);

       
    codeEnvelope (noise_level[1], res,
                  &h_envChan[1]->sbrCodeNoiseFloor,
                  h_envChan[1]->encEnvData.domain_vec_noise, 0,
                  (frame_info[1]->nEnvelopes > 1 ? 2 : 1), 0,
                  sbrBitstreamData->HeaderActive);

    break;

  case SBR_COUPLING:

    
    coupleNoiseFloor(noiseFloor[0],noiseFloor[1]);

    
    sbrNoiseFloorLevelsQuantisation (noise_level[0],noiseFloor[0], 0);

       
    codeEnvelope (noise_level[0], res,
                  &h_envChan[0]->sbrCodeNoiseFloor,
                  h_envChan[0]->encEnvData.domain_vec_noise, 1,
                  (frame_info[0]->nEnvelopes > 1 ? 2 : 1), 0,
                  sbrBitstreamData->HeaderActive);


    
    sbrNoiseFloorLevelsQuantisation (noise_level[1],noiseFloor[1], 1);

       
    codeEnvelope (noise_level[1], res,
                  &h_envChan[1]->sbrCodeNoiseFloor,
                  h_envChan[1]->encEnvData.domain_vec_noise, 1,
                  (frame_info[1]->nEnvelopes > 1 ? 2 : 1), 1,
                  sbrBitstreamData->HeaderActive);



    break;

  case SBR_SWITCH_LRC:

    
    sbrNoiseFloorLevelsQuantisation (noise_level[0],noiseFloor[0], 0);

    
    sbrNoiseFloorLevelsQuantisation (noise_level[1],noiseFloor[1], 0);

    
    coupleNoiseFloor(noiseFloor[0],noiseFloor[1]);

    
    sbrNoiseFloorLevelsQuantisation (noise_level_coupling[0],noiseFloor[0], 0);

    
    sbrNoiseFloorLevelsQuantisation (noise_level_coupling[1],noiseFloor[1], 1);
    break;


#endif /* (MAX_CHANNELS>1) */

  default:
    assert(0);

  }


  
  switch (stereoMode) {

  case SBR_MONO:

     
    sbrHeaderData->coupling = 0;
    h_envChan[0]->encEnvData.balance = 0;

     
    codeEnvelope (sfb_nrg[0], frame_info[0]->freqRes,
                  &h_envChan[0]->sbrCodeEnvelope,
                  h_envChan[0]->encEnvData.domain_vec,
                  sbrHeaderData->coupling,
                  frame_info[0]->nEnvelopes, 0,
                  sbrBitstreamData->HeaderActive);
    break;

#if (MAX_CHANNELS>1)

  case SBR_LEFT_RIGHT:

     
    sbrHeaderData->coupling = 0;

     
    h_envChan[0]->encEnvData.balance = 0;
    h_envChan[1]->encEnvData.balance = 0;


     
    codeEnvelope (sfb_nrg[0], frame_info[0]->freqRes,
                  &h_envChan[0]->sbrCodeEnvelope,
                  h_envChan[0]->encEnvData.domain_vec,
                  sbrHeaderData->coupling,
                  frame_info[0]->nEnvelopes, 0,
                  sbrBitstreamData->HeaderActive);

     
    codeEnvelope (sfb_nrg[1], frame_info[1]->freqRes,
                  &h_envChan[1]->sbrCodeEnvelope,
                  h_envChan[1]->encEnvData.domain_vec,
                  sbrHeaderData->coupling,
                  frame_info[1]->nEnvelopes, 0,
                  sbrBitstreamData->HeaderActive);
    break;

  case SBR_COUPLING:

     
    sbrHeaderData->coupling = 1;
    h_envChan[0]->encEnvData.balance = 0;
    h_envChan[1]->encEnvData.balance = 1;

     
    codeEnvelope (sfb_nrg[0], frame_info[0]->freqRes,
                  &h_envChan[0]->sbrCodeEnvelope,
                  h_envChan[0]->encEnvData.domain_vec,
                  sbrHeaderData->coupling,
                  frame_info[0]->nEnvelopes, 0,
                  sbrBitstreamData->HeaderActive);

     
    codeEnvelope (sfb_nrg[1], frame_info[1]->freqRes,
                  &h_envChan[1]->sbrCodeEnvelope,
                  h_envChan[1]->encEnvData.domain_vec,
                  sbrHeaderData->coupling,
                  frame_info[1]->nEnvelopes, 1,
                  sbrBitstreamData->HeaderActive);
    break;

  case SBR_SWITCH_LRC:
    {
      int payloadbitsLR;
      int payloadbitsCOUPLING;

      int sfbNrgPrevTemp[MAX_CHANNELS][MAX_FREQ_COEFFS];
      int noisePrevTemp[MAX_CHANNELS][MAX_NUM_NOISE_COEFFS];
      int upDateNrgTemp[MAX_CHANNELS];
      int upDateNoiseTemp[MAX_CHANNELS];
      int domainVecTemp[MAX_CHANNELS][MAX_ENVELOPES];
      int domainVecNoiseTemp[MAX_CHANNELS][MAX_ENVELOPES];

      int tempFlagRight = 0;
      int tempFlagLeft = 0;


       /* counting previous operations */


       /* pointers for h_envChan[ch]->sbrCodeEnvelope,
                                   h_envChan[ch]->sbrCodeNoiseFloor,
                                   sfbNrgPrevTemp[ch],
                                   noisePrevTemp[ch],
                                   upDateNrgTemp[ch],
                                   upDateNoiseTemp[ch]
                   */
      
      for(ch = 0; ch < nChannels;ch++){

            
        memcpy (sfbNrgPrevTemp[ch], h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev,
              MAX_FREQ_COEFFS * sizeof (int));

            
        memcpy (noisePrevTemp[ch], h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev,
              MAX_NUM_NOISE_COEFFS * sizeof (int));

        
        upDateNrgTemp[ch] = h_envChan[ch]->sbrCodeEnvelope.upDate;
        upDateNoiseTemp[ch] = h_envChan[ch]->sbrCodeNoiseFloor.upDate;

         
        if(sbrHeaderData->prev_coupling){

          
          h_envChan[ch]->sbrCodeEnvelope.upDate = 0;
          h_envChan[ch]->sbrCodeNoiseFloor.upDate = 0;
        }
      } /* ch */


       
      codeEnvelope (sfb_nrg[0], frame_info[0]->freqRes,
                    &h_envChan[0]->sbrCodeEnvelope,
                    h_envChan[0]->encEnvData.domain_vec, 0,
                    frame_info[0]->nEnvelopes, 0,
                    sbrBitstreamData->HeaderActive);

       
      codeEnvelope (sfb_nrg[1], frame_info[1]->freqRes,
                    &h_envChan[1]->sbrCodeEnvelope,
                    h_envChan[1]->encEnvData.domain_vec, 0,
                    frame_info[1]->nEnvelopes, 0,
                    sbrBitstreamData->HeaderActive);

      c = 0;

       /* pointers for h_envChan[0]->encEnvData.noScfBands[i]
                                   h_envChan[0]->encEnvData.ienvelope[i][j]
                                   h_envChan[1]->encEnvData.ienvelope[i][j]
                                   sfb_nrg[0][c]
                                   sfb_nrg[1][c]
                   */
      
      for (i = 0; i < nEnvelopes[0]; i++) {

        
        for (j = 0; j < h_envChan[0]->encEnvData.noScfBands[i]; j++)
          {
            
            h_envChan[0]->encEnvData.ienvelope[i][j] = sfb_nrg[0][c];
            h_envChan[1]->encEnvData.ienvelope[i][j] = sfb_nrg[1][c];

            c++;
          }
      }



         
      codeEnvelope (noise_level[0], res,
                    &h_envChan[0]->sbrCodeNoiseFloor,
                    h_envChan[0]->encEnvData.domain_vec_noise, 0,
                    (frame_info[0]->nEnvelopes > 1 ? 2 : 1), 0,
                    sbrBitstreamData->HeaderActive);


       /* pointers for h_envChan[0]->encEnvData.sbr_noise_levels[i],
                                   noise_level[0][i]
                   */
      
      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {

        
        h_envChan[0]->encEnvData.sbr_noise_levels[i] = noise_level[0][i];
      }


         
      codeEnvelope (noise_level[1], res,
                    &h_envChan[1]->sbrCodeNoiseFloor,
                    h_envChan[1]->encEnvData.domain_vec_noise, 0,
                    (frame_info[1]->nEnvelopes > 1 ? 2 : 1), 0,
                    sbrBitstreamData->HeaderActive);

       /* pointers for h_envChan[1]->encEnvData.sbr_noise_levels[i],
                                   noise_level[1][i]
                   */
      
      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {

        
        h_envChan[1]->encEnvData.sbr_noise_levels[i] = noise_level[1][i];
      }


       
      sbrHeaderData->coupling = 0;
      h_envChan[0]->encEnvData.balance = 0;
      h_envChan[1]->encEnvData.balance = 0;

       
      payloadbitsLR = CountSbrChannelPairElement (sbrHeaderData,
                                                  sbrBitstreamData,
                                                  &h_envChan[0]->encEnvData,
                                                  &h_envChan[1]->encEnvData,
                                                  hCmonData);


       /* pointers for h_envChan[ch]->sbrCodeEnvelope,
                                   h_envChan[ch]->sbrCodeNoiseFloor,
                                   h_envChan[ch]->encEnvData.domain_vec,
                                   h_envChan[ch]->encEnvData.domain_vec_noise,
                                   upDateNrgTemp[ch],
                                   upDateNoiseTemp[ch],
                                   domainVecTemp[ch],
                                   domainVecNoiseTemp[ch]
                   */
      
      for(ch = 0; ch < nChannels;ch++){
        int   itmp;

         /* pointers for h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev[i],
                                     sfbNrgPrevTemp[ch][i]
                     */
        
        for(i=0;i<MAX_FREQ_COEFFS;i++){


           
          itmp =  h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev[i];
          h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev[i]=sfbNrgPrevTemp[ch][i];
          sfbNrgPrevTemp[ch][i]=itmp;
        }

         /* pointers for h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev[i],
                                     noisePrevTemp[ch][i]
                     */
        
        for(i=0;i<MAX_NUM_NOISE_COEFFS;i++){


           
          itmp =  h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev[i];
          h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev[i]=noisePrevTemp[ch][i];
          noisePrevTemp[ch][i]=itmp;
       }

         
        itmp  = h_envChan[ch]->sbrCodeEnvelope.upDate;
        h_envChan[ch]->sbrCodeEnvelope.upDate=upDateNrgTemp[ch];
        upDateNrgTemp[ch] = itmp;

         
        itmp =  h_envChan[ch]->sbrCodeNoiseFloor.upDate;
        h_envChan[ch]->sbrCodeNoiseFloor.upDate=upDateNoiseTemp[ch];
        upDateNoiseTemp[ch]=itmp;

            
        memcpy(domainVecTemp[ch],h_envChan[ch]->encEnvData.domain_vec,sizeof(int)*MAX_ENVELOPES);

            
        memcpy(domainVecNoiseTemp[ch],h_envChan[ch]->encEnvData.domain_vec_noise,sizeof(int)*MAX_ENVELOPES);


         
        if(!sbrHeaderData->prev_coupling){

          
          h_envChan[ch]->sbrCodeEnvelope.upDate = 0;
          h_envChan[ch]->sbrCodeNoiseFloor.upDate = 0;
        }
      } /* ch */



       
      codeEnvelope (sfb_nrg_coupling[0], frame_info[0]->freqRes,
                    &h_envChan[0]->sbrCodeEnvelope,
                    h_envChan[0]->encEnvData.domain_vec, 1,
                    frame_info[0]->nEnvelopes, 0,
                    sbrBitstreamData->HeaderActive);

       
      codeEnvelope (sfb_nrg_coupling[1], frame_info[1]->freqRes,
                    &h_envChan[1]->sbrCodeEnvelope,
                    h_envChan[1]->encEnvData.domain_vec, 1,
                    frame_info[1]->nEnvelopes, 1,
                    sbrBitstreamData->HeaderActive);


      c = 0;

       /* pointers for h_envChan[0]->encEnvData.noScfBands[i]
                                   h_envChan[0]->encEnvData.ienvelope[i][j]
                                   h_envChan[1]->encEnvData.ienvelope[i][j]
                                   sfb_nrg_coupling[0][c]
                                   sfb_nrg_coupling[1][c]
                   */
      
      for (i = 0; i < nEnvelopes[0]; i++) {

        
        for (j = 0; j < h_envChan[0]->encEnvData.noScfBands[i]; j++) {

          
          h_envChan[0]->encEnvData.ienvelope[i][j] = sfb_nrg_coupling[0][c];
          h_envChan[1]->encEnvData.ienvelope[i][j] = sfb_nrg_coupling[1][c];

          c++;
        }
      }



         
      codeEnvelope (noise_level_coupling[0], res,
                    &h_envChan[0]->sbrCodeNoiseFloor,
                    h_envChan[0]->encEnvData.domain_vec_noise, 1,
                    (frame_info[0]->nEnvelopes > 1 ? 2 : 1), 0,
                     sbrBitstreamData->HeaderActive);

       /* pointers for h_envChan[0]->encEnvData.sbr_noise_levels[i],
                                   noise_level_coupling[0][i]
                   */
      
      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {

        
        h_envChan[0]->encEnvData.sbr_noise_levels[i] = noise_level_coupling[0][i];
      }


         
      codeEnvelope (noise_level_coupling[1], res,
                    &h_envChan[1]->sbrCodeNoiseFloor,
                    h_envChan[1]->encEnvData.domain_vec_noise, 1,
                    (frame_info[1]->nEnvelopes > 1 ? 2 : 1), 1,
                    sbrBitstreamData->HeaderActive);

       /* pointers for h_envChan[1]->encEnvData.sbr_noise_levels[i],
                                   noise_level_coupling[1][i]
                   */
      
      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {

        
        h_envChan[1]->encEnvData.sbr_noise_levels[i] = noise_level_coupling[1][i];
      }

       
      sbrHeaderData->coupling = 1;

       
      h_envChan[0]->encEnvData.balance  = 0;
      h_envChan[1]->encEnvData.balance  = 1;

       
      tempFlagLeft  = h_envChan[0]->encEnvData.addHarmonicFlag;
      tempFlagRight = h_envChan[1]->encEnvData.addHarmonicFlag;

        
      payloadbitsCOUPLING =
        CountSbrChannelPairElement (sbrHeaderData,
                                    sbrBitstreamData,
                                    &h_envChan[0]->encEnvData,
                                    &h_envChan[1]->encEnvData,
                                    hCmonData);

       
      h_envChan[0]->encEnvData.addHarmonicFlag = tempFlagLeft;
      h_envChan[1]->encEnvData.addHarmonicFlag = tempFlagRight;

       
      if (payloadbitsCOUPLING < payloadbitsLR) {

          
          for(ch = 0; ch < nChannels;ch++){

                
            memcpy (sfb_nrg[ch], sfb_nrg_coupling[ch],
                  MAX_NUM_ENVELOPE_VALUES * sizeof (int));

                
            memcpy (noise_level[ch], noise_level_coupling[ch],
                  MAX_NUM_NOISE_VALUES * sizeof (int));
          }

           
          sbrHeaderData->coupling = 1;
          h_envChan[0]->encEnvData.balance  = 0;
          h_envChan[1]->encEnvData.balance  = 1;
      }
      else{
           /* pointer for h_envChan[ch]->sbrCodeEnvelope.upDate,
                                      h_envChan[ch]->sbrCodeNoiseFloor.upDate,
                                      upDateNrgTemp[ch],
                                      upDateNoiseTemp[ch]
                       */
          
          for(ch = 0; ch < nChannels;ch++){

                
            memcpy (h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev,
                    sfbNrgPrevTemp[ch], MAX_FREQ_COEFFS * sizeof (int));

            
            h_envChan[ch]->sbrCodeEnvelope.upDate = upDateNrgTemp[ch];

                
            memcpy (h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev,
                    noisePrevTemp[ch], MAX_NUM_NOISE_COEFFS * sizeof (int));

                
            memcpy(h_envChan[ch]->encEnvData.domain_vec,domainVecTemp[ch],sizeof(int)*MAX_ENVELOPES);

                
            memcpy(h_envChan[ch]->encEnvData.domain_vec_noise,domainVecNoiseTemp[ch],sizeof(int)*MAX_ENVELOPES);

            
            h_envChan[ch]->sbrCodeNoiseFloor.upDate = upDateNoiseTemp[ch];
          }

           
          sbrHeaderData->coupling = 0;
          h_envChan[0]->encEnvData.balance  = 0;
          h_envChan[1]->encEnvData.balance  = 0;
        }
    }
    break;

#endif /* (MAX_CHANNELS>1) */

  default:
    assert(0);

  } /* end switch(stereoMode) */

   
  if (nChannels == 1) {

      
    if (h_envChan[0]->encEnvData.domain_vec[0] == TIME) {

       
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac++;
    }
    else {

       
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
    }
  }
  else {

       
    if (h_envChan[0]->encEnvData.domain_vec[0] == TIME ||
        h_envChan[1]->encEnvData.domain_vec[0] == TIME) {

       
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac++;
      h_envChan[1]->sbrCodeEnvelope.dF_edge_incr_fac++;
    }
    else {

       
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
      h_envChan[1]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
    }
  }


   /* pointers for h_envChan[ch]->encEnvData.noScfBands[i],
                               h_envChan[ch]->encEnvData.ienvelope[i][j],
                               h_envChan[ch]->encEnvData.sbr_noise_levels[i],
                               sfb_nrg[ch][c]
                               noise_level[ch][i]
               */
  
  for(ch = 0; ch < nChannels;ch++){
    c = 0;

    
    for (i = 0; i < nEnvelopes[ch]; i++) {

      
      for (j = 0; j < h_envChan[ch]->encEnvData.noScfBands[i]; j++) {

        
        h_envChan[ch]->encEnvData.ienvelope[i][j] = sfb_nrg[ch][c];

        c++;
      }
    }

    
    for (i = 0; i < MAX_NUM_NOISE_VALUES; i++){

      
      h_envChan[ch]->encEnvData.sbr_noise_levels[i] = noise_level[ch][i];

    }

  }/* ch */

   
  if (nChannels == 2) {

      
    WriteEnvChannelPairElement(sbrHeaderData,
                               sbrBitstreamData,
                               &h_envChan[0]->encEnvData,
                               &h_envChan[1]->encEnvData,
                               hCmonData);
  }
  else {

      
    WriteEnvSingleChannelElement(sbrHeaderData,
                                 sbrBitstreamData,
                                 &h_envChan[0]->encEnvData,
                                 hPsEnc,
                                 hCmonData);
  }


   /* pointers for h_envChan[ch]->sbrExtractEnvelope.YBufferWriteOffset,
                               h_envChan[ch]->sbrExtractEnvelope,
                               h_envChan[ch]->sbrExtractEnvelope.YBuffer[]
               */
  
  for(ch = 0; ch < nChannels;ch++){

    
    for (i = 0; i < h_envChan[ch]->sbrExtractEnvelope.YBufferWriteOffset; i++) {
      float *temp;

      temp = h_envChan[ch]->sbrExtractEnvelope.YBuffer[i];

      
      h_envChan[ch]->sbrExtractEnvelope.YBuffer[i] = h_envChan[ch]->sbrExtractEnvelope.YBuffer[i + h_envChan[ch]->sbrExtractEnvelope.no_cols/2];
      h_envChan[ch]->sbrExtractEnvelope.YBuffer[i + h_envChan[ch]->sbrExtractEnvelope.no_cols/2] = temp;
    }
  }/* ch */

   
  sbrHeaderData->prev_coupling = sbrHeaderData->coupling;

  
}



/***************************************************************************/
/*!

  \brief  creates an envelope extractor handle

  \return error status

****************************************************************************/
int
CreateExtractSbrEnvelope (SBRRam_t *sbrram, int chan,
                          HANDLE_SBR_EXTRACT_ENVELOPE  hSbrCut,
                          int start_index)
{
  int i;
  int YBufferLength, rBufferLength;

  

      
  memset(hSbrCut,0,sizeof(SBR_EXTRACT_ENVELOPE));

   
  hSbrCut->YBufferWriteOffset = 32;

   
  YBufferLength = hSbrCut->YBufferWriteOffset + 32;

  
  rBufferLength = 32;

   
  hSbrCut->pre_transient_info[0] = 0;
  hSbrCut->pre_transient_info[1] = 0;

   
  hSbrCut->no_cols = 32;
  hSbrCut->no_rows = 64;
  hSbrCut->start_index = start_index;

   
  hSbrCut->time_slots = 16;
  hSbrCut->time_step = 2;

  assert(rBufferLength  ==   QMF_TIME_SLOTS);
  assert(YBufferLength  == 2*QMF_TIME_SLOTS);


  
  YBufferLength /= 2;

    
  hSbrCut->YBufferWriteOffset /= 2;

   /* hSbrCut->YBuffer[] */
   
  for (i = 0; i < YBufferLength; i++) {
      
    hSbrCut->YBuffer[i] = &sbrram->sbr_envYBuffer[chan*YBufferLength*64 + i*QMF_CHANNELS];
  }

   /* hSbrCut->rBuffer[]
                  hSbrCut->iBuffer[]
               */
   
  for (i = 0; i < rBufferLength; i++) {
      
    hSbrCut->rBuffer[i] = &sbrram->sbr_envRBuffer[chan*QMF_TIME_SLOTS*QMF_CHANNELS + i*QMF_CHANNELS];

      
    hSbrCut->iBuffer[i] = &sbrram->sbr_envIBuffer[chan*QMF_TIME_SLOTS*QMF_CHANNELS + i*QMF_CHANNELS];
  }


      
  memset(hSbrCut->envelopeCompensation,0,sizeof(char)*MAX_FREQ_COEFFS);


  

  return 0;
}




/***************************************************************************/
/*!

  \brief  deinitializes an envelope extractor handle

  \return void

****************************************************************************/
void
deleteExtractSbrEnvelope (HANDLE_SBR_EXTRACT_ENVELOPE hSbrCut)
{

  

  
  if (hSbrCut) {


  }

  
}



