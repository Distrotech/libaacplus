/*
  General tonality correction detector module.
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


/*!
  Auto correlation coefficients.
*/
typedef struct {
  float  r00r;
  float  r11r;
  float  r01r;
  float  r01i;
  float  r02r;
  float  r02i;
  float  r12r;
  float  r12i;
  float  r22r;
  float  det;
} ACORR_COEFS;



#define LPC_ORDER 2



/**************************************************************************/
/*!
  \brief     Calculates the second order auto correlation.

  Written by Per Ekstrand CTS AB.

  \return    none.

*/
/**************************************************************************/
static void
calcAutoCorrSecondOrder (ACORR_COEFS *ac,
                         float **realBuf,
                         float **imagBuf,
                         int bd
                               )
{
  int   j, jminus1, jminus2;
  float rel = 1.0f / (1.0f + RELAXATION);

  

   /* counting previous operations */

      
  memset(ac, 0, sizeof(ACORR_COEFS));

   /* pointer for realBuf[][],
                              imagBuf[][]
               */
  
  for ( j = 0; j < 14 - 1; j++ ) {

    jminus1 = j - 1;
    jminus2 = jminus1 - 1;

    
    ac->r00r += realBuf[j][bd] * realBuf[j][bd] +
                imagBuf[j][bd] * imagBuf[j][bd];


    
    ac->r11r += realBuf[jminus1][bd] * realBuf[jminus1][bd] +
                imagBuf[jminus1][bd] * imagBuf[jminus1][bd];

    
    ac->r01r += realBuf[j][bd] * realBuf[jminus1][bd] +
                imagBuf[j][bd] * imagBuf[jminus1][bd];

      
    ac->r01i += imagBuf[j][bd] * realBuf[jminus1][bd] -
                realBuf[j][bd] * imagBuf[jminus1][bd];

    
    ac->r02r += realBuf[j][bd] * realBuf[jminus2][bd] +
                imagBuf[j][bd] * imagBuf[jminus2][bd];

      
    ac->r02i += imagBuf[j][bd] * realBuf[jminus2][bd] -
                realBuf[j][bd] * imagBuf[jminus2][bd];

  }

   
  ac->r22r = ac->r11r + realBuf[-2][bd] * realBuf[-2][bd] +
                        imagBuf[-2][bd] * imagBuf[-2][bd];

   
  ac->r12r = ac->r01r + realBuf[-1][bd] * realBuf[-2][bd] +
                        imagBuf[-1][bd] * imagBuf[-2][bd];

    
  ac->r12i = ac->r01i + imagBuf[-1][bd] * realBuf[-2][bd] -
                        realBuf[-1][bd] * imagBuf[-2][bd];

  jminus1 = j - 1;
  jminus2 = jminus1 - 1;

  
  ac->r00r += realBuf[j][bd] * realBuf[j][bd] +
              imagBuf[j][bd] * imagBuf[j][bd];

  
  ac->r11r += realBuf[jminus1][bd] * realBuf[jminus1][bd] +
              imagBuf[jminus1][bd] * imagBuf[jminus1][bd];

  
  ac->r01r += realBuf[j][bd] * realBuf[jminus1][bd] +
              imagBuf[j][bd] * imagBuf[jminus1][bd];

    
  ac->r01i += imagBuf[j][bd] * realBuf[jminus1][bd] -
              realBuf[j][bd] * imagBuf[jminus1][bd];

  
  ac->r02r += realBuf[j][bd] * realBuf[jminus2][bd] +
              imagBuf[j][bd] * imagBuf[jminus2][bd];

  
  ac->r02i += imagBuf[j][bd] * realBuf[jminus2][bd] -
              realBuf[j][bd] * imagBuf[jminus2][bd];


    
  ac->det = ac->r11r * ac->r22r - rel * (ac->r12r * ac->r12r + ac->r12i * ac->r12i);

    /* move all register variables to the structure ac->... */

  
}


/**************************************************************************/
/*!
  \brief Calculates the tonal to noise ration

  \return none.

*/
/**************************************************************************/
void
CalculateTonalityQuotas(HANDLE_SBR_TON_CORR_EST hTonCorr,
                        float **sourceBufferReal,
                        float **sourceBufferImag,
                        int usb
                        )
{
  int    i, k, r, timeIndex;
  float  alphar[2], alphai[2], r01r, r02r, r11r, r12r, r22r, r01i, r02i,
          r12i, det, r00r;

  ACORR_COEFS ac;


  int     startIndexMatrix  = hTonCorr->startIndexMatrix;
  int     totNoEst          = hTonCorr->numberOfEstimates;
  int     noEstPerFrame     = hTonCorr->numberOfEstimatesPerFrame;
  int     move              = hTonCorr->move;
  int     noQmfChannels     = hTonCorr->noQmfChannels;
  float*  nrgVector         = hTonCorr->nrgVector;
  float** quotaMatrix       = hTonCorr->quotaMatrix;

  

  /*
   * Buffering of the quotaMatrix and the quotaMatrixTransp.
   *********************************************************/
  for(i =  0 ; i < move; i++){


        
    memcpy(quotaMatrix[i],
           quotaMatrix[i + noEstPerFrame],
           noQmfChannels * sizeof(float));
  }


      
  memmove(nrgVector,nrgVector+noEstPerFrame,move*sizeof(float));

      
  memset(nrgVector+startIndexMatrix,0,(totNoEst-startIndexMatrix)*sizeof(float));


  
  for (r = 0; r < usb; r++) {

     
    k = 2;
    timeIndex = startIndexMatrix;

     /* pointer for quotaMatrix[][]
                                nrgVector[]
                 */
     
    while(k <= 32 - 14) {

      
      r01r = r02r = r11r = r12r = r22r = r00r = 0;
      r01i = r02i = r12i = 0;

      
      calcAutoCorrSecondOrder (&ac, &sourceBufferReal[k], &sourceBufferImag[k], r);

      r00r = (float) ac.r00r;
      r11r = (float) ac.r11r;
      r12r = (float) ac.r12r;
      r12i = (float) ac.r12i;
      r22r = (float) ac.r22r;
      r01r = (float) ac.r01r;
      r01i = (float) ac.r01i;
      r02r = (float) ac.r02r;
      r02i = (float) ac.r02i;
      det  = (float) ac.det;



      
      if (ac.det == 0 ) {

        
        alphar[1] = alphai[1] = 0;
      } else {

          
        alphar[1] = ( r01r * r12r - r01i * r12i - r02r * r11r ) / det;

           
        alphai[1] = ( r01i * r12r + r01r * r12i - r02i * r11r ) / det;
      }

      
      if ( ac.r11r == 0 ) {

        
        alphar[0] =  alphai[0] = 0;
      } else {

           
        alphar[0] = - ( r01r + alphar[1] * r12r + alphai[1] * r12i ) / r11r;

          
        alphai[0] = - ( r01i + alphai[1] * r12r - alphar[1] * r12i ) / r11r;
      }


      
      if(r00r){
        float tmp;

          
        tmp =  - ( alphar[0]*r01r + alphai[0]*r01i + alphar[1]*r02r + alphai[1]*r02i) / (r00r);

         
        quotaMatrix[timeIndex][r] = tmp / ((float)1.0 - tmp +(float)RELAXATION);
      }
      else {

        
        quotaMatrix[timeIndex][r] = 0;
      }

      
      nrgVector[timeIndex] += r00r;

      k += 16;
      timeIndex++;
    }
  }

  
}



/**************************************************************************/
/*!
  \brief Extracts the parameters required in the decoder.

  \return none.

*/
/**************************************************************************/
void
TonCorrParamExtr(HANDLE_SBR_TON_CORR_EST hTonCorr,
                 INVF_MODE* infVec,
                 float* noiseLevels,
                 int* missingHarmonicFlag,
                 unsigned char * missingHarmonicsIndex,
                 char * envelopeCompensation,
                 const SBR_FRAME_INFO *frameInfo,
                 int* transientInfo,
                 unsigned char* freqBandTable,
                 int nSfb,
                 XPOS_MODE xposType
                 )
{
  int band;
  int transientFlag = transientInfo[1] ;
  int transientPos  = transientInfo[0];
  int transientFrame, transientFrameInvfEst;
  INVF_MODE* infVecPtr;

  

   /* counting previous operations */

  
  transientFrame = 0;

   
  if(hTonCorr->transientNextFrame){

     
    transientFrame = 1;
    hTonCorr->transientNextFrame = 0;

    
    if(transientFlag){

        
      if(transientPos + hTonCorr->transientPosOffset >= frameInfo->borders[frameInfo->nEnvelopes]){

         
        hTonCorr->transientNextFrame = 1;
      }
    }
  }
  else{

    
    if(transientFlag){

        
      if(transientPos + hTonCorr->transientPosOffset < frameInfo->borders[frameInfo->nEnvelopes]){

         
        transientFrame = 1;
        hTonCorr->transientNextFrame = 0;
      }
      else{

         
        hTonCorr->transientNextFrame = 1;
      }
    }
  }


   
  transientFrameInvfEst = hTonCorr->transientNextFrame;


   
  if (hTonCorr->switchInverseFilt)
  {
       
    qmfInverseFilteringDetector (&hTonCorr->sbrInvFilt,
                                 hTonCorr->quotaMatrix,
                                 hTonCorr->nrgVector,
                                 hTonCorr->indexVector,
                                 hTonCorr->frameStartIndexInvfEst,
                                 hTonCorr->numberOfEstimatesPerFrame + hTonCorr->frameStartIndexInvfEst,
                                 transientFrameInvfEst,
                                 infVec);
  }


   
  if (xposType == XPOS_LC ){

      
    SbrMissingHarmonicsDetectorQmf(&hTonCorr->sbrMissingHarmonicsDetector,
                                   hTonCorr->quotaMatrix,
                                   hTonCorr->indexVector,
                                   frameInfo,
                                   transientInfo,
                                   missingHarmonicFlag,
                                   missingHarmonicsIndex,
                                   freqBandTable,
                                   nSfb,
                                   envelopeCompensation);
  }
  else{

    
    *missingHarmonicFlag = 0;

        
    memset(missingHarmonicsIndex,0,nSfb*sizeof(int));
  }


   
  infVecPtr = hTonCorr->sbrInvFilt.prevInvfMode;

    
  sbrNoiseFloorEstimateQmf (&hTonCorr->sbrNoiseFloorEstimate,
                            frameInfo,
                            noiseLevels,
                            hTonCorr->quotaMatrix,
                            hTonCorr->indexVector,
                            *missingHarmonicFlag,
                            hTonCorr->frameStartIndex,
                            hTonCorr->numberOfEstimatesPerFrame,
                            hTonCorr->numberOfEstimates,
                            transientFrame,
                            infVecPtr);


   /* hTonCorr->sbrInvFilt.prevInvfMode[]
                  infVec[]
               */
   
  for(band = 0 ; band < hTonCorr->sbrInvFilt.noDetectorBands; band++){

    
    hTonCorr->sbrInvFilt.prevInvfMode[band] = infVec[band];
  }

  
}


/**************************************************************************/
/*!
  \brief     Searches for the closest match

  \return   closest entry.

*/
/**************************************************************************/
static int
findClosestEntry(int goalSb,
                 unsigned char *v_k_master,
                 int numMaster,
                 int direction)
{
  int index;

  

   
  if( goalSb <= v_k_master[0] )
  {
    
    return v_k_master[0];
  }

    
  if( goalSb >= v_k_master[numMaster] )
  {
    
    return v_k_master[numMaster];
  }

  
  if(direction) {

    
    index = 0;

     /* v_k_master[index] */
    
    while( v_k_master[index] < goalSb ) {

      
      index++;
    }
  } else {

    
    index = numMaster;

     /* v_k_master[index] */
    
    while( v_k_master[index] > goalSb ) {

      
      index--;
    }
  }

  

  return v_k_master[index];
}


/**************************************************************************/
/*!
  \brief     resets the patch

  \return   errorCode
*/
/**************************************************************************/
static int
resetPatch(HANDLE_SBR_TON_CORR_EST hTonCorr,
           int xposctrl,
           int highBandStartSb,
           int channelOffset,
           unsigned char *v_k_master,
           int numMaster,
           int fs,
           int noChannels)
{
  int patch,k,i;
  int targetStopBand;

  PATCH_PARAM  *patchParam = hTonCorr->patchParam;

  int sbGuard = hTonCorr->guard;
  int sourceStartBand;
  int patchDistance;
  int numBandsInPatch;

  int lsb = v_k_master[0];
  int usb = v_k_master[numMaster];
  int xoverOffset = highBandStartSb - v_k_master[0];

  int goalSb;

  

      /* counting previous operations */


   
  if (xposctrl == 1) {

    
    lsb += xoverOffset;

    
    xoverOffset = 0;
  }

    
  goalSb = (int)( 2 * noChannels * 16000.0f / fs  + 0.5f );

  
  goalSb = findClosestEntry(goalSb, v_k_master, numMaster, 1);

   
  sourceStartBand = hTonCorr->shiftStartSb + xoverOffset;

  
  targetStopBand = lsb + xoverOffset;

  
  patch = 0;

   /* patchParam[patch] */
  
  while(targetStopBand < usb) {

     
    if (patch >= MAX_NUM_PATCHES)
    {
      
      return(1);
    }

    
    patchParam[patch].guardStartBand = targetStopBand;

    
    targetStopBand += sbGuard;

    
    patchParam[patch].targetStartBand = targetStopBand;

    
    numBandsInPatch = goalSb - targetStopBand;

     
    if ( numBandsInPatch >= lsb - sourceStartBand ) {
      
      patchDistance   = targetStopBand - sourceStartBand;

      
      patchDistance   = patchDistance & ~1;

      
      numBandsInPatch = lsb - (targetStopBand - patchDistance);

       
      numBandsInPatch = findClosestEntry(targetStopBand + numBandsInPatch, v_k_master, numMaster, 0) -
                        targetStopBand;
    }


    
    patchDistance   = numBandsInPatch + targetStopBand - lsb;

     
    patchDistance   = (patchDistance + 1) & ~1;

    
    if (numBandsInPatch <= 0) {

      
      patch--;
    } else {

       
      patchParam[patch].sourceStartBand = targetStopBand - patchDistance;

      
      patchParam[patch].targetBandOffs  = patchDistance;
      patchParam[patch].numBandsInPatch = numBandsInPatch;

       
      patchParam[patch].sourceStopBand  = patchParam[patch].sourceStartBand + numBandsInPatch;

      
      targetStopBand += patchParam[patch].numBandsInPatch;
    }

     
    sourceStartBand = hTonCorr->shiftStartSb;

      
    if( abs(targetStopBand - goalSb) < 3) {
      
      goalSb = usb;
    }

    
    patch++;

  }

  
  patch--;


    
  if ( patchParam[patch].numBandsInPatch < 3 && patch > 0 ) {
    
    patch--;

    
    targetStopBand = patchParam[patch].targetStartBand + patchParam[patch].numBandsInPatch;
  }

    
  hTonCorr->noOfPatches = patch + 1;

   /* hTonCorr->indexVector[] */
   
  for(k = 0; k < hTonCorr->patchParam[0].guardStartBand; k++)
  {
    
    hTonCorr->indexVector[k] = k;
  }

   /* hTonCorr->patchParam[] */
   
  for(i = 0; i < hTonCorr->noOfPatches; i++)
  {
    int sourceStart    = hTonCorr->patchParam[i].sourceStartBand;
    int targetStart    = hTonCorr->patchParam[i].targetStartBand;
    int numberOfBands  = hTonCorr->patchParam[i].numBandsInPatch;
    int startGuardBand = hTonCorr->patchParam[i].guardStartBand;

     /* counting previous operations */

     /* hTonCorr->indexVector[] */
     
    for(k = 0; k < (targetStart- startGuardBand); k++)
    {
       
      hTonCorr->indexVector[startGuardBand+k] = -1;
    }

     /* hTonCorr->indexVector[] */
    
    for(k = 0; k < numberOfBands; k++)
    {
       
      hTonCorr->indexVector[targetStart+k] = sourceStart+k;
    }
  }

  

  return (0);
}



/**************************************************************************/
/*!
  \brief  Creates an instance of the tonality correction parameter module.

  \return   errorCode
*/
/**************************************************************************/
int
CreateTonCorrParamExtr (SBRRam_t *sbrram,
                        int chan,
                        HANDLE_SBR_TON_CORR_EST hTonCorr,
                        int fs,
                        int usb,
                        int noQmfChannels,
                        int xposCtrl,
                        int highBandStartSb,
                        int channelOffset,
                        unsigned char *v_k_master,
                        int numMaster,
                        int ana_max_level,
                        unsigned char *freqBandTable[2],
                        int* nSfb,
                        int noiseBands,
                        int noiseFloorOffset,
                        unsigned int useSpeechConfig
                        )
{
  int i;

  

      
  memset(hTonCorr,0,sizeof(SBR_TON_CORR_EST));

   
  hTonCorr->numberOfEstimates         = NO_OF_ESTIMATES;
  hTonCorr->numberOfEstimatesPerFrame = 2;
  hTonCorr->frameStartIndexInvfEst    = 0;
  hTonCorr->transientPosOffset        = 4;
  hTonCorr->move                      = 2;
  hTonCorr->startIndexMatrix          = 2;


   
  hTonCorr->frameStartIndex    = 0;
  hTonCorr->prevTransientFlag  = 0;
  hTonCorr->transientNextFrame = 0;

   
  hTonCorr->noQmfChannels = noQmfChannels;


    
   /* hTonCorr->quotaMatrix[]
                  sbr_quotaMatrix[]
               */
   
  for(i=0;i<hTonCorr->numberOfEstimates;i++) {
    
    hTonCorr->quotaMatrix[i] = &(sbrram->sbr_quotaMatrix[chan* NO_OF_ESTIMATES*QMF_CHANNELS + i*noQmfChannels]);

        
    memset(hTonCorr->quotaMatrix[i] ,0, sizeof(float)*QMF_CHANNELS);
  }

   /* Reset the patch.*/
   
  hTonCorr->guard = 0;
  hTonCorr->shiftStartSb = 1;

   
  if(resetPatch(hTonCorr,
                xposCtrl,
                highBandStartSb,
                channelOffset,
                v_k_master,
                numMaster,
                fs,
                noQmfChannels))
  {
    
    return(1);
  }

    
  if(CreateSbrNoiseFloorEstimate (&hTonCorr->sbrNoiseFloorEstimate,
                                   ana_max_level,
                                   freqBandTable[LO],
                                   nSfb[LO],
                                   noiseBands,
                                   noiseFloorOffset,
                                   //timeSlots,
                                   useSpeechConfig))
  {
    
    return(1);
  }


    
  if(createInvFiltDetector(&hTonCorr->sbrInvFilt,
                            hTonCorr->sbrNoiseFloorEstimate.freqBandTableQmf,
                            hTonCorr->sbrNoiseFloorEstimate.noNoiseBands,
                            hTonCorr->numberOfEstimatesPerFrame,
                            useSpeechConfig))
  {
    
    return(1);
  }


    
  if(CreateSbrMissingHarmonicsDetector (sbrram,
                                        chan,
                                        &hTonCorr->sbrMissingHarmonicsDetector,
                                        fs,
                                        freqBandTable[HI],
                                        nSfb[HI],
                                        noQmfChannels,
                                        hTonCorr->numberOfEstimates,
                                        hTonCorr->move,
                                        hTonCorr->numberOfEstimatesPerFrame))
  {
    
    return(1);
  }



  

  return (0);
}



/**************************************************************************/
/*!
  \brief  Deletes the tonality correction parameter module.

  \return   none
*/
/**************************************************************************/
void
DeleteTonCorrParamExtr (HANDLE_SBR_TON_CORR_EST hTonCorr)
{

  

  
  if (hTonCorr) {

     
   deleteInvFiltDetector (&hTonCorr->sbrInvFilt);

  }

  
}
