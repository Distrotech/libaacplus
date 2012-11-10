/*
  Sbr QMF inverse filtering detector
*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */



#define MAX_NUM_REGIONS 10



#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif






static const float quantStepsSbr[4]  = {1, 10, 14, 19};
static const float quantStepsOrig[4] = {0,  3,  7, 10};
static const float nrgBorders[4] =     {25.0f, 30.0f, 35.0f, 40.0f};

static const DETECTOR_PARAMETERS detectorParamsAAC = {
    quantStepsSbr,
    quantStepsOrig,
    nrgBorders,
    4,
    4,
    4,
    {
      {INVF_MID_LEVEL,   INVF_LOW_LEVEL,  INVF_OFF,        INVF_OFF, INVF_OFF},
      {INVF_MID_LEVEL,   INVF_LOW_LEVEL,  INVF_OFF,        INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_MID_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF}
    },
    {
      {INVF_LOW_LEVEL,   INVF_LOW_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_LOW_LEVEL,   INVF_LOW_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_MID_LEVEL,  INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF}
    },
    {-4, -3, -2, -1, 0}
};

static const float hysteresis = 1.0f;


static const DETECTOR_PARAMETERS detectorParamsAACSpeech = {
    quantStepsSbr,
    quantStepsOrig,
    nrgBorders,
    4,
    4,
    4,
    {
      {INVF_MID_LEVEL,   INVF_MID_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_MID_LEVEL,   INVF_MID_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_MID_LEVEL,  INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF}
    },
    {
      {INVF_MID_LEVEL,   INVF_MID_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_MID_LEVEL,   INVF_MID_LEVEL,  INVF_LOW_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_MID_LEVEL,  INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF},
      {INVF_HIGH_LEVEL,  INVF_HIGH_LEVEL, INVF_MID_LEVEL,  INVF_OFF, INVF_OFF}
    },
    {-4, -3, -2, -1, 0}
};


typedef const float FIR_FILTER[5];

static const FIR_FILTER fir_0 = { 1.0f };
static const FIR_FILTER fir_1 = { 0.3333333f, 0.6666666f };
static const FIR_FILTER fir_2 = { 0.125f, 0.375f, 0.5f };
static const FIR_FILTER fir_3 = { 0.0585786f, 0.2f, 0.3414214f, 0.4f };
static const FIR_FILTER fir_4 = { 0.0318305f, 0.1151638f, 0.2181695f, 0.3015028f, 0.3333333f };

static const FIR_FILTER *fir_table[5] = {
  &fir_0,
  &fir_1,
  &fir_2,
  &fir_3,
  &fir_4
};



/**************************************************************************/
/*!
  \brief     Calculates the values used for the detector.

  \return    none

*/
/**************************************************************************/
static void
calculateDetectorValues(float **quotaMatrixOrig,
                        char * indexVector,
                        float* nrgVector,
                        DETECTOR_VALUES *detectorValues,
                        int startChannel,
                        int stopChannel,
                        int startIndex,
                        int stopIndex
                        )
{
  int i,j;

  float origQuota, sbrQuota;

  const float* filter = *fir_table[INVF_SMOOTHING_LENGTH];

  float quotaVecOrig[64], quotaVecSbr[64];

  

   /* counting previous operation */

      
  memset(quotaVecOrig,0,64*sizeof(float));

      
  memset(quotaVecSbr ,0,64*sizeof(float));


   
  detectorValues->avgNrg = 0;

   /* quotaVecOrig[],
                  indexVector[],
                  quotaVecSbr[]
               */
  
  for(j = startIndex ; j < stopIndex ; j++){

     /* quotaMatrixOrig[][] */
    
    for(i = startChannel; i < stopChannel; i++) {

       
      quotaVecOrig[i] += (quotaMatrixOrig[j][i]);

       
      if(indexVector[i] != -1)
      {
          
        quotaVecSbr[i] += (quotaMatrixOrig[j][indexVector[i]]);
      }
    }
      
    detectorValues->avgNrg += nrgVector[j];
  }

     
  detectorValues->avgNrg /= (stopIndex-startIndex);


   /* quotaVecOrig[]
                  quotaVecSbr[]
               */
  
  for(i = startChannel; i < stopChannel; i++) {

     
    quotaVecOrig[i] /= (stopIndex-startIndex);
    quotaVecSbr[i]  /= (stopIndex-startIndex);
  }

  
  origQuota = 0.0f;
  sbrQuota  = 0.0f;

   /* quotaVecOrig[],
                  quotaVecSbr[]
               */
  
  for(i = startChannel; i < stopChannel; i++) {

    
    origQuota += quotaVecOrig[i];
    sbrQuota  += quotaVecSbr[i];
  }

   
  origQuota /= (stopChannel - startChannel);
  sbrQuota  /= (stopChannel - startChannel);


      
  memmove(detectorValues->origQuotaMean, detectorValues->origQuotaMean + 1, INVF_SMOOTHING_LENGTH*sizeof(float));

      
  memmove(detectorValues->sbrQuotaMean, detectorValues->sbrQuotaMean + 1, INVF_SMOOTHING_LENGTH*sizeof(float));


   
  detectorValues->origQuotaMean[INVF_SMOOTHING_LENGTH]          = origQuota;
  detectorValues->sbrQuotaMean[INVF_SMOOTHING_LENGTH]           = sbrQuota;



   
  detectorValues->origQuotaMeanFilt = 0;
  detectorValues->sbrQuotaMeanFilt = 0;


   /* detectorValues->origQuotaMean[]
                  detectorValues->sbrQuotaMean[]
                  filter[]
               */
  for(i=0;i<INVF_SMOOTHING_LENGTH+1;i++) {
    
    detectorValues->origQuotaMeanFilt += detectorValues->origQuotaMean[i]*filter[i];
    detectorValues->sbrQuotaMeanFilt  += detectorValues->sbrQuotaMean[i]*filter[i];
  }
    /* store previous calculations */

  
}





/**************************************************************************/
/*!
  \brief     Returns the region in which the input value belongs.

  \return    region.

*/
/**************************************************************************/
static int
findRegion(float currVal,
           const float* borders,
           const int numBorders,
           int prevRegion
           )
{
  int i;

  

   
  if(currVal < borders[0])
  {
    
    return 0;
  }

   /* borders[] */
  
  for(i = 1; i < numBorders; i++){

      
    if( currVal >= borders[i-1] && currVal < borders[i])
    {
      
      return i;
    }
  }

   
  if(currVal > borders[numBorders-1])
  {
    
    return numBorders;
  }

  
  return 0;
}





/**************************************************************************/
/*!
  \brief     Makes a decision based on the quota vector.

  \return     decision on which invf mode to use

*/
/**************************************************************************/
static INVF_MODE
decisionAlgorithm(const DETECTOR_PARAMETERS* detectorParams,
                  DETECTOR_VALUES detectorValues,
                  int transientFlag,
                  INVF_MODE prevInvfMode,
                  int* prevRegionSbr,
                  int* prevRegionOrig
                  )
{
  int invFiltLevel, regionSbr, regionOrig, regionNrg;

  const float *quantStepsSbr  = detectorParams->quantStepsSbr;
  const float *quantStepsOrig = detectorParams->quantStepsOrig;
  const float *nrgBorders     = detectorParams->nrgBorders;
  const int numRegionsSbr     = detectorParams->numRegionsSbr;
  const int numRegionsOrig    = detectorParams->numRegionsOrig;
  const int numRegionsNrg     = detectorParams->numRegionsNrg;

  float quantStepsSbrTmp[MAX_NUM_REGIONS];
  float quantStepsOrigTmp[MAX_NUM_REGIONS];

  float origQuotaMeanFilt;
  float sbrQuotaMeanFilt;
  float nrg;

  

     /* counting previous operations */

    
  origQuotaMeanFilt = (float) (ILOG2*3.0*log(detectorValues.origQuotaMeanFilt+EPS));
  sbrQuotaMeanFilt  = (float) (ILOG2*3.0*log(detectorValues.sbrQuotaMeanFilt+EPS));
  nrg               = (float) (ILOG2*1.5*log(detectorValues.avgNrg+EPS));




      
  memcpy(quantStepsSbrTmp,quantStepsSbr,numRegionsSbr*sizeof(float));

      
  memcpy(quantStepsOrigTmp,quantStepsOrig,numRegionsOrig*sizeof(float));


   /* quantStepsSbrTmp[*prevRegionSbr]
                  quantStepsOrigTmp[*prevRegionOrig]
               */

   
  if(*prevRegionSbr < numRegionsSbr)
  {
     
    quantStepsSbrTmp[*prevRegionSbr] = quantStepsSbr[*prevRegionSbr] + hysteresis;
  }

  
  if(*prevRegionSbr > 0)
  {
     
    quantStepsSbrTmp[*prevRegionSbr - 1] = quantStepsSbr[*prevRegionSbr - 1] - hysteresis;
  }


  if(*prevRegionOrig < numRegionsOrig)
  {
     
    quantStepsOrigTmp[*prevRegionOrig] = quantStepsOrig[*prevRegionOrig] + hysteresis;
  }

  
  if(*prevRegionOrig > 0)
  {
     
    quantStepsOrigTmp[*prevRegionOrig - 1] = quantStepsOrig[*prevRegionOrig - 1] - hysteresis;
  }

  
  regionSbr  = findRegion(sbrQuotaMeanFilt, quantStepsSbrTmp, numRegionsSbr, *prevRegionSbr);

  
  regionOrig = findRegion(origQuotaMeanFilt, quantStepsOrigTmp, numRegionsOrig, *prevRegionOrig);

  
  regionNrg  = findRegion(nrg,nrgBorders,numRegionsNrg, 0);


  
  *prevRegionSbr = regionSbr;
  *prevRegionOrig = regionOrig;

   
  if(transientFlag == 1){

     
    invFiltLevel = detectorParams->regionSpaceTransient[regionSbr][regionOrig];
  }
  else{
     
    invFiltLevel = detectorParams->regionSpace[regionSbr][regionOrig];
  }

     
  invFiltLevel = max(invFiltLevel + detectorParams->EnergyCompFactor[regionNrg],0);

  

  return (INVF_MODE) (invFiltLevel);
}




/**************************************************************************/
/*!
  \brief     Estiamtion of the inverse filtering level required
             in the decoder.

 \return    none.

*/
/**************************************************************************/
void
qmfInverseFilteringDetector (HANDLE_SBR_INV_FILT_EST hInvFilt,
                             float ** quotaMatrix,
                             float *nrgVector,
                             char* indexVector,
                             int startIndex,
                             int stopIndex,
                             int transientFlag,
                             INVF_MODE* infVec
                             )
{
  int band;


  

   /* hInvFilt->freqBandTableInvFilt[band]
                  hInvFilt->detectorValues[band]
                  hInvFilt->prevInvfMode[band]
                  hInvFilt->prevRegionSbr[band]
                  hInvFilt->prevRegionOrig[band]
                  infVec[band]
               */
   
  for(band = 0 ; band < hInvFilt->noDetectorBands; band++){
    int startChannel = hInvFilt->freqBandTableInvFilt[band];
    int stopChannel  = hInvFilt->freqBandTableInvFilt[band+1];

     /* counting previous operations */

      
    calculateDetectorValues(quotaMatrix,
                            indexVector,
                            nrgVector,
                            &hInvFilt->detectorValues[band],
                            startChannel,
                            stopChannel,
                            startIndex,
                            stopIndex
                            );


       
    infVec[band]= decisionAlgorithm(hInvFilt->detectorParams,
                                    hInvFilt->detectorValues[band],
                                    transientFlag,
                                    hInvFilt->prevInvfMode[band],
                                    &hInvFilt->prevRegionSbr[band],
                                    &hInvFilt->prevRegionOrig[band]);

  }

  
}


/**************************************************************************/
/*!
  \brief     Creates an instance of the inverse filtering level estimator.

  \return   errorCode

*/
/**************************************************************************/
int
createInvFiltDetector (HANDLE_SBR_INV_FILT_EST hInvFilt,
                       int* freqBandTableDetector,
                       int numDetectorBands,
                       int numberOfEstimatesPerFrame,
                       unsigned int useSpeechConfig
                       )
{
  int i;

  


      
  memset( hInvFilt,0,sizeof(SBR_INV_FILT_EST));

  
  if (useSpeechConfig) {
     
    hInvFilt->detectorParams=&detectorParamsAACSpeech;
  }
  else {
     
    hInvFilt->detectorParams=&detectorParamsAAC;
  }

   
  hInvFilt->noDetectorBandsMax = numDetectorBands;


   /* hInvFilt->detectorValues[]
                  hInvFilt->prevInvfMode[]
                  hInvFilt->prevRegionOrig[]
                  hInvFilt->prevRegionSbr[]
               */
  
  for(i=0;i<hInvFilt->noDetectorBandsMax;i++){
        
    memset(&hInvFilt->detectorValues[i],0,sizeof(DETECTOR_VALUES));

    
    hInvFilt->prevInvfMode[i]   = INVF_OFF;
    hInvFilt->prevRegionOrig[i] = 0;
    hInvFilt->prevRegionSbr[i]  = 0;
  }


   
  resetInvFiltDetector(hInvFilt,
                       freqBandTableDetector,
                       hInvFilt->noDetectorBandsMax);


  

  return (0);
}


/**************************************************************************/
/*!
  \brief     resets sbr inverse filtering structure.

  \return   errorCode

*/
/**************************************************************************/
int
resetInvFiltDetector(HANDLE_SBR_INV_FILT_EST hInvFilt,
                     int* freqBandTableDetector,
                     int numDetectorBands)
{

  

   

      
  memcpy(hInvFilt->freqBandTableInvFilt,freqBandTableDetector,(numDetectorBands+1)*sizeof(int));

   
  hInvFilt->noDetectorBands = numDetectorBands;

  

  return (0);
}


/**************************************************************************/
/*!
  \brief    deletes sbr inverse filtering structure.

  \return   none

*/
/**************************************************************************/
void
deleteInvFiltDetector (HANDLE_SBR_INV_FILT_EST hs)
{
  

  /*
    nothing to do
  */

  
}
