/*
  Missing harmonics detection
*/


#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define DELTA_TIME          9
#define MAX_COMP            2
#define TONALITY_QUOTA      0.1f
#define DIFF_QUOTA          0.75f
#define THR_DIFF           25.0f
#define THR_DIFF_GUIDE      1.26f
#define THR_TONE           15.0f
#define I_THR_TONE         (1.0f/15.0f)
#define THR_TONE_GUIDE      1.26f
#define THR_SFM_SBR         0.3f
#define THR_SFM_ORIG        0.1f
#define DECAY_GUIDE_ORIG    0.3f
#define DECAY_GUIDE_DIFF    0.5f



/**************************************************************************/
/*!
  \brief     Calculates the difference.

  \return    none.

*/
/**************************************************************************/
static void diff(float* pTonalityOrig,          /*!< */
                 float* pDiffMapped2Scfb,
                 const unsigned char* pFreqBandTable,
                 int nScfb,
                 char * indexVector)
{

  int i, ll, lu, k;
  float maxValOrig, maxValSbr;

  

   /* pointers for pFreqBandTable[],
                               pDiffMapped2Scfb[]
               */
  
  for(i=0; i < nScfb; i++){

    
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i+1];

    
    maxValOrig = 0;
    maxValSbr = 0;

     /* pointers for pTonalityOrig[],
                                 indexVector[]
                 */
    
    for(k=ll;k<lu;k++){

       
      if(pTonalityOrig[k] > maxValOrig)
      {
        
        maxValOrig = pTonalityOrig[k];
      }

        
      if(pTonalityOrig[indexVector[k]] > maxValSbr)
      {
        
        maxValSbr = pTonalityOrig[indexVector[k]];
      }
    }

     
    if(maxValSbr >= 1)
    {
       
      pDiffMapped2Scfb[i] = maxValOrig/maxValSbr;
    }
    else
    {
      
      pDiffMapped2Scfb[i] = maxValOrig;
    }
  }

  
}


/**************************************************************************/
/*!
  \brief     Calculates a flatness measure.
  \return    none.

*/
/**************************************************************************/
static void calculateFlatnessMeasure(float* pQuotaBuffer,
                                     char* indexVector,
                                     float* pSfmOrigVec,
                                     float* pSfmSbrVec,
                                     const unsigned char* pFreqBandTable,
                                     int nSfb)
{
  int i,j;

  

   /* pointers for pFreqBandTable[],
                               pSfmOrigVec[],
                               pSfmSbrVec[],
                               pSfmOrigVec[]
               */
  
  for(i=0;i<nSfb;i++)
  {
    int ll = pFreqBandTable[i];
    int lu = pFreqBandTable[i+1];

     /* counting previous operations */

    
    pSfmOrigVec[i] = 1;
    pSfmSbrVec[i]  = 1;

     
    if(lu -ll > 1){
      float amOrig,amTransp,gmOrig,gmTransp,sfmOrig,sfmTransp;

      
      amOrig = amTransp = 0;
      gmOrig = gmTransp = 1;

       /* pointers for pQuotaBuffer[],
                                   indexVector[]
                    */
      
      for(j= ll;j<lu;j++){

        
        sfmOrig   = pQuotaBuffer[j];

         
        sfmTransp = pQuotaBuffer[indexVector[j]];

         
        amOrig   += sfmOrig;
        gmOrig   *= sfmOrig;
        amTransp += sfmTransp;
        gmTransp *= sfmTransp;
      }

       /* lu-ll */ 
      amOrig   /= lu-ll;
      amTransp /= lu-ll;

       /* 1.0f/(lu-ll) */ 
      gmOrig   = (float) pow(gmOrig, 1.0f/(lu-ll));
      gmTransp = (float) pow(gmTransp,1.0f/(lu-ll));

      
      if(amOrig)
      {
         
        pSfmOrigVec[i] = gmOrig/amOrig;
      }

      
      if(amTransp)
      {
         
        pSfmSbrVec[i] = gmTransp/amTransp;
      }
    }
  }

  
}






/**************************************************************************/
/*!
  \brief     Calculates the input to the detection.
  \return    none.

*/
/**************************************************************************/
static void calculateDetectorInput(float** pQuotaBuffer,
                                   char* indexVector,
                                   float** tonalityDiff,
                                   float** pSfmOrig,
                                   float** pSfmSbr,
                                   const unsigned char* freqBandTable,
                                   int nSfb,
                                   int noEstPerFrame,
                                   int move,
                                   int noQmfBands)
{
  int est;

  


   
   for(est =  0 ; est < move; est++){

        
    memcpy(tonalityDiff[est],
           tonalityDiff[est + noEstPerFrame],
           noQmfBands * sizeof(float));

        
    memcpy(pSfmOrig[est],
           pSfmOrig[est + noEstPerFrame],
           noQmfBands * sizeof(float));

        
    memcpy(pSfmSbr[est],
           pSfmSbr[est + noEstPerFrame],
           noQmfBands * sizeof(float));

  }


   /* pointers for pQuotaBuffer[],
                               tonalityDiff[],
                               pSfmOrig[],
                               pSfmSbr[]
               */
  
  for(est = 0; est < noEstPerFrame; est++){

    
    diff(pQuotaBuffer[est+move],
         tonalityDiff[est+move],
         freqBandTable,
         nSfb,
         indexVector);

    
    calculateFlatnessMeasure(pQuotaBuffer[est+ move],
                             indexVector,
                             pSfmOrig[est + move],
                             pSfmSbr[est + move],
                             freqBandTable,
                             nSfb);
  }

  
}


/**************************************************************************/
/*!
  \brief     Checks if it is allowed to detect a tone.

  \return    newDetectionAllowed flag.

*/
/**************************************************************************/
static int isDetectionOfNewToneAllowed(const SBR_FRAME_INFO *pFrameInfo,
                                       int prevTransientFrame,
                                       int prevTransientPos,
                                       int prevTransientFlag,
                                       int transientPosOffset,
                                       int transientFlag,
                                       int transientPos,
                                       HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMissingHarmonicsDetector)
{
  int transientFrame, newDetectionAllowed;

  

  
  transientFrame = 0;

  
  if(transientFlag){

      
    if(transientPos + transientPosOffset < pFrameInfo->borders[pFrameInfo->nEnvelopes])
    {
      
      transientFrame = 1;
    }
  }
  else{
     
    if(prevTransientFlag && !prevTransientFrame){
      
      transientFrame = 1;
    }
  }


  
  newDetectionAllowed = 0;

  
  if(transientFrame){
    
    newDetectionAllowed = 1;
  }
  else {
        
    if(prevTransientFrame &&
       abs(pFrameInfo->borders[0] - (prevTransientPos + transientPosOffset -
                                     h_sbrMissingHarmonicsDetector->timeSlots)) < DELTA_TIME)
    {
      
      newDetectionAllowed = 1;
    }
  }

   
  h_sbrMissingHarmonicsDetector->previousTransientFlag  = transientFlag;
  h_sbrMissingHarmonicsDetector->previousTransientFrame = transientFrame;
  h_sbrMissingHarmonicsDetector->previousTransientPos   = transientPos;

  

  return (newDetectionAllowed);
}




/**************************************************************************/
/*!
  \brief     Cleans up after a transient.
  \return    none.

*/
/**************************************************************************/
static void transientCleanUp(float** quotaBuffer,
                             float** pDiffVecScfb,
                             int nSfb,
                             unsigned char ** detectionVectors,
                             const unsigned char* pFreqBandTable,
                             GUIDE_VECTORS guideVectors,
                             int newDetectionAllowed,
                             int start,
                             int stop)
{
  int i,j,li, ui,est;

  unsigned char pHarmVec[MAX_FREQ_COEFFS];

  

      
  memset(pHarmVec,0,MAX_FREQ_COEFFS*sizeof(unsigned char));

   /* pHarmVec[] */
  
  for(est = start; est < stop; est++){

     /* detectionVectors[][] */
    
    for(i=0;i<nSfb-1;i++){

       
      pHarmVec[i] = pHarmVec[i] || detectionVectors[est][i];
    }
  }


   /* pFreqBandTable[]
                  pHarmVec[]
                  guideVectors.guideVectorDetected[]
                  guideVectors.guideVectorOrig[]
                  guideVectors.guideVectorDiff[]
                */
  
  for(i=0;i<nSfb-1;i++)
  {
    
    li = pFreqBandTable[i];
    ui = pFreqBandTable[i+1];


     
    if(pHarmVec[i] && pHarmVec[i+1]){
      float maxVal1, maxVal2;
      int maxPos1, maxPos2;

      
      li = pFreqBandTable[i];
      ui = pFreqBandTable[i+1];


      
      maxPos1 = li;

       
      maxVal1 = quotaBuffer[start][li];

       /* quotaBuffer[][] */
      
      for(j = li; j<ui; j++){

         
        if(quotaBuffer[start][j] > maxVal1){

          
          maxVal1 = quotaBuffer[start][j];
          maxPos1 = j;
        }
      }

      
      for(est = start + 1; est < stop; est++){

         /* quotaBuffer[][] */
        
        for(j = li; j<ui; j++){

           
          if(quotaBuffer[est][j] > maxVal1){

            
            maxVal1 = quotaBuffer[est][j];
            maxPos1 = j;
          }
        }
      }


      
      li = pFreqBandTable[i+1];
      ui = pFreqBandTable[i+2];


      
      maxPos2 = li;

       
      maxVal2 = quotaBuffer[start][li];

       /* quotaBuffer[][] */
      
      for(j = li; j<ui; j++){

         
        if((float)quotaBuffer[start][j] > maxVal2){

          
          maxVal2 = quotaBuffer[start][j];
          maxPos2 = j;
        }
      }

      
      for(est = start + 1; est < stop; est++){

         /* quotaBuffer[][] */
        
        for(j = li; j<ui; j++){

           
          if((float)quotaBuffer[est][j] > maxVal2){

            
            maxVal2 = quotaBuffer[est][j];
            maxPos2 = j;
          }
        }
      }


       
      if(maxPos2-maxPos1 < 2){

         
        if(maxVal1 > maxVal2){

          
          guideVectors.guideVectorDetected[i+1] = 0;
          guideVectors.guideVectorOrig[i+1] = 0;
          guideVectors.guideVectorDiff[i+1] = 0;

           /* detectionVectors[][] */
          
          for(est = start; est<stop; est++){
            
            detectionVectors[est][i+1] = 0;
          }
        }
        else{

          
          guideVectors.guideVectorDetected[i] = 0;
          guideVectors.guideVectorOrig[i] = 0;
          guideVectors.guideVectorDiff[i] = 0;

           /* detectionVectors[][] */
          
          for(est = start; est<stop; est++){
            
            detectionVectors[est][i] = 0;
          }
        }
      }
    }
  }
  
}

/*************************************************************************/
/*!
  \brief     Do detection.
  \return    none.

*/
/**************************************************************************/
static void detection(float* quotaBuffer,
                      char* indexVector,
                      float* pDiffVecScfb,
                      int nSfb,
                      unsigned char* pHarmVec,
                      const unsigned char* pFreqBandTable,
                      float* sfmOrig,
                      float* sfmSbr,
                      GUIDE_VECTORS guideVectors,
                      GUIDE_VECTORS newGuideVectors,
                      int newDetectionAllowed)
{

  int i,j,ll, lu;
  float thresTemp,thresOrig;

  

   /* pointers for guideVectors.guideVectorDiff[],
                               guideVectors.guideVectorOrig[],
                               pDiffVecScfb[],
                               pHarmVec[],
                               newGuideVectors.guideVectorDiff[]
               */
  
  for(i=0;i<nSfb;i++)
  {
    
    if (guideVectors.guideVectorDiff[i]) {

         
      thresTemp = max(DECAY_GUIDE_DIFF*guideVectors.guideVectorDiff[i],THR_DIFF_GUIDE);
    } else {

      
      thresTemp = THR_DIFF;
    }

      
    thresTemp = min(thresTemp, THR_DIFF);

     
    if(pDiffVecScfb[i] > thresTemp){

      
      pHarmVec[i] = 1;

      
      newGuideVectors.guideVectorDiff[i] = pDiffVecScfb[i];
    }
    else{

      
      if(guideVectors.guideVectorDiff[i]){

        
        guideVectors.guideVectorOrig[i] = THR_TONE_GUIDE ;
      }
    }
  }


   /* pointers for guideVectors.guideVectorOrig[],
                               pHarmVec[],
                               newGuideVectors.guideVectorDiff[],
                               pFreqBandTable[]
               */
  
  for(i=0;i<nSfb;i++){

    
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i+1];

       
    thresOrig   = max(guideVectors.guideVectorOrig[i]*DECAY_GUIDE_ORIG,THR_TONE_GUIDE );

      
    thresOrig   = min(thresOrig,THR_TONE);

    
    if(guideVectors.guideVectorOrig[i]){

       /* pointers for quotaBuffer[] */
      
      for(j= ll;j<lu;j++){

         
        if(quotaBuffer[j] > thresOrig){

          
          pHarmVec[i] = 1;
          newGuideVectors.guideVectorOrig[i] = quotaBuffer[j];
        }
      }
    }
  }


  
  thresOrig   = THR_TONE;

   /* pointers for newGuideVectors.guideVectorDiff[],
                               pHarmVec[],
                               pFreqBandTable[],
                               pDiffVecScfb[],
                               sfmSbr[],
                               sfmOrig[]
               */
  
  for(i=0;i<nSfb;i++){

    
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i+1];

     
    if(lu -ll > 1){

       /* pointers for quotaBuffer[] */
      
      for(j= ll;j<lu;j++){

          
        if(quotaBuffer[j] > thresOrig && (sfmSbr[i] > THR_SFM_SBR && sfmOrig[i] < THR_SFM_ORIG)){

          
          pHarmVec[i] = 1;
          newGuideVectors.guideVectorOrig[i] = quotaBuffer[j];
        }
      }
    }
    else{
       
      if(i < nSfb -1){

        
        ll = pFreqBandTable[i];

         /* pointers for quotaBuffer[] */

        
        if(i>0){

            
          if(quotaBuffer[ll] > THR_TONE &&
             (pDiffVecScfb[+1] < I_THR_TONE ||
              pDiffVecScfb[i-1] < I_THR_TONE)
             ){

              
              pHarmVec[i] = 1;
              newGuideVectors.guideVectorOrig[i] = quotaBuffer[ll];
          }
        }
        else{
            
          if(quotaBuffer[ll] > THR_TONE &&
             pDiffVecScfb[i+1] < I_THR_TONE){

              
              pHarmVec[i] = 1;
              newGuideVectors.guideVectorOrig[i] = quotaBuffer[ll];
          }
        }
      }
    }
  }

  
}








/**************************************************************************/
/*!
  \brief     Do detection for every estimate
  \return    none.

*/
/**************************************************************************/
static void detectionWithPrediction(float** quotaBuffer,
                                    char* indexVector,
                                    float** pDiffVecScfb,
                                    int nSfb,
                                    const unsigned char* pFreqBandTable,
                                    float** sfmOrig,
                                    float** sfmSbr,
                                    unsigned char ** detectionVectors,
                                    unsigned char* prevFrameSfbHarm,
                                    GUIDE_VECTORS* guideVectors,
                                    int noEstPerFrame,
                                    int totNoEst,
                                    int newDetectionAllowed,
                                    unsigned char* pAddHarmonicsScaleFactorBands)
{
  int est = 0,i;
  int start;

  

   /* counting previous operation */

      
  memset(pAddHarmonicsScaleFactorBands,0,nSfb*sizeof(unsigned char));

  
  if(newDetectionAllowed){

     
    if(totNoEst > 1){

      
      start = noEstPerFrame;

          
      memcpy(guideVectors[noEstPerFrame].guideVectorDiff,guideVectors[0].guideVectorDiff,nSfb*sizeof(float));

          
      memcpy(guideVectors[noEstPerFrame].guideVectorOrig,guideVectors[0].guideVectorOrig,nSfb*sizeof(float));

          
      memset(guideVectors[noEstPerFrame-1].guideVectorDetected,0,nSfb*sizeof(unsigned char));
    }
    else{

      
      start = 0;
    }
  }
  else{

    
    start = 0;
  }


   /* guideVectors[],
                  detectionVectors[]
                  quotaBuffer[]
                  pDiffVecScfb[]
                  sfmOrig[]
                  sfmSbr[]
               */
  
  for(est = start; est < totNoEst; est++){

    
    if(est > 0){
          
      memcpy(guideVectors[est].guideVectorDetected,detectionVectors[est-1],nSfb*sizeof(unsigned char));
    }

        
    memset(detectionVectors[est], 0, nSfb*sizeof(unsigned char));

     
    if(est < totNoEst-1){
          
      memset(guideVectors[est+1].guideVectorDiff,0,nSfb*sizeof(float));

          
      memset(guideVectors[est+1].guideVectorOrig,0,nSfb*sizeof(float));

          
      memset(guideVectors[est+1].guideVectorDetected,0,nSfb*sizeof(unsigned char));

      
      detection(quotaBuffer[est],
                indexVector,
                pDiffVecScfb[est],
                nSfb,
                detectionVectors[est],
                pFreqBandTable,
                sfmOrig[est],
                sfmSbr[est],
                guideVectors[est],
                guideVectors[est+1],
                newDetectionAllowed);
    }
    else{
          
      memset(guideVectors[est].guideVectorDiff,0,nSfb*sizeof(float));

          
      memset(guideVectors[est].guideVectorOrig,0,nSfb*sizeof(float));

          
      memset(guideVectors[est].guideVectorDetected,0,nSfb*sizeof(unsigned char));

      
      detection(quotaBuffer[est],
                indexVector,
                pDiffVecScfb[est],
                nSfb,
                detectionVectors[est],
                pFreqBandTable,
                sfmOrig[est],
                sfmSbr[est],
                guideVectors[est],
                guideVectors[est],
                newDetectionAllowed);
    }


  }



  
  if(newDetectionAllowed){

     
    if(totNoEst > 1){

       
      transientCleanUp(quotaBuffer,
                       pDiffVecScfb,
                       nSfb,
                       detectionVectors,
                       pFreqBandTable,
                       guideVectors[noEstPerFrame],
                       newDetectionAllowed,
                       start,
                       totNoEst);
    }
    else{

       
       transientCleanUp(quotaBuffer,
                       pDiffVecScfb,
                       nSfb,
                       detectionVectors,
                       pFreqBandTable,
                       guideVectors[0],
                       newDetectionAllowed,
                       start,
                       totNoEst);
    }
  }



   /* pAddHarmonicsScaleFactorBands[] */
  
  for(i = 0; i< nSfb; i++){

     /* detectionVectors[][] */
    
    for(est = start; est < totNoEst; est++){

       
      pAddHarmonicsScaleFactorBands[i] = pAddHarmonicsScaleFactorBands[i] || detectionVectors[est][i];
    }
  }


  
  if(!newDetectionAllowed){
     /* pAddHarmonicsScaleFactorBands[]
                    prevFrameSfbHarm[]
                 */
    
    for(i=0;i<nSfb;i++){

       
      if(pAddHarmonicsScaleFactorBands[i] - prevFrameSfbHarm[i] > 0)
      {

        
        pAddHarmonicsScaleFactorBands[i] = 0;
      }
    }
  }

  
}


/**************************************************************************/
/*!
  \brief     Calculates a compensation vector.
  \return    none.

*/
/**************************************************************************/
static void calculateCompVector(unsigned char* pAddHarmonicsScaleFactorBands, /*!<  */
                                float** tonality,
                                char* envelopeCompensation,
                                int nSfb,
                                const unsigned char* freqBandTable,
                                float** diff,
                                int totNoEst,
                                char *prevEnvelopeCompensation,
                                int newDetectionAllowed)
{
  int i,j,l,ll,lu,maxPosF,maxPosT;
  float maxVal;
  int compValue;

  



      
  memset(envelopeCompensation,0,nSfb*sizeof(char));

   /* pAddHarmonicsScaleFactorBands[]
                  freqBandTable[]
                  envelopeCompensation[]
               */
  
  for(i=0 ; i < nSfb; i++){

    
    if(pAddHarmonicsScaleFactorBands[i]){ /* A missing sine was detected*/

      
      ll = freqBandTable[i];
      lu = freqBandTable[i+1];

      
      maxPosF = 0;
      maxPosT = 0;
      maxVal = 0;

      
      for(j=0;j<totNoEst;j++){

         /* tonality[][] */
        
        for(l=ll; l<lu; l++){

           
          if(tonality[j][l] > maxVal)
          {
            
            maxVal = tonality[j][l];
            maxPosF = l;
            maxPosT = j;
          }
        }
      }


        
      if(maxPosF == ll && i){

            
        compValue = (int) (fabs(ILOG2*log(diff[maxPosT][i - 1]+EPS)) + 0.5f);

         
        if (compValue > MAX_COMP)
        {
          
          compValue = MAX_COMP;
        }

        
        if(!pAddHarmonicsScaleFactorBands[i-1]) {

             
          if(tonality[maxPosT][maxPosF -1] > TONALITY_QUOTA*tonality[maxPosT][maxPosF]){

             
            envelopeCompensation[i-1] = -1*compValue;
          }
        }
      }

       
      if(maxPosF == lu-1 && i+1 < nSfb){

            
        compValue = (int) (fabs(ILOG2*log(diff[maxPosT][i + 1]+EPS)) + 0.5f);

         
        if (compValue > MAX_COMP)
        {
          
          compValue = MAX_COMP;
        }

        
        if(!pAddHarmonicsScaleFactorBands[i+1]) {

             
          if(tonality[maxPosT][maxPosF+1] > TONALITY_QUOTA*tonality[maxPosT][maxPosF]){

            
            envelopeCompensation[i+1] = compValue;
          }
        }
      }

        
      if(i && i < nSfb - 1){

            
        compValue = (int) (fabs(ILOG2*log(diff[maxPosT][i -1]+EPS)) + 0.5f);

         
        if (compValue > MAX_COMP)
        {
          
          compValue = MAX_COMP;
        }

           
        if((float)1.0f/diff[maxPosT][i-1] > DIFF_QUOTA*diff[maxPosT][i]){

           
          envelopeCompensation[i-1] = -1*compValue;
        }

            
        compValue = (int) (fabs(ILOG2*log(diff[maxPosT][i + 1]+EPS)) + 0.5f);

         
        if (compValue > MAX_COMP)
        {
          
          compValue = MAX_COMP;
        }

           
        if((float)1.0f/diff[maxPosT][i+1] > DIFF_QUOTA*diff[maxPosT][i]){

          
          envelopeCompensation[i+1] = compValue;
        }
      }
    }
  }



  
  if(!newDetectionAllowed){

     /* envelopeCompensation[]
                    prevEnvelopeCompensation[]
                 */
    
    for(i=0;i<nSfb;i++){

       
      if(envelopeCompensation[i] != 0 && prevEnvelopeCompensation[i] == 0)
      {
        
        envelopeCompensation[i] = 0;
      }
    }
  }

  
}



/**************************************************************************/
/*!
  \brief     Detects missing harmonics in the QMF
  \return    none.

*/
/**************************************************************************/
void
SbrMissingHarmonicsDetectorQmf(HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMHDet,
                               float ** pQuotaBuffer,
                               char* indexVector,
                               const SBR_FRAME_INFO *pFrameInfo,
                               const int* pTranInfo,
                               int* pAddHarmonicsFlag,
                               unsigned char* pAddHarmonicsScaleFactorBands,
                               const unsigned char* freqBandTable,
                               int nSfb,
                               char* envelopeCompensation)
{

  int i;
  int transientFlag = pTranInfo[1];
  int transientPos  = pTranInfo[0];
  int newDetectionAllowed;


  unsigned char ** detectionVectors  = h_sbrMHDet->detectionVectors;
  int move                = h_sbrMHDet->move;
  int noEstPerFrame       = h_sbrMHDet->noEstPerFrame;
  int totNoEst            = h_sbrMHDet->totNoEst;
  float** sfmSbr          = h_sbrMHDet->sfmSbr;
  float** sfmOrig         = h_sbrMHDet->sfmOrig;
  float** tonalityDiff    = h_sbrMHDet->tonalityDiff;
  int prevTransientFlag   = h_sbrMHDet->previousTransientFlag;
  int prevTransientFrame  = h_sbrMHDet->previousTransientFrame;
  int transientPosOffset  = h_sbrMHDet->transientPosOffset;
  int prevTransientPos    = h_sbrMHDet->previousTransientPos;
  GUIDE_VECTORS* guideVectors = h_sbrMHDet->guideVectors;

  int noQmfBands =  freqBandTable[nSfb] - freqBandTable[0];

  

     /* counting previous operation */

  
  newDetectionAllowed = isDetectionOfNewToneAllowed(pFrameInfo,
                                                    prevTransientFrame,
                                                    prevTransientPos,
                                                    prevTransientFlag,
                                                    transientPosOffset,
                                                    transientFlag,
                                                    transientPos,
                                                    h_sbrMHDet);


  
  calculateDetectorInput(pQuotaBuffer,
                         indexVector,
                         tonalityDiff,
                         sfmOrig,
                         sfmSbr,
                         freqBandTable,
                         nSfb,
                         noEstPerFrame,
                         move,
                         noQmfBands);


   
  detectionWithPrediction(pQuotaBuffer,
                          indexVector,
                          tonalityDiff,
                          nSfb,
                          freqBandTable,
                          sfmOrig,
                          sfmSbr,
                          detectionVectors,
                          h_sbrMHDet->guideScfb,
                          guideVectors,
                          noEstPerFrame,
                          totNoEst,
                          newDetectionAllowed,
                          pAddHarmonicsScaleFactorBands);



   
  calculateCompVector(pAddHarmonicsScaleFactorBands,
                      pQuotaBuffer,
                      envelopeCompensation,
                      nSfb,
                      freqBandTable,
                      tonalityDiff,
                      totNoEst,
                      h_sbrMHDet->prevEnvelopeCompensation,
                      newDetectionAllowed);



  
  *pAddHarmonicsFlag = 0;

   /* pAddHarmonicsScaleFactorBands[i] */
  
  for(i=0;i<nSfb;i++){

    
    if(pAddHarmonicsScaleFactorBands[i]){

      
      *pAddHarmonicsFlag = 1;
      break;
    }
  }


      
  memcpy(h_sbrMHDet->prevEnvelopeCompensation, envelopeCompensation, nSfb*sizeof(char));

      
  memcpy(h_sbrMHDet->guideScfb, pAddHarmonicsScaleFactorBands, nSfb*sizeof(unsigned char));

      
  memcpy(guideVectors[0].guideVectorDetected,pAddHarmonicsScaleFactorBands,nSfb*sizeof(unsigned char));

   
  if(totNoEst > noEstPerFrame){
        
    memcpy(guideVectors[0].guideVectorDiff,guideVectors[noEstPerFrame].guideVectorDiff,nSfb*sizeof(float));

        
    memcpy(guideVectors[0].guideVectorOrig,guideVectors[noEstPerFrame].guideVectorOrig,nSfb*sizeof(float));
  }
  else {
        
    memcpy(guideVectors[0].guideVectorDiff,guideVectors[noEstPerFrame-1].guideVectorDiff,nSfb*sizeof(float));

        
    memcpy(guideVectors[0].guideVectorOrig,guideVectors[noEstPerFrame-1].guideVectorOrig,nSfb*sizeof(float));
  }

   /* guideVectors[0].guideVectorDiff[]
                  guideVectors[0].guideVectorOrig[]
                  pAddHarmonicsScaleFactorBands[]
               */
  
  for(i=0;i<nSfb;i++){

     
    if((guideVectors[0].guideVectorDiff[i] ||
        guideVectors[0].guideVectorOrig[i]) &&
      !pAddHarmonicsScaleFactorBands[i]){

      
      guideVectors[0].guideVectorDiff[i] = 0;
      guideVectors[0].guideVectorOrig[i] = 0;
    }
  }

  
}

/**************************************************************************/
/*!
  \brief     Creates an instance of the missing harmonics detector.
  \return    errorCode, noError if OK.

*/
/**************************************************************************/
int
CreateSbrMissingHarmonicsDetector (SBRRam_t *sbrram,
                                   int chan,
                                   HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMHDet,
                                   int sampleFreq,
                                   unsigned char* freqBandTable,
                                   int nSfb,
                                   int qmfNoChannels,
                                   int totNoEst,
                                   int move,
                                   int noEstPerFrame)
{
  int i;
  float* ptr;

  HANDLE_SBR_MISSING_HARMONICS_DETECTOR hs =hSbrMHDet;

  

   /* counting previous operation */

  assert(totNoEst == NO_OF_ESTIMATES);

      
  memset(hs,0,sizeof( SBR_MISSING_HARMONICS_DETECTOR));


   
  hs->transientPosOffset =  4;
  hs->timeSlots          = 16;


   
  hs->qmfNoChannels = qmfNoChannels;
  hs->sampleFreq = sampleFreq;
  hs->nSfb = nSfb;

   
  hs->totNoEst = totNoEst;
  hs->move = move;
  hs->noEstPerFrame = noEstPerFrame;

    
  ptr = &sbrram->sbr_toncorrBuff[chan*5*NO_OF_ESTIMATES*MAX_FREQ_COEFFS];

   /* hs->tonalityDiff[i]
                  hs->sfmOrig[i]
                  hs->sfmSbr[i]
                  hs->guideVectors[i]
                  hs->detectionVectors[i]
               */
  
  for(i=0; i < totNoEst; i++) {

     
    hs->tonalityDiff[i] = ptr; ptr += MAX_FREQ_COEFFS;

        
    memset(hs->tonalityDiff[i],0,sizeof(float)*MAX_FREQ_COEFFS);

     
    hs->sfmOrig[i] = ptr; ptr += MAX_FREQ_COEFFS;

        
    memset(hs->sfmOrig[i],0,sizeof(float)*MAX_FREQ_COEFFS);

     
    hs->sfmSbr[i]  = ptr; ptr += MAX_FREQ_COEFFS;

        
    memset(hs->sfmSbr[i],0,sizeof(float)*MAX_FREQ_COEFFS);

     
    hs->guideVectors[i].guideVectorDiff = ptr; ptr += MAX_FREQ_COEFFS;

        
    memset(hs->guideVectors[i].guideVectorDiff,0,sizeof(float)*MAX_FREQ_COEFFS);

     
    hs->guideVectors[i].guideVectorOrig = ptr; ptr += MAX_FREQ_COEFFS;

        
    memset(hs->guideVectors[i].guideVectorOrig,0,sizeof(float)*MAX_FREQ_COEFFS);

       
    hs->detectionVectors[i] = &(sbrram->sbr_detectionVectors[chan*NO_OF_ESTIMATES*MAX_FREQ_COEFFS + i*MAX_FREQ_COEFFS]);

        
    memset(hs->detectionVectors[i],0,sizeof(unsigned char)*MAX_FREQ_COEFFS);

       
    hs->guideVectors[i].guideVectorDetected = &(sbrram->sbr_guideVectorDetected[chan*NO_OF_ESTIMATES*MAX_FREQ_COEFFS + i*MAX_FREQ_COEFFS]);

        
    memset(hs->guideVectors[i].guideVectorDetected,0,sizeof(unsigned char)*MAX_FREQ_COEFFS);

  }

    
  hs->prevEnvelopeCompensation = &(sbrram->sbr_prevEnvelopeCompensation[chan*MAX_FREQ_COEFFS]);

      
  memset( hs->prevEnvelopeCompensation,0, sizeof(char)*MAX_FREQ_COEFFS);

    
  hs->guideScfb = &(sbrram->sbr_guideScfb[chan*MAX_FREQ_COEFFS]);

      
  memset( hs->guideScfb,0, sizeof(unsigned char)*MAX_FREQ_COEFFS);

   
  hs->previousTransientFlag = 0;
  hs->previousTransientFrame = 0;
  hs->previousTransientPos = 0;

  assert (ptr-&sbrram->sbr_toncorrBuff[0] <= 5*MAX_CHANNELS*NO_OF_ESTIMATES*MAX_FREQ_COEFFS);

  

  return 0;
}

