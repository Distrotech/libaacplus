/*
  Noise floor estimation
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static const float smoothFilter[4]  = {0.05857864376269f, 0.2f, 0.34142135623731f, 0.4f};

#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif


/**************************************************************************/
/*!
  \brief     The function applies smoothing to the noise levels.
  \return    none

*/
/**************************************************************************/
static void
smoothingOfNoiseLevels (float *NoiseLevels,
                        int nEnvelopes,
                        int noNoiseBands,
                        float prevNoiseLevels[NF_SMOOTHING_LENGTH][MAX_NUM_NOISE_VALUES],
                        const float* smoothFilter)
{
  int i,band,env;

  

  
  for(env = 0; env < nEnvelopes; env++){

    
    for (i = 1; i < NF_SMOOTHING_LENGTH; i++){
          
      memcpy(prevNoiseLevels[i - 1],prevNoiseLevels[i],noNoiseBands*sizeof(float));
    }

        
    memcpy(prevNoiseLevels[NF_SMOOTHING_LENGTH - 1],NoiseLevels+env*noNoiseBands,noNoiseBands*sizeof(float));

      /* pointer for NoiseLevels[] */
    
    for (band = 0; band < noNoiseBands; band++){

      
      NoiseLevels[band+ env*noNoiseBands] = 0;

       /* pointer for smoothFilter[],
                                  prevNoiseLevels[][]
                   */
      
      for (i = 0; i < NF_SMOOTHING_LENGTH; i++){

        
        NoiseLevels[band+ env*noNoiseBands] += smoothFilter[i]*prevNoiseLevels[i][band];
      }
       /* NoiseLevels[band+ env*noNoiseBands] */
    }
  }

  
}




/**************************************************************************/
/*!
  \brief     Does the noise floor level estiamtion.
  \return    none

*/
/**************************************************************************/
static void
qmfBasedNoiseFloorDetection(float *noiseLevel,
                            float** quotaMatrixOrig,
                            char* indexVector,
                            int startIndex,
                            int stopIndex,
                            int startChannel,
                            int stopChannel,
                            float ana_max_level,
                            float noiseFloorOffset,
                            int missingHarmonicFlag,
                            float weightFac,
                            INVF_MODE diffThres,
                            INVF_MODE inverseFilteringLevel)
{
  int l,k;
  float meanOrig=0, meanSbr=0, diff;
  float tonalityOrig, tonalitySbr;

  

   /* counting previous operations */

   
  if(missingHarmonicFlag == 1){

     /* pointers for indexVector[l] */
    
    for(l = startChannel; l < stopChannel;l++){

      
      tonalityOrig = 0;

       /* pointers for quotaMatrixOrig[][] */
      
      for(k = startIndex ; k < stopIndex; k++){

        
        tonalityOrig += quotaMatrixOrig[k][l];
      }

       
      tonalityOrig /= (stopIndex-startIndex);

      
      tonalitySbr = 0;

      
      for(k = startIndex ; k < stopIndex; k++){

         
        tonalitySbr += quotaMatrixOrig[k][indexVector[l]];
      }

      /*  already calculated */ 
      tonalitySbr /= (stopIndex-startIndex);


       
      if(tonalityOrig > meanOrig)
      {
        
        meanOrig = tonalityOrig;
      }

       
      if(tonalitySbr > meanSbr)
      {
        
        meanSbr = tonalitySbr;
      }
    }
  }
  else{

     /* pointers for indexVector[l] */
    
    for(l = startChannel; l < stopChannel;l++){

      
      tonalityOrig = 0;

       /* pointers for quotaMatrixOrig[][] */
      
      for(k = startIndex ; k < stopIndex; k++){

        
        tonalityOrig += quotaMatrixOrig[k][l];
      }

       
      tonalityOrig /= (stopIndex-startIndex);

      
      tonalitySbr = 0;

      
      for(k = startIndex ; k < stopIndex; k++){

         
        tonalitySbr += quotaMatrixOrig[k][indexVector[l]];
      }

      /*  already calculated */ 
      tonalitySbr /= (stopIndex-startIndex);

      
      meanOrig += tonalityOrig;
      meanSbr += tonalitySbr;
    }

     
    meanOrig /= (stopChannel - startChannel);
    meanSbr /= (stopChannel - startChannel);
  }


    
  if(meanOrig < (float)0.000976562f && meanSbr < (float)0.000976562f){

    
    meanOrig = 101.5936673f;
    meanSbr  = 101.5936673f;
  }

   
  if(meanOrig < (float)1.0f)
  {
    
    meanOrig = 1.0f;
  }

   
  if(meanSbr < (float)1.0f)
  {
    
    meanSbr  = 1.0f;
  }



   
  if(missingHarmonicFlag == 1){

    
    diff = 1.0f;
  }
  else{

        
    diff = max ((float)1.0f, weightFac*meanSbr/meanOrig);
  }

    
  if(inverseFilteringLevel == INVF_MID_LEVEL ||
     inverseFilteringLevel == INVF_LOW_LEVEL ||
     inverseFilteringLevel == INVF_OFF){

    
    diff = 1.0f;
  }

   
  if (inverseFilteringLevel <= diffThres) {

    
    diff = 1.0f;
  }


   /* value still in register */
  *noiseLevel =  diff/meanOrig;

   /* value still in register */
  *noiseLevel *= noiseFloorOffset;


    
  *noiseLevel = min (*noiseLevel, (float)ana_max_level);

  
}









/**************************************************************************/
/*!
  \brief     Does the noise floor level estiamtion..
  \return    none

*/
/**************************************************************************/
void
sbrNoiseFloorEstimateQmf (HANDLE_SBR_NOISE_FLOOR_ESTIMATE h_sbrNoiseFloorEstimate,
                          const SBR_FRAME_INFO *frame_info,
                          float *noiseLevels,
                          float** quotaMatrixOrig,
                          char* indexVector,
                          int missingHarmonicsFlag,
                          int startIndex,
                          int numberOfEstiamtesPerFrame,
                          int totalNumberOfEstimates,
                          int transientFlag,
                          INVF_MODE* pInvFiltLevels
                          )

{

  int nNoiseEnvelopes, startPos[2], stopPos[2], env, band;

  int noNoiseBands      = h_sbrNoiseFloorEstimate->noNoiseBands;
  int *freqBandTable    = h_sbrNoiseFloorEstimate->freqBandTableQmf;

  

     /* counting previous operations */

   
  nNoiseEnvelopes = frame_info->nNoiseEnvelopes;

   
  if(nNoiseEnvelopes == 1){

    
    startPos[0] = startIndex;

     
    stopPos[0]  = startIndex + 2;
  }
  else{

    
    startPos[0] = startIndex;

     
    stopPos[0]  = startIndex + 1;
    startPos[1] = startIndex + 1;
    stopPos[1]  = startIndex + 2;
  }


   /* pointers for startPos[],
                               stopPos[]
               */
  
  for(env = 0; env < nNoiseEnvelopes; env++){

     /* pointers for h_sbrNoiseFloorEstimate->noiseFloorOffset[],
                                 noiseLevels[],
                                 freqBandTable[],
                                 pInvFiltLevels[]
                 */
    
    for(band = 0; band < noNoiseBands; band++){


         
      qmfBasedNoiseFloorDetection(&noiseLevels[band + env*noNoiseBands],
                                  quotaMatrixOrig,
                                  indexVector,
                                  startPos[env],
                                  stopPos[env],
                                  freqBandTable[band],
                                  freqBandTable[band+1],
                                  h_sbrNoiseFloorEstimate->ana_max_level,
                                  h_sbrNoiseFloorEstimate->noiseFloorOffset[band],
                                  missingHarmonicsFlag,
                                  h_sbrNoiseFloorEstimate->weightFac,
                                  h_sbrNoiseFloorEstimate->diffThres,
                                  pInvFiltLevels[band]);

    }

  }


   
  smoothingOfNoiseLevels(noiseLevels,
                         nNoiseEnvelopes,
                         h_sbrNoiseFloorEstimate->noNoiseBands,
                         h_sbrNoiseFloorEstimate->prevNoiseLevels,
                         h_sbrNoiseFloorEstimate->smoothFilter);




  
  for(env = 0; env < nNoiseEnvelopes; env++){

     /* pointers for noiseLevels[] */
    
    for(band = 0; band < noNoiseBands; band++){

       /* env*noNoiseBands */    
      noiseLevels[band + env*noNoiseBands] = NOISE_FLOOR_OFFSET - (float) (ILOG2*log(noiseLevels[band + env*noNoiseBands]));
    }
  }

  
}



/**************************************************************************/
/*!
  \brief
  \return    errorCode,

*/
/**************************************************************************/
static int
downSampleLoRes(int *v_result,
                int num_result,
                const unsigned char *freqBandTableRef,
                int num_Ref)
{
  int step;
  int i,j;
  int org_length,result_length;
  int v_index[MAX_FREQ_COEFFS/2];

  

  
  org_length=num_Ref;
  result_length=num_result;

  
  v_index[0]=0;	/* Always use left border */

  i=0;

   /* v_index[] */
  
  while(org_length > 0)
    {
      i++;

      
      step=org_length/result_length;

      
      org_length=org_length - step;

      
      result_length--;

       
      v_index[i]=v_index[i-1]+step;
    }

   
  if(i != num_result )
  {
    
    return (1);/* error */
  }

   /* v_result[j]
                  v_index[j]
               */
  
  for(j=0;j<=i;j++)
    {
       
      v_result[j]=freqBandTableRef[v_index[j]];
    }

  

  return (0);
}



/**************************************************************************/
/*!
  \brief    Creates an instance of the noise floor level estimation module.
  \return    errorCode

*/
/**************************************************************************/
int
CreateSbrNoiseFloorEstimate (HANDLE_SBR_NOISE_FLOOR_ESTIMATE  h_sbrNoiseFloorEstimate,
                             int ana_max_level,
                             const unsigned char *freqBandTable,
                             int nSfb,
                             int noiseBands,
                             int noiseFloorOffset,
                             unsigned int useSpeechConfig
                            )
{
  int i;


  

      
  memset(h_sbrNoiseFloorEstimate,0,sizeof(SBR_NOISE_FLOOR_ESTIMATE));



   
  h_sbrNoiseFloorEstimate->smoothFilter = smoothFilter;

  
  if (useSpeechConfig) {

     
    h_sbrNoiseFloorEstimate->weightFac = 1.0;
    h_sbrNoiseFloorEstimate->diffThres = INVF_LOW_LEVEL;
  }
  else {

     
    h_sbrNoiseFloorEstimate->weightFac = 0.25;
    h_sbrNoiseFloorEstimate->diffThres = INVF_MID_LEVEL;
  }

     
  h_sbrNoiseFloorEstimate->ana_max_level          = (float) pow(2,ana_max_level/3.0f);

   
  h_sbrNoiseFloorEstimate->noiseBands             = noiseBands;


   
  if(resetSbrNoiseFloorEstimate(h_sbrNoiseFloorEstimate,freqBandTable,nSfb))
  {
    
    return(1);
  }


   /* h_sbrNoiseFloorEstimate->noiseFloorOffset[] */
   
  for(i=0;i<h_sbrNoiseFloorEstimate->noNoiseBands;i++)
  {
      
    h_sbrNoiseFloorEstimate->noiseFloorOffset[i] = (float) pow(2,noiseFloorOffset/3.0f);
  }

  

  return 0;
}



/**************************************************************************/
/*!
  \brief     Resets the current instance of the noise floor estiamtion
          module.
  \return    errorCode

*/
/**************************************************************************/
int
resetSbrNoiseFloorEstimate (HANDLE_SBR_NOISE_FLOOR_ESTIMATE h_sbrNoiseFloorEstimate,
                            const unsigned char *freqBandTable,
                            int nSfb)
{
  int k2,kx;

  

   
  k2=freqBandTable[nSfb];
  kx=freqBandTable[0];

   
  if(h_sbrNoiseFloorEstimate->noiseBands == 0){
     
    h_sbrNoiseFloorEstimate->noNoiseBands = 1;
  }
  else{
         
    h_sbrNoiseFloorEstimate->noNoiseBands =(int) ( (h_sbrNoiseFloorEstimate->noiseBands * log( (float)k2/kx) * ILOG2)+0.5);

    
    if( h_sbrNoiseFloorEstimate->noNoiseBands==0)
    {
      
      h_sbrNoiseFloorEstimate->noNoiseBands=1;
    }
  }

    /* counting post-operations */

  

  return(downSampleLoRes(h_sbrNoiseFloorEstimate->freqBandTableQmf,
                         h_sbrNoiseFloorEstimate->noNoiseBands,
                         freqBandTable,nSfb));
}





/**************************************************************************/
/*!
  \brief     Deletes the current instancce of the noise floor level
  estimation module.
  \return    none

*/
/**************************************************************************/
void
deleteSbrNoiseFloorEstimate (HANDLE_SBR_NOISE_FLOOR_ESTIMATE h_sbrNoiseFloorEstimate)
{

  

  
  if (h_sbrNoiseFloorEstimate) {
    /*
      nothing to do
    */
  }

  
}
