/*
  Temporal Noise Shaping
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define RATIO_MULT        0.25f

#define TNS_MODIFY_BEGIN           2600 /* Hz */
#define RATIO_PATCH_LOWER_BORDER   380  /* Hz */
#define TNS_PARCOR_THRESH          0.1f


static void CalcGaussWindow(float *win, const int winSize, const int samplingRate,
                            const int blockType, const float timeResolution );

static void CalcWeightedSpectrum(const float spectrum[],
                                 float weightedSpectrum[],
                                 float* sfbEnergy,
                                 const int* sfbOffset, int lpcStartLine,
                                 int lpcStopLine, int lpcStartBand,int lpcStopBand);

static void AutoCorrelation(const float input[], float corr[],
                            int samples, int corrCoeff);

static float AutoToParcor(const float input[],
                          float reflCoeff[], int numOfCoeff, float workBuffer[]);


static float CalcTnsFilter(const float* signal, const float window[], int numOfLines,
                           int tnsOrder, float parcor[]);

static void Parcor2Index(const float parcor[], int index[], int order,
                         int bitsPerCoeff);

static void Index2Parcor(const int index[], float parcor[], int order,
                         int bitsPerCoeff);



static void AnalysisFilterLattice(const float signal[], int numOfLines,
                                  const float parCoeff[], int order,
                                  float output[]);


/*****************************************************************************

    functionname: FreqToBandWithRounding
    description:  Returns index of nearest band border
    returns:
    input:        frequency, sampling frequency, total number of bands,
                  table of band borders
    output:
    globals:

*****************************************************************************/
static int FreqToBandWithRounding(int freq, int fs, int numOfBands,
                                  const int *bandStartOffset)
{
  int lineNumber, band;

  

  /*  assert(freq >= 0);  */

     
  lineNumber = (freq*bandStartOffset[numOfBands]*4/fs+1)/2;

  /* freq > fs/2 */
   
  if (lineNumber >= bandStartOffset[numOfBands])
  {
    
    return numOfBands;
  }

  /* find band the line number lies in */
   /* bandStartOffset[] */
  
  for (band=0; band<numOfBands; band++) {
     
    if (bandStartOffset[band+1]>lineNumber) break;
  }

  /* round to nearest band border */
   
  if (lineNumber - bandStartOffset[band] >
      bandStartOffset[band+1] - lineNumber )
    {
      
      band++;
    }

  

  return band;
};


/*****************************************************************************

    functionname: InitTnsConfiguration
    description:  fill TNS_CONFIG structure with sensible content
    returns:
    input:        bitrate, samplerate, number of channels,
                  TNS Config struct (modified),
                  psy config struct,
                  tns active flag
    output:

*****************************************************************************/
int InitTnsConfiguration(   int bitRate,
                            long sampleRate,
                            int channels,
                            TNS_CONFIG *tC,
                            PSY_CONFIGURATION_LONG pc,
                            int active) 
{
  

   
  tC->maxOrder     = TNS_MAX_ORDER;
  tC->tnsStartFreq = 1275;
  tC->coefRes      = 4;
  
       
  if ( GetTnsParam(&tC->confTab, bitRate/channels, channels, LONG_WINDOW) )
    return 1;


    
  CalcGaussWindow(tC->acfWindow, 
                  tC->maxOrder+1,
                  sampleRate, 
                  LONG_WINDOW,
                  tC->confTab.tnsTimeResolution);

    
  GetTnsMaxBands(sampleRate, LONG_WINDOW, &tC->tnsMaxSfb);

   
  tC->tnsActive = 1;

  
  if (active==0)  {
    
    tC->tnsActive=0;
  }

  /*now calc band and line borders */

     
  tC->tnsStopBand=min(pc.sfbCnt, tC->tnsMaxSfb);

   
  tC->tnsStopLine = pc.sfbOffset[tC->tnsStopBand];

    
  tC->tnsStartBand=FreqToBandWithRounding(tC->tnsStartFreq, sampleRate,
                                          pc.sfbCnt, pc.sfbOffset);

    
  tC->tnsModifyBeginCb = FreqToBandWithRounding(TNS_MODIFY_BEGIN,
                                                sampleRate,
                                                pc.sfbCnt,
                                                pc.sfbOffset);

    
  tC->tnsRatioPatchLowestCb = FreqToBandWithRounding(RATIO_PATCH_LOWER_BORDER,
                                                     sampleRate,
                                                     pc.sfbCnt,
                                                     pc.sfbOffset);


   
  tC->tnsStartLine = pc.sfbOffset[tC->tnsStartBand];

    
  tC->lpcStopBand=FreqToBandWithRounding(tC->confTab.lpcStopFreq, sampleRate,
                                         pc.sfbCnt, pc.sfbOffset);

     
  tC->lpcStopBand=min(tC->lpcStopBand, pc.sfbActive);

   
  tC->lpcStopLine = pc.sfbOffset[tC->lpcStopBand];

    
  tC->lpcStartBand=FreqToBandWithRounding(tC->confTab.lpcStartFreq, sampleRate,
                                           pc.sfbCnt, pc.sfbOffset);

   
  tC->lpcStartLine = pc.sfbOffset[tC->lpcStartBand];

   
  tC->threshold =tC->confTab.threshOn;
  
  

  return 0;
}

/*****************************************************************************

    functionname: InitTnsConfigurationShort
    description:  fill TNS_CONFIG structure with sensible content
    returns:
    input:        bitrate, samplerate, number of channels,
                  TNS Config struct (modified),
                  psy config struct,
                  tns active flag
    output:

*****************************************************************************/
int InitTnsConfigurationShort(   int bitRate,
                                 long sampleRate,
                                 int channels,
                                 TNS_CONFIG *tC,
                                 PSY_CONFIGURATION_SHORT pc,
                                 int active) 
{
  

   
  tC->maxOrder     = TNS_MAX_ORDER_SHORT;
  tC->tnsStartFreq = 2750;
  tC->coefRes      = 3;
  
       
  if ( GetTnsParam(&tC->confTab, bitRate/channels, channels, SHORT_WINDOW) )
    return 1;


    
  CalcGaussWindow(tC->acfWindow, 
                  tC->maxOrder+1,
                  sampleRate, 
                  SHORT_WINDOW,
                  tC->confTab.tnsTimeResolution);

    
  GetTnsMaxBands(sampleRate, SHORT_WINDOW, &tC->tnsMaxSfb);

   
  tC->tnsActive = 1;

  
  if (active==0)  {
    
    tC->tnsActive=0;
  }

  /*now calc band and line borders */

     
  tC->tnsStopBand=min(pc.sfbCnt, tC->tnsMaxSfb);

   
  tC->tnsStopLine = pc.sfbOffset[tC->tnsStopBand];

    
  tC->tnsStartBand=FreqToBandWithRounding(tC->tnsStartFreq, sampleRate,
                                          pc.sfbCnt, pc.sfbOffset);

    
  tC->tnsModifyBeginCb = FreqToBandWithRounding(TNS_MODIFY_BEGIN,
                                                sampleRate,
                                                pc.sfbCnt,
                                                pc.sfbOffset);

    
  tC->tnsRatioPatchLowestCb = FreqToBandWithRounding(RATIO_PATCH_LOWER_BORDER,
                                                     sampleRate,
                                                     pc.sfbCnt,
                                                     pc.sfbOffset);


   
  tC->tnsStartLine = pc.sfbOffset[tC->tnsStartBand];

    
  tC->lpcStopBand=FreqToBandWithRounding(tC->confTab.lpcStopFreq, sampleRate,
                                         pc.sfbCnt, pc.sfbOffset);

     
  tC->lpcStopBand=min(tC->lpcStopBand, pc.sfbActive);

   
  tC->lpcStopLine = pc.sfbOffset[tC->lpcStopBand];

    
  tC->lpcStartBand=FreqToBandWithRounding(tC->confTab.lpcStartFreq, sampleRate,
                                          pc.sfbCnt, pc.sfbOffset);

   
  tC->lpcStartLine = pc.sfbOffset[tC->lpcStartBand];

   
  tC->threshold =tC->confTab.threshOn;
  
  

  return 0;
};


/*****************************************************************************
    functionname: TnsDetect
    description:  do decision, if TNS shall be used or not
    returns:
    input:        tns data structure (modified),
                  tns config structure,
                  pointer to scratch space,
                  scalefactor size and table,
                  spectrum,
                  subblock num, blocktype,
                  sfb-wise energy.

*****************************************************************************/
int TnsDetect(TNS_DATA* tnsData,
              TNS_CONFIG tC,
              float* pScratchTns,
              const int sfbOffset[],
              float* spectrum,
              int subBlockNumber,
              int blockType,
              float * sfbEnergy)
{
  float predictionGain;
  float* pWeightedSpectrum = pScratchTns + subBlockNumber*FRAME_LEN_SHORT;

  

    /* counting previous operation */

  
  
  if (tC.tnsActive) {

    
    CalcWeightedSpectrum(spectrum,
                         pWeightedSpectrum,
                         sfbEnergy,
                         sfbOffset,
                         tC.lpcStartLine,
                         tC.lpcStopLine,
                         tC.lpcStartBand,
                         tC.lpcStopBand);

     
    if (blockType!=SHORT_WINDOW) {

           
        predictionGain = CalcTnsFilter( &pWeightedSpectrum[tC.lpcStartLine],
                                        tC.acfWindow,
                                        tC.lpcStopLine-tC.lpcStartLine,
                                        tC.maxOrder,
                                        tnsData->dataRaw.Long.subBlockInfo.parcor);

           
        tnsData->dataRaw.Long.subBlockInfo.tnsActive =
            (predictionGain > tC.threshold)?1:0;

         
        tnsData->dataRaw.Long.subBlockInfo.predictionGain = predictionGain;
    }
    else{

           
        predictionGain = CalcTnsFilter( &pWeightedSpectrum[tC.lpcStartLine],
                                        tC.acfWindow,
                                        tC.lpcStopLine-tC.lpcStartLine,
                                        tC.maxOrder,
                                        tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor);

           
        tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].tnsActive =
            (predictionGain > tC.threshold)?1:0;

         
        tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].predictionGain = predictionGain;
    }

  }
  else{

     
    if (blockType!=SHORT_WINDOW){

         
        tnsData->dataRaw.Long.subBlockInfo.tnsActive = 0;
        tnsData->dataRaw.Long.subBlockInfo.predictionGain = 0.0f;
    }
    else {

         
        tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].tnsActive = 0;
        tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].predictionGain = 0.0f;
    }
  }

  

  return 0;
}


/*****************************************************************************
    functionname: TnsSync
    description:

*****************************************************************************/
void TnsSync(TNS_DATA *tnsDataDest,
             const TNS_DATA *tnsDataSrc,
             const TNS_CONFIG tC,
             const int subBlockNumber,
             const int blockType)
{
   TNS_SUBBLOCK_INFO *sbInfoDest;
   const TNS_SUBBLOCK_INFO *sbInfoSrc;
   int i;

   

    
   if (blockType!=SHORT_WINDOW) {

       
      sbInfoDest = &tnsDataDest->dataRaw.Long.subBlockInfo;
      sbInfoSrc  = &tnsDataSrc->dataRaw.Long.subBlockInfo;
   }
   else {

       
      sbInfoDest = &tnsDataDest->dataRaw.Short.subBlockInfo[subBlockNumber];
      sbInfoSrc  = &tnsDataSrc->dataRaw.Short.subBlockInfo[subBlockNumber];
   }

       
   if (fabs(sbInfoDest->predictionGain - sbInfoSrc->predictionGain) <
       ((float)0.03f * sbInfoDest->predictionGain)) {

       
      sbInfoDest->tnsActive = sbInfoSrc->tnsActive;

       /* sbInfoDest->parcor[]
                      sbInfoSrc->parcor[]
                   */
      
      for ( i=0; i< tC.maxOrder; i++) {

        
        sbInfoDest->parcor[i] = sbInfoSrc->parcor[i];
      }
   }

   
}


/*****************************************************************************
    functionname: TnsEncodeData
    description:

*****************************************************************************/
int TnsEncodeData(TNS_INFO* tnsInfo,
              TNS_DATA* tnsData,
              int numOfSfb,
              TNS_CONFIG tC,
              int lowPassLine,
              float* spectrum,
              int subBlockNumber,
              int blockType)
{
  int i;

  

   
  if (blockType!=SHORT_WINDOW) {

     /* tnsInfo->numOfFilters[subBlockNumber]
                    tnsInfo->coef[subBlockNumber*TRANS_FAC]
                    tnsInfo->order[subBlockNumber]
                    tnsInfo->tnsActive[subBlockNumber]
                    tnsInfo->coefRes[subBlockNumber]
                    tnsInfo->length[subBlockNumber]
                 */

     
    if (tnsData->dataRaw.Long.subBlockInfo.tnsActive == 0) {

       
      tnsInfo->tnsActive[subBlockNumber] = 0;

      
      return 0;
    }
    else {

       
      Parcor2Index(tnsData->dataRaw.Long.subBlockInfo.parcor,
                   &tnsInfo->coef[0],
                   tC.maxOrder,
                   tC.coefRes);

       
      Index2Parcor(&tnsInfo->coef[0],
                   tnsData->dataRaw.Long.subBlockInfo.parcor,
                   tC.maxOrder,
                   tC.coefRes);

       /* tnsData->dataRaw.Long.subBlockInfo.parcor[i] */
      
      for (i=tC.maxOrder-1; i>=0; i--)  {

          
        if (tnsData->dataRaw.Long.subBlockInfo.parcor[i]>TNS_PARCOR_THRESH  ||
            tnsData->dataRaw.Long.subBlockInfo.parcor[i]<-TNS_PARCOR_THRESH)
          break;
      }
        
      tnsInfo->order[subBlockNumber]=i+1;

       
      tnsInfo->tnsActive[subBlockNumber]=1;

       /* tnsInfo->tnsActive[] */
      
      for (i=subBlockNumber+1; i<TRANS_FAC; i++) {

        
        tnsInfo->tnsActive[i]=0;
      }

       
      tnsInfo->coefRes[subBlockNumber]=tC.coefRes;

        
      tnsInfo->length[subBlockNumber]= numOfSfb-tC.tnsStartBand;


           
      AnalysisFilterLattice(&(spectrum[tC.tnsStartLine]),
                            min(tC.tnsStopLine,lowPassLine)-tC.tnsStartLine,
                            tnsData->dataRaw.Long.subBlockInfo.parcor,
                            tnsInfo->order[subBlockNumber],
                            &(spectrum[tC.tnsStartLine]));

    }
  }     /* if (blockType!=SHORT_WINDOW) */
  else /*short block*/ {

     /* tnsData->dataRaw.Short.subBlockInfo[subBlockNumber]
                    tnsInfo->tnsActive[subBlockNumber]
                    tnsInfo->coef[subBlockNumber*TRANS_FAC]
                    tnsInfo->order[subBlockNumber]
                    tnsInfo->coefRes[subBlockNumber]
                    tnsInfo->length[subBlockNumber]
                 */

     
    if (tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].tnsActive == 0) {

       
      tnsInfo->tnsActive[subBlockNumber] = 0;

      
      return(0);
    }
    else {

       
      Parcor2Index(tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor,
                   &tnsInfo->coef[subBlockNumber*TNS_MAX_ORDER_SHORT],
                   tC.maxOrder,
                   tC.coefRes);

       
      Index2Parcor(&tnsInfo->coef[subBlockNumber*TNS_MAX_ORDER_SHORT],
                   tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor,
                   tC.maxOrder,
                   tC.coefRes);

       /* tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor[] */
      
      for (i=tC.maxOrder-1; i>=0; i--)  {

           
        if (tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor[i]>TNS_PARCOR_THRESH  ||
            tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor[i]<-TNS_PARCOR_THRESH)
          break;
      }

        
      tnsInfo->order[subBlockNumber]=i+1;

       
      tnsInfo->tnsActive[subBlockNumber]=1;

       
      tnsInfo->coefRes[subBlockNumber]=tC.coefRes;

        
      tnsInfo->length[subBlockNumber]= numOfSfb-tC.tnsStartBand;


         
      AnalysisFilterLattice(&(spectrum[tC.tnsStartLine]), tC.tnsStopLine-tC.tnsStartLine,
                 tnsData->dataRaw.Short.subBlockInfo[subBlockNumber].parcor,
                 tnsInfo->order[subBlockNumber],
                 &(spectrum[tC.tnsStartLine]));
      
    }
  }

  

  return 0;
}







/*****************************************************************************

    functionname: CalcGaussWindow
    edscription:  calculates Gauss window for acf windowing depending on desired
                  temporal resolution, transform size and sampling rate
    returns:      -
    input:        window size, fs, bitRate, no. of transform lines, time res.
    output:       window coefficients (right half)

*****************************************************************************/
static void CalcGaussWindow(float  *win,
                            const int winSize,
                            const int samplingRate,
                            const int blockType,
                            const float timeResolution )
{
  int     i;
  float gaussExp = 3.14159265358979323f * samplingRate * 0.001f * (float)timeResolution / (blockType != SHORT_WINDOW ? 1024.0f:128.0f);
   
  

       /* .. != .. ? */  

    
    gaussExp = -0.5f * gaussExp * gaussExp;

     /* win[] */
    
    for(i=0; i<winSize; i++) {

         
      win[i] = (float) exp( gaussExp * (i+0.5) * (i+0.5) );
    }

  
}




/*****************************************************************************

    functionname: CalcWeightedSpectrum
    description:  calculate weighted spectrum for LPC calculation
    returns:      -
    input:        input spectrum, ptr. to weighted spectrum, no. of lines,
                  sfb energies
    output:       weighted spectrum coefficients

*****************************************************************************/
static void CalcWeightedSpectrum(const float    spectrum[],
                                 float          weightedSpectrum[],
                                 float        *sfbEnergy,
                                 const int     *sfbOffset,
                                 int            lpcStartLine,
                                 int            lpcStopLine,
                                 int            lpcStartBand,
                                 int            lpcStopBand)
{
    int     i, sfb;
    float   tmp;
    float   tnsSfbMean[MAX_SFB];    /* length [lpcStopBand-lpcStartBand] should be sufficient here */

    

    /* calc 1/sqrt(en) */
     /* tnsSfbMean[]
                    sfbEnergy[]
                 */
    
    for( sfb = lpcStartBand; sfb < lpcStopBand; sfb++){

         
      tnsSfbMean[sfb] = (float) ( 1.0 / sqrt(sfbEnergy[sfb] + 1e-30f) );
    }

    /* spread normalized values from sfbs to lines */
    
    sfb = lpcStartBand;

     /* tnsSfbMean[] */
    
    tmp = tnsSfbMean[sfb];

     /* sfbOffset[]
                    weightedSpectrum[]
                 */
    
    for ( i=lpcStartLine; i<lpcStopLine; i++){

         
        if (sfbOffset[sfb+1]==i){

            
            sfb++;

             
            if (sfb+1 < lpcStopBand){

                
                tmp = tnsSfbMean[sfb];
            }
        }

        
        weightedSpectrum[i] = tmp;
    }

    /*filter down*/
     /* weightedSpectrum[] */
     
    for (i=lpcStopLine-2; i>=lpcStartLine; i--){
        
      weightedSpectrum[i] = (weightedSpectrum[i] + weightedSpectrum[i+1]) * 0.5f;
    }

    /* filter up */
     /* weightedSpectrum[] */
     
    for (i=lpcStartLine+1; i<lpcStopLine; i++){

        
      weightedSpectrum[i] = (weightedSpectrum[i] + weightedSpectrum[i-1]) * 0.5f;
    }

    /* weight and normalize */
     /* weightedSpectrum[] */
    
    for (i=lpcStartLine; i<lpcStopLine; i++){

       
      weightedSpectrum[i] = weightedSpectrum[i] * spectrum[i];
    }

    
}



/*****************************************************************************

    functionname: CalcTnsFilter
    description:  LPC calculation for one TNS filter
    returns:      prediction gain
    input:        signal spectrum, acf window, no. of spectral lines,
                  max. TNS order, ptr. to reflection ocefficients
    output:       reflection coefficients

*****************************************************************************/
static float CalcTnsFilter(const float *signal,
                            const float window[],
                            int numOfLines,
                            int tnsOrder,
                            float parcor[])
{
    float   autoCorrelation[TNS_MAX_ORDER+1];
    float   parcorWorkBuffer[2*TNS_MAX_ORDER];
    float  predictionGain;
    int     i;

    

    assert(tnsOrder <= TNS_MAX_ORDER);

     
    AutoCorrelation(signal, autoCorrelation, numOfLines, tnsOrder+1);

    
    if(window) {

         /* autoCorrelation[]
                        window[]
                     */
        
        for(i=0; i<tnsOrder+1; i++) {

           
          autoCorrelation[i] = autoCorrelation[i] * window[i];
        }
    }

    
    predictionGain = AutoToParcor(autoCorrelation, parcor, tnsOrder, parcorWorkBuffer);

    

    return(predictionGain);
}

/*****************************************************************************

    functionname: AutoCorrelation
    description:  calc. autocorrelation (acf)
    returns:      -
    input:        input values, no. of input values, no. of acf values
    output:       acf values

*****************************************************************************/
static void AutoCorrelation(const float input[],
                            float       corr[],
                            int         samples,
                            int         corrCoeff) {
    int         i, j;
    float       accu;
    int         scf;

    
  
    /* 
      get next power of 2 of samples 
    */
    
    for(scf=0;(1<<scf) < samples;scf++);
   
    
    accu = 0.0;

     /* pointer for input[] */
    
    for(j=0; j<samples; j++) {

      
      accu += input[j] * input[j];
    }

    
    corr[0] = accu;

    
    samples--;

     /* pointer for corr[] */
    
    for(i=1; i<corrCoeff; i++) {

        
        accu = 0.0;

         /* pointer for input[], input[j+i] */
        
        for(j=0; j<samples; j++) {

          
          accu += input[j] * input[j+i];
        }

        
        corr[i] = accu;

        
        samples--;
    }

    
}    


/*****************************************************************************

    functionname: AutoToParcor
    description:  conversion autocorrelation to reflection coefficients
    returns:      prediction gain
    input:        <order+1> input values, no. of output values (=order),
                  ptr. to workbuffer (required size: 2*order)
    output:       <order> reflection coefficients

*****************************************************************************/
static float AutoToParcor(const float input[], float reflCoeff[], int numOfCoeff,
                          float workBuffer[]) {
  int i, j;
  float  *pWorkBuffer; /* temp pointer */
  float predictionGain = 1.0f;

  

   /* counting previous operation */

   /* reflCoeff[] */
  
  for(i=0;i<numOfCoeff;i++)
  {
    
    reflCoeff[i] = 0;
  }

  
  if(input[0] == 0.0)
  {
    
    return(predictionGain);
  }

   /* workBuffer[]
                  workBuffer[i+numOfCoeff]
                  input[]
               */
  
  for(i=0; i<numOfCoeff; i++) {

    
    workBuffer[i] = input[i];
    workBuffer[i+numOfCoeff] = input[i+1];
  }

   /* workBuffer[i+numOfCoeff]
                  reflCoeff[i]
               */
  
  for(i=0; i<numOfCoeff; i++) {
    float refc, tmp;

    
    tmp = workBuffer[numOfCoeff + i];

    
    if(tmp < 0.0)
    {
      
      tmp = -tmp;
    }

     
    if(workBuffer[0] < tmp)
      break;

    
    if(workBuffer[0] == 0.0f) {

      
      refc = 0.0f;
    } else {

      
      refc = tmp / workBuffer[0];
    }

    
    if(workBuffer[numOfCoeff + i] > 0.0)
    {
      
      refc = -refc;
    }

    
    reflCoeff[i] = refc;

    
    pWorkBuffer = &(workBuffer[numOfCoeff]);

     /* workBuffer[j-i] */
    
    for(j=i; j<numOfCoeff; j++) {
      float accu1, accu2;

       
      accu1 = pWorkBuffer[j]  + refc * workBuffer[j-i];
      accu2 = workBuffer[j-i] + refc * pWorkBuffer[j];
      pWorkBuffer[j] = accu1;
      workBuffer[j-i] = accu2;
    }
  }

   
  predictionGain = (input[0] + 1e-30f) / (workBuffer[0] + 1e-30f);

  

  return(predictionGain);
}



static int Search3(float parcor)
{
  int index=0;
  int i;

  

   /* counting previous operation */

   /* tnsCoeff3Borders[] */
  
  for(i=0;i<8;i++){

     
    if(parcor > tnsCoeff3Borders[i])
    {
      
      index=i;
    }
  }

   /* counting post-operation */
  

  return(index-4);
}

static int Search4(float parcor)
{
  int index=0;
  int i;

  

   /* counting previous operation */

   /* tnsCoeff4Borders[] */
  
  for(i=0;i<16;i++){

     
    if(parcor > tnsCoeff4Borders[i])
    {
      
      index=i;
    }
  }

   /* counting post-operation */
  

  return(index-8);
}



/*****************************************************************************

    functionname: Parcor2Index

*****************************************************************************/
static void Parcor2Index(const float parcor[], int index[], int order,
                         int bitsPerCoeff) {
  int i;

  

   /* index[]
                  parcor[]
               */
  
  for(i=0; i<order; i++) {

     
    if(bitsPerCoeff == 3)
    {
       
      index[i] = Search3(parcor[i]);
    }
    else
    {
       
      index[i] = Search4(parcor[i]);
    }
  }

  
}


/*****************************************************************************

    functionname: Index2Parcor
    description:  inverse quantization for reflection coefficients
    returns:      -
    input:        quantized values, ptr. to reflection coefficients,
                  no. of coefficients, resolution
    output:       reflection coefficients

*****************************************************************************/
static void Index2Parcor(const int index[], float parcor[], int order,
                         int bitsPerCoeff) {
  int i;

  

   /* parcor[]
                  index[]
               */
  
  for(i=0; i<order; i++) {

      /* .. == .. ? */  
    parcor[i] = bitsPerCoeff == 4 ? tnsCoeff4[index[i]+8] : tnsCoeff3[index[i]+4];
  }

  
}



/*****************************************************************************

    functionname: FIRLattice

*****************************************************************************/
static float FIRLattice(int order, 
                        float x,
                        float *state_par,
                        const float *coef_par)
{
   int i;
   float accu, tmp, tmpSave;

   

   
   tmpSave = x;

    /* coef_par[]
                   state_par[]
                */
   
   for(i=0; i<order-1; i++) {

     
     tmp      = coef_par[i] * x;

     
     tmp     += state_par[i];

     
     accu     = coef_par[i] * state_par[i];

     
     accu    += x;

     
     x        = accu;
     state_par[i] = tmpSave;
     tmpSave  = tmp;
  }

  /* last stage: only need half operations */

  
  accu  = state_par[order-1] * coef_par[order-1];

  
  accu += x;

  
  x     = accu;
  state_par[order-1] = tmpSave;

  

  return x;
}



/*****************************************************************************

    functionname: AnalysisFilterLattice
    description:  TNS analysis filter
    returns:      -
    input:        input signal values, no. of lines, PARC coefficients,
                  filter order, ptr. to output values, filtering direction
    output:       filtered signal values

*****************************************************************************/
static void AnalysisFilterLattice(const float signal[], int numOfLines,
                                  const float parCoeff[], int order,
                                  float output[])
{

  float state_par[TNS_MAX_ORDER] = {0.0f};
  int j;

  

   /* counting previous operation */

   /* output[]
                  signal[]
               */
  
  for(j=0; j<numOfLines; j++) {

     
    output[j] = FIRLattice(order,signal[j],state_par,parCoeff);
  }

  
}


/*****************************************************************************

    functionname: ApplyTnsMultTableToRatios
    description:  change thresholds according to tns
    returns:      -
    input:        modifier table,
                  TNS subblock info,
                  thresholds (modified),
                  number of bands

*****************************************************************************/
void ApplyTnsMultTableToRatios(int startCb,
                               int stopCb,
                               float *thresholds) 
{
  int i;

  

   /* thresholds[]
               */
  
  for(i=startCb; i<stopCb; i++) {
    
     
    thresholds[i] = thresholds[i]*RATIO_MULT;
  }
  
  
}
