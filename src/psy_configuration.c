/*
  Psychoaccoustic configuration
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */



typedef struct{
  long  sampleRate;
  const unsigned char *paramLong;
  const unsigned char *paramShort;
}SFB_INFO_TAB;

static const float ABS_LEV = 20.0f;
static const float ABS_LOW = 16887.8f*NORM_PCM_ENERGY; /* maximum peak sine - 96 db*/
static const float BARC_THR_QUIET[] = {15.0f, 10.0f,  7.0f,  2.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f, 3.0f,  5.0f, 10.0f, 20.0f, 30.0f};


static const SFB_INFO_TAB sfbInfoTab[] ={
  {11025, sfb_11025_long_1024, sfb_11025_short_128},
  {12000, sfb_12000_long_1024, sfb_12000_short_128},
  {16000, sfb_16000_long_1024, sfb_16000_short_128},
  {22050, sfb_22050_long_1024, sfb_22050_short_128},
  {24000, sfb_24000_long_1024, sfb_24000_short_128}
};




#define max_barc   24.0f /* maximum barc-value */
#define c_ratio     0.001258925f /* pow(10.0f, -(29.0f/10.0f)) */

#define maskLow                 (30.0f*0.1f) /* in 0.1dB/bark */
#define maskHigh                (15.0f*0.1f) /* in 0.1*dB/bark */
#define maskLowSprEnLong        (30.0f*0.1f) /* in 0.1dB/bark */
#define maskHighSprEnLong       (20.0f*0.1f) /* in 0.1dB/bark */
#define maskHighSprEnLongLowBr  (15.0f*0.1f) /* in 0.1dB/bark */
#define maskLowSprEnShort       (20.0f*0.1f) /* in 0.1dB/bark */
#define maskHighSprEnShort      (15.0f*0.1f) /* in 0.1dB/bark */


static int initSfbTable(long sampleRate,int blockType,int *sfbOffset,int *sfbCnt)
{
  const unsigned char *sfbParam = 0;
  unsigned int i;
  int  specStartOffset,specLines = 0;

  

   /* counting previous operation */

  /*
    select table
  */
   /* sfbInfoTab[] */
  
  for(i = 0; i < sizeof(sfbInfoTab)/sizeof(SFB_INFO_TAB); i++){

     
    if(sfbInfoTab[i].sampleRate == sampleRate){

      
      switch(blockType){
      case LONG_WINDOW:
      case START_WINDOW:
      case STOP_WINDOW:
        
        sfbParam = sfbInfoTab[i].paramLong;
        specLines = FRAME_LEN_LONG;
        break;

      case SHORT_WINDOW:
        
        sfbParam = sfbInfoTab[i].paramShort;
        specLines = FRAME_LEN_SHORT;
        break;
      }
      break;
    }
  }

  
  if(sfbParam==0)
  {
    
    return(1);
  }

  /*
    calc sfb offsets
  */
  
  *sfbCnt=0;
  specStartOffset = 0;
  
  
  do{
     
    sfbOffset[*sfbCnt] = specStartOffset;

     
    specStartOffset += sfbParam[*sfbCnt];

     
    (*sfbCnt)++;
  } while(specStartOffset< specLines);
  
  assert(specStartOffset == specLines);

   
  sfbOffset[*sfbCnt] = specStartOffset;

  

  return 0;
}


/*****************************************************************************

    functionname: atan_approx
    description:  Calculates atan , val >=0
    returns:      approx of atan(val), error is less then 0.5 %
    input:
    output:

*****************************************************************************/
static float atan_approx(float val)
{
  

   
  if(val < (float)1.0)
  {
       /* counting post-operations */
    
    return(val/((float)1.0f+(float)0.280872f*val*val));
  }
  else
  {
       /* counting post-operations */
    
    return((float)1.57079633f-val/((float)0.280872f +val*val));
  }

}


/*****************************************************************************

    functionname: BarcLineValue
    description:  Calculates barc value for one frequency line
    returns:      barc value of line
    input:        number of lines in transform, index of line to check, Fs
    output:

*****************************************************************************/
static float BarcLineValue(int noOfLines, int fftLine, long samplingFreq) {

  float center_freq, temp, bvalFFTLine;

  

  /*
    center frequency of fft line
  */
   
  center_freq = (float) fftLine * ((float)samplingFreq * (float)0.5f)/(float)noOfLines;

   
  temp = (float) atan_approx((float)1.3333333e-4f * center_freq);

    
  bvalFFTLine = (float)13.3f * atan_approx((float)0.00076f * center_freq) + (float)3.5f * temp * temp;

  

  return(bvalFFTLine);

}



static void initThrQuiet(int   numPb,
                         int   *pbOffset,
                         float *pbBarcVal,
                         float *pbThresholdQuiet) {
  int i;
  float barcThrQuiet;

  

   /* pbBarcVal[]
                  pbThresholdQuiet[]
                  pbOffset[]
               */
  
  for(i=0; i<numPb; i++) {

    /* because one pb equals a sfb in psych it's better to use
       the minimum of the bark-wise threshold in quiet */
    int bv1, bv2;

    
    if (i>0)
    {
        
       bv1 = (int)(pbBarcVal[i] + pbBarcVal[i-1])>>1;
    }
    else
    {
       
       bv1 = (int)(pbBarcVal[i])>>1;
    }

     
    if (i<numPb-1)
    {
        
       bv2 = (int)(pbBarcVal[i] + pbBarcVal[i+1])>>1;
    }
    else
    {
       
       bv2 = (int)(pbBarcVal[i]);
    }

      
    bv1=min(bv1,(int)max_barc);
    bv2=min(bv2,(int)max_barc);


       
    barcThrQuiet = min(BARC_THR_QUIET[bv1], BARC_THR_QUIET[bv2]);

       
    pbThresholdQuiet[i] = (float) pow(10.0f,(barcThrQuiet - ABS_LEV) * (float)0.1f)*ABS_LOW*(float)(pbOffset[i+1] - pbOffset[i]);

  }

  
}



static void initSpreading(int numPb,
                          float *pbBarcValue,
                          float *pbMaskLoFactor,
                          float *pbMaskHiFactor,
                          float *pbMaskLoFactorSprEn,
                          float *pbMaskHiFactorSprEn,
                          const long bitrate,
                          const int blockType)
{
  int i;
  float maskLowSprEn, maskHighSprEn;

  

   
  if (blockType != SHORT_WINDOW) {

    
    maskLowSprEn = maskLowSprEnLong;

      
    maskHighSprEn = (bitrate>22000)?maskHighSprEnLong:maskHighSprEnLongLowBr;
  }
  else {

    
    maskLowSprEn = maskLowSprEnShort;
    maskHighSprEn = maskHighSprEnShort;
  }

   /* pbBarcValue[i]
                  pbMaskHiFactor[i]
                  pbMaskLoFactor[i-1]
                  pbMaskHiFactorSprEn[i]
                  pbMaskLoFactorSprEn[i-1]
                  pbMaskLoFactor[numPb-1]
                  pbMaskLoFactorSprEn[numPb-1]
               */
  
  for(i=0; i<numPb; i++) {

    
    if (i > 0) {
       float dbVal;

        
       dbVal = maskHigh * (pbBarcValue[i]-pbBarcValue[i-1]);

         
       pbMaskHiFactor[i] = (float) pow(10.0f, -dbVal);

        
       dbVal = maskLow * (pbBarcValue[i]-pbBarcValue[i-1]);

         
       pbMaskLoFactor[i-1] = (float) pow(10.0f, -dbVal);

        
       dbVal = maskHighSprEn * (pbBarcValue[i]-pbBarcValue[i-1]);

         
       pbMaskHiFactorSprEn[i] = (float) pow(10.0f, -dbVal);

        
       dbVal = maskLowSprEn * (pbBarcValue[i]-pbBarcValue[i-1]);

         
       pbMaskLoFactorSprEn[i-1] = (float) pow(10.0f, -dbVal);
    }
    else {
       
       pbMaskHiFactor[i] = 0.0f;
       pbMaskLoFactor[numPb-1] = 0.0f;

       
       pbMaskHiFactorSprEn[i] = 0.0f;
       pbMaskLoFactorSprEn[numPb-1] = 0.0f;
    }

  }

  
}



static void initBarcValues(int numPb,
                           int *pbOffset,
                           int numLines,
                           long samplingFrequency,
                           float *pbBval)
{
  int   i;
  float pbBval0,pbBval1;

  

  
  pbBval0=0.0;

   /* pbOffset[] */
  
  for(i=0; i<numPb; i++){
    
    pbBval1   = BarcLineValue(numLines, pbOffset[i+1], samplingFrequency);

      
    pbBval[i] = (pbBval0+pbBval1)*(float)0.5f;

    
    pbBval0=pbBval1;
  }

  
}



static void initMinSnr(const long    bitrate,
                       const long    samplerate,
                       const int     numLines,
                       const int     *sfbOffset,
                       const float *pbBarcVal,
                       const int     sfbActive,
                       float       *sfbMinSnr)
{
   int sfb;
   float barcFactor;
   float barcWidth;
   float pePerWindow, pePart;
   float snr;
   float pbVal0,pbVal1;

   

   /* relative number of active barks */
     
   barcFactor = (float)1.0/min(pbBarcVal[sfbActive-1]/max_barc,(float)1.0);

     
   pePerWindow = bits2pe((float)bitrate / (float)samplerate * (float)numLines);

   
   pbVal0=(float)0.0;

    /* pbBarcVal[]
                   sfbOffset[]
                   sfbMinSnr[]
                */
   
   for (sfb = 0; sfb < sfbActive; sfb++) {

       
      pbVal1 = (float)2.0*pbBarcVal[sfb]-pbVal0;

      
      barcWidth=pbVal1-pbVal0;

      
      pbVal0=pbVal1;

      /* allow at least 2.4% of pe for each active barc */
      
      pePart = pePerWindow * (float)0.024f * barcFactor;

      /* adapt to sfb bands */
      
      pePart *= barcWidth;

      /* pe -> snr calculation */
       
      pePart /= (float)(sfbOffset[sfb+1] - sfbOffset[sfb]);

       
      snr = (float) pow(2.0f, pePart) - 1.5f;

        
      snr = 1.0f / max(snr, 1.0f);

      /* upper limit is -1 dB */
        
      snr = min(snr, 0.8f);

      /* lower limit is -25 dB */
        
      snr = max(snr, 0.003f);

      
      sfbMinSnr[sfb] = snr;
   }

   
}


int InitPsyConfiguration(long  bitrate,
                         long  samplerate,
                         int   bandwidth,
                         PSY_CONFIGURATION_LONG *psyConf) 
{
  int sfb;
  float sfbBarcVal[MAX_SFB_LONG];

  

  /*
    init sfb table
  */
    
  if(initSfbTable(samplerate,LONG_WINDOW,psyConf->sfbOffset,&(psyConf->sfbCnt)))
  {
    
    return(1);
  }

  /*
    calculate barc values for each pb
  */
   
  initBarcValues(psyConf->sfbCnt,
                 psyConf->sfbOffset,
                 psyConf->sfbOffset[psyConf->sfbCnt],
                 samplerate,
                 sfbBarcVal);


  /*
    init thresholds in quiet
  */
   
  initThrQuiet(psyConf->sfbCnt,
               psyConf->sfbOffset,
               sfbBarcVal,
               psyConf->sfbThresholdQuiet);


  /*
    calculate spreading function
  */
   
  initSpreading(psyConf->sfbCnt,
                sfbBarcVal,
                psyConf->sfbMaskLowFactor,
                psyConf->sfbMaskHighFactor,
                psyConf->sfbMaskLowFactorSprEn,
                psyConf->sfbMaskHighFactorSprEn,
                bitrate,
                LONG_WINDOW);


  /*
    init ratio
  */
   
  psyConf->ratio = c_ratio;

   
  psyConf->maxAllowedIncreaseFactor = 2.0f;
  psyConf->minRemainingThresholdFactor = 0.01f;

   
  psyConf->clipEnergy = 1.0e9f*NORM_PCM_ENERGY;

     
  psyConf->lowpassLine = (int)((2*bandwidth*FRAME_LEN_LONG)/samplerate);
  
   /* psyConf->sfbOffset[] */
   
  for (sfb = 0; sfb < psyConf->sfbCnt; sfb++){
     
    if (psyConf->sfbOffset[sfb] >= psyConf->lowpassLine)
       break;
   }
   
  psyConf->sfbActive  = sfb;

  /*
    calculate minSnr
  */
   
  initMinSnr(bitrate,
             samplerate,
             psyConf->sfbOffset[psyConf->sfbCnt],
             psyConf->sfbOffset,
             sfbBarcVal,
             psyConf->sfbActive,
             psyConf->sfbMinSnr);

  

  return 0;
}


int InitPsyConfigurationShort(long  bitrate,
                              long  samplerate,
                              int   bandwidth,
                              PSY_CONFIGURATION_SHORT *psyConf) 
{
  int sfb;
  float sfbBarcVal[MAX_SFB_SHORT];

  

  /*
    init sfb table
  */
    
  if(initSfbTable(samplerate,SHORT_WINDOW,psyConf->sfbOffset,&(psyConf->sfbCnt)))
  {
    
    return(1);
  }

  /*
    calculate barc values for each pb
  */
   
  initBarcValues(psyConf->sfbCnt,
                 psyConf->sfbOffset,
                 psyConf->sfbOffset[psyConf->sfbCnt],
                 samplerate,
                 sfbBarcVal);


  /*
    init thresholds in quiet
  */
   
  initThrQuiet(psyConf->sfbCnt,
               psyConf->sfbOffset,
               sfbBarcVal,
               psyConf->sfbThresholdQuiet);


  /*
    calculate spreading function
  */
   
  initSpreading(psyConf->sfbCnt,
                sfbBarcVal,
                psyConf->sfbMaskLowFactor,
                psyConf->sfbMaskHighFactor,
                psyConf->sfbMaskLowFactorSprEn,
                psyConf->sfbMaskHighFactorSprEn,
                bitrate,
                SHORT_WINDOW);


  /*
    init ratio
  */
   
  psyConf->ratio = c_ratio;

   
  psyConf->maxAllowedIncreaseFactor = 2.0f;
  psyConf->minRemainingThresholdFactor = 0.01f;

    
  psyConf->clipEnergy = 1.0e9f*NORM_PCM_ENERGY / (TRANS_FAC * TRANS_FAC);

     
  psyConf->lowpassLine = (int)((2*bandwidth*FRAME_LEN_SHORT)/samplerate);
  
   /* psyConf->sfbOffset[] */
   
  for (sfb = 0; sfb < psyConf->sfbCnt; sfb++){
     
    if (psyConf->sfbOffset[sfb] >= psyConf->lowpassLine)
       break;
   }
   
  psyConf->sfbActive  = sfb;

  /*
    calculate minSnr
  */
   
  initMinSnr(bitrate,
             samplerate,
             psyConf->sfbOffset[psyConf->sfbCnt],
             psyConf->sfbOffset,
             sfbBarcVal,
             psyConf->sfbActive,
             psyConf->sfbMinSnr);

  

  return 0;
}

