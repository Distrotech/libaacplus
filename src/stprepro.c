/*
  stereo pre-processing
*/
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif



/***************************************************************************/
/*!
 
  \brief  create and initialize a handle for stereo preprocessing
 
  \return an error state
 
****************************************************************************/
int InitStereoPreProcessing(HANDLE_STEREO_PREPRO hStPrePro, /*! handle (modified) */
                            int nChannels,    /*! number of channels */ 
                            int bitRate,      /*! the bit rate */
                            int sampleRate,   /*! the sample rate */
                            float usedScfRatio /*! the amount of scalefactors used (0-1.0) */
                            )
{
  float bpf = bitRate*1024.0f/sampleRate;
  float tmp;
  
  

    /* counting previous operation */

      
  memset(hStPrePro,0,sizeof(struct STEREO_PREPRO));
  
   
  if(nChannels == 2) {

     
    (hStPrePro)->stereoAttenuationFlag = 1;


       
    (hStPrePro)->normPeFac    = 230.0f * usedScfRatio / bpf;

         
    (hStPrePro)->ImpactFactor = max (1, 400000.0f / (float) (bitRate - sampleRate*sampleRate/72000.0f) );

       
    (hStPrePro)->stereoAttenuationInc = 22050.0f / sampleRate * 400.0f / bpf;
    (hStPrePro)->stereoAttenuationDec = 22050.0f / sampleRate * 200.0f / bpf;

     
    (hStPrePro)->ConstAtt     = 0.0f;
    (hStPrePro)->stereoAttMax = 12.0f;

    /* energy ratio thresholds (dB) */
     
    (hStPrePro)->SMMin = 0.0f;
    (hStPrePro)->SMMax = 15.0f;
    (hStPrePro)->LRMin = 10.0f;
    (hStPrePro)->LRMax = 30.0f;

    /* pe thresholds */
     
    (hStPrePro)->PeCrit =  1200.0f; 
    (hStPrePro)->PeMin  =  700.0f;
    (hStPrePro)->PeImpactMax = 100.0f;  

    /* init start values */
     
    (hStPrePro)->avrgFreqEnergyL  = 0.0f;
    (hStPrePro)->avrgFreqEnergyR  = 0.0f;
    (hStPrePro)->avrgFreqEnergyS  = 0.0f;
    (hStPrePro)->avrgFreqEnergyM  = 0.0f;
    (hStPrePro)->smoothedPeSumSum = 7000.0f; /* typical start value */
    (hStPrePro)->avgStoM          = -10.0f;  /* typical start value */
    (hStPrePro)->lastLtoR         = 0.0f;
    (hStPrePro)->lastNrgLR        = 0.0f;

     
    tmp = 1.0f - (bpf / 2600.0f);

     
    tmp = max (tmp, 0.0f);
    
      
    (hStPrePro)->stereoAttenuation =  tmp * (hStPrePro)->stereoAttMax;
  }

  

  return 0;
}


/***************************************************************************/
/*!
 
  \brief  do an appropriate attenuation on the side channel of a stereo
          signal
 
  \return nothing
 
****************************************************************************/
void ApplyStereoPreProcess(HANDLE_STEREO_PREPRO hStPrePro, /*!st.-preproc handle */
                           int                 nChannels, /*! total number of channels */              
                           ELEMENT_INFO        *elemInfo,
                           float *timeData,     /*! lr time data (modified) */
                           int granuleLen) /*! num. samples to be processed */
{
  /* inplace operation on inData ! */

  float SMRatio, StoM;
  float LRRatio, LtoR, deltaLtoR, deltaNrg;
  float EnImpact, PeImpact, PeNorm;
  float Att, AttAimed;
  float maxInc, maxDec, swiftfactor;
  float DELTA=0.1f;
  
  float fac = hStPrePro->stereoAttFac;
  float mPart, upper, div;
  float lFac,rFac;

  int i;

  

    /* counting previous operations */

   
  if (!hStPrePro->stereoAttenuationFlag) {
    
    return;
  }

  
  /* calc L/R ratio */
    
  mPart = 2.0f * hStPrePro->avrgFreqEnergyM * (1.0f - fac*fac);

     
  upper = hStPrePro->avrgFreqEnergyL * (1.0f+fac) + hStPrePro->avrgFreqEnergyR * (1.0f-fac) - mPart;
  div   = hStPrePro->avrgFreqEnergyR * (1.0f+fac) + hStPrePro->avrgFreqEnergyL * (1.0f-fac) - mPart;
  
   
  if (div == 0.0f || upper == 0.0f) {

     
    LtoR = hStPrePro->LRMax;
  }
  else {

     
    LRRatio = (float) fabs ( upper / div ) ;

      
    LtoR = (float) fabs(10.0 * log10(LRRatio));
  }
  

  /* calc delta energy to previous frame */
    
  deltaNrg = ( hStPrePro->avrgFreqEnergyL + hStPrePro->avrgFreqEnergyR + 1.0f) / 
             ( hStPrePro->lastNrgLR + 1.0f );

    
  deltaNrg = (float) (fabs(10.0 * log10(deltaNrg)));
  
  

  /* Smooth S/M over time */
    
  SMRatio = (hStPrePro->avrgFreqEnergyS + 1.0f) / (hStPrePro->avrgFreqEnergyM + 1.0f);

   
  StoM = (float) (10.0 * log10(SMRatio));

     
  hStPrePro->avgStoM = DELTA * StoM + (1-DELTA) * hStPrePro->avgStoM;
  
  

  
  EnImpact = 1.0f;
  
    
  if (hStPrePro->avgStoM > hStPrePro->SMMin) {

      
    if (hStPrePro->avgStoM > hStPrePro->SMMax) {

      
      EnImpact = 0.0f;
    }
    else {

       
      EnImpact = (hStPrePro->SMMax - hStPrePro->avgStoM) / (hStPrePro->SMMax - hStPrePro->SMMin);
    }
  }
  
    
  if (LtoR > hStPrePro->LRMin) {

      
    if ( LtoR > hStPrePro->LRMax) {

      
      EnImpact = 0.0f;
    }
    else {

        
      EnImpact *= (hStPrePro->LRMax - LtoR) / (hStPrePro->LRMax - hStPrePro->LRMin);
    }
  }
  
  
  
  PeImpact = 0.0f;
  
   
  PeNorm = hStPrePro->smoothedPeSumSum * hStPrePro->normPeFac;
  
    
  if ( PeNorm > hStPrePro->PeMin )  {

      
    PeImpact= ((PeNorm - hStPrePro->PeMin) / (hStPrePro->PeCrit - hStPrePro->PeMin));
  }
  
    
  if (PeImpact > hStPrePro->PeImpactMax) {

    
    PeImpact = hStPrePro->PeImpactMax;
  }
  
  
   
  AttAimed = EnImpact * PeImpact * hStPrePro->ImpactFactor;
  
  
    
  if (AttAimed > hStPrePro->stereoAttMax) {

    
    AttAimed = hStPrePro->stereoAttMax;
  }
  
  /* only accept changes if they are large enough */

      
  if ( fabs(AttAimed - hStPrePro->stereoAttenuation) < 1.0f && AttAimed != 0.0f) {

    
    AttAimed = hStPrePro->stereoAttenuation;
  }
  
  
  Att = AttAimed;


      /* max() */   
  swiftfactor = (6.0f + hStPrePro->stereoAttenuation) / (10.0f + LtoR) * max(1.0f, 0.2f * deltaNrg );
  
     
  deltaLtoR = max(3.0f, LtoR - hStPrePro->lastLtoR );
  
   
  maxDec = deltaLtoR * deltaLtoR / 9.0f * swiftfactor;
  
    
  maxDec = min( maxDec, 5.0f );
  
   
  maxDec *= hStPrePro->stereoAttenuationDec;
  
  
     
  if (maxDec > hStPrePro->stereoAttenuation * 0.8f) {
    
    
    maxDec = hStPrePro->stereoAttenuation * 0.8f;
  }
  
     
  deltaLtoR = max(3.0f, hStPrePro->lastLtoR - LtoR );
  
   
  maxInc = deltaLtoR * deltaLtoR / 9.0f * swiftfactor;
  
    
  maxInc = min( maxInc, 5.0f );
  
   
  maxInc *= hStPrePro->stereoAttenuationInc;
  
  
     
  if (maxDec > hStPrePro->stereoAttenuation * 0.8f) {
    
    
    maxDec = hStPrePro->stereoAttenuation * 0.8f;
  }
  
     
  deltaLtoR = max(3.0f, hStPrePro->lastLtoR - LtoR );
  
   
  maxInc = deltaLtoR * deltaLtoR / 9.0f * swiftfactor;
  
    
  maxInc = min( maxInc, 5.0f );
  
   
  maxInc *= hStPrePro->stereoAttenuationInc;
  
  
    
  if (Att > hStPrePro->stereoAttenuation + maxInc) {

    
    Att = hStPrePro->stereoAttenuation + maxInc;
  }

    
  if (Att < hStPrePro->stereoAttenuation - maxDec) {

    
    Att = hStPrePro->stereoAttenuation - maxDec;
  }
  
    
  if (hStPrePro->ConstAtt == 0) hStPrePro->stereoAttenuation = Att;
  else                          hStPrePro->stereoAttenuation = hStPrePro->ConstAtt;
  

  /* perform attenuation of Side Channel */
  
     
  hStPrePro->stereoAttFac = (float) pow(10.0f, 0.05f*(-hStPrePro->stereoAttenuation ));
  
    
  lFac = 0.5f * (1.0f + hStPrePro->stereoAttFac);
  rFac = 0.5f * (1.0f - hStPrePro->stereoAttFac);
  
   /* pointer for timeData[nChannels * i + elemInfo->ChannelIndex[0]], 
                              timeData[nChannels * i + elemInfo->ChannelIndex[1]]
               */
  
  for (i = 0; i < granuleLen; i++){
    float L = lFac * timeData[nChannels * i+elemInfo->ChannelIndex[0]] + rFac * timeData[nChannels * i + elemInfo->ChannelIndex[1]];
    float R = rFac * timeData[nChannels * i+elemInfo->ChannelIndex[0]] + lFac * timeData[nChannels * i + elemInfo->ChannelIndex[1]];

      /* counting operations above */
    
    
    timeData[nChannels * i + elemInfo->ChannelIndex[0]]= L;
    timeData[nChannels * i + elemInfo->ChannelIndex[1]]= R;
  }
  
   
  hStPrePro->lastLtoR = LtoR;

    
  hStPrePro->lastNrgLR = hStPrePro->avrgFreqEnergyL + hStPrePro->avrgFreqEnergyR;
  
  
}


/***************************************************************************/
/*!
 
  \brief  calc attenuation parameters - this has to take place after
          the 'QCMain' call because PeSum of the merged l/r-m/s signal
          is relevant
 
  \return nothing
 
****************************************************************************/
void UpdateStereoPreProcess(PSY_OUT_CHANNEL psyOutChannel[MAX_CHANNELS],
                            QC_OUT_ELEMENT* qcOutElement,   /*! provides access to PE */
                            HANDLE_STEREO_PREPRO hStPrePro,
                            float weightPeFac               /*! ratio of ms PE vs. lr PE */
                            )
{
  

   
  if (hStPrePro->stereoAttenuationFlag)
  { 
    float DELTA = 0.1f;

     
    hStPrePro->avrgFreqEnergyL = psyOutChannel[0].sfbEnSumLR;
    hStPrePro->avrgFreqEnergyR = psyOutChannel[1].sfbEnSumLR;
    hStPrePro->avrgFreqEnergyM = psyOutChannel[0].sfbEnSumMS;
    hStPrePro->avrgFreqEnergyS = psyOutChannel[1].sfbEnSumMS;


       
    hStPrePro->smoothedPeSumSum = 
      DELTA * qcOutElement->pe * weightPeFac + (1 - DELTA) * hStPrePro->smoothedPeSumSum;


  }

  
}
