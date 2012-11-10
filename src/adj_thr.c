/*
  adjust thresholds
*/
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define redExp       0.25f
#define invRedExp    (1.0f/redExp)
#define minSnrLimit  0.8f /* 1dB */

/* values for avoid hole flag */
enum _avoid_hole_state {
  NO_AH              =0,
  AH_INACTIVE        =1,
  AH_ACTIVE          =2
};


typedef struct {
   PE_CHANNEL_DATA peChannelData[MAX_CHANNELS];
   float pe;
   float constPart;
   float nActiveLines;
   float offset;
} PE_DATA;


/* convert from bits to pe */
float bits2pe(const float bits) {

   
    /* counting post-operation */
   

   return (bits * 1.18f);
}

/* loudness calculation (threshold to the power of redExp) */
static void calcThreshExp(float thrExp[MAX_CHANNELS][MAX_GROUPED_SFB],
                          PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                          const int nChannels)
{
   int ch, sfb,sfbGrp;

   

   
   for (ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];

     /* counting operation above */

     
    for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup) {

       /* pointer for thrExp[][],
                                  psyOutChan->sfbThreshold[]
                   */
       
      for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

         
        thrExp[ch][sfbGrp+sfb] = (float) pow(psyOutChan->sfbThreshold[sfbGrp+sfb], redExp);
      }
    }
   }

   
}

/* reduce minSnr requirements for bands with relative low energies */
static void adaptMinSnr(PSY_OUT_CHANNEL     psyOutChannel[MAX_CHANNELS],
                        MINSNR_ADAPT_PARAM *msaParam,
                        const int           nChannels)
{
  int ch, sfb, sfbOffs, nSfb;
  float avgEn, dbRatio, minSnrRed;

  

  
  for (ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL* psyOutChan = &psyOutChannel[ch];

     /* counting operation above */

    /* calc average energy per scalefactor band */
    
    avgEn = 0.0f;
    nSfb = 0;

     
    for (sfbOffs=0; 
         sfbOffs<psyOutChan->sfbCnt;
         sfbOffs+=psyOutChan->sfbPerGroup) {

       /* pointer for psyOutChan->sfbEnergy[] */
       
      for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

        
        avgEn += psyOutChan->sfbEnergy[sfbOffs+sfb];

        
        nSfb++;
      }
    }

    
    if (nSfb > 0) {

      
      avgEn /= nSfb;
    }

    /* reduce minSnr requirement by minSnr^minSnrRed dependent on avgEn/sfbEn */
     
    for (sfbOffs=0; 
         sfbOffs<psyOutChan->sfbCnt;
         sfbOffs+=psyOutChan->sfbPerGroup) {

       /* pointer for psyOutChan->sfbEnergy[],
                                  psyOutChan->sfbMinSnr[]
                   */
       
      for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

           
        if (msaParam->startRatio*psyOutChan->sfbEnergy[sfbOffs+sfb] < avgEn) {

             
          dbRatio = (float) (10.0*log10((FLT_MIN+avgEn) /
                                        (FLT_MIN+psyOutChan->sfbEnergy[sfbOffs+sfb])));

            
          minSnrRed = msaParam->redOffs + msaParam->redRatioFac * dbRatio;

             
          minSnrRed = max(minSnrRed, msaParam->maxRed);

           
          psyOutChan->sfbMinSnr[sfbOffs+sfb] =
            (float)pow(psyOutChan->sfbMinSnr[sfbOffs+sfb], minSnrRed);

            
          psyOutChan->sfbMinSnr[sfbOffs+sfb] =
            min(minSnrLimit, psyOutChan->sfbMinSnr[sfbOffs+sfb]);
        }
      }
    }

  }

  
}

/* determine bands where avoid hole is not necessary resp. possible */
static void initAvoidHoleFlag(int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                              PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                              PSY_OUT_ELEMENT* psyOutElement,
                              const int nChannels,
                              AH_PARAM *ahParam)
{
  int ch, sfb,sfbGrp;
  float sfbEn;
  float scaleSprEn;

  
  
  /* decrease spreaded energy by 3dB for long blocks, resp. 2dB for shorts
     (avoid more holes in long blocks) */
  
  for (ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
    
     /* counting operation above */
    
      
    if (psyOutChan->windowSequence != SHORT_WINDOW) {
      
      
      scaleSprEn = 0.5f;
    }
    else {
      
      
      scaleSprEn = 0.63f;
    }
    
     
    for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){
      
       /* pointer for psyOutChan->sfbSpreadedEnergy[] */
      
      for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

         
        psyOutChan->sfbSpreadedEnergy[sfbGrp+sfb] *= scaleSprEn;
      }
    }
  }


  /* increase minSnr for local peaks, decrease it for valleys */
  
   
  if (ahParam->modifyMinSnr) {
    
    
    for(ch=0; ch<nChannels; ch++) {
      PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
      
       /* counting operation above */
      
       
      for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){
        
         /* pointers for psyOutChan->sfbEnergy[sfbGrp],
                        psyOutChan->sfbEnergy[sfbGrp+sfb],
                        psyOutChan->sfbMinSnr[]
                     */
         
        for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {
          float sfbEnm1, sfbEnp1, avgEn;
          
           
          if (sfb > 0)
            sfbEnm1 = psyOutChan->sfbEnergy[sfbGrp+sfb-1];
          else
            sfbEnm1 = psyOutChan->sfbEnergy[sfbGrp];
          
             
          if (sfb < psyOutChan->maxSfbPerGroup-1)
            sfbEnp1 = psyOutChan->sfbEnergy[sfbGrp+sfb+1];
          else
            sfbEnp1 = psyOutChan->sfbEnergy[sfbGrp+sfb];
          
           
          avgEn = (sfbEnm1 + sfbEnp1)/(float)2.0f;
          
          
          sfbEn = psyOutChan->sfbEnergy[sfbGrp+sfb];
          
          /* peak ? */
           
          if (sfbEn > avgEn) {
            float tmpMinSnr = max((float)0.8f*avgEn/sfbEn,
                                  (float)0.316f);
            
                 /* counting previous operation */
            
               /* if() */    /* max() */
            if (psyOutChan->windowSequence!=SHORT_WINDOW)
              tmpMinSnr = max(tmpMinSnr, (float)0.316f);
            else
              tmpMinSnr = max(tmpMinSnr, (float)0.5f);
            
              
            psyOutChan->sfbMinSnr[sfbGrp+sfb] =
              min(psyOutChan->sfbMinSnr[sfbGrp+sfb], tmpMinSnr);
          }
          
          /* valley ? */
             
          if (((float)2.0f*sfbEn < avgEn) && (sfbEn > (float)0.0f)) {
            float tmpMinSnr = avgEn/((float)2.0f*sfbEn) *
              psyOutChan->sfbMinSnr[sfbGrp+sfb];
            
              /* counting previous operation */
            
              
            tmpMinSnr = min((float)0.8f, tmpMinSnr);
            
               
            psyOutChan->sfbMinSnr[sfbGrp+sfb] = min(tmpMinSnr,
                                                    psyOutChan->sfbMinSnr[sfbGrp+sfb]*(float)3.16f);
          }
        }
      }
    }
  }
  
   
  if (nChannels == 2) {
    PSY_OUT_CHANNEL *psyOutChanM = &psyOutChannel[0];
    PSY_OUT_CHANNEL *psyOutChanS = &psyOutChannel[1];
    
     /* counting previous operation */
    
     /* pointers for psyOutElement->toolsInfo.msMask[],
                    psyOutChanM->sfbEnergy[],
                    psyOutChanS->sfbEnergy[],
                    psyOutChanM->sfbMinSnr[],
                    psyOutChanS->sfbMinSnr[],
                    psyOutChanM->sfbSpreadedEnergy[],
                    psyOutChanS->sfbSpreadedEnergy[]
                 */
     
    for (sfb=0; sfb<psyOutChanM->sfbCnt; sfb++) {
      
      
      if (psyOutElement->toolsInfo.msMask[sfb]) {
        float sfbEnM = psyOutChanM->sfbEnergy[sfb];
        float sfbEnS = psyOutChanS->sfbEnergy[sfb];
        float maxSfbEn = max(sfbEnM, sfbEnS);
        float maxThr = 0.25f * psyOutChanM->sfbMinSnr[sfb] * maxSfbEn;
        
           /* max() */   /* counting previous operations */
        
             /* min() */    /* max() */
        psyOutChanM->sfbMinSnr[sfb] = (float)max(psyOutChanM->sfbMinSnr[sfb],
                                                 min ( FLT_MAX, ((double)maxThr/(FLT_MIN+(double)sfbEnM))));
        
         
        if (psyOutChanM->sfbMinSnr[sfb] <= 1.0f) {
          
            
          psyOutChanM->sfbMinSnr[sfb] = min(psyOutChanM->sfbMinSnr[sfb], 0.8f);
        }
        
             /* min() */    /* max() */       
        psyOutChanS->sfbMinSnr[sfb] = (float)max(psyOutChanS->sfbMinSnr[sfb],
                                                 min ( FLT_MAX, ((double)maxThr/(FLT_MIN+(double)sfbEnS))));
        
         
        if (psyOutChanS->sfbMinSnr[sfb] <= 1.0f) {
          
            
          psyOutChanS->sfbMinSnr[sfb] = min(psyOutChanS->sfbMinSnr[sfb], 0.8f);
        }
        
         
        if (sfbEnM > psyOutChanM->sfbSpreadedEnergy[sfb]) {
          
           
          psyOutChanS->sfbSpreadedEnergy[sfb] = 0.9f * sfbEnS;
        }
        
         
        if (sfbEnS > psyOutChanS->sfbSpreadedEnergy[sfb]) {
          
           
          psyOutChanM->sfbSpreadedEnergy[sfb] = 0.9f * sfbEnM;
        }
      }
    }
  }
  
  
  /* init ahFlag (0: no ah necessary, 1: ah possible, 2: ah active */
  
  for(ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
    
     /* counting previous operation */
    
     
    for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){
      
       /* pointers for psyOutChan->sfbEnergy[],
                      psyOutChan->sfbMinSnr[],
                      psyOutChan->sfbSpreadedEnergy[],
                      ahFlag[][]
                   */
       
      for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {
        
          
        if (psyOutChan->sfbSpreadedEnergy[sfbGrp+sfb] > psyOutChan->sfbEnergy[sfbGrp+sfb] ||
            psyOutChan->sfbMinSnr[sfbGrp+sfb] > (float)1.0) {
          
          
          ahFlag[ch][sfbGrp+sfb] = NO_AH;
        }
        else {
          
          
          ahFlag[ch][sfbGrp+sfb] = AH_INACTIVE;
        }
      }
      
       /* pointer for ahFlag[][] */
       
      for (sfb=psyOutChan->maxSfbPerGroup; sfb<psyOutChan->sfbPerGroup; sfb++) {
        
        
        ahFlag[ch][sfbGrp+sfb] = NO_AH;
      }
    }
  }
  
  
}




/* constants that do not change during successive pe calculations */
static void preparePe(PE_DATA *peData,
                      PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                      float sfbFormFactor[MAX_CHANNELS][MAX_GROUPED_SFB],
                      const int nChannels,
                      const float peOffset)
{
  int ch;
  
  
  
  
  for(ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
    
     /* counting previous operation */
    
      
    prepareSfbPe(&peData->peChannelData[ch],
                 psyOutChan->sfbEnergy,
                 psyOutChan->sfbThreshold,
                 sfbFormFactor[ch],
                 psyOutChan->sfbOffsets,
                 psyOutChan->sfbCnt,
                 psyOutChan->sfbPerGroup,
                 psyOutChan->maxSfbPerGroup);
  }
  
   
  peData->offset = peOffset;
  
  
}



/* calculate pe for both channels */
static void calcPe(PE_DATA *peData,
                   PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                   const int nChannels)
{
  int ch;

  
  
   
  peData->pe = peData->offset;
  peData->constPart = 0.0f;
  peData->nActiveLines = 0.0f;
  
  
  for(ch=0; ch<nChannels; ch++) {
    PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
    PE_CHANNEL_DATA *peChanData = &peData->peChannelData[ch];
    
      /* counting previous operations */
    
      
    calcSfbPe(&peData->peChannelData[ch],
              psyOutChan->sfbEnergy,
              psyOutChan->sfbThreshold,
              psyOutChan->sfbCnt,
              psyOutChan->sfbPerGroup,
              psyOutChan->maxSfbPerGroup);
    
      
    peData->pe += peChanData->pe;
    peData->constPart += peChanData->constPart;
    peData->nActiveLines += peChanData->nActiveLines;
    
     
    psyOutChannel[ch].pe = peData->pe;                 /* update pe for stereo preprocessing */
  } 
  
  
}



/* sum the pe data only for bands where avoid hole is inactive */
static void calcPeNoAH(float *pe,
                       float *constPart,
                       float  *nActiveLines,
                       PE_DATA *peData,
                       int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                       PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                       const int nChannels)
{
   int ch, sfb,sfbGrp;

   

   
   *pe = 0.0f;
   *constPart = 0.0f;
   *nActiveLines = 0;

   
   for(ch=0; ch<nChannels; ch++) {
      PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];
      PE_CHANNEL_DATA *peChanData = &peData->peChannelData[ch];

       /* counting previous operations */

       
      for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){

         /* pointers for ahFlag[ch][sfbGrp+sfb],
                                     peChanData->sfbPe[sfbGrp+sfb],
                                     peChanData->sfbConstPart[sfbGrp+sfb],
                                     peChanData->sfbNActiveLines[sfbGrp+sfb]
                     */
         
        for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

           
          if(ahFlag[ch][sfbGrp+sfb] < AH_ACTIVE) {

             
            *pe += peChanData->sfbPe[sfbGrp+sfb];
            *constPart += peChanData->sfbConstPart[sfbGrp+sfb];
            *nActiveLines += peChanData->sfbNActiveLines[sfbGrp+sfb];
          }
        }
      }
   }

   
}


/* apply reduction formula */
static void reduceThresholds(PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                             int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                             float thrExp[MAX_CHANNELS][MAX_GROUPED_SFB],
                             const int nChannels,
                             const float redVal)
{
   int ch, sfb,sfbGrp;
   float sfbEn, sfbThr,sfbThrReduced;

   

   
   for(ch=0; ch<nChannels; ch++) {
      PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];

       /* counting previous operation */

       
      for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){

         /* pointers for psyOutChan->sfbMinSnr[sfbGrp+sfb],
                                     psyOutChan->sfbEnergy[sfbGrp+sfb],
                                     psyOutChan->sfbThreshold[sfbGrp+sfb],
                                     thrExp[ch][sfbGrp+sfb],
                                     ahFlag[ch][sfbGrp+sfb]
                     */
         
        for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

         
         sfbEn  = psyOutChan->sfbEnergy[sfbGrp+sfb];
         sfbThr = psyOutChan->sfbThreshold[sfbGrp+sfb];

          
         if (sfbEn > sfbThr) {
            /* threshold reduction formula */

             
            sfbThrReduced = (float) pow(thrExp[ch][sfbGrp+sfb]+redVal, invRedExp);

            /* avoid holes */
               
            if ((sfbThrReduced > psyOutChan->sfbMinSnr[sfbGrp+sfb] * sfbEn) && (ahFlag[ch][sfbGrp+sfb] != NO_AH)){

                
              sfbThrReduced = max(psyOutChan->sfbMinSnr[sfbGrp+sfb] * sfbEn, sfbThr);

              
              ahFlag[ch][sfbGrp+sfb] = AH_ACTIVE;
            }

            
            psyOutChan->sfbThreshold[sfbGrp+sfb] = sfbThrReduced;
         }
        }
      }
   }

   
}


static void correctThresh(PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                          int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                          PE_DATA *peData,
                          float thrExp[MAX_CHANNELS][MAX_GROUPED_SFB],
                          const float redVal,
                          const int nChannels,
                          const float deltaPe)
{
   int ch, sfb,sfbGrp;
   PSY_OUT_CHANNEL *psyOutChan;
   PE_CHANNEL_DATA *peChanData;
   float deltaSfbPe;
   float thrFactor;
   float sfbPeFactors[MAX_CHANNELS][MAX_GROUPED_SFB], normFactor;
   float sfbEn, sfbThr, sfbThrReduced;

   

   /* for each sfb calc relative factors for pe changes */
   
   normFactor = FLT_MIN; /* to avoid division by zero */

   
   for(ch=0; ch<nChannels; ch++) {
      psyOutChan = &psyOutChannel[ch];
      peChanData = &peData->peChannelData[ch];

        /* counting previous operations */

       
      for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){

         /* pointers for ahFlag[ch][sfbGrp+sfb],
                                     thrExp[ch][sfbGrp+sfb],
                                     sfbPeFactors[ch][sfbGrp+sfb],
                                     peChanData->sfbNActiveLines[sfbGrp+sfb]
                     */
         
        for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

           
         if ((ahFlag[ch][sfbGrp+sfb] < AH_ACTIVE) || (deltaPe > 0)) {

              
            sfbPeFactors[ch][sfbGrp+sfb] = peChanData->sfbNActiveLines[sfbGrp+sfb] /
                                              (thrExp[ch][sfbGrp+sfb] + redVal);

            
            normFactor += sfbPeFactors[ch][sfbGrp+sfb];
         }
         else {
            
            sfbPeFactors[ch][sfbGrp+sfb] = 0.0f;
         }
        }
      }
   }

   
   normFactor = 1.0f/normFactor;

   
   for(ch=0; ch<nChannels; ch++) {
      psyOutChan = &psyOutChannel[ch];
      peChanData = &peData->peChannelData[ch];

        /* counting previous operations */

       
      for(sfbGrp = 0;sfbGrp < psyOutChan->sfbCnt;sfbGrp+= psyOutChan->sfbPerGroup){

         /* pointers for ahFlag[ch][sfbGrp+sfb],
                                     sfbPeFactors[ch][sfbGrp+sfb],
                                     psyOutChan->sfbMinSnr[sfbGrp+sfb],
                                     peChanData->sfbNActiveLines[sfbGrp+sfb],
                                     psyOutChan->sfbEnergy[sfbGrp+sfb],
                                     psyOutChan->sfbThreshold[sfbGrp+sfb]
                     */
         
        for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {

         /* pe difference for this sfb */
         
         deltaSfbPe = sfbPeFactors[ch][sfbGrp+sfb] * normFactor * deltaPe;

          
         if (peChanData->sfbNActiveLines[sfbGrp+sfb] > (float)0.5f) {

            /* new threshold */
            
            sfbEn  = psyOutChan->sfbEnergy[sfbGrp+sfb];
            sfbThr = psyOutChan->sfbThreshold[sfbGrp+sfb];

                
            thrFactor = min(-deltaSfbPe/peChanData->sfbNActiveLines[sfbGrp+sfb],20.f);

            
            thrFactor = (float) pow(2.0f,thrFactor);

            
            sfbThrReduced = sfbThr * thrFactor;
            /* avoid hole */

              
            if ((sfbThrReduced > psyOutChan->sfbMinSnr[sfbGrp+sfb] * sfbEn) &&
                (ahFlag[ch][sfbGrp+sfb] == AH_INACTIVE)) {

                 
               sfbThrReduced = max(psyOutChan->sfbMinSnr[sfbGrp+sfb] * sfbEn, sfbThr);

               
               ahFlag[ch][sfbGrp+sfb] = AH_ACTIVE;
            }

            
            psyOutChan->sfbThreshold[sfbGrp+sfb] = sfbThrReduced;
         }
        }
      }
   }

   
}

static void reduceMinSnr(PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                         PE_DATA *peData, 
                         int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                         const int nChannels,
                         const float desiredPe)
{
  int ch, sfb, sfbSubWin;
  float deltaPe;

  

  /* start at highest freq down to 0 */
  
  sfbSubWin = psyOutChannel[0].maxSfbPerGroup;
  
   /* pointer for ahFlag[ch][sfb],
                  psyOutChannel[ch].sfbMinSnr[sfb],
                  psyOutChannel[ch].sfbThreshold[sfb],
                  psyOutChannel[ch].sfbEnergy[sfb],
                  peData->peChannelData[ch].sfbNLines[sfb],
                  peData->peChannelData[ch].sfbPe[sfb],
                  peData->peChannelData[ch].pe
               */
   
  while (peData->pe > desiredPe && sfbSubWin > 0) {
      /* while() condition */
    
    
    sfbSubWin--;
    
    /* loop over all subwindows */
    
    for (sfb=sfbSubWin; sfb<psyOutChannel[0].sfbCnt;
         sfb+=psyOutChannel[0].sfbPerGroup) {
      
      /* loop over all channels */
      
      for (ch=0; ch<nChannels; ch++) {
        
         
        if (ahFlag[ch][sfb] != NO_AH &&
            psyOutChannel[ch].sfbMinSnr[sfb] < minSnrLimit) {
          
          
          psyOutChannel[ch].sfbMinSnr[sfb] = minSnrLimit;
          
           
          psyOutChannel[ch].sfbThreshold[sfb] =
            psyOutChannel[ch].sfbEnergy[sfb] * psyOutChannel[ch].sfbMinSnr[sfb];
          
           
          deltaPe = peData->peChannelData[ch].sfbNLines[sfb] * 1.5f -
            peData->peChannelData[ch].sfbPe[sfb];
          
            
          peData->pe += deltaPe;
          
           
          peData->peChannelData[ch].pe += deltaPe;
        }
      }
      
        
      if (peData->pe <= desiredPe)
        break;
    }
  }
  
  
}


static void allowMoreHoles(PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS], 
                           PSY_OUT_ELEMENT *psyOutElement,
                           PE_DATA *peData, 
                           int ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB],
                           const AH_PARAM *ahParam,
                           const int nChannels,
                           const float desiredPe)
{
  int ch, sfb;
  float actPe;

  

   
  actPe = peData->pe;

    
  if (nChannels==2 &&
      psyOutChannel[0].windowSequence==psyOutChannel[1].windowSequence) {
    PSY_OUT_CHANNEL *psyOutChanL = &psyOutChannel[0];
    PSY_OUT_CHANNEL *psyOutChanR = &psyOutChannel[1];

     /* counting previous operations */

     /* pointers for ahFlag[0][sfb],
                                  ahFlag[1][sfb],
                                  psyOutElement->toolsInfo.msMask[sfb],
                                  psyOutChanL->sfbMinSnr[sfb],
                                  psyOutChanR->sfbMinSnr[sfb],
                                  psyOutChanL->sfbEnergy[sfb],
                                  psyOutChanR->sfbEnergy[sfb],
                                  psyOutChanL->sfbThreshold[sfb],
                                  psyOutChanR->sfbThreshold[sfb],
                                  peData->peChannelData[0].sfbPe[sfb],
                                  peData->peChannelData[1].sfbPe[sfb]
                  */
     
    for (sfb=0; sfb<psyOutChanL->sfbCnt; sfb++) {

      
      if (psyOutElement->toolsInfo.msMask[sfb]) {

           
        if (ahFlag[1][sfb] != NO_AH &&
            0.4f*psyOutChanL->sfbMinSnr[sfb]*psyOutChanL->sfbEnergy[sfb] >
            psyOutChanR->sfbEnergy[sfb]) {

          
          ahFlag[1][sfb] = NO_AH;

           
          psyOutChanR->sfbThreshold[sfb] = 2.0f * psyOutChanR->sfbEnergy[sfb];

          
          actPe -= peData->peChannelData[1].sfbPe[sfb];
        }
        else {
             
          if (ahFlag[0][sfb] != NO_AH &&
            0.4f*psyOutChanR->sfbMinSnr[sfb]*psyOutChanR->sfbEnergy[sfb] >
            psyOutChanL->sfbEnergy[sfb]) {

          
          ahFlag[0][sfb] = NO_AH;

           
          psyOutChanL->sfbThreshold[sfb] = 2.0f * psyOutChanL->sfbEnergy[sfb];

          
          actPe -= peData->peChannelData[0].sfbPe[sfb];
          }
        }
         
        if (actPe < desiredPe)
          break;
      }
    }
  }

   
  if (actPe > desiredPe) {
    int startSfb[2];
    float avgEn, minEn;
    int ahCnt;
    int enIdx;
    float en[4];
    int minSfb, maxSfb;
    int done;

     /* pointers for psyOutChannel[ch].windowSequence,
                                 startSfb[ch]
                 */
    
    for (ch=0; ch<nChannels; ch++) {

         
      if (psyOutChannel[ch].windowSequence != SHORT_WINDOW)
        startSfb[ch] = ahParam->startSfbL;
      else
        startSfb[ch] = ahParam->startSfbS;
    }

    
    avgEn = 0.0f;
    minEn = FLT_MAX;
    ahCnt = 0;

     /* pointers for startSfb[ch] */
    
    for (ch=0; ch<nChannels; ch++) {
      PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];

       /* counting previous operation */

       /* pointers for ahFlag[ch][sfb],
                                   psyOutChan->sfbMinSnr[sfb],
                                   psyOutChan->sfbEnergy[sfb],
                                   psyOutChan->sfbThreshold[sfb],
                   */
       
      for (sfb=startSfb[ch]; sfb<psyOutChan->sfbCnt; sfb++){

          
        if ((ahFlag[ch][sfb]!=NO_AH) &&
            (psyOutChan->sfbEnergy[sfb] > psyOutChan->sfbThreshold[sfb])){

            
          minEn = min(minEn, psyOutChan->sfbEnergy[sfb]);

          
          avgEn += psyOutChan->sfbEnergy[sfb];

          
          ahCnt++;
        }
      }
    }
       /* min() */  
    avgEn = min ( FLT_MAX , avgEn /(ahCnt+FLT_MIN));

    /* calc some energy borders between minEn and avgEn */
     /* pointers for en[enIdx] */
    
    for (enIdx=0; enIdx<4; enIdx++) {
          
      en[enIdx] = minEn * (float)pow(avgEn/(minEn+FLT_MIN), (2*enIdx+1)/7.0f);
    }

    
    maxSfb = psyOutChannel[0].sfbCnt - 1;

    
    minSfb = startSfb[0];

     
    if (nChannels==2) {

        
      maxSfb = max(maxSfb, psyOutChannel[1].sfbCnt-1);

        
      minSfb = min(minSfb, startSfb[1]);
    }

    
    sfb = maxSfb;
    enIdx = 0;
    done = 0;

    
    while (!done) {

       /* pointer for en[enIdx] */

      
      for (ch=0; ch<nChannels; ch++) {
        PSY_OUT_CHANNEL *psyOutChan = &psyOutChannel[ch];

         /* counting previous operation */

         /* pointers for startSfb[ch]
                                     ahFlag[ch][sfb],
                                     psyOutChan->sfbEnergy[sfb],
                                     psyOutChan->sfbThreshold[sfb]
                     */

           
        if (sfb>=startSfb[ch] && sfb < psyOutChan->sfbCnt) {

            
          if (ahFlag[ch][sfb] != NO_AH && psyOutChan->sfbEnergy[sfb] < en[enIdx]){
            
            
            ahFlag[ch][sfb] = NO_AH;
            
             
            psyOutChan->sfbThreshold[sfb] = 2.0f * psyOutChan->sfbEnergy[sfb];
            
            
            actPe -= peData->peChannelData[ch].sfbPe[sfb];
          }

           
          if (actPe < desiredPe) {
            
            
            done = 1;
            break;
          }
        }
      }
      
      sfb--;

       
      if (sfb < minSfb) {
        
        sfb = maxSfb;

        
        enIdx++;

         
        if (enIdx >= 4) {

          
          done = 1;
        }
      }
    }
  }

  
}


static void adaptThresholdsToPe(PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
                                PSY_OUT_ELEMENT* psyOutElement,
                                PE_DATA *peData,
                                const int nChannels,
                                const float desiredPe,
                                AH_PARAM *ahParam,
                                MINSNR_ADAPT_PARAM *msaParam)
{
   float noRedPe, redPe, redPeNoAH;
   float constPart, constPartNoAH;
   float nActiveLines, nActiveLinesNoAH;
   float desiredPeNoAH;
   float avgThrExp, redVal;
   int   ahFlag[MAX_CHANNELS][MAX_GROUPED_SFB];
   float thrExp[MAX_CHANNELS][MAX_GROUPED_SFB];
   int iter;

   

   
   calcThreshExp(thrExp, psyOutChannel, nChannels);

   
   adaptMinSnr(psyOutChannel, msaParam, nChannels);

   
   initAvoidHoleFlag(ahFlag, psyOutChannel, psyOutElement, nChannels, ahParam);

    
   noRedPe = peData->pe;
   constPart = peData->constPart;
   nActiveLines = peData->nActiveLines;

      
   avgThrExp = (float)pow(2.0f, (constPart - noRedPe)/(invRedExp*nActiveLines));

   /*  calculated above */   
   redVal = (float)pow(2.0f, (constPart - desiredPe)/(invRedExp*nActiveLines)) -
            avgThrExp;

    
   redVal = max(0.0f, redVal);

   
   reduceThresholds(psyOutChannel, ahFlag, thrExp, nChannels, redVal);

   
   calcPe(peData, psyOutChannel, nChannels);

    
   redPe = peData->pe;

   
   iter = 0;

   
   do {

      
     calcPeNoAH(&redPeNoAH, &constPartNoAH, &nActiveLinesNoAH,
                peData, ahFlag, psyOutChannel, nChannels);

       
     desiredPeNoAH = max(desiredPe - (redPe - redPeNoAH), 0);

     
     if (nActiveLinesNoAH > 0) {

          
       avgThrExp = (float)pow(2.0f, (constPartNoAH - redPeNoAH) / (invRedExp*nActiveLinesNoAH));

       /*  calculated above */   
       redVal   += (float)pow(2.0f, (constPartNoAH - desiredPeNoAH) / (invRedExp*nActiveLinesNoAH))
                   - avgThrExp;

        
       redVal = max(0.0f, redVal);

       
       reduceThresholds(psyOutChannel, ahFlag, thrExp, nChannels, redVal);
     }

     
     calcPe(peData, psyOutChannel, nChannels);

      
     redPe = peData->pe;

     
     iter++;

         /* while() condition */
   } while ((fabs(redPe - desiredPe) > (0.05f)*desiredPe) && (iter < 2));


     
   if (redPe < 1.15f*desiredPe) {

      
     correctThresh(psyOutChannel, ahFlag, peData, thrExp, redVal,
                   nChannels, desiredPe - redPe);
   }
   else {

      
     reduceMinSnr(psyOutChannel, peData, ahFlag, nChannels, 1.05f*desiredPe);

      
     allowMoreHoles(psyOutChannel, psyOutElement, peData, ahFlag,
                    ahParam, nChannels, 1.05f*desiredPe);
   }

   
}


static float calcBitSave(float fillLevel,
                         const float clipLow,
                         const float clipHigh,
                         const float minBitSave,
                         const float maxBitSave)
{
   float bitsave;

   

     
   fillLevel = max(fillLevel, clipLow);
   fillLevel = min(fillLevel, clipHigh);

     
   bitsave = maxBitSave - ((maxBitSave-minBitSave) / (clipHigh-clipLow)) *
                          (fillLevel-clipLow);

   

   return (bitsave);
}

static float calcBitSpend(float fillLevel,
                          const float clipLow,
                          const float clipHigh,
                          const float minBitSpend,
                          const float maxBitSpend)
{
   float bitspend;

   

     
   fillLevel = max(fillLevel, clipLow);
   fillLevel = min(fillLevel, clipHigh);

     
   bitspend = minBitSpend + ((maxBitSpend-minBitSpend) / (clipHigh-clipLow)) *
                            (fillLevel-clipLow);

   

   return (bitspend);
}


static void adjustPeMinMax(const float currPe,
                           float      *peMin,
                           float      *peMax)
{
  float minFacHi = 0.3f, maxFacHi = 1.0f, minFacLo = 0.14f, maxFacLo = 0.07f;
  float diff;
  float minDiff = currPe * (float)0.1666666667f;

  

    /* counting previous operations */

   
  if (currPe > *peMax) {

     
     diff = (currPe-*peMax) ;

      
     *peMin += diff * minFacHi;
     *peMax += diff * maxFacHi;

  } else {

     
    if (currPe < *peMin) {

     
     diff = (*peMin-currPe) ;

      
     *peMin -= diff * minFacLo;
     *peMax -= diff * maxFacLo;
  } else {

       
     *peMin += (currPe - *peMin) * minFacHi;

       
     *peMax -= (*peMax - currPe) * maxFacLo;
  }
  }

   
  if ((*peMax - *peMin) < minDiff) {
     float partLo, partHi;

       
     partLo = max((float)0.0f, currPe - *peMin);

       
     partHi = max((float)0.0f, *peMax - currPe);

        
     *peMax = currPe + partHi/(partLo+partHi) * minDiff;
     *peMin = currPe - partLo/(partLo+partHi) * minDiff;

       
     *peMin = max((float)0.0f, *peMin);
  }

  
}



static float bitresCalcBitFac( const int       bitresBits,
                                const int       maxBitresBits,
                                const float    pe,
                                const int       windowSequence,
                                const int       avgBits,
                                const float    maxBitFac,
                                ADJ_THR_STATE   *AdjThr,
                                ATS_ELEMENT     *adjThrChan)
{
   BRES_PARAM *bresParam;
   float pex;
   float fillLevel;
   float bitSave, bitSpend, bitresFac;

   

   
   fillLevel = (float)bitresBits / (float)maxBitresBits;


      
   if (windowSequence != SHORT_WINDOW)
      bresParam = &(AdjThr->bresParamLong);
   else
      bresParam = &(AdjThr->bresParamShort);

      
   pex = max(pe, adjThrChan->peMin);

      
   pex = min(pex,adjThrChan->peMax);

    
   bitSave = calcBitSave(fillLevel,
                         bresParam->clipSaveLow, bresParam->clipSaveHigh,
                         bresParam->minBitSave, bresParam->maxBitSave);

    
   bitSpend = calcBitSpend(fillLevel,
                           bresParam->clipSpendLow, bresParam->clipSpendHigh,
                           bresParam->minBitSpend, bresParam->maxBitSpend);

      
   bitresFac = (float)1.0f - bitSave +
               ((bitSpend + bitSave) / (adjThrChan->peMax - adjThrChan->peMin)) *
               (pex - adjThrChan->peMin);

   /* limit bitresFac for small bitreservoir */
      
   bitresFac = min(bitresFac, (float)1.0f - (float)0.3f + (float)bitresBits/(float)avgBits);

   /* limit bitresFac for high bitrates */
     
   bitresFac = min(bitresFac, maxBitFac);

     
   adjustPeMinMax(pe, &adjThrChan->peMin, &adjThrChan->peMax);

   

   return bitresFac;
}



void AdjThrInit(ADJ_THR_STATE  *hAdjThr,
                const float    meanPe,
                int            chBitrate)
{
  ATS_ELEMENT* atsElem = &hAdjThr->adjThrStateElem;
  MINSNR_ADAPT_PARAM *msaParam = &atsElem->minSnrAdaptParam;

  

    /* counting previous operations */

  /* common for all elements: */
  /* parameters for bitres control */
   
  hAdjThr->bresParamLong.clipSaveLow   =  0.2f;
  hAdjThr->bresParamLong.clipSaveHigh  =  0.95f;
  hAdjThr->bresParamLong.minBitSave    = -0.05f;
  hAdjThr->bresParamLong.maxBitSave    =  0.3f;
  hAdjThr->bresParamLong.clipSpendLow  =  0.2f;
  hAdjThr->bresParamLong.clipSpendHigh =  0.95f;
  hAdjThr->bresParamLong.minBitSpend   = -0.10f;
  hAdjThr->bresParamLong.maxBitSpend   =  0.4f;
  
   
  hAdjThr->bresParamShort.clipSaveLow   =  0.2f;
  hAdjThr->bresParamShort.clipSaveHigh  =  0.75f;
  hAdjThr->bresParamShort.minBitSave    = 0.0f;
  hAdjThr->bresParamShort.maxBitSave    =  0.2f;
  hAdjThr->bresParamShort.clipSpendLow  =  0.2f;
  hAdjThr->bresParamShort.clipSpendHigh =  0.75f;
  hAdjThr->bresParamShort.minBitSpend   = -0.05f;
  hAdjThr->bresParamShort.maxBitSpend   =  0.5f ;
  
    
  atsElem->peMin = (float)0.8f * meanPe;
  atsElem->peMax = (float)1.2f * meanPe;
  
   
  atsElem->peOffset = 0.0f;
  
   
  if (chBitrate < 32000) {
        
    atsElem->peOffset = max((float)50.0f, (float)(100.0f)-(float)(100.0f/32000)*(float)chBitrate);
  }
  
  /* avoid hole parameters */
   
  if (chBitrate > 20000) {
     
    atsElem->ahParam.modifyMinSnr = TRUE;
    atsElem->ahParam.startSfbL = 15;
    atsElem->ahParam.startSfbS = 3;
  }
  else {
     
    atsElem->ahParam.modifyMinSnr = FALSE;
    atsElem->ahParam.startSfbL = 0;
    atsElem->ahParam.startSfbS = 0;
  }
  
  /* minSnr adaptation */
  /* maximum reduction of minSnr goes down to minSnr^maxRed */
   
  msaParam->maxRed = (float)0.25f;
  
  /* start adaptation of minSnr for avgEn/sfbEn > startRatio */
   
  msaParam->startRatio = (float)1.e1f;
  
  /* maximum minSnr reduction to minSnr^maxRed is reached for
     avgEn/sfbEn >= maxRatio */
   
  msaParam->maxRatio = (float)1.e3f;
  
  /* helper variables to interpolate minSnr reduction for
     avgEn/sfbEn between startRatio and maxRatio */
       
  msaParam->redRatioFac = (1.0f - msaParam->maxRed) /
    (10.0f*(float)log10(msaParam->startRatio/msaParam->maxRatio));
  
      
  msaParam->redOffs = 1.0f - msaParam->redRatioFac *
    10.0f * (float) log10(msaParam->startRatio);
  
  /* pe correction */
   
  atsElem->peLast = (float)0.0f;
  atsElem->dynBitsLast = 0;
  atsElem->peCorrectionFactor = (float)1.0f;

  
}


static void calcPeCorrection(float *correctionFac,
                             const float peAct,
                             const float peLast, 
                             const int bitsLast) 
{
  

       
  if ((bitsLast > 0) && (peAct < (float)1.5f*peLast) && (peAct > (float)0.7f*peLast) &&
      ((float)1.2f*bits2pe((float)bitsLast) > peLast) && 
      ((float)0.65f*bits2pe((float)bitsLast) < peLast))
  {
    float newFac = peLast / bits2pe((float)bitsLast);

      /* counting previous operation */

    /* dead zone */
     
    if (newFac < (float)1.0f) {

         
      newFac = min((float)1.1f*newFac, (float)1.0f);

        
      newFac = max(newFac, (float)0.85f);
    }
    else {
         
      newFac = max((float)0.9f*newFac, (float)1.0f);

        
      newFac = min(newFac, (float)1.15f);
    }
      
    if (((newFac > (float)1.0f) && (*correctionFac < (float)1.0f)) ||
        ((newFac < (float)1.0f) && (*correctionFac > (float)1.0f))) {

      
      *correctionFac = (float)1.0f;
    }

    /* faster adaptation towards 1.0, slower in the other direction */
       /* if() */   
    if ((*correctionFac < (float)1.0f && newFac < *correctionFac) ||
        (*correctionFac > (float)1.0f && newFac > *correctionFac))
      *correctionFac = (float)0.85f * (*correctionFac) + (float)0.15f * newFac;
    else
      *correctionFac = (float)0.7f * (*correctionFac) + (float)0.3f * newFac;

      
    *correctionFac = min(*correctionFac, (float)1.15f);
    *correctionFac = max(*correctionFac, (float)0.85f);
  }
  else {
    
    *correctionFac = (float)1.0f;
  }

  
}


void AdjustThresholds(ADJ_THR_STATE     *adjThrState,
                      ATS_ELEMENT       *AdjThrStateElement,
                      PSY_OUT_CHANNEL   psyOutChannel[MAX_CHANNELS],
                      PSY_OUT_ELEMENT   *psyOutElement,
                      float             *chBitDistribution,
                      float            sfbFormFactor[MAX_CHANNELS][MAX_GROUPED_SFB],
                      const int         nChannels,
                      QC_OUT_ELEMENT    *qcOE,
                      const int         avgBits,
                      const int         bitresBits,
                      const int         maxBitresBits,
                      const float       maxBitFac,
                      const int         sideInfoBits)
{
   float noRedPe, grantedPe, grantedPeCorr;
   int curWindowSequence;
   PE_DATA peData;
   float bitFactor;
   int ch;

   

     
   preparePe(&peData, psyOutChannel, sfbFormFactor, nChannels, AdjThrStateElement->peOffset);

    
   calcPe(&peData, psyOutChannel, nChannels);

   
   noRedPe = peData.pe;

   
   curWindowSequence = LONG_WINDOW;

    
   if (nChannels==2) {

      
     if ((psyOutChannel[0].windowSequence == SHORT_WINDOW) ||
         (psyOutChannel[1].windowSequence == SHORT_WINDOW)) {

       
       curWindowSequence = SHORT_WINDOW;
     }
   }
   else {

     
     curWindowSequence = psyOutChannel[0].windowSequence;
   }


     
   bitFactor = bitresCalcBitFac(bitresBits, maxBitresBits,
                                noRedPe+5.0f*sideInfoBits,
                                curWindowSequence, avgBits, maxBitFac,
                                adjThrState,
                                AdjThrStateElement);

    
   grantedPe = bitFactor * bits2pe((float)avgBits);

       /* min() */  
   calcPeCorrection(&(AdjThrStateElement->peCorrectionFactor), 
                    min(grantedPe, noRedPe),
                    AdjThrStateElement->peLast, 
                    AdjThrStateElement->dynBitsLast);

    
   grantedPeCorr = grantedPe * AdjThrStateElement->peCorrectionFactor;

    
   if (grantedPeCorr < noRedPe) {

        
      adaptThresholdsToPe(psyOutChannel,
                          psyOutElement,
                          &peData,
                          nChannels,
                          grantedPeCorr,
                          &AdjThrStateElement->ahParam,
                          &AdjThrStateElement->minSnrAdaptParam);
   }

    /* pointers for chBitDistribution[],
                                peData.peChannelData[]
                */
   
   for (ch=0; ch<nChannels; ch++) {

     
     if (peData.pe) {

          
       chBitDistribution[ch] = 0.2f + (float)(1.0f-nChannels*0.2f) * (peData.peChannelData[ch].pe/peData.pe);
     } else {

       
       chBitDistribution[ch] = 0.2f;
     }
   }

    
   qcOE->pe= noRedPe;

    
   AdjThrStateElement->peLast = grantedPe;

   
}

void AdjThrUpdate(ATS_ELEMENT *AdjThrStateElement,
                  const int dynBitsUsed)
{
   

    
   AdjThrStateElement->dynBitsLast = dynBitsUsed;

   
}
