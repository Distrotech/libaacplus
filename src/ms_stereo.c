/*
  MS stereo processing
*/
#include <math.h> /* for atan() */
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

void MsStereoProcessing(float       *sfbEnergyLeft,
                        float       *sfbEnergyRight,
                        const float *sfbEnergyMid,
                        const float *sfbEnergySide,
                        float       *mdctSpectrumLeft,
                        float       *mdctSpectrumRight,
                        float       *sfbThresholdLeft,
                        float       *sfbThresholdRight,
                        float       *sfbSpreadedEnLeft,
                        float       *sfbSpreadedEnRight,
                        int         *msDigest,
                        int         *msMask,
                        const int   sfbCnt,
                        const int   sfbPerGroup,
                        const int   maxSfbPerGroup,
                        const int   *sfbOffset,
                        float       *weightMsLrPeRatio) 
{
  int sfb,sfboffs, j, cnt = 0;
  int msMaskTrueSomewhere = 0;
  int msMaskFalseSomewhere = 0;
  float sumMsLrPeRatio = 0;

  

   /* counting previous operations */

   /* pointers for sfbThresholdLeft[sfb+sfboffs]
                                sfbThresholdRight[sfb+sfboffs]
                                sfbEnergyLeft[sfb+sfboffs]
                                sfbEnergyRight[sfb+sfboffs]
                                sfbEnergyMid[sfb+sfboffs]
                                sfbEnergySide[sfb+sfboffs]
                                sfbSpreadedEnLeft[sfb+sfboffs]
                                sfbSpreadedEnRight[sfb+sfboffs]
                                sfbOffset[sfb+sfboffs]
                                msMask[sfb+sfboffs]
               */
  
  for(sfb=0; sfb<sfbCnt; sfb+=sfbPerGroup) {

    
    for(sfboffs=0;sfboffs<maxSfbPerGroup;sfboffs++) {
      float pnlr,pnms,minThreshold;
      int useMS;
      
        
      minThreshold=min(sfbThresholdLeft[sfb+sfboffs], sfbThresholdRight[sfb+sfboffs]);
      
         /* max() */  
      pnlr = (sfbThresholdLeft[sfb+sfboffs]/
              max(sfbEnergyLeft[sfb+sfboffs],sfbThresholdLeft[sfb+sfboffs]))*
             (sfbThresholdRight[sfb+sfboffs]/
              max(sfbEnergyRight[sfb+sfboffs],sfbThresholdRight[sfb+sfboffs]));
      
         /* max() */  
      pnms=  (minThreshold/max(sfbEnergyMid[sfb+sfboffs],minThreshold))*
             (minThreshold/max(sfbEnergySide[sfb+sfboffs],minThreshold));
      
       
      sumMsLrPeRatio += (pnlr + 1.0e-9f) / (pnms + 1.0e-9f) ;

      
      cnt++;

      
      useMS = (pnms >= pnlr);
      
      
      if(useMS){

        
        msMask[sfb+sfboffs] = 1;
        msMaskTrueSomewhere = 1;

         /* pointers for mdctSpectrumLeft[],
                                     mdctSpectrumRight[]
                     */
        
        for(j=sfbOffset[sfb+sfboffs]; j<sfbOffset[sfb+sfboffs+1]; j++) {
          float tmp = mdctSpectrumLeft[j];

           /* counting operation above */

            
          mdctSpectrumLeft[j]  = 0.5f * (mdctSpectrumLeft[j] + mdctSpectrumRight[j]) ;

            
          mdctSpectrumRight[j] = 0.5f * (tmp - mdctSpectrumRight[j]) ;
        }

        
        sfbThresholdLeft[sfb+sfboffs] = sfbThresholdRight[sfb+sfboffs] = minThreshold;

        
        sfbEnergyLeft[sfb+sfboffs] = sfbEnergyMid[sfb+sfboffs];
        sfbEnergyRight[sfb+sfboffs] = sfbEnergySide[sfb+sfboffs];

           /* min() */  
        sfbSpreadedEnLeft[sfb+sfboffs] = sfbSpreadedEnRight[sfb+sfboffs] =
          min(sfbSpreadedEnLeft[sfb+sfboffs], sfbSpreadedEnRight[sfb+sfboffs]) * 0.5f;

      }else {
        
        
        msMask[sfb+sfboffs] = 0;
        msMaskFalseSomewhere = 1;
      }
    }
  }

   
  if(msMaskTrueSomewhere == 1) {

     
    if(msMaskFalseSomewhere == 1) {

      
      *msDigest = SI_MS_MASK_SOME;
    } else {

      
      *msDigest = SI_MS_MASK_ALL;
    }
  } else {

    
    *msDigest = SI_MS_MASK_NONE;
  }
  
    
  cnt = max (1,cnt);

      
  *weightMsLrPeRatio = (float) (0.28 * atan( 0.37*(sumMsLrPeRatio/cnt-6.5) ) + 1.25);
  
  
}
