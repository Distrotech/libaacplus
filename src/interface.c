/*
  Interface psychoaccoustic/quantizer
*/
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

void BuildInterface(float                  *groupedMdctSpectrum,
                    SFB_THRESHOLD           *groupedSfbThreshold,
                    SFB_ENERGY              *groupedSfbEnergy,
                    SFB_ENERGY              *groupedSfbSpreadedEnergy,
                    const SFB_ENERGY_SUM    sfbEnergySumLR,
                    const SFB_ENERGY_SUM    sfbEnergySumMS,
                    const int               windowSequence,
                    const int               windowShape,
                    const int               groupedSfbCnt,
                    const int              *groupedSfbOffset,
                    const int               maxSfbPerGroup,
                    const float          *groupedSfbMinSnr,
                    const int               noOfGroups,
                    const int              *groupLen,
                    PSY_OUT_CHANNEL        *psyOutCh) // output
{
  int j;
  int grp;
  int mask;

  

  /*
    copy values to psyOut
  */

     
  psyOutCh->maxSfbPerGroup = maxSfbPerGroup;
  psyOutCh->sfbCnt         = groupedSfbCnt;
  psyOutCh->sfbPerGroup    = groupedSfbCnt / noOfGroups;
  psyOutCh->windowSequence = windowSequence;
  psyOutCh->windowShape    = windowShape;
  psyOutCh->mdctSpectrum   = groupedMdctSpectrum;
  psyOutCh->sfbEnergy      = groupedSfbEnergy->Long;
  psyOutCh->sfbThreshold   = groupedSfbThreshold->Long;
  psyOutCh->sfbSpreadedEnergy = groupedSfbSpreadedEnergy->Long;
  

   /* psyOutCh->sfbOffsets[]
                  groupedSfbOffset[]
               */
  
  for(j=0; j<groupedSfbCnt+1; j++)
  {
    
    psyOutCh->sfbOffsets[j]=groupedSfbOffset[j];
  }

   /* psyOutCh->sfbMinSnr[]
                  groupedSfbMinSnr[]
               */
  for(j=0;j<groupedSfbCnt; j++){

    
    psyOutCh->sfbMinSnr[j] = groupedSfbMinSnr[j];
  }

  /* generate grouping mask */
  
  mask = 0;

   /* groupLen[grp] */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
    
    mask <<= 1;

    
    for (j=1; j<groupLen[grp]; j++) {

       
      mask<<=1; 
      mask|=1;
    }
  }
  /* copy energy sums to psy Out for stereo preprocessing */
  
   
  if (windowSequence != SHORT_WINDOW) {

     
    psyOutCh->sfbEnSumLR =  sfbEnergySumLR.Long;
    psyOutCh->sfbEnSumMS =  sfbEnergySumMS.Long;
  }
  else {
    int i;

     
    psyOutCh->sfbEnSumMS=0;
    psyOutCh->sfbEnSumLR=0;

     /* sfbEnergySumLR.Short[]
                    sfbEnergySumMS.Short[]
                 */
    
    for (i=0;i< TRANS_FAC; i++) {

      
      psyOutCh->sfbEnSumLR += sfbEnergySumLR.Short[i];
      psyOutCh->sfbEnSumMS += sfbEnergySumMS.Short[i];
    }
     
  }

   
  psyOutCh->groupingMask = mask;

  
}
