/*
  Short block grouping
*/
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

/*
* this routine does not work in-place
*/

void
groupShortData(float         *mdctSpectrum,
               float         *tmpSpectrum,
               SFB_THRESHOLD *sfbThreshold,
               SFB_ENERGY    *sfbEnergy,
               SFB_ENERGY    *sfbEnergyMS,
               SFB_ENERGY    *sfbSpreadedEnergy,
               const int      sfbCnt,
               const int     *sfbOffset,
               const float *sfbMinSnr,
               int           *groupedSfbOffset,
               int           *maxSfbPerGroup,
               float       *groupedSfbMinSnr,
               const int     noOfGroups,
               const int     *groupLen)
{
  int i,j;
  int line;
  int sfb;
  int grp;
  int wnd;
  int offset;
  int highestSfb;

  



  /* for short: regroup and  */
  /* cumulate energies und thresholds group-wise . */
 
  
    /* calculate sfbCnt */
    
    highestSfb = 0;

    
    for (wnd = 0; wnd < TRANS_FAC; wnd++)
    {
       /* sfbOffset[] */
      
      for (sfb = sfbCnt-1; sfb >= highestSfb; sfb--)
      {
         /* mdctSpectrum[] */
        
        for (line = sfbOffset[sfb+1]-1; line >= sfbOffset[sfb]; line--)
        {
          
          if (mdctSpectrum[wnd*FRAME_LEN_SHORT+line] != 0.0) break;  // this band is not completely zero
        }
         
        if (line >= sfbOffset[sfb]) break; // this band was not completely zero
      }
        
      highestSfb = max(highestSfb, sfb);
    }
     
    highestSfb = highestSfb > 0 ? highestSfb : 0;

     
    *maxSfbPerGroup = highestSfb+1;

  /* calculate sfbOffset */
  
  i = 0;
  offset = 0;

   /* groupedSfbOffset[]
                  groupLen[]
               */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
     /* sfbOffset[] */
    
    for (sfb = 0; sfb < sfbCnt; sfb++)
    {
        
      groupedSfbOffset[i++] = offset + sfbOffset[sfb] * groupLen[grp];
    }

    
    offset += groupLen[grp] * FRAME_LEN_SHORT;
  }
  
  groupedSfbOffset[i++] = FRAME_LEN_LONG;

   /* calculate minSnr */

  
  i = 0;
  offset = 0;

   /* groupedSfbMinSnr[] */
  for (grp = 0; grp < noOfGroups; grp++)
  {
     /* sfbMinSnr[] */
    
    for (sfb = 0; sfb < sfbCnt; sfb++)
    {
      
      groupedSfbMinSnr[i++] = sfbMinSnr[sfb];
    }

    
    offset += groupLen[grp] * FRAME_LEN_SHORT;
  }



  /* sum up sfbThresholds */
  
  wnd = 0;
  i = 0;

   /* groupLen[]
                  sfbThreshold->Long[]
               */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
     /* sfbThreshold->Short[][] */
    
    for (sfb = 0; sfb < sfbCnt; sfb++)
    {
      float thresh = sfbThreshold->Short[wnd][sfb];

       /* counting previous operation */

      
      for (j=1; j<groupLen[grp]; j++)
      {
        
        thresh += sfbThreshold->Short[wnd+j][sfb];
      }

      
      sfbThreshold->Long[i++] = thresh;
    }
    wnd += groupLen[grp];
  }

  /* sum up sfbEnergies left/right */
  
  wnd = 0;
  i = 0;

   /* groupLen[]
                  sfbEnergy->Long[]
               */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
     /* sfbEnergy->Short[][] */
    
    for (sfb = 0; sfb < sfbCnt; sfb++)
    {
      float energy = sfbEnergy->Short[wnd][sfb];

       /* counting previous operation */

      
      for (j=1; j<groupLen[grp]; j++)
      {
        
        energy += sfbEnergy->Short[wnd+j][sfb];
      }
      
      sfbEnergy->Long[i++] = energy;
    }
    wnd += groupLen[grp];
  }

  /* sum up sfbEnergies mid/side */
  
  wnd = 0;
  i = 0;

   /* groupLen[]
                  sfbEnergy->Long[]
               */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
     /* sfbEnergy->Short[][] */
    
    for (sfb = 0; sfb < sfbCnt; sfb++)
    {
      float energy = sfbEnergyMS->Short[wnd][sfb];

       /* counting previous operation */

      
      for (j=1; j<groupLen[grp]; j++)
      {
        
        energy += sfbEnergyMS->Short[wnd+j][sfb];
      }
      
      sfbEnergyMS->Long[i++] = energy;
    }
    wnd += groupLen[grp];
  }

  /* sum up sfbSpreadedEnergies */
  
  wnd = 0;
  i = 0;

   /* groupLen[]
                  sfbEnergy->Long[]
               */
  
  for (grp = 0; grp < noOfGroups; grp++)
  {
      /* sfbEnergy->Short[][] */
     
     for (sfb = 0; sfb < sfbCnt; sfb++)
     {
        float energy = sfbSpreadedEnergy->Short[wnd][sfb];

         /* counting previous operation */

        
        for (j=1; j<groupLen[grp]; j++)
        {
           
           energy += sfbSpreadedEnergy->Short[wnd+j][sfb];
        }
        
        sfbSpreadedEnergy->Long[i++] = energy;
     }
     wnd += groupLen[grp];
  }

  /* re-group spectrum */
    
    wnd = 0;
    i = 0;

     /* groupLen[]
                    tmpSpectrum[]
                 */
    
    for (grp = 0; grp < noOfGroups; grp++)
    {
       /* sfbOffset[] */
      
      for (sfb = 0; sfb < sfbCnt; sfb++)
      {
        
        for (j = 0; j < groupLen[grp]; j++)
        {
           /* mdctSpectrum[] */
          
          for (line = sfbOffset[sfb]; line < sfbOffset[sfb+1]; line++)
          {
            
            tmpSpectrum[i++] = mdctSpectrum[(wnd+j)*FRAME_LEN_SHORT+line];
          }
        }
      }
     wnd += groupLen[grp];
  }

   /* mdctSpectrum[]
                  tmpSpectrum[]
               */
  
  for(i=0;i<FRAME_LEN_LONG;i++)
  {
    
    mdctSpectrum[i]=tmpSpectrum[i];
  }

  
}
