/*
  Perceptual entropie module
*/
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define LOG2_1    1.442695041f
#define C1        3.0f        /* log(8.0)/log(2) */
#define C2        1.3219281f  /* log(2.5)/log(2) */
#define C3        0.5593573f  /* 1-C2/C1 */


/* constants that do not change during successive pe calculations */
void prepareSfbPe(PE_CHANNEL_DATA *peChanData,
                  const float *sfbEnergy,
                  const float *sfbThreshold,
                  const float *sfbFormFactor,
                  const int     *sfbOffset,
                  const int     sfbCnt,
                  const int     sfbPerGroup,
                  const int     maxSfbPerGroup)
{
   int sfbGrp,sfb;
   int sfbWidth;
   float avgFormFactor;

   

   
   for(sfbGrp = 0;sfbGrp < sfbCnt;sfbGrp+=sfbPerGroup){

     /* pointers for sfbEnergy[],
                                 sfbThreshold[],
                                 sfbOffset[],
                                 sfbFormFactor[],
                                 peChanData->sfbNLines[],
                                 peChanData->sfbLdEnergy[]
                 */
    
    for (sfb=0; sfb<maxSfbPerGroup; sfb++) {

       
      if (sfbEnergy[sfbGrp+sfb] > sfbThreshold[sfbGrp+sfb]) {

         
         sfbWidth = sfbOffset[sfbGrp+sfb+1] - sfbOffset[sfbGrp+sfb];

         /* estimate number of active lines */
          
         avgFormFactor = (float) pow(sfbEnergy[sfbGrp+sfb]/(float)sfbWidth, 0.25f);

          
         peChanData->sfbNLines[sfbGrp+sfb] =

         sfbFormFactor[sfbGrp+sfb]/avgFormFactor;

          /* ld(sfbEn) */
           
         peChanData->sfbLdEnergy[sfbGrp+sfb] = (float) (log(sfbEnergy[sfbGrp+sfb]) * LOG2_1);
      }
      else {

         
         peChanData->sfbNLines[sfbGrp+sfb] = 0.0f;
         peChanData->sfbLdEnergy[sfbGrp+sfb] = 0.0f;
      }
    }
   }

   
}


void calcSfbPe(PE_CHANNEL_DATA *peChanData,
               const float *sfbEnergy,
               const float *sfbThreshold,
               const int    sfbCnt,
               const int    sfbPerGroup,
               const int    maxSfbPerGroup)
{
   int sfbGrp,sfb;
   float nLines;
   float ldThr, ldRatio;

   

    
   peChanData->pe = 0.0f;
   peChanData->constPart = 0.0f;
   peChanData->nActiveLines = 0.0f;

   
   for(sfbGrp = 0;sfbGrp < sfbCnt;sfbGrp+=sfbPerGroup){

     /* pointers for sfbEnergy[]
                                 sfbThreshold[]
                                 peChanData->sfbLdEnergy[]
                                 peChanData->sfbNLines[]
                                 peChanData->sfbPe[]
                                 peChanData->sfbConstPart[]
                                 peChanData->sfbLdEnergy[]
                                 peChanData->sfbNActiveLines[]
                 */
    
    for (sfb=0; sfb<maxSfbPerGroup; sfb++) {

       
      if (sfbEnergy[sfbGrp+sfb] > sfbThreshold[sfbGrp+sfb]) {

          
         ldThr = (float)log(sfbThreshold[sfbGrp+sfb]) * LOG2_1;

         
         ldRatio = peChanData->sfbLdEnergy[sfbGrp+sfb] - ldThr;

         
         nLines = peChanData->sfbNLines[sfbGrp+sfb];

          
         if (ldRatio >= C1) {

             
            peChanData->sfbPe[sfbGrp+sfb] = nLines * ldRatio;

             
            peChanData->sfbConstPart[sfbGrp+sfb] = nLines*peChanData->sfbLdEnergy[sfbGrp+sfb];
         }
         else {

              
            peChanData->sfbPe[sfbGrp+sfb] = nLines * (C2 + C3 * ldRatio);

              
            peChanData->sfbConstPart[sfbGrp+sfb] = nLines *
               (C2 + C3 * peChanData->sfbLdEnergy[sfbGrp+sfb]);

            
            nLines = nLines * C3;
         }

         
         peChanData->sfbNActiveLines[sfbGrp+sfb] = nLines;
      }
      else {

         
         peChanData->sfbPe[sfbGrp+sfb] = 0.0f;
         peChanData->sfbConstPart[sfbGrp+sfb] = 0.0f;
         peChanData->sfbNActiveLines[sfbGrp+sfb] = 0.0;
      }

        
      peChanData->pe += peChanData->sfbPe[sfbGrp+sfb];
      peChanData->constPart += peChanData->sfbConstPart[sfbGrp+sfb];
      peChanData->nActiveLines += peChanData->sfbNActiveLines[sfbGrp+sfb];
    }
   }

   
}
