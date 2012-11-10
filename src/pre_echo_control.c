/*
  Pre echo control
*/
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

void InitPreEchoControl(float *pbThresholdNm1,
                        int numPb,
                        float *pbThresholdQuiet)
{
  int pb;

  

   /* pbThresholdNm1[]
                  pbThresholdQuiet[]
               */
  
  for(pb=0;pb<numPb;pb++)
  {
    
    pbThresholdNm1[pb]=pbThresholdQuiet[pb];
  }

  
}



void PreEchoControl(float *pbThresholdNm1,
                    int     numPb,
                    float  maxAllowedIncreaseFactor,
                    float  minRemainingThresholdFactor,
                    float  *pbThreshold)
{
  int i;
  float tmpThreshold1, tmpThreshold2;

  

     /* pbThresholdNm1[]
                    pbThreshold[]
                 */
    
    for(i = 0; i < numPb; i++) {
      
      
      tmpThreshold1 = maxAllowedIncreaseFactor * (pbThresholdNm1[i]);
      tmpThreshold2 = minRemainingThresholdFactor * pbThreshold[i];
      
      /* copy thresholds to internal memory */
      
      pbThresholdNm1[i] = pbThreshold[i];
      
       
      if(pbThreshold[i] > tmpThreshold1) {

        
        pbThreshold[i] = tmpThreshold1;
      }

       
      if(tmpThreshold2 > pbThreshold[i]) {

        
        pbThreshold[i] = tmpThreshold2;
      }
      
    }
  
  
}
