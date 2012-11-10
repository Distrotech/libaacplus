/*
  Spreading of energy and weighted tonality
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */



void SpreadingMax(const int    pbCnt,
                  const float *maskLowFactor,
                  const float *maskHighFactor,
                  float       *pbSpreadedEnergy)
{
   int i;

   

   /* slope to higher frequencies */
    /* pointers for pbSpreadedEnergy[],
                                maskHighFactor[]
                */
   
   for (i=1; i<pbCnt; i++) {

         
      pbSpreadedEnergy[i] = max(pbSpreadedEnergy[i],
                                maskHighFactor[i] * pbSpreadedEnergy[i-1]);
   }

   /* slope to lower frequencies */
    /* pointers for pbSpreadedEnergy[],
                                maskLowFactor[]
                */
   
   for (i=pbCnt-2; i>=0; i--) {

         
      pbSpreadedEnergy[i] = max(pbSpreadedEnergy[i],
                                maskLowFactor[i] * pbSpreadedEnergy[i+1]);
   }

   

}
