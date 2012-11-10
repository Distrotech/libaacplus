/*
  Scalefactor estimation
*/
#include <memory.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define MAX_SCF_DELTA  60

#define C1   -69.33295f     /* -16/3*log(MAX_QUANT+0.5-logCon)/log(2) */
#define C2     5.77078f     /* 4/log(2) */

#define LOG2_1 1.442695041f /* 1/log(2) */

#define PE_C1  3.0f         /* log(8.0)/log(2) */
#define PE_C2  1.3219281f   /* log(2.5)/log(2) */
#define PE_C3  0.5593573f   /* 1-C2/C1 */


static void
CalcFormFactorChannel(float *sfbFormFactor,
                      float *sfbNRelevantLines,
                      PSY_OUT_CHANNEL *psyOutChan)
{
  int i, j;
  int sfbOffs, sfb,sfbWidth;

  

      
  memset(sfbNRelevantLines,0,sizeof(float)*psyOutChan->sfbCnt);

 
   
  for (sfbOffs=0; sfbOffs<psyOutChan->sfbCnt; sfbOffs+=psyOutChan->sfbPerGroup){

     /* pointer for sfbFormFactor[],
                                psyOutChan->sfbEnergy[]
                                psyOutChan->sfbThreshold[]
                                psyOutChan->sfbOffsets[]
                                sfbNRelevantLines[]
                 */
     
    for (sfb=0; sfb<psyOutChan->maxSfbPerGroup; sfb++) {
      i = sfbOffs+sfb;

      
      sfbFormFactor[i] = FLT_MIN;

       
      if (psyOutChan->sfbEnergy[i] > psyOutChan->sfbThreshold[i]) {
        float avgFormFactor;

        /* calc sum of sqrt(spec) */
         /* pointer for psyOutChan->mdctSpectrum[j] */
        
        for (j=psyOutChan->sfbOffsets[i]; j<psyOutChan->sfbOffsets[i+1]; j++) {

            
          sfbFormFactor[i] += (float) sqrt(fabs(psyOutChan->mdctSpectrum[j]));
        }
         /* sfbFormFactor[i] */

        /* estimate number of relevant lines in this sfb */
        
        sfbWidth = psyOutChan->sfbOffsets[i+1]-psyOutChan->sfbOffsets[i];

         
        avgFormFactor = (float) pow(psyOutChan->sfbEnergy[i]/(float)sfbWidth, 0.25);

         
        sfbNRelevantLines[i] = sfbFormFactor[i] / avgFormFactor;
      }
    }
  }

  
}


static int improveScf(float *spec, 
                      float *expSpec,
                      short *quantSpec,
                      short *quantSpecTmp,
                      int    sfbWidth, 
                      float thresh, 
                      short    scf,
                      short    minScf,
                      float *dist, 
                      short    *minScfCalculated)
{
   float sfbDist;
   short scfBest = scf;
   int j;
   
   

    /* counting previous operations */

   /* calc real distortion */
   
   sfbDist = calcSfbDist(spec, expSpec, quantSpec, sfbWidth, scf);

   
   *minScfCalculated = scf;

   /* nmr > 1.25 -> try to improve nmr */
     
   if (sfbDist > (1.25 * thresh)) {
      short scfEstimated = scf;
      float sfbDistBest = sfbDist;
      int cnt;

       /* counting previous operations */

      /* improve by bigger scf ? */
      
      cnt = 0;

      
      while ((sfbDist > (1.25 * thresh)) && (cnt++ < 3)) {
            /* while() condition */

         
         scf++;

         
         sfbDist = calcSfbDist(spec, expSpec, quantSpecTmp, sfbWidth, scf);

          
         if (sfbDist < sfbDistBest) {

            
            scfBest = scf;
            sfbDistBest = sfbDist;

             /* pointers for quantSpec[],
                                         quantSpecTmp[]
                         */
            
            for (j=0; j<sfbWidth; j++) {

              
              quantSpec[j] = quantSpecTmp[j];
            }
         }
      }

      /* improve by smaller scf ? */
      
      cnt = 0;
      scf = scfEstimated;
      sfbDist = sfbDistBest;

      
      while ((sfbDist > (1.25 * thresh)) && (cnt++ < 1) && (scf > minScf)) {
            /* while() condition */

         
         scf--;

         
         sfbDist = calcSfbDist(spec, expSpec, quantSpecTmp, sfbWidth, scf);

          
         if (sfbDist < sfbDistBest) {

            
            scfBest = scf;
            sfbDistBest = sfbDist;

             /* pointers for quantSpec[],
                                         quantSpecTmp[]
                         */
            
            for (j=0; j<sfbWidth; j++) {

              
              quantSpec[j] = quantSpecTmp[j];
            }
         }
         
         *minScfCalculated = scf;
      }
      
      *dist = sfbDistBest;
   }
   else { /* nmr <= 1.25 -> try to find bigger scf to use less bits */
      float sfbDistBest = sfbDist; 
      float sfbDistAllowed = min(sfbDist * 1.25f,thresh);
      int cnt;

          /* min() */  /* counting previous operations */

      
      for (cnt=0; cnt<3; cnt++) {

         
         scf++;

         
         sfbDist = calcSfbDist(spec, expSpec, quantSpecTmp, sfbWidth, scf);

          
         if (sfbDist < sfbDistAllowed) {

            
           *minScfCalculated = scfBest+1;

           
           scfBest = scf;
           sfbDistBest = sfbDist;

            /* pointers for quantSpec[],
                                        quantSpecTmp[]
                        */
           
           for (j=0; j<sfbWidth; j++) {

             
             quantSpec[j] = quantSpecTmp[j];
           }
         }
      }
      
      *dist = sfbDistBest;
   }

   /* sign for quantSpec */
    /* pointers for spec[],
                                quantSpec[]
                */
   
   for (j=0; j<sfbWidth; j++) {

      
      if (spec[j] < 0) {

          
         quantSpec[j] = -quantSpec[j];
      }
   }

   

   /* return best scalefactor */
   return scfBest;
}

static int countSingleScfBits(int scf, int scfLeft, int scfRight)
{
  int scfBits;

  

    
  scfBits = bitCountScalefactorDelta(scfLeft-scf) +
            bitCountScalefactorDelta(scf-scfRight);

  

  return scfBits;
}

static float calcSingleSpecPe(int scf, float sfbConstPePart, float nLines)
{
  float specPe = 0.0f;
  float ldRatio;

  

   
  ldRatio = sfbConstPePart - (float)0.375f * (float)scf;

   
  if (ldRatio >= PE_C1) {

    
    specPe = (float)0.7f * nLines * ldRatio;
  }
  else {

     
    specPe = (float)0.7f * nLines * (PE_C2 + PE_C3*ldRatio);
  }

  
 
  return specPe;
}



static int countScfBitsDiff(short *scfOld, short *scfNew, 
                            int sfbCnt, int startSfb, int stopSfb)
{
  int scfBitsDiff = 0;
  int sfb = 0, sfbLast;
  int sfbPrev, sfbNext;

  

  /* search for first relevant sfb */
  
  sfbLast = startSfb;

   /* pointers for scfOld[sfbLast],
                               scfNew[sfbLast]
               */
  
  while ((sfbLast<stopSfb) && (scfOld[sfbLast]==SHRT_MIN)) {
      /* while() condition */

    sfbLast++;
  }

  /* search for previous relevant sfb and count diff */
  
  sfbPrev = startSfb - 1;

   /* pointers for scfOld[sfbPrev],
                               scfNew[sfbPrev]
               */
  
  while ((sfbPrev>=0) && (scfOld[sfbPrev]==SHRT_MIN)) {
      /* while() condition */

    sfbPrev--;
  }

  
  if (sfbPrev>=0) {

      
    scfBitsDiff += bitCountScalefactorDelta(scfNew[sfbPrev]-scfNew[sfbLast]) -
                   bitCountScalefactorDelta(scfOld[sfbPrev]-scfOld[sfbLast]);
  }

  /* now loop through all sfbs and count diffs of relevant sfbs */
   /* pointers for scfOld[sfb],
                               scfNew[sfb]
               */
  
  for (sfb=sfbLast+1; sfb<stopSfb; sfb++) {

     
    if (scfOld[sfb]!=SHRT_MIN) {

        
      scfBitsDiff += bitCountScalefactorDelta(scfNew[sfbLast]-scfNew[sfb]) - 
                     bitCountScalefactorDelta(scfOld[sfbLast]-scfOld[sfb]);

      
      sfbLast = sfb;
    }
  }

  /* search for next relevant sfb and count diff */
  
  sfbNext = stopSfb;

  
  while ((sfbNext<sfbCnt) && (scfOld[sfbNext]==SHRT_MIN)) {
      /* while() condition */

    sfbNext++;
  }

   
  if (sfbNext<sfbCnt) {

      
    scfBitsDiff += bitCountScalefactorDelta(scfNew[sfbLast]-scfNew[sfbNext]) -
                   bitCountScalefactorDelta(scfOld[sfbLast]-scfOld[sfbNext]);
  }

  

  return scfBitsDiff;
}

static float calcSpecPeDiff(PSY_OUT_CHANNEL *psyOutChan,
                            short *scfOld,
                            short *scfNew,
                            float *sfbConstPePart,
                            float *sfbFormFactor,
                            float *sfbNRelevantLines,
                            int startSfb, 
                            int stopSfb)
{
  float specPeDiff = 0.0f;
  int sfb;
  int sfbWidth;

  

   /* counting previous operation */

  /* loop through all sfbs and count pe difference */
   /* pointers for scfOld[sfb],
                               scfNew[sfb],
                               sfbConstPePart[sfb],
                               sfbFormFactor[sfb],
                               psyOutChan->sfbOffsets[sfb+1],
                               psyOutChan->sfbEnergy[sfb]
               */
  
  for (sfb=startSfb; sfb<stopSfb; sfb++) {

     
    if (scfOld[sfb]!=SHRT_MIN) {
      float ldRatioOld, ldRatioNew, pOld, pNew;

      
      sfbWidth = psyOutChan->sfbOffsets[sfb+1] - psyOutChan->sfbOffsets[sfb];

       
      if (sfbConstPePart[sfb]==FLT_MIN) {

           
        sfbConstPePart[sfb] = (float) log(psyOutChan->sfbEnergy[sfb] * (float)6.75f /
                                          sfbFormFactor[sfb]) * LOG2_1 ;
      }

       
      ldRatioOld = sfbConstPePart[sfb] - 0.375f * scfOld[sfb];
      ldRatioNew = sfbConstPePart[sfb] - 0.375f * scfNew[sfb];

       
      if (ldRatioOld >= PE_C1) {

        
        pOld = ldRatioOld;
      }
      else {

         
        pOld = PE_C2 + PE_C3 * ldRatioOld;
      }

       
      if (ldRatioNew >= PE_C1) {

        
        pNew = ldRatioNew;
      }
      else {

         
        pNew = PE_C2 + PE_C3 * ldRatioNew;
      }

        
      specPeDiff += (float)0.7f * sfbNRelevantLines[sfb] * (pNew - pOld);
    }
  }

  

  return specPeDiff;
}



static void assimilateSingleScf(PSY_OUT_CHANNEL *psyOutChan,
                                float *expSpec,
                                short *quantSpec,
                                short *quantSpecTmp,
                                short *scf, 
                                short *minScf,
                                float *sfbDist, 
                                float *sfbConstPePart,
                                float *sfbFormFactor,
                                float *sfbNRelevantLines,
                                short *minScfCalculated,
                                int restartOnSuccess)
{
  int sfbLast, sfbAct, sfbNext;
  short scfAct, *scfLast, *scfNext, scfMin, scfMax;
  int sfbWidth, sfbOffs;
  float en;
  float sfbPeOld, sfbPeNew;
  float sfbDistNew;
  int j;
  int success = 0;
  float deltaPe = 0.0f, deltaPeNew, deltaPeTmp;
  short prevScfLast[MAX_GROUPED_SFB], prevScfNext[MAX_GROUPED_SFB];
  float deltaPeLast[MAX_GROUPED_SFB];
  int updateMinScfCalculated;

  

   /* counting previous operations */

   /* pointers for prevScfLast[],
                               prevScfNext[],
                               deltaPeLast[]
               */
   
  for(j=0;j<psyOutChan->sfbCnt;j++){

    
    prevScfLast[j]=SHRT_MAX;
    prevScfNext[j]=SHRT_MAX;
    deltaPeLast[j]=FLT_MAX;
  }

    
 
 
  
  sfbLast = -1;
  sfbAct  = -1;
  sfbNext = -1;
  scfLast = 0;
  scfNext = 0;
  scfMin  = SHRT_MAX;
  scfMax  = SHRT_MAX;

   /* pointers for scf[sfbNext],
                                scf[sfbAct],
                                minScf[sfbAct],
                                minScfCalculated[sfbAct],
                                prevScfLast[sfbAct],
                                prevScfNext[sfbAct],
                                deltaPeLast[sfbAct],
                                sfbDist[sfbAct],
                                sfbConstPePart[sfbAct],
                                sfbFormFactor[sfbAct],
                                sfbNRelevantLines[sfbAct],
                                psyOutChan->sfbOffsets[sfbAct],
                                psyOutChan->sfbEnergy[sfbAct]
                */
   
  do {
    /* search for new relevant sfb */
    
    sfbNext++;

     
    while ((sfbNext < psyOutChan->sfbCnt) && (scf[sfbNext] == SHRT_MIN)) {
        /* while() condition */

      sfbNext++;
    }

       
    if ((sfbLast>=0) && (sfbAct>=0) && (sfbNext<psyOutChan->sfbCnt)) {
      /* relevant scfs to the left and to the right */

      
      scfAct  = scf[sfbAct];

      
      scfLast = scf + sfbLast;
      scfNext = scf + sfbNext;

        
      scfMin  = min(*scfLast, *scfNext);

        
      scfMax  = max(*scfLast, *scfNext); 
    }
    else {

         
      if ((sfbLast==-1) && (sfbAct>=0) && (sfbNext<psyOutChan->sfbCnt)) {
      /* first relevant scf */

      
      scfAct  = scf[sfbAct];

      
      scfLast = &scfAct;

      
      scfNext = scf + sfbNext;

      
      scfMin  = *scfNext;

      
      scfMax  = *scfNext;
    }
    else {

         
      if ((sfbLast>=0) && (sfbAct>=0) && (sfbNext==psyOutChan->sfbCnt)) {
      /* last relevant scf */

      
      scfAct  = scf[sfbAct];

      
      scfLast = scf + sfbLast;

      
      scfNext = &scfAct;

      
      scfMin  = *scfLast;

      
      scfMax  = *scfLast;
    }
    }
    }

     
    if (sfbAct>=0) {
     
        
      scfMin = max(scfMin, minScf[sfbAct]);
    }

       
    if ((sfbAct >= 0) && 
        (sfbLast>=0 || sfbNext<psyOutChan->sfbCnt) && 
        (scfAct > scfMin) && 
        (scfAct <= scfMin+MAX_SCF_DELTA) &&
        (scfAct >= scfMax-MAX_SCF_DELTA) &&
        (*scfLast != prevScfLast[sfbAct] || 
         *scfNext != prevScfNext[sfbAct] || 
         deltaPe < deltaPeLast[sfbAct])) {
      /* bigger than neighbouring scf found, try to use smaller scf */
      
      success = 0;

      
      sfbWidth = psyOutChan->sfbOffsets[sfbAct+1] - 
                 psyOutChan->sfbOffsets[sfbAct];

      
      sfbOffs = psyOutChan->sfbOffsets[sfbAct];

      /* estimate required bits for actual scf */
      
      en = psyOutChan->sfbEnergy[sfbAct];

       
      if (sfbConstPePart[sfbAct] == FLT_MIN) {

           
        sfbConstPePart[sfbAct] = (float)log(en*(float)6.75f/sfbFormFactor[sfbAct]) * LOG2_1 ;
      }
                                  
        
      sfbPeOld = calcSingleSpecPe(scfAct, sfbConstPePart[sfbAct], 
                                  sfbNRelevantLines[sfbAct]) +
                 countSingleScfBits(scfAct, *scfLast, *scfNext);

      
      deltaPeNew = deltaPe;
      updateMinScfCalculated = 1;

      
      do {
        /* estimate required bits for smaller scf */

        
        scfAct--;

        /* check only if the same check was not done before */
          
        if (scfAct < minScfCalculated[sfbAct] && scfAct>=scfMax-MAX_SCF_DELTA){
          /* estimate required bits for new scf */

            
          sfbPeNew = calcSingleSpecPe(scfAct, sfbConstPePart[sfbAct], 
                                      sfbNRelevantLines[sfbAct]) +
                     (float)countSingleScfBits(scfAct, *scfLast, *scfNext);

          /* use new scf if no increase in pe and 
             quantization error is smaller */
          
          deltaPeTmp = deltaPe + sfbPeNew - sfbPeOld;

           
          if (deltaPeTmp < (float)10.0f) {
            /* distortion of new scf */

              
            sfbDistNew = calcSfbDist(psyOutChan->mdctSpectrum+sfbOffs,
                                     expSpec+sfbOffs,
                                     quantSpecTmp+sfbOffs,
                                     sfbWidth,
                                     scfAct);

             
            if (sfbDistNew < sfbDist[sfbAct]) {
              /* success, replace scf by new one */

              
              scf[sfbAct] = scfAct;
              sfbDist[sfbAct] = sfbDistNew;

               /* pointers for psyOutChan->mdctSpectrum[],
                                           quantSpec[],
                                           quantSpecTmp[]
                           */
              
              for (j=sfbOffs; j<sfbOffs+sfbWidth; j++) {

                
                quantSpec[j] = quantSpecTmp[j];

                /* sign for quant spec */
                
                if (psyOutChan->mdctSpectrum[j] < 0.0f) {
                   
                  quantSpec[j] = -quantSpec[j];
                }
              }
              
              deltaPeNew = deltaPeTmp;
              success = 1;
            }

            /* mark as already checked */
            
            if (updateMinScfCalculated) {
              
              minScfCalculated[sfbAct] = scfAct;
            }
          }
          else {
            /* from this scf value on not all new values have been checked */
            
            updateMinScfCalculated = 0;
          }
        }
      } while (scfAct > scfMin);

      
      deltaPe = deltaPeNew;

      /* save parameters to avoid multiple computations of the same sfb */
      
      prevScfLast[sfbAct] = *scfLast;
      prevScfNext[sfbAct] = *scfNext;
      deltaPeLast[sfbAct] = deltaPe;
    }

     
    if (success && restartOnSuccess) {
      /* start again at first sfb */

      
      sfbLast = -1;
      sfbAct  = -1;
      sfbNext = -1;
      scfLast = 0;
      scfNext = 0;
      scfMin  = SHRT_MAX;
      scfMax  = SHRT_MAX;
      success = 0;
    }
    else {
      /* shift sfbs for next band */

      
      sfbLast = sfbAct;
      sfbAct  = sfbNext;
    }
  } while (sfbNext < psyOutChan->sfbCnt);

  
}


static void assimilateMultipleScf(PSY_OUT_CHANNEL *psyOutChan,
                                  float *expSpec,
                                  short *quantSpec,
                                  short *quantSpecTmp,
                                  short *scf, 
                                  short *minScf,
                                  float *sfbDist, 
                                  float *sfbConstPePart,
                                  float *sfbFormFactor,
                                  float *sfbNRelevantLines)
{
  int sfb, startSfb, stopSfb;
  short scfTmp[MAX_GROUPED_SFB], scfMin, scfMax, scfAct;
  int possibleRegionFound;
  int sfbWidth, sfbOffs, j;
  float sfbDistNew[MAX_GROUPED_SFB], distOldSum, distNewSum;
  int deltaScfBits;
  float deltaSpecPe;
  float deltaPe = 0.0f, deltaPeNew;
  int sfbCnt = psyOutChan->sfbCnt;

  

   /* counting previous operations */

  /* calc min and max scalfactors */
  
  scfMin = SHRT_MAX;
  scfMax = SHRT_MIN;

   /* pointer for scf[sfb] */
  
  for (sfb=0; sfb<sfbCnt; sfb++) {

     
    if (scf[sfb]!=SHRT_MIN) {

        
      scfMin = min(scfMin, scf[sfb]);

        
      scfMax = max(scfMax, scf[sfb]);
    }
  }

    
  if (scfMax != SHRT_MIN && scfMax <= scfMin+MAX_SCF_DELTA) {

    
    scfAct = scfMax;

     /* pointers for scf[sfb],
                                 minScf[sfb],
                                 scfTmp[sfb],
                                 psyOutChan->sfbOffsets[sfb],
                                 sfbDistNew[sfb],
                                 psyOutChan->sfbThreshold[sfb]
                 */
    
    do {
      /* try smaller scf */

      
      scfAct--;

          
      memcpy(scfTmp,scf,MAX_GROUPED_SFB*sizeof(short));

      
      stopSfb = 0;

      
      do {
        /* search for region where all scfs are bigger than scfAct */

        
        sfb = stopSfb;

        
        while (sfb<sfbCnt && (scf[sfb]==SHRT_MIN || scf[sfb] <= scfAct)) {
            /* while() condition */

          
          sfb++;
        }

        
        startSfb = sfb;

        
        sfb++;

        
        while (sfb<sfbCnt && (scf[sfb]==SHRT_MIN || scf[sfb] > scfAct)) {
            /* while() condition */

          
          sfb++;
        }

        
        stopSfb = sfb;

        /* check if in all sfb of a valid region scfAct >= minScf[sfb] */
        
        possibleRegionFound = 0;

         
        if (startSfb < sfbCnt) {

          
          possibleRegionFound = 1;

          
          for (sfb=startSfb; sfb<stopSfb; sfb++) {

             
            if (scf[sfb] != SHRT_MIN) {

               
              if (scfAct < minScf[sfb]) {

                
                possibleRegionFound = 0;
                break;
              }
            }

          }
        }

        
        if (possibleRegionFound) { /* region found */

          /* replace scfs in region by scfAct */
          
          for (sfb=startSfb; sfb<stopSfb; sfb++) {

             
            if (scfTmp[sfb] != SHRT_MIN) {

              
              scfTmp[sfb] = scfAct;
            }
          }

          /* estimate change in bit demand for new scfs */
          
          deltaScfBits = countScfBitsDiff(scf,scfTmp,sfbCnt,startSfb,stopSfb);

          
          deltaSpecPe = calcSpecPeDiff(psyOutChan, scf, scfTmp, sfbConstPePart,
                                       sfbFormFactor, sfbNRelevantLines, 
                                       startSfb, stopSfb);

          
          deltaPeNew = deltaPe + (float)deltaScfBits + deltaSpecPe;

          /* new bit demand small enough ? */
           
          if (deltaPeNew < (float)10.0f) {

            /* quantize and calc sum of new distortion */
            
            distOldSum = distNewSum = 0.0f;

            
            for (sfb=startSfb; sfb<stopSfb; sfb++) {

               
              if (scfTmp[sfb] != SHRT_MIN) {

                
                distOldSum += sfbDist[sfb];

                
                sfbWidth = psyOutChan->sfbOffsets[sfb+1] - 
                           psyOutChan->sfbOffsets[sfb];

                
                sfbOffs = psyOutChan->sfbOffsets[sfb];

                   
                sfbDistNew[sfb] = calcSfbDist(psyOutChan->mdctSpectrum+sfbOffs, 
                                              expSpec+sfbOffs,
                                              quantSpecTmp+sfbOffs,
                                              sfbWidth, 
                                              scfAct);

                 
                if (sfbDistNew[sfb] > psyOutChan->sfbThreshold[sfb]) {
                  /* no improvement, skip further dist. calculations */

                  
                  distNewSum = (float)2.0f * distOldSum;
                  break;
                }

                
                distNewSum += sfbDistNew[sfb];
              }
            }

            /* distortion smaller ? -> use new scalefactors */
             
            if (distNewSum < distOldSum) {

              
              deltaPe = deltaPeNew;

              
              for (sfb=startSfb; sfb<stopSfb; sfb++) {

                 
                if (scf[sfb] != SHRT_MIN) {

                  
                  sfbWidth = psyOutChan->sfbOffsets[sfb+1] - 
                             psyOutChan->sfbOffsets[sfb];

                  
                  sfbOffs = psyOutChan->sfbOffsets[sfb];
                  scf[sfb] = scfAct;
                  sfbDist[sfb] = sfbDistNew[sfb];

                   /* pointers for psyOutChan->mdctSpectrum[],
                                               quantSpec[],
                                               quantSpecTmp[]
                               */
                  
                  for (j=sfbOffs; j<sfbOffs+sfbWidth; j++) {

                    
                    quantSpec[j] = quantSpecTmp[j];

                    /* sign for quant spec */
                    
                    if (psyOutChan->mdctSpectrum[j] < 0.0f) {

                       
                      quantSpec[j] = -quantSpec[j];
                    }
                  }
                }
              }
            }

          }
        }


      } while (stopSfb <= sfbCnt);

    } while (scfAct > scfMin);
  }

  
}




static void
estimateScaleFactorsChannel(AACRam_t *aacram,
                            PSY_OUT_CHANNEL *psyOutChan,
                            short *scf,
                            int *globalGain,
                            float *sfbFormFactor,
                            float *sfbNRelevantLines,
                            short *quantSpec)
{
  int   i, j;
  float thresh, energy, energyPart, thresholdPart;
  float scfFloat;
  short scfInt, minScf, maxScf;
  float maxSpec;
  float sfbDist[MAX_GROUPED_SFB];
  short minSfMaxQuant[MAX_GROUPED_SFB];
  short minScfCalculated[MAX_GROUPED_SFB];

  



   /* pointer for quantSpec[],
                  expSpec[]
               */
  
  for (i=0; i<FRAME_LEN_LONG; i++) {
    
    
    aacram->expSpec[i] = 0.0f;
    quantSpec[i] = 0;
  }

   /* pointer for psyOutChan->sfbThreshold[i],
                              psyOutChan->sfbEnergy[i],
                              psyOutChan->sfbOffsets[i],
                              scf[i],
                              minSfMaxQuant[i],
                              sfbFormFactor[i],
                              sfbDist[i]
               */
   
  for(i=0; i<psyOutChan->sfbCnt; i++) {

    
    thresh = psyOutChan->sfbThreshold[i];
    energy = psyOutChan->sfbEnergy[i];

    
    maxSpec = 0.0f;

    /* maximum of spectrum */
     /* pointer for psyOutChan->mdctSpectrum[j] */
    
    for(j=psyOutChan->sfbOffsets[i]; j<psyOutChan->sfbOffsets[i+1]; j++ ){

        /* max() */   /* for -mdctSpectrum[j] */
      maxSpec = max(maxSpec, psyOutChan->mdctSpectrum[j] > 0.0 ? psyOutChan->mdctSpectrum[j]:-psyOutChan->mdctSpectrum[j] );
    }

    /* scfs without energy or with thresh>energy are marked with SHRT_MIN */
    
    scf[i] = SHRT_MIN;
    minSfMaxQuant[i]=SHRT_MIN; 

      
    if( (maxSpec > 0.0) && (energy > thresh) ){

         
        energyPart = (float) log10(sfbFormFactor[i]);

        /* influence of allowed distortion */
          
        thresholdPart = (float) log10(6.75*thresh+FLT_MIN);

         
        scfFloat = 8.8585f *(thresholdPart - energyPart);    /* scf calc */

        
        scfInt = (int)floor(scfFloat);           /* integer scalefactor */

        /* avoid quantized values bigger than MAX_QUANT */
            
        minSfMaxQuant[i] = (int)ceil(C1 + C2*log(maxSpec));

          
        scfInt = max(scfInt, minSfMaxQuant[i]);


        /* find better scalefactor with analysis by synthesis */
         
        calcExpSpec(aacram->expSpec+psyOutChan->sfbOffsets[i],
                    psyOutChan->mdctSpectrum+psyOutChan->sfbOffsets[i],
                    psyOutChan->sfbOffsets[i+1]-psyOutChan->sfbOffsets[i]);
        
         
        scfInt = improveScf(psyOutChan->mdctSpectrum+psyOutChan->sfbOffsets[i],
                            aacram->expSpec+psyOutChan->sfbOffsets[i],
                            quantSpec+psyOutChan->sfbOffsets[i],
                            aacram->quantSpecTmp+psyOutChan->sfbOffsets[i],
                            psyOutChan->sfbOffsets[i+1]-psyOutChan->sfbOffsets[i],
                            thresh, scfInt, minSfMaxQuant[i], 
                            &sfbDist[i], &minScfCalculated[i]);
        
        
        scf[i] = scfInt;
    }
  }

  
  {
     /* try to decrease scf differences */
    float sfbConstPePart[MAX_GROUPED_SFB];

     /* pointer for sfbConstPePart[] */
     
    for(i=0;i<psyOutChan->sfbCnt;i++) {
      
      sfbConstPePart[i]=FLT_MIN;
    }
     
    
    assimilateSingleScf(psyOutChan, aacram->expSpec, quantSpec, aacram->quantSpecTmp, scf,
                         minSfMaxQuant, sfbDist, sfbConstPePart,
                         sfbFormFactor, sfbNRelevantLines, minScfCalculated, 1);

    
    assimilateMultipleScf(psyOutChan, aacram->expSpec, quantSpec, aacram->quantSpecTmp, scf,
                          minSfMaxQuant, sfbDist, sfbConstPePart,
                          sfbFormFactor, sfbNRelevantLines);
  }

  /* get max scalefac for global gain */
  
  maxScf = SHRT_MIN;
  minScf = SHRT_MAX;

   /* pointer for scf[] */
   
  for(i=0; i<psyOutChan->sfbCnt; i++) {

     
    if( maxScf < scf[i] ) {

      
      maxScf = scf[i];
    }

      
    if ((scf[i] != SHRT_MIN) && (minScf > scf[i])) {

      
      minScf = scf[i];
    }
  }

  /* limit scf delta */
   /* pointer for scf[] */
   
  for(i=0; i<psyOutChan->sfbCnt; i++) {

     
    if ((scf[i] != SHRT_MIN) && (minScf+MAX_SCF_DELTA) < scf[i]) {

       
      scf[i] = minScf + MAX_SCF_DELTA;

      /* changed bands need to be quantized again */
        /* ??? */
      sfbDist[i] =
        calcSfbDist(psyOutChan->mdctSpectrum+psyOutChan->sfbOffsets[i],
                    aacram->expSpec+psyOutChan->sfbOffsets[i],
                    quantSpec+psyOutChan->sfbOffsets[i],
                    psyOutChan->sfbOffsets[i+1]-psyOutChan->sfbOffsets[i],
                    scf[i]);
    }
  }

  /* new maxScf if any scf has been limited */
    
  maxScf = min( (minScf+MAX_SCF_DELTA), maxScf);

  /* calc loop scalefactors, if spec is not all zero (i.e. maxScf == -99) */
   
  if( maxScf > SHRT_MIN ) {

    
    *globalGain = maxScf;

     /* pointer for scf[],
                                psyOutChan->sfbOffsets[i]
                 */
     
    for(i=0; i<psyOutChan->sfbCnt; i++) {

      /* intermediate bands with SHRT_MIN get the value of previous band */
       
      if( scf[i] == SHRT_MIN ) {

        
        scf[i] = 0;

        /* set band explicitely to zero */
         /* pointer for psyOutChan->mdctSpectrum[j] */
        
        for (j=psyOutChan->sfbOffsets[i]; j<psyOutChan->sfbOffsets[i+1]; j++) {

          
          psyOutChan->mdctSpectrum[j] = 0.0f;
        }
      }
      else {
         
        scf[i] = maxScf - scf[i];
      }
    }
  }
  else{

    
    *globalGain = 0;

    /* set spectrum explicitely to zero */
     /* pointer for scf[],
                                 psyOutChan->sfbOffsets[i]
                 */
     
    for(i=0; i<psyOutChan->sfbCnt; i++) {

      
      scf[i]=0;

       /* pointer for psyOutChan->mdctSpectrum[j] */
      
      for (j=psyOutChan->sfbOffsets[i]; j<psyOutChan->sfbOffsets[i+1]; j++) {

        
        psyOutChan->mdctSpectrum[j] = 0.0f;
      }
    }
  }

  
}


void
CalcFormFactor(float sfbFormFactor[MAX_CHANNELS][MAX_GROUPED_SFB],
               float sfbNRelevantLines[MAX_CHANNELS][MAX_GROUPED_SFB],
               PSY_OUT_CHANNEL psyOutChannel[MAX_CHANNELS],
               const int nChannels)
{
   int j;

   

   
   for (j=0; j<nChannels; j++) {

      
      CalcFormFactorChannel(sfbFormFactor[j], sfbNRelevantLines[j],&psyOutChannel[j]);
   }

   
}


void
EstimateScaleFactors(AACRam_t *aacram,
                     PSY_OUT_CHANNEL psyOutChannel[MAX_CHANNELS],
                     QC_OUT_CHANNEL  qcOutChannel[MAX_CHANNELS],
                     float  sfbFormFactor[MAX_CHANNELS][MAX_GROUPED_SFB],
                     float sfbNRelevantLines[MAX_CHANNELS][MAX_GROUPED_SFB],
                     const    int nChannels)
{
   int j;

   

   
   for (j=0; j<nChannels; j++) {

      
      estimateScaleFactorsChannel(aacram,
                                  &psyOutChannel[j],
                                  qcOutChannel[j].scf,
                                  &(qcOutChannel[j].globalGain),
                                  sfbFormFactor[j],
                                  sfbNRelevantLines[j],
                                  qcOutChannel[j].quantSpec);
   }

  
}
