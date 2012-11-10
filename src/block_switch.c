/*
  Block Switching
*/
#include <stdio.h>
#include <math.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


static float
IIRFilter(const float in, const float coeff[], float states[]);

static float
SrchMaxWithIndex(const float *in, int *index, int n);

static int
CalcWindowEnergy(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                 float *timeSignal,
                 int chIncrement,
                 int windowLen);

/*
  IIR high pass coeffs
*/
static const float hiPassCoeff[BLOCK_SWITCHING_IIR_LEN]=
{
  -0.5095f, 0.7548f
};

#define accWindowNrgFac         0.3f
#define oneMinusAccWindowNrgFac 0.7f
#define invAttackRatioHighBr    0.1f
#define invAttackRatioLowBr     0.056f;
#define minAttackNrg (1e+6*NORM_PCM_ENERGY)


int InitBlockSwitching(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                       const int bitRate, const int nChannels)
{
  
  
  /* select attackRatio */
     
  if ((nChannels==1 && bitRate > 24000) || (nChannels>1 && bitRate/nChannels > 16000)) {
     
    blockSwitchingControl->invAttackRatio = invAttackRatioHighBr;
  }
  else  {
     
    blockSwitchingControl->invAttackRatio = invAttackRatioLowBr;
  }
  
  
  return TRUE;
}




static const int suggestedGroupingTable[TRANS_FAC][MAX_NO_OF_GROUPS] =
{
    /* Attack in Window 0 */ {1,  3,  3,  1},
    /* Attack in Window 1 */ {1,  1,  3,  3},
    /* Attack in Window 2 */ {2,  1,  3,  2},
    /* Attack in Window 3 */ {3,  1,  3,  1},
    /* Attack in Window 4 */ {3,  1,  1,  3},
    /* Attack in Window 5 */ {3,  2,  1,  2},
    /* Attack in Window 6 */ {3,  3,  1,  1},
    /* Attack in Window 7 */ {3,  3,  1,  1}
};


int BlockSwitching(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                   float *timeSignal,
                   int chIncrement)
{
    int i,w;
    float enM1, enMax;

    

     /* blockSwitchingControl->groupLen[] */
    
    for (i = 0; i < TRANS_FAC; i++)
    {
        
        blockSwitchingControl->groupLen[i] = 0;
    }


       
    blockSwitchingControl->maxWindowNrg = SrchMaxWithIndex( &blockSwitchingControl->windowNrg[0][BLOCK_SWITCH_WINDOWS-1],
                                                            &blockSwitchingControl->attackIndex,
                                                            BLOCK_SWITCH_WINDOWS);

     
    blockSwitchingControl->attackIndex = blockSwitchingControl->lastAttackIndex;

     
    blockSwitchingControl->noOfGroups = MAX_NO_OF_GROUPS;

    
     /* blockSwitchingControl->groupLen[]
                    suggestedGroupingTable[][]
                 */
    
    for (i = 0; i < MAX_NO_OF_GROUPS; i++)
    {
        
        blockSwitchingControl->groupLen[i] = suggestedGroupingTable[blockSwitchingControl->attackIndex][i];
    }


     /* blockSwitchingControl->windowNrg[0][]
                    blockSwitchingControl->windowNrg[1][]
                    blockSwitchingControl->windowNrgF[0][]
                    blockSwitchingControl->windowNrgF[1][]
                 */
    
    for (w = 0; w < BLOCK_SWITCH_WINDOWS;w++) {

        
        blockSwitchingControl->windowNrg[0][w] = blockSwitchingControl->windowNrg[1][w];
        blockSwitchingControl->windowNrgF[0][w] = blockSwitchingControl->windowNrgF[1][w];
    }


    
    CalcWindowEnergy(blockSwitchingControl, timeSignal,chIncrement,BLOCK_SWITCH_WINDOW_LEN);


     
    blockSwitchingControl->attack = FALSE;

    
    enMax = 0.0;
    enM1 = blockSwitchingControl->windowNrgF[0][BLOCK_SWITCH_WINDOWS-1];

     /* blockSwitchingControl->windowNrgF[1][] */
    
    for (w=0; w<BLOCK_SWITCH_WINDOWS; w++) {
      
         
      blockSwitchingControl->accWindowNrg = (oneMinusAccWindowNrgFac * blockSwitchingControl->accWindowNrg) + 
                                            (accWindowNrgFac * enM1);

      
         
      if (blockSwitchingControl->windowNrgF[1][w] * blockSwitchingControl->invAttackRatio > blockSwitchingControl->accWindowNrg ) {

         
        blockSwitchingControl->attack = TRUE;
        blockSwitchingControl->lastAttackIndex = w;
      }

      
      enM1 = blockSwitchingControl->windowNrgF[1][w];

        
      enMax = max(enMax, enM1);
    }

     
    if (enMax < minAttackNrg) {

       
      blockSwitchingControl->attack = FALSE;
    }

      
    if((!blockSwitchingControl->attack) &&
       (blockSwitchingControl->lastattack))
    {
          
        if (blockSwitchingControl->attackIndex == TRANS_FAC-1)
        {
           
          blockSwitchingControl->attack = TRUE;
        }

         
        blockSwitchingControl->lastattack = FALSE;
    }
    else
    {
         
        blockSwitchingControl->lastattack = blockSwitchingControl->attack;
    }


   
  blockSwitchingControl->windowSequence =  blockSwitchingControl->nextwindowSequence;

   
  if(blockSwitchingControl->attack){

     
    blockSwitchingControl->nextwindowSequence = SHORT_WINDOW;
  }
  else{

     
    blockSwitchingControl->nextwindowSequence=LONG_WINDOW;
  }

    
  if(blockSwitchingControl->nextwindowSequence == SHORT_WINDOW){

     
    if(blockSwitchingControl->windowSequence == LONG_WINDOW)
    {
       
      blockSwitchingControl->windowSequence=START_WINDOW;
    }

      
    if(blockSwitchingControl->windowSequence == STOP_WINDOW) {

       
      blockSwitchingControl->windowSequence = SHORT_WINDOW;

       
      blockSwitchingControl->noOfGroups = 3;
      blockSwitchingControl->groupLen[0] = 3;
      blockSwitchingControl->groupLen[1] = 3;
      blockSwitchingControl->groupLen[2] = 2;
    }
  }

    
  if(blockSwitchingControl->nextwindowSequence == LONG_WINDOW){

      
    if(blockSwitchingControl->windowSequence == SHORT_WINDOW)
    {
       
      blockSwitchingControl->nextwindowSequence = STOP_WINDOW;
    }
  }

  

  return TRUE;
}



static float SrchMaxWithIndex(const float in[], int *index, int n)
{
  float max;
  int i,idx;
  
  
  
  /* Search maximum value in array and return index and value */
  
  max = 0.0f;
  idx = 0;
  
   /* in[] */
  
  for (i = 0; i < n; i++) {

     
    if (in[i+1] > max) {
      
      max = in[i+1];
      idx = i;
    }
  }
  
  *index = idx;
  
  
  
  return max;
}



static int CalcWindowEnergy(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                            float *timeSignal,
                            int chIncrement,
                            int windowLen)
{
    int w, i;
    float accuUE,accuFE;
    float tempUnfiltered, tempFiltered;

    

     /* pointers for blockSwitchingControl->windowNrg[1][w],
                                 blockSwitchingControl->windowNrgF[1][w]
                 */
    
    for (w=0; w < BLOCK_SWITCH_WINDOWS; w++) {

        
        accuUE = 0.0;
        accuFE = 0.0;

         /* pointer for timeSignal[] */
        
        for(i=0;i<windowLen;i++) {

          
          tempUnfiltered = timeSignal[(windowLen * w + i) * chIncrement] ;

           
          tempFiltered   = IIRFilter(tempUnfiltered,hiPassCoeff,blockSwitchingControl->iirStates);
          
          
          accuUE += (tempUnfiltered * tempUnfiltered);
          accuFE += (tempFiltered   * tempFiltered);
        }

        
        blockSwitchingControl->windowNrg[1][w] = accuUE;
        blockSwitchingControl->windowNrgF[1][w] =accuFE;
  }

  

  return(TRUE);
}


static float IIRFilter(const float in, const float coeff[], float states[])
{
  float accu1, accu2;
  float out;

  

  
  accu1 = coeff[1]*in;
  accu2 = coeff[1]*states[0];

  
  accu1 = accu1 - accu2;

  
  accu2 = coeff[0]*states[1];

  
  accu1 = accu1 - accu2;

  
  out = accu1;
  states[0] = in;
  states[1] = out;

  

  return out;
}


static const int synchronizedBlockTypeTable[4][4] =
{
  /*                 LONG_WINDOW   START_WINDOW   SHORT_WINDOW   STOP_WINDOW */
  /* LONG_WINDOW  */{LONG_WINDOW,  START_WINDOW,  SHORT_WINDOW,  STOP_WINDOW},
  /* START_WINDOW */{START_WINDOW, START_WINDOW,  SHORT_WINDOW,  SHORT_WINDOW},
  /* SHORT_WINDOW */{SHORT_WINDOW, SHORT_WINDOW,  SHORT_WINDOW,  SHORT_WINDOW},
  /* STOP_WINDOW  */{STOP_WINDOW,  SHORT_WINDOW,  SHORT_WINDOW,  STOP_WINDOW}
};


int SyncBlockSwitching(BLOCK_SWITCHING_CONTROL *blockSwitchingControlLeft,
                       BLOCK_SWITCHING_CONTROL *blockSwitchingControlRight,
                       const int nChannels)
{
  int i;
  int patchType = LONG_WINDOW;
  
  
  
   /* counting previous operations */
  
   
  if( nChannels == 1) { 
    /* Mono */
      
    if (blockSwitchingControlLeft->windowSequence!=SHORT_WINDOW){
      
       
      blockSwitchingControlLeft->noOfGroups = 1;
      blockSwitchingControlLeft->groupLen[0] = 1;
      
       /* blockSwitchingControlLeft->groupLen[] */
      
      for (i = 1; i < TRANS_FAC; i++) {
        
        blockSwitchingControlLeft->groupLen[i] = 0;
      }
    }
  }
  else {
    /* Stereo */
          
     
    patchType = synchronizedBlockTypeTable[patchType][blockSwitchingControlLeft->windowSequence];
    patchType = synchronizedBlockTypeTable[patchType][blockSwitchingControlRight->windowSequence];
    
     
    blockSwitchingControlLeft->windowSequence = patchType;
    blockSwitchingControlRight->windowSequence = patchType;
    
     
    if(patchType != SHORT_WINDOW) {     /* Long Blocks */

       
      blockSwitchingControlLeft->noOfGroups   = 1;
      blockSwitchingControlRight->noOfGroups  = 1;
      blockSwitchingControlLeft->groupLen[0]  = 1;
      blockSwitchingControlRight->groupLen[0] = 1;
      
       /* blockSwitchingControlLeft->groupLen[]
                      blockSwitchingControlRight->groupLen[]
                   */
      
      for (i = 1; i < TRANS_FAC; i++) {
        
        blockSwitchingControlLeft->groupLen[i] = 0;
        blockSwitchingControlRight->groupLen[i] = 0;
      }
    }
    else  {     /* Short Blocks */

        
      if(blockSwitchingControlLeft->maxWindowNrg > blockSwitchingControlRight->maxWindowNrg) {  
         
        blockSwitchingControlRight->noOfGroups = blockSwitchingControlLeft->noOfGroups;
        
         /* blockSwitchingControlRight->groupLen[]
                        blockSwitchingControlLeft->groupLen[]
                     */
        
        for (i = 0; i < TRANS_FAC; i++) {
          
          blockSwitchingControlRight->groupLen[i] = blockSwitchingControlLeft->groupLen[i];
        }
      }
      else  {
         
        blockSwitchingControlLeft->noOfGroups = blockSwitchingControlRight->noOfGroups;
        
         /* blockSwitchingControlRight->groupLen[]
                        blockSwitchingControlLeft->groupLen[]
                     */
        
        for (i = 0; i < TRANS_FAC; i++) {
          
          blockSwitchingControlLeft->groupLen[i] = blockSwitchingControlRight->groupLen[i];
        }
      }

    }

  } /*endif Mono or Stereo */
  
  
  
  return TRUE;
}
