/*
  Channel Mapping
*/
#include <string.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


static int initElement (ELEMENT_INFO* elInfo, ELEMENT_TYPE elType) {

  int error=0;

  

   /* counting previous operations */

   
  elInfo->elType=elType;

  
  switch(elInfo->elType) {

  case ID_SCE:
     
    elInfo->nChannelsInEl=1;

     
    elInfo->ChannelIndex[0]=0;

     
    elInfo->instanceTag=0;
    break;

  case ID_CPE:

     
    elInfo->nChannelsInEl=2;

     
    elInfo->ChannelIndex[0]=0;
    elInfo->ChannelIndex[1]=1;

     
    elInfo->instanceTag=0;
    break;

  default:
    
    error=1;
  }

  

  return error;
}

int InitElementInfo (int nChannels, ELEMENT_INFO* elInfo) {

  int error=0;

  

   /* counting previous operation */


  
  switch(nChannels) {

  case 1: 
    
    initElement(elInfo, ID_SCE);
    break;

  case 2:
    
    initElement(elInfo, ID_CPE);
    break;

  default:
    
    error=1;

    
    return error;
  }

  

  return 0;
}


int InitElementBits(ELEMENT_BITS*   elementBits,
                    ELEMENT_INFO    elInfo,
                    int bitrateTot,
                    int averageBitsTot,
                    int staticBitsTot)
{
  int error=0;

  

   /* counting previous operation */

   
  switch(elInfo.nChannelsInEl) {

  case 1:
     
    elementBits->chBitrate = bitrateTot;

     
    elementBits->averageBits = (averageBitsTot - staticBitsTot);

    
    elementBits->maxBits = MAX_CHANNEL_BITS;

     
    elementBits->maxBitResBits = MAX_CHANNEL_BITS - averageBitsTot;

      
    elementBits->maxBitResBits -= (elementBits[0].maxBitResBits%8);

    
    elementBits->bitResLevel = elementBits[0].maxBitResBits;
    elementBits->relativeBits  = 1.0f;
    break;

  case 2:
     
    elementBits->chBitrate   = (int)(bitrateTot*0.5f);

     
    elementBits->averageBits = (averageBitsTot - staticBitsTot);

     
    elementBits->maxBits     = 2*MAX_CHANNEL_BITS;

      
    elementBits->maxBitResBits = 2*MAX_CHANNEL_BITS - averageBitsTot;

      
    elementBits->maxBitResBits -= (elementBits[0].maxBitResBits%8);

    
    elementBits->bitResLevel = elementBits[0].maxBitResBits;
    elementBits->relativeBits = 1.0f;
    break;

  default:
    
    error=1;
  }

  

  return error;
}

