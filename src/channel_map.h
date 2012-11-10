/*
  Channel Mapping
*/
#ifndef _CHANNEL_MAP_H
#define _CHANNEL_MAP_H

#include "aacplusenc.h"


int InitElementInfo (int nChannels,
                     ELEMENT_INFO *elInfo);

int InitElementBits(ELEMENT_BITS*   elementBits,
                    ELEMENT_INFO elInfo,
                    int bitrateTot,
                    int averageBitsTot,
                    int staticBitsTot);

#endif /* CHANNEL_MAP_H */
