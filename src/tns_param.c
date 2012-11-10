/*
  Temporal Noise Shaping parameters
*/
#include <stdio.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


/*****************************************************************************

    functionname: GetTnsParam
    description:  Get threshold calculation parameter set that best matches
                  the bit rate
    returns:      the parameter set
    input:        blockType, bitRate
    output:       the parameter set

*****************************************************************************/
int GetTnsParam(TNS_CONFIG_TABULATED *tnsConfigTab, int bitRate, int channels, int blockType) {

  unsigned int i;

  

   /* tnsInfoTab[] */
  
   
  if (tnsConfigTab == NULL)
    return 1;

  
  tnsConfigTab->threshOn = -1.0f;

  
  for(i = 0; i < sizeof(tnsInfoTab)/sizeof(TNS_INFO_TAB); i++) {

      
    if((bitRate >= tnsInfoTab[i].bitRateFrom) && bitRate <= tnsInfoTab[i].bitRateTo) {

      
      switch(blockType) {

      case LONG_WINDOW :
        
        switch(channels) {
        case 1 :
          
          *tnsConfigTab=*tnsInfoTab[i].paramMono_Long;
          break;
        case 2 :
          
          *tnsConfigTab=*tnsInfoTab[i].paramStereo_Long;
          break;
        }
      break;

      case SHORT_WINDOW :
        
        switch(channels) {
        case 1 :
          
          
          *tnsConfigTab=*tnsInfoTab[i].paramMono_Short;
          break;
        case 2 :
          
          *tnsConfigTab=*tnsInfoTab[i].paramStereo_Short;
          break;
        }

        break;
      }
    }
  }

   
  if (tnsConfigTab->threshOn == -1.0f) {
    return 1;
  }

  
  return 0;
}

/*****************************************************************************

    functionname: GetTnsMaxBands
    description:  sets tnsMaxSfbLong, tnsMaxSfbShort according to sampling rate
    returns:
    input:        samplingRate, profile, granuleLen
    output:       tnsMaxSfbLong, tnsMaxSfbShort

*****************************************************************************/
void GetTnsMaxBands(int samplingRate, int blockType, int* tnsMaxSfb){
  unsigned int i;

  

   
  *tnsMaxSfb=-1;

   /* tnsMaxBandsTab[] */
  
  for(i=0;i<sizeof(tnsMaxBandsTab)/sizeof(TNS_MAX_TAB_ENTRY);i++){

     
    if(samplingRate == tnsMaxBandsTab[i].samplingRate){

          
        *tnsMaxSfb=(blockType==2)?tnsMaxBandsTab[i].maxBandShort:
                                  tnsMaxBandsTab[i].maxBandLong;
        break;
    }
  }

  
}
