/*
  Static bit counter $Revision
*/
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static int countMsMaskBits(int   sfbCnt,
                           int   sfbPerGroup,
                           int   maxSfbPerGroup,
                           struct TOOLSINFO *toolsInfo)
{

  int msBits=0,sfbOff,sfb;

  

   /* counting previous operation */

   
  switch(toolsInfo->msDigest)
  {
    case MS_NONE:
    case MS_ALL:
    break;

    case MS_SOME:

      
      for(sfbOff = 0; sfbOff < sfbCnt; sfbOff+=sfbPerGroup)
      {
        
        for(sfb=0; sfb<maxSfbPerGroup; sfb++)
        {
        
        msBits++;
        }
      }
    break;
  }

  

  return(msBits);

}


static int tnsCount(TNS_INFO *tnsInfo, int blockType)
{

  int i,k;
  int tnsPresent;
  int numOfWindows;
  int count;
  int coefBits;

  

  
  count = 0;

    
  numOfWindows=(blockType==2?8:1);

  
  tnsPresent=0;

   /* tnsInfo->numOfFilters[] */
  
  for (i=0; i<numOfWindows; i++) {

     
    if (tnsInfo->tnsActive[i]==1) {

      
      tnsPresent=1;
    }
  }

  
  if (tnsPresent==0) {
    /* count+=1; */
  }
  else{ /* there is data to be written*/
    /*count += 1; */

     /* tnsInfo->tnsActive[]
                    tnsInfo->coefRes[]
                 */
    
    for (i=0; i<numOfWindows; i++) {

        /* .. == .. ? */ 
      count +=(blockType==SHORT_WINDOW?1:2);

      
      if (tnsInfo->tnsActive[i]) {

        
        count += 1;

          /* .. == .. ? */ 
        count +=(blockType==SHORT_WINDOW?4:6);
        count +=(blockType==SHORT_WINDOW?3:5);
        
        
        if (tnsInfo->order[i]){

          
          count +=1; /*direction*/
          count +=1; /*coef_compression */
          
           
          if(tnsInfo->coefRes[i] == 4) {

            
            coefBits=3;
            
             /* tnsInfo->coef[]*/
            
            for(k=0; k<tnsInfo->order[i]; k++) {
              
                
              if (tnsInfo->coef[i*TNS_MAX_ORDER_SHORT+k]> 3 ||
                  tnsInfo->coef[i*TNS_MAX_ORDER_SHORT+k]< -4) {
                
                
                coefBits = 4;
                break;
              }
            }
          }
          else {
            
            
            coefBits = 2;
            
             /* tnsInfo->coef[]*/
            
            for(k=0; k<tnsInfo->order[i]; k++) {
              
                
              if (tnsInfo->coef[i*TNS_MAX_ORDER_SHORT+k]> 1 ||
                  tnsInfo->coef[i*TNS_MAX_ORDER_SHORT+k]< -2) {
                
                
                coefBits = 3;
                break;
              }
            }
          }
          
          
          for (k=0; k<tnsInfo->order[i]; k++ ) {
            
            
            count +=coefBits;
          }
        }
      }
    }
  }

  

  return count;
}


static int countTnsBits(TNS_INFO *tnsInfo,int blockType)
{
  
  
  

  return(tnsCount(tnsInfo, blockType));
}





int countStaticBitdemand(PSY_OUT_CHANNEL psyOutChannel[MAX_CHANNELS],
                         PSY_OUT_ELEMENT *psyOutElement,
                         int channels)
{

  int statBits=0;
  int ch;

  

   /* counting previous operation */

  
  switch(channels) {

  case 1:

    
    statBits+=SI_ID_BITS+SI_SCE_BITS+SI_ICS_BITS;

      
    statBits+=countTnsBits(&(psyOutChannel[0].tnsInfo),
                           psyOutChannel[0].windowSequence);

    
    switch(psyOutChannel[0].windowSequence){
      case LONG_WINDOW:
      case START_WINDOW:
      case STOP_WINDOW:

        
        statBits+=SI_ICS_INFO_BITS_LONG;
      break;
      case SHORT_WINDOW:

      
      statBits+=SI_ICS_INFO_BITS_SHORT;
      break;
    }
    break;

  case 2:
    
    statBits+=SI_ID_BITS+SI_CPE_BITS+2*SI_ICS_BITS;
    
    
    statBits+=SI_CPE_MS_MASK_BITS;
    
       
    statBits+=countMsMaskBits(psyOutChannel[0].sfbCnt,
                              psyOutChannel[0].sfbPerGroup,
                              psyOutChannel[0].maxSfbPerGroup,
                              &psyOutElement->toolsInfo);
    
     /* psyOutChannel[] */
    
    switch(psyOutChannel[0].windowSequence) {
      
    case LONG_WINDOW:
    case START_WINDOW:
    case STOP_WINDOW:
      
      
      statBits+=SI_ICS_INFO_BITS_LONG;
      break;
      
    case SHORT_WINDOW:
      
      
      statBits+=SI_ICS_INFO_BITS_SHORT;
      break;
    }
    
     /* psyOutChannel[ch] */
    
    for(ch=0;ch<2;ch++) {
      
        
      statBits+=countTnsBits(&(psyOutChannel[ch].tnsInfo),
                             psyOutChannel[ch].windowSequence);
    }
    
    break;
  }

  
  
  return statBits;
}

