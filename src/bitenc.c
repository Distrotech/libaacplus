/*
  Bitstream encoder
*/
#include <stdio.h>
#include <stdlib.h>

#include "aacplusenc.h"


 /* the 3GPP instrumenting tool */

#define globalGainOffset 100
#define icsReservedBit     0


/*****************************************************************************

    functionname: encodeSpectralData
    description:  encode spectral data
    returns:      none
    input:
    output:

*****************************************************************************/
static int encodeSpectralData(int                 *sfbOffset,
                              SECTION_DATA        *sectionData,
                              short               *quantSpectrum,
                              HANDLE_BIT_BUF      hBitStream)
{
  int i,sfb;
  int dbgVal;

  

    /* GetBitsAvail() */
  dbgVal = GetBitsAvail(hBitStream);

   /* sectionData->section[] */
   
  for(i=0;i<sectionData->noOfSections;i++)
  {
  /*
   huffencode spectral data for this section
   */
     /* sfbOffset[] */
    
    for(sfb=sectionData->section[i].sfbStart;
        sfb<sectionData->section[i].sfbStart+sectionData->section[i].sfbCnt;
        sfb++)
    {

        
      codeValues(quantSpectrum+sfbOffset[sfb],
                 sfbOffset[sfb+1]-sfbOffset[sfb],
                 sectionData->section[i].codeBook,
                 hBitStream);

    }
  }

    /* counting post-operation */

  

  return(GetBitsAvail(hBitStream)-dbgVal);

}

/*****************************************************************************

    functionname:encodeGlobalGain
    description: encodes Global Gain (common scale factor)
    returns:     none
    input:
    output:

*****************************************************************************/
static void encodeGlobalGain(int globalGain,
                             int logNorm,
                             int scalefac,
                             HANDLE_BIT_BUF hBitStream)
{
  

    
  WriteBits(hBitStream,globalGain - scalefac + globalGainOffset-4*logNorm,8);

  
}


/*****************************************************************************

    functionname:encodeIcsInfo
    description: encodes Ics Info
    returns:     none
    input:
    output:

*****************************************************************************/

static void encodeIcsInfo(int blockType,
                          int windowShape,
                          int groupingMask,
                          SECTION_DATA *sectionData,
                          HANDLE_BIT_BUF  hBitStream)
{
  

  
  WriteBits(hBitStream,icsReservedBit,1);

  
  WriteBits(hBitStream,blockType,2);

  
  WriteBits(hBitStream,windowShape,1);

  
  switch(blockType){
  case LONG_WINDOW:
  case START_WINDOW:
  case STOP_WINDOW:
     
    WriteBits(hBitStream,sectionData->maxSfbPerGroup,6);

    /* No predictor data present */
    
    WriteBits(hBitStream, 0, 1);
    break;

  case SHORT_WINDOW:
     
    WriteBits(hBitStream,sectionData->maxSfbPerGroup,4);

    /*
     Write grouping bits
    */
    
    WriteBits(hBitStream,groupingMask,TRANS_FAC-1);
    break;
  }

  
}

/*****************************************************************************

    functionname: encodeSectionData
    description:  encode section data (common Huffman codebooks for adjacent
                  SFB's)
    returns:      none
    input:
    output:

*****************************************************************************/
static int encodeSectionData(SECTION_DATA *sectionData,
                             HANDLE_BIT_BUF hBitStream)
{
  int sectEscapeVal=0,sectLenBits=0;
  int sectLen;
  int i;
  int dbgVal=GetBitsAvail(hBitStream);

  

    /* GetBitsAvail() */  /* counting previous operations */

   
  switch(sectionData->blockType)
  {
  case LONG_WINDOW:
  case START_WINDOW:
  case STOP_WINDOW:

    
    sectEscapeVal = SECT_ESC_VAL_LONG;
    sectLenBits   = SECT_BITS_LONG;
    break;

  case SHORT_WINDOW:

    
    sectEscapeVal = SECT_ESC_VAL_SHORT;
    sectLenBits   = SECT_BITS_SHORT;
    break;
  }

   /* sectionData->section[] */
   
  for(i=0;i<sectionData->noOfSections;i++)
  {
    
    WriteBits(hBitStream,sectionData->section[i].codeBook,4);

    
    sectLen = sectionData->section[i].sfbCnt;

    while(sectLen >= sectEscapeVal)
    {
       /* while() condition */

      
      WriteBits(hBitStream,sectEscapeVal,sectLenBits);

      
      sectLen-=sectEscapeVal;
    }

    
    WriteBits(hBitStream,sectLen,sectLenBits);
  }

    /* counting post-operation */

  

  return(GetBitsAvail(hBitStream)-dbgVal);
}

/*****************************************************************************

    functionname: encodeScaleFactorData
    description:  encode DPCM coded scale factors
    returns:      none
    input:
    output:

*****************************************************************************/
static int encodeScaleFactorData(unsigned short        *maxValueInSfb,
                                 SECTION_DATA        *sectionData,
                                 short                 *scalefac,
                                 HANDLE_BIT_BUF      hBitStream)
{
  int i,j,lastValScf,deltaScf;
  int dbgVal = GetBitsAvail(hBitStream);

  

    /* counting previous operation */

   
  lastValScf=scalefac[sectionData->firstScf];

   /* sectionData->section[] */
   
  for(i=0;i<sectionData->noOfSections;i++){

     
    if(sectionData->section[i].codeBook != CODE_BOOK_ZERO_NO){

       /* maxValueInSfb[]
                      scalefac[]
                   */
      
      for(j=sectionData->section[i].sfbStart;j<sectionData->section[i].sfbStart+sectionData->section[i].sfbCnt;j++){

        /*
          check if we can repeat the last value to save bits
        */
        
        if(maxValueInSfb[j] == 0)
        {
          
          deltaScf = 0;
        }
        else{
          
          deltaScf = -(scalefac[j]-lastValScf);

          
          lastValScf = scalefac[j];
        }
         
        if(codeScalefactorDelta(deltaScf,hBitStream)){
#ifdef DEBUG
          fprintf(stderr,"\nLAV for scf violated\n");
#endif

          
          return(1);
        }

      }
    }

  }

    /* counting post-operation */

  

  return(GetBitsAvail(hBitStream)-dbgVal);
}

/*****************************************************************************

    functionname:encodeMsInfo
    description: encodes MS-Stereo Info
    returns:     none
    input:
    output:

*****************************************************************************/
static void encodeMSInfo(int            sfbCnt,
                         int            grpSfb,
                         int            maxSfb,
                         int            msDigest,
                         int           *jsFlags,
                         HANDLE_BIT_BUF hBitStream)
{
  int sfb, sfbOff;

  

  
  switch(msDigest)
  {
  case MS_NONE:
    
    WriteBits(hBitStream,SI_MS_MASK_NONE,2);
    break;

  case MS_ALL:
    
    WriteBits(hBitStream,SI_MS_MASK_ALL,2);
    break;

  case MS_SOME:
    
    WriteBits(hBitStream,SI_MS_MASK_SOME,2);

    
    for(sfbOff = 0; sfbOff < sfbCnt; sfbOff+=grpSfb)
    {
      
      for(sfb=0; sfb<maxSfb; sfb++)
      {
         
        if(jsFlags[sfbOff+sfb] & MS_ON){

          
          WriteBits(hBitStream,1,1);

        }
        else{

          
          WriteBits(hBitStream,0,1);

        }
      }
    }

    break;

  }

  
}

/*****************************************************************************

    functionname: encodeTnsData
    description:  encode TNS data (filter order, coeffs, ..)
    returns:      none
    input:
    output:

*****************************************************************************/
static void encodeTnsData(TNS_INFO tnsInfo,
                          int blockType,
                          HANDLE_BIT_BUF hBitStream) 
{
  int i,k;
  int tnsPresent;
  int numOfWindows;
  int coefBits;

  

    
  numOfWindows=(blockType==2?TRANS_FAC:1);

  
  tnsPresent=0;

   /* tnsInfo.numOfFilters[] */
  
  for (i=0; i<numOfWindows; i++) {
    
     
    if (tnsInfo.tnsActive[i]==1) {

      
      tnsPresent=1;
    }
  }

  
  if (tnsPresent==0) {

    
    WriteBits(hBitStream,0,1);
  }
  else { /* there is data to be written*/
    
    
    WriteBits(hBitStream,1,1); /*data_present */

     /* tnsInfo.tnsActive[]
                    tnsInfo.coefRes[]
                 */
    
    for (i=0; i<numOfWindows; i++) {

        /* .. == .. ? */ 
      WriteBits(hBitStream,tnsInfo.tnsActive[i],(blockType==2?1:2));
      
      
      if (tnsInfo.tnsActive[i]) {
        
         /* tnsInfo.length[]
                        tnsInfo.order[] */

          /* .. == .. ? */ 
        WriteBits(hBitStream,(tnsInfo.coefRes[i]==4?1:0),1);
        
          /* .. == .. ? */ 
        WriteBits(hBitStream,tnsInfo.length[i],(blockType==2?4:6));
        
          /* .. == .. ? */ 
        WriteBits(hBitStream,tnsInfo.order[i],(blockType==2?3:5));
        
        
        if (tnsInfo.order[i]){
          
          
          WriteBits(hBitStream,FILTER_DIRECTION,1);
          
           
          if(tnsInfo.coefRes[i] == 4) {
            
            
            coefBits = 3;
            
             /* tnsInfo.coef[][] */
            
            for(k=0; k<tnsInfo.order[i]; k++) {
              
                
              if (tnsInfo.coef[i*TNS_MAX_ORDER_SHORT+k]> 3 ||
                  tnsInfo.coef[i*TNS_MAX_ORDER_SHORT+k]< -4) {
                
                
                coefBits = 4;
                break;
              }
            }
          }
          else {
            
            
            coefBits = 2;
            
             /* tnsInfo.coef[][] */
                           
            for(k=0; k<tnsInfo.order[i]; k++) {
              
                
              if (tnsInfo.coef[i*TNS_MAX_ORDER_SHORT+k]> 1 ||
                  tnsInfo.coef[i*TNS_MAX_ORDER_SHORT+k]< -2) {
                
                
                coefBits = 3;
                break;
              }
            }
          }
          
           
          WriteBits(hBitStream,-(coefBits - tnsInfo.coefRes[i]),1); /*coef_compres*/
          
           /* tnsInfo.coef[] */
           /* rmask[] */ 
          for (k=0; k<tnsInfo.order[i]; k++ ) {
            static const int rmask[] = {0,1,3,7,15};
            
             
            WriteBits(hBitStream,tnsInfo.coef[i*TNS_MAX_ORDER_SHORT+k] & rmask[coefBits],coefBits);
          }
        }
      }
    }
  }
  
  
}
/*****************************************************************************

    functionname: encodeGainControlData
    description:  unsupported
    returns:      none
    input:
    output:

*****************************************************************************/
static void encodeGainControlData(HANDLE_BIT_BUF hBitStream)
{
  

  
  WriteBits(hBitStream,0,1);

  
}

/*****************************************************************************

    functionname: encodePulseData
    description:  not supported yet (dummy)
    returns:      none
    input:
    output:

*****************************************************************************/
static void encodePulseData(HANDLE_BIT_BUF hBitStream)
{
  

  
  WriteBits(hBitStream,0,1);

  
}



/*****************************************************************************

    functionname: WriteIndividualChannelStream
    description:  management of write process of individual channel stream
    returns:      none
    input:
    output:

*****************************************************************************/
static void
writeIndividualChannelStream(int commonWindow,
                             int windowShape,
                             int groupingMask,
                             int *sfbOffset,
                             short scf[],
                             unsigned short *maxValueInSfb,
                             int globalGain,
                             short quantSpec[],
                             SECTION_DATA *sectionData,
                             HANDLE_BIT_BUF hBitStream,
                             TNS_INFO tnsInfo)
{
  int logNorm = -1;

  

   /* counting previous operation */

   
  encodeGlobalGain(globalGain, logNorm,scf[sectionData->firstScf], hBitStream);

  
  if(!commonWindow)
  {
     
    encodeIcsInfo(sectionData->blockType, windowShape, groupingMask, sectionData, hBitStream);
  }

     
  if(encodeSectionData(sectionData, hBitStream) != sectionData->sideInfoBits)
  {
#ifdef DEBUG
    fprintf(stderr,"ERROR writing sectionData\n");
#endif
  }

     
  if(encodeScaleFactorData( maxValueInSfb,
                            sectionData,
                            scf,
                            hBitStream) != sectionData->scalefacBits)
  {
#ifdef DEBUG
    fprintf(stderr,"ERROR writing scalefacData\n");
#endif
  }



  
  encodePulseData(hBitStream);

   
  encodeTnsData(tnsInfo, sectionData->blockType, hBitStream);

  
  encodeGainControlData(hBitStream);

     
  if(encodeSpectralData(sfbOffset,
                     sectionData,
                     quantSpec,
                     hBitStream) != sectionData->huffmanBits)
  {
#ifdef DEBUG
   fprintf(stderr,"ERROR writing spectralData\n");
#endif
  }


  
}

/*****************************************************************************

    functionname: writeSingleChannelElement
    description:  write single channel element to bitstream
    returns:      none
    input:
    output:

*****************************************************************************/
static int writeSingleChannelElement(int instanceTag,
                                     int *sfbOffset,
                                     QC_OUT_CHANNEL* qcOutChannel,
                                     HANDLE_BIT_BUF hBitStream,
                                     TNS_INFO tnsInfo)
{
  

  
  WriteBits(hBitStream,ID_SCE,3);

  
  WriteBits(hBitStream,instanceTag,4);

    
  writeIndividualChannelStream(0,
                               qcOutChannel->windowShape,
                               qcOutChannel->groupingMask,
                               sfbOffset,
                               qcOutChannel->scf,
                               qcOutChannel->maxValueInSfb,
                               qcOutChannel->globalGain,
                               qcOutChannel->quantSpec,
                               &(qcOutChannel->sectionData),
                               hBitStream,
                               tnsInfo
                               );

  
  return 0;
}



/*****************************************************************************

    functionname: writeChannelPairElement
    description:
    returns:      none
    input:
    output:

*****************************************************************************/
static int writeChannelPairElement(int instanceTag,
                                   int msDigest,
                                   int msFlags[MAX_GROUPED_SFB],
                                   int *sfbOffset[2],
                                   QC_OUT_CHANNEL qcOutChannel[2],
                                   HANDLE_BIT_BUF hBitStream,
                                   TNS_INFO tnsInfo[2]
                                   )
{
  

  
  WriteBits(hBitStream,ID_CPE,3);

  
  WriteBits(hBitStream,instanceTag,4);

  
  WriteBits(hBitStream,1,1); /* common window */

   
  encodeIcsInfo(qcOutChannel[0].sectionData.blockType,
                qcOutChannel[0].windowShape,
                qcOutChannel[0].groupingMask,
                &(qcOutChannel[0].sectionData),
                hBitStream);
  
  
  encodeMSInfo(qcOutChannel[0].sectionData.sfbCnt,
               qcOutChannel[0].sectionData.sfbPerGroup,
               qcOutChannel[0].sectionData.maxSfbPerGroup,
               msDigest,
               msFlags,
               hBitStream);

   
  writeIndividualChannelStream(1,
                               qcOutChannel[0].windowShape,
                               qcOutChannel[0].groupingMask,
                               sfbOffset[0],
                               qcOutChannel[0].scf,
                               qcOutChannel[0].maxValueInSfb,
                               qcOutChannel[0].globalGain,
                               qcOutChannel[0].quantSpec,
                               &(qcOutChannel[0].sectionData),
                               hBitStream,
                               tnsInfo[0]);
  
   
  writeIndividualChannelStream(1,
                               qcOutChannel[1].windowShape,
                               qcOutChannel[1].groupingMask,
                               sfbOffset[1],
                               qcOutChannel[1].scf,
                               qcOutChannel[1].maxValueInSfb,
                               qcOutChannel[1].globalGain,
                               qcOutChannel[1].quantSpec,
                               &(qcOutChannel[1].sectionData),
                               hBitStream,
                               tnsInfo[1]);

  

  return 0;
}



/*****************************************************************************

    functionname: writeFillElement
    description:  write fill elements to bitstream
    returns:      none
    input:
    output:

*****************************************************************************/
static void writeFillElement( const unsigned char *ancBytes,
                             int totFillBits,
                             HANDLE_BIT_BUF hBitStream)
{
  int i, cnt, esc_count;
  
  
  
  /*
    Write fill Element(s):
    amount of a fill element can be 7+X*8 Bits, X element of [0..270]
  */
  
  
  while(totFillBits>=(3+4)) {

       
    cnt = min((totFillBits-(3+4))/8, ((1<<4)-1));
    
    
    WriteBits(hBitStream,ID_FIL,3);
    
    
    WriteBits(hBitStream,cnt,4);
    
    
    totFillBits-=(3+4);
    
     
    if(cnt == (1<<4)-1){
      
         
      esc_count = min( (totFillBits/8)- ((1<<4)-1), (1<<8)-1);
      
      
      WriteBits(hBitStream,esc_count,8);
      
      
      totFillBits-=8;
      
      
      cnt+=esc_count-1;
    }
    
    
    for(i=0;i<cnt;i++){
      
       
      if(ancBytes)
        WriteBits(hBitStream,*ancBytes++,8);
      else
        WriteBits(hBitStream,0,8);
      
      
      totFillBits-=8;
    }
    
  }
  
}


/*****************************************************************************

    functionname: WriteBitStreamData
    description:  main function of write process
    returns:
    input:
    output:

*****************************************************************************/
int WriteBitstreamData (HANDLE_BIT_BUF hBitStream,
                    ELEMENT_INFO elInfo,
                    QC_OUT *qcOut,
                    PSY_OUT* psyOut,
                    int* globUsedBits,
                    const unsigned char *ancBytes
                    ) /* returns error code */
{
  int bitMarkUp, elementUsedBits, frameBits;

  

    /* GetBitsAvail() */
  bitMarkUp = GetBitsAvail(hBitStream);

  
  *globUsedBits = 0;

  {
    int* sfbOffset[2];
    TNS_INFO tnsInfo[2];
    
    
    elementUsedBits =0;

     
    switch (elInfo.elType) {

    case ID_SCE:      /* single channel */
      
       
      sfbOffset[0] =
        psyOut->psyOutChannel[elInfo.ChannelIndex[0]].sfbOffsets;
      tnsInfo[0]=psyOut->psyOutChannel[elInfo.ChannelIndex[0]].tnsInfo;
      
       
      writeSingleChannelElement(  elInfo.instanceTag,
                                  sfbOffset[0],
                                  &qcOut->qcChannel[elInfo.ChannelIndex[0]],
                                  hBitStream,
                                  tnsInfo[0] );
      break;
      
    case ID_CPE:     /* channel pair */
      {
        int msDigest = psyOut->psyOutElement.toolsInfo.msDigest;
        int* msFlags = psyOut->psyOutElement.toolsInfo.msMask;
        
          /* counting previous operation */
        
         
        sfbOffset[0] = psyOut->psyOutChannel[elInfo.ChannelIndex[0]].sfbOffsets;
        sfbOffset[1] = psyOut->psyOutChannel[elInfo.ChannelIndex[1]].sfbOffsets;
        
        
        tnsInfo[0] = psyOut->psyOutChannel[elInfo.ChannelIndex[0]].tnsInfo;
        tnsInfo[1] = psyOut->psyOutChannel[elInfo.ChannelIndex[1]].tnsInfo;
        
          
        writeChannelPairElement(elInfo.instanceTag,
                                msDigest,
                                msFlags,
                                sfbOffset,
                                &qcOut->qcChannel[elInfo.ChannelIndex[0]],
                                hBitStream,
                                tnsInfo);
      }
      break;
      
    default:
      
      return 1;

    }   /* switch */

     
    elementUsedBits -= bitMarkUp;

      /* GetBitsAvail() */
    bitMarkUp = GetBitsAvail(hBitStream);

     
    frameBits = elementUsedBits + bitMarkUp;
    
  }

  /*
    ancillary data
  */
   
  writeFillElement(ancBytes,
                   qcOut->totAncBitsUsed, 
                   hBitStream);
  

   
  writeFillElement(NULL,
                   qcOut->totFillBits, 
                   hBitStream);

  
  WriteBits(hBitStream,ID_END,3);

  /* byte alignement */

     
  WriteBits(hBitStream,0,(8-(hBitStream->cntBits % 8)) % 8);
  
   
  *globUsedBits -= bitMarkUp;

    /* GetBitsAvail() */
  bitMarkUp = GetBitsAvail(hBitStream);

   
  *globUsedBits+=bitMarkUp;

  
  frameBits+= *globUsedBits;

    
  if(frameBits != qcOut->totStaticBitsUsed+qcOut->totDynBitsUsed + qcOut->totAncBitsUsed + 
     + qcOut->totFillBits+qcOut->alignBits){
#ifdef DEBUG
    fprintf(stderr,"\nNo of written frame bits != expected Bits !!!\n");
#endif
    
    return -1;
  }

  
  return 0;
}
