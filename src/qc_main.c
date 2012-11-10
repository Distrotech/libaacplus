/*
  Quantizing & coding main module
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aacplusenc.h"

#include <math.h>

 /* the 3GPP instrumenting tool */

typedef enum{
    FRAME_LEN_BYTES_MODULO =  1,
    FRAME_LEN_BYTES_INT    =  2
}FRAME_LEN_RESULT_MODE;

static const int maxFillElemBits = 7 + 270*8;

/* forward declarations */

static int calcMaxValueInSfb(int sfbCnt,
                             int maxSfbPerGroup,
                             int sfbPerGroup,
                             int sfbOffset[MAX_GROUPED_SFB],
                             short quantSpectrum[FRAME_LEN_LONG],
                             unsigned short maxValue[MAX_GROUPED_SFB]);


/*****************************************************************************

    functionname: calcFrameLen
    description:
    returns:
    input:
    output:

*****************************************************************************/
static int calcFrameLen(int bitRate,
                        int sampleRate,
                        FRAME_LEN_RESULT_MODE mode)
{

   int result;

   

   
   result = ((FRAME_LEN_LONG)>>3)*(bitRate);

   
   switch(mode) {
     case FRAME_LEN_BYTES_MODULO:
         
         result %= sampleRate;
     break;
     case FRAME_LEN_BYTES_INT:
         
         result /= sampleRate;
     break;
   }

   

   return(result);
}

/*****************************************************************************

    functionname:framePadding
    description: Calculates if padding is needed for actual frame
    returns:
    input:
    output:

*****************************************************************************/
static int framePadding(int bitRate,
                        int sampleRate,
                        int *paddingRest)
{
  int paddingOn;
  int difference;

  

  
  paddingOn = 0;

  
  difference = calcFrameLen( bitRate,
                             sampleRate,
                             FRAME_LEN_BYTES_MODULO );

   
  *paddingRest-=difference;

  
  if (*paddingRest <= 0 ) {

    
    paddingOn = 1;

     
    *paddingRest += sampleRate;
  }

  

  return( paddingOn );
}


/*********************************************************************************

         functionname: QCOutNew
         description:
         return:

**********************************************************************************/

int QCOutNew(AACRam_t *aacram, QC_OUT *hQC, int nChannels)
{
  int error=0;
  int i;

  

   /* counting previous operation */

   /* hQC->qcChannel[]
                  quantSpec[]
                  maxValueInSfb[]
                  scf[]
               */
  
  for (i=0; i<nChannels; i++) {
    
    hQC->qcChannel[i].quantSpec = &aacram->quantSpec[i*FRAME_LEN_LONG];

    
    hQC->qcChannel[i].maxValueInSfb = &aacram->maxValueInSfb[i*MAX_GROUPED_SFB];

    
    hQC->qcChannel[i].scf = &aacram->scf[i*MAX_GROUPED_SFB];
  }
 
  
  
  if (error){
    
    QCOutDelete(hQC);

    
    hQC = 0;
    }

   /* counting post-operation */

  

  return (hQC == 0);
}


/*********************************************************************************

         functionname: QCOutDelete
         description:
         return:

**********************************************************************************/
void QCOutDelete(QC_OUT* hQC)
{
  

  /* 
    nothing to do
  */
 
  
  hQC=NULL;

  
}

/*********************************************************************************

         functionname: QCNew
         description:
         return:

**********************************************************************************/
int QCNew(QC_STATE *hQC)
{
  
  

      
  memset(hQC,0,sizeof(QC_STATE));

  

  return (0);
}

/*********************************************************************************

         functionname: QCDelete
         description:
         return:

**********************************************************************************/
void QCDelete(QC_STATE *hQC)
{
 
  

  /* 
    nothing to do
  */
  
  hQC=NULL;

  
}

/*********************************************************************************

         functionname: QCInit
         description:
         return:

**********************************************************************************/
int QCInit(AACRam_t *aacram, QC_STATE *hQC,
           struct QC_INIT *init)
{
  

   
  hQC->nChannels       = init->elInfo->nChannelsInEl;
  hQC->maxBitsTot      = init->maxBits;

    
  hQC->bitResTot       = init->bitRes - init->averageBits;

   
  hQC->averageBitsTot  = init->averageBits;
  hQC->maxBitFac       = init->maxBitFac;
  hQC->padding.paddingRest = init->padding.paddingRest;

   
  hQC->globStatBits    = 3;                                  /* for ID_END */

   
  InitElementBits(&hQC->elementBits,
                  *init->elInfo,
                  init->bitrate,
                  init->averageBits,
                  hQC->globStatBits);

   
  AdjThrInit(&hQC->adjThr,
             init->meanPe,
             hQC->elementBits.chBitrate);

  
  BCInit(aacram);

  

  return 0;
}


/*********************************************************************************

         functionname: QCMain
         description:
         return:

**********************************************************************************/
int QCMain(AACRam_t *aacram,
           QC_STATE* hQC,
           int nChannels,
           ELEMENT_BITS* elBits,
           ATS_ELEMENT* adjThrStateElement,
           PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],  /* may be modified in-place */
           PSY_OUT_ELEMENT* psyOutElement,
           QC_OUT_CHANNEL  qcOutChannel[MAX_CHANNELS],    /* out                      */
           QC_OUT_ELEMENT* qcOutElement,
           int ancillaryDataBytes)      
{
  int ch;
  float sfbFormFactor[MAX_CHANNELS][MAX_GROUPED_SFB];
  float sfbNRelevantLines[MAX_CHANNELS][MAX_GROUPED_SFB];
  int maxChDynBits[MAX_CHANNELS];
  float chBitDistribution[MAX_CHANNELS];  
  int loopCnt=0;

  

   
  if (elBits->bitResLevel < 0) {
#ifdef DEBUG
     fprintf(stderr, "\ntoo little bits in bitres\n");
#endif
     return -1;
  }

    
  if (elBits->bitResLevel > elBits->maxBitResBits) {
#ifdef DEBUG
     fprintf(stderr, "\ntoo many bits in bitres\n");
#endif
     return -1;
  }

    
  qcOutElement->staticBitsUsed = countStaticBitdemand(psyOutChannel,
                                                      psyOutElement,
                                                      nChannels);

  
  
  if (ancillaryDataBytes) {

        
    qcOutElement->ancBitsUsed = 7+8*(ancillaryDataBytes + (ancillaryDataBytes >=15));
  } else {

     
    qcOutElement->ancBitsUsed = 0;
  }

  
  CalcFormFactor(sfbFormFactor,sfbNRelevantLines, psyOutChannel, nChannels);

     
  AdjustThresholds(&hQC->adjThr,
                  adjThrStateElement,
                  psyOutChannel,
                  psyOutElement,
                  chBitDistribution,
                  sfbFormFactor,
                  nChannels,
                  qcOutElement,
                  elBits->averageBits-qcOutElement->staticBitsUsed - qcOutElement->ancBitsUsed,
                  elBits->bitResLevel,
                  elBits->maxBits,
                  hQC->maxBitFac,
                  qcOutElement->staticBitsUsed+qcOutElement->ancBitsUsed);

  
  EstimateScaleFactors(aacram,
                       psyOutChannel,
                       qcOutChannel,
                       sfbFormFactor,
                       sfbNRelevantLines,
                       nChannels);


  /* condition to prevent empty bitreservoir */
   /* pointers for maxChDynBits[],
                               chBitDistribution[]
               */
  
  for (ch = 0; ch < nChannels; ch++) {

          
      maxChDynBits[ch] = (int)floor(chBitDistribution[ch] * (float)
                                    (elBits->averageBits + elBits->bitResLevel - 7 -
                                     qcOutElement->staticBitsUsed - qcOutElement->ancBitsUsed));
                                 
      /* -7 bec. of align bits */

  }

   
  qcOutElement->dynBitsUsed = 0;

   /* pointers for psyOutChannel[],
                               qcOutChannel[],
                               maxChDynBits[]
               */
  
  for (ch = 0; ch < nChannels; ch++)
  {
    /* now loop until bitstream constraints (chDynBits < maxChDynBits)
       are fulfilled */
    int chDynBits;
    int constraintsFulfilled;
    int iter=0;

     /* counting previous operation */

    do
    {
        
        constraintsFulfilled = 1;

        
        if (iter>0) {

          
          QuantizeSpectrum(psyOutChannel[ch].sfbCnt,
                           psyOutChannel[ch].maxSfbPerGroup,
                           psyOutChannel[ch].sfbPerGroup,
                           psyOutChannel[ch].sfbOffsets,
                           psyOutChannel[ch].mdctSpectrum,
                           qcOutChannel[ch].globalGain,
                           qcOutChannel[ch].scf,
                           qcOutChannel[ch].quantSpec);
        }

          
        if (calcMaxValueInSfb(psyOutChannel[ch].sfbCnt,
                                psyOutChannel[ch].maxSfbPerGroup,
                                psyOutChannel[ch].sfbPerGroup,
                                psyOutChannel[ch].sfbOffsets,
                                qcOutChannel[ch].quantSpec,
                                qcOutChannel[ch].maxValueInSfb) > MAX_QUANT)
        {
#ifdef DEBUG
            fprintf(stderr,"\nMAXQUANT violated\n");
#endif
            
            constraintsFulfilled=0;
        }

         
        chDynBits =
            dynBitCount(aacram,
                        qcOutChannel[ch].quantSpec,
                        qcOutChannel[ch].maxValueInSfb,
                        qcOutChannel[ch].scf,
                        psyOutChannel[ch].windowSequence,
                        psyOutChannel[ch].sfbCnt,
                        psyOutChannel[ch].maxSfbPerGroup,
                        psyOutChannel[ch].sfbPerGroup,
                        psyOutChannel[ch].sfbOffsets,
                        &qcOutChannel[ch].sectionData);

         
        if (chDynBits >= maxChDynBits[ch])
        {
            
            constraintsFulfilled = 0;
#ifdef DEBUG
            fprintf(stderr, "\nWARNING: too many dynBits, extra quantization necessary\n");
#endif
        }

        
        if (!constraintsFulfilled)
        {
           
          qcOutChannel[ch].globalGain++;
        }

        
        iter++;
        loopCnt+=100;


    } while(!constraintsFulfilled);

      
    qcOutElement->dynBitsUsed += chDynBits;
    

    
    qcOutChannel[ch].groupingMask = psyOutChannel[ch].groupingMask;
    qcOutChannel[ch].windowShape  = psyOutChannel[ch].windowShape;
  }

  /* save dynBitsUsed for correction of bits2pe relation */
   
  AdjThrUpdate(adjThrStateElement, qcOutElement->dynBitsUsed);

  {
    int bitResSpace = elBits->maxBitResBits - elBits->bitResLevel;
    int deltaBitRes = elBits->averageBits - (qcOutElement->staticBitsUsed
                                             + qcOutElement->dynBitsUsed + qcOutElement->ancBitsUsed);

      /* counting previous operations */

       
    qcOutElement->fillBits = max(0,(deltaBitRes - bitResSpace));
  }

  

  return 0; /* OK */
}


/*********************************************************************************

         functionname: calcMaxValueInSfb
         description:
         return:

**********************************************************************************/

static int calcMaxValueInSfb(int sfbCnt,
                             int maxSfbPerGroup,
                             int sfbPerGroup,
                             int sfbOffset[MAX_GROUPED_SFB],
                             short quantSpectrum[FRAME_LEN_LONG],
                             unsigned short maxValue[MAX_GROUPED_SFB])
{
  int sfbOffs,sfb;
  int maxValueAll = 0;

  

   /* counting previous operation */

  
  for(sfbOffs=0;sfbOffs<sfbCnt;sfbOffs+=sfbPerGroup)
  {

   /* pointers for sfbOffset[],
                               maxValue[]
               */
  
  for (sfb = 0; sfb < maxSfbPerGroup; sfb++)
  {
    int line;
    int maxThisSfb = 0;

     /* counting previous operation */
 
     /* pointers for quantSpectrum[] */
    
    for (line = sfbOffset[sfbOffs+sfb]; line < sfbOffset[sfbOffs+sfb+1]; line++)
    {
        
      if (abs(quantSpectrum[line]) > maxThisSfb)
      {
        
        maxThisSfb = abs(quantSpectrum[line]);
      }
    }

    
    maxValue[sfbOffs+sfb] = maxThisSfb;

     
    if (maxThisSfb > maxValueAll)
    {
      
      maxValueAll = maxThisSfb;
    }
  }
  }

  

  return maxValueAll;
}


/*********************************************************************************

         functionname: updateBitres
         description:
         return:

**********************************************************************************/
void UpdateBitres(QC_STATE* qcKernel,
                  QC_OUT*   qcOut)
                  
{
  ELEMENT_BITS* elBits;

  
 
   
  qcKernel->bitResTot=0;

   
  elBits = &qcKernel->elementBits;
  
   
  if (elBits->averageBits > 0) {
    
    /* constant bitrate */
      
    elBits->bitResLevel += 
      elBits->averageBits - (qcOut->qcElement.staticBitsUsed
                             + qcOut->qcElement.dynBitsUsed
                             + qcOut->qcElement.ancBitsUsed
                             + qcOut->qcElement.fillBits);
    
      
    qcKernel->bitResTot += elBits->bitResLevel;
  }
  else {
    
    /* variable bitrate */
     
    elBits->bitResLevel = elBits->maxBits;
    qcKernel->bitResTot = qcKernel->maxBitsTot;
  }

  
}

/*********************************************************************************

         functionname: FinalizeBitConsumption
         description:
         return:

**********************************************************************************/
int FinalizeBitConsumption( QC_STATE *qcKernel,
                            QC_OUT* qcOut)
{
  int nFullFillElem, diffBits;
  int totFillBits = 0;

  

   /* prev. instructions */

   
  qcOut->totStaticBitsUsed = qcKernel->globStatBits;
  qcOut->totDynBitsUsed = 0;
  qcOut->totAncBitsUsed = 0;
  qcOut->totFillBits=0;

   /* qcOut->qcElement[] */
  {
      
    qcOut->totStaticBitsUsed += qcOut->qcElement.staticBitsUsed;
    qcOut->totDynBitsUsed    += qcOut->qcElement.dynBitsUsed;
    qcOut->totAncBitsUsed    += qcOut->qcElement.ancBitsUsed;
    qcOut->totFillBits       += qcOut->qcElement.fillBits;

    
    if (qcOut->qcElement.fillBits) {

       
      totFillBits += qcOut->qcElement.fillBits;
    }
  }

    
  nFullFillElem = (qcOut->totFillBits-1) / maxFillElemBits;

  
  if (nFullFillElem) {

       
    qcOut->totFillBits -= nFullFillElem * maxFillElemBits;
  }

  /* check fill elements */
   
  if (qcOut->totFillBits > 0) {

    /* minimum Fillelement contains 7 (TAG + byte cnt) bits */
      
    qcOut->totFillBits = max(7, qcOut->totFillBits);

    /* fill element size equals n*8 + 7 */
      /* .. % 8 */ 
    qcOut->totFillBits += ((8-(qcOut->totFillBits-7)%8)%8);
  }

   
  qcOut->totFillBits += nFullFillElem * maxFillElemBits;

  /* now distribute extra fillbits and alignbits over channel elements */
     /* .. % 8 */
  qcOut->alignBits = 7 - (qcOut->totDynBitsUsed + qcOut->totStaticBitsUsed + qcOut->totAncBitsUsed + 
                          + qcOut->totFillBits -1)%8;

     
  if( ((qcOut->alignBits + qcOut->totFillBits - totFillBits)==8) && (qcOut->totFillBits>8) ) {
      
    qcOut->totFillBits -= 8;
  }

    
  diffBits = (qcOut->alignBits + qcOut->totFillBits) - totFillBits;
  
  
  if (diffBits) {

    
    if(diffBits<0) {
#ifdef DEBUG
      fprintf(stderr,"\nbuffer fullness underflow \n");
#endif
      
      return -1;
    }
    else {

        
      qcOut->qcElement.fillBits += diffBits;

    }
  }

  
    
  if ((qcOut->totDynBitsUsed + qcOut->totStaticBitsUsed + qcOut->totAncBitsUsed + 
       qcOut->totFillBits +  qcOut->alignBits) > qcKernel->maxBitsTot) {
#ifdef DEBUG
    fprintf(stderr,"\ntoo many bits used in Q&C\n");
#endif

    
    return -1;
  }

  
  return 0;
}


/*********************************************************************************

         functionname: AdjustBitrate
         description:  adjusts framelength via padding on a frame to frame basis,
                       to achieve a bitrate that demands a non byte aligned
                       framelength
         return:       errorcode

**********************************************************************************/
int AdjustBitrate(QC_STATE* hQC,
                  int       bitRate,    /* total bitrate */
                  int       sampleRate  /* output sampling rate */
                  )
{
  int paddingOn;
  int frameLen;
  int codeBits;
  int codeBitsLast;

  

  /* Do we need a extra padding byte? */
    
  paddingOn = framePadding(bitRate,
                           sampleRate,
                           &hQC->padding.paddingRest);

   
  frameLen = paddingOn + calcFrameLen(bitRate,
                                      sampleRate,
                                      FRAME_LEN_BYTES_INT);

  
  frameLen <<= 3;

   
  codeBitsLast = hQC->averageBitsTot - hQC->globStatBits;

   
  codeBits     = frameLen - hQC->globStatBits;


  /* calculate bits for every channel element */
   
  if (codeBits != codeBitsLast) {
    int totalBits = 0;

     /* counting previous operation */

     
    hQC->elementBits.averageBits =  (int)(hQC->elementBits.relativeBits * codeBits);
    
    
    totalBits += hQC->elementBits.averageBits;

     
    hQC->elementBits.averageBits += codeBits - totalBits;
  }

   
  hQC->averageBitsTot = frameLen;

  

  return 0;
}

