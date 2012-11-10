/*
  fast aac coder interface library functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */




/*-----------------------------------------------------------------------------

     functionname: AacInitDefaultConfig
     description:  gives reasonable default configuration
     returns:      ---

 ------------------------------------------------------------------------------*/
void AacInitDefaultConfig(aacplusEncConfiguration *config)
{
  

  /* make the pre initialization of the structs flexible */
      
  memset(config, 0, sizeof(aacplusEncConfiguration));

  /* default configurations */
   
  config->bitRate         = 48000;
  config->bandWidth       = 0;
  config->nSamplesPerFrame = AACENC_BLOCKSIZE;

  
}






/*---------------------------------------------------------------------------

    functionname: AacEncOpen
    description:  allocate and initialize a new encoder instance
    returns:      0 if success

  ---------------------------------------------------------------------------*/

int  
AacEncOpen ( struct AAC_ENCODER      *hAacEnc,       /* pointer to an encoder handle, initialized on return */
             aacplusEncConfiguration     *config          /* pre-initialized config struct */
             )
{
  int error = 0;
  int profile = 1;
  ELEMENT_INFO* elInfo = NULL;
  
   /* counting previous operations */
  
  
  if (!error) {
    /* sanity checks on config structure */
       
    error = (&config == 0            || hAacEnc == 0             ||
             config->nChannelsIn  < 1 || config->nChannelsIn  > MAX_CHANNELS ||
             config->nChannelsOut < 1 || config->nChannelsOut > MAX_CHANNELS ||
             config->nChannelsIn  < config->nChannelsOut            ||
             (config->bitRate!=0 && (config->bitRate / config->nChannelsOut < 8000      ||
                                    config->bitRate / config->nChannelsOut > 160000)));
  }
  
  /* check sample rate */
  
  if (!error)  {

    
    switch (config->sampleRate) {
    case  8000: case 11025: case 12000:
    case 16000: case 22050: case 24000:
    case 32000: case 44100: case 48000:
      break;

    default:
      
      error = 1; break;
    }
  }

  
  /* check if bit rate is not too high for sample rate */
  
  if (!error) {

      
    if (config->bitRate > ((float)(MAX_CHANNEL_BITS-744)/FRAME_LEN_LONG*
                          config->sampleRate*config->nChannelsOut))
      {
        
        error=1;
      }
  }
  
  
  
  if (!error) {
     
    hAacEnc->config = config;
  }
  
  
  
  if (!error) {
      
    error = InitElementInfo (config->nChannelsOut,
                             &hAacEnc->elInfo);
  }
  
  
  if (!error) {
     
    elInfo = &hAacEnc->elInfo;
  }
  
  /* allocate the Psy aud Psy Out structure */
  
  if (!error) {

        
    error = (PsyNew(hAacEnc->aac_ram, hAacEnc->sbr_ram, &hAacEnc->psyKernel, elInfo->nChannelsInEl) ||
             PsyOutNew(&hAacEnc->psyOut));
  }
  
  
  if (!error) {
    int tnsMask=3;
    
     /* counting previous operation */
    
     
    hAacEnc->bandwidth90dB = (int)hAacEnc->config->bandWidth;
    
      
    error = psyMainInit(&hAacEnc->psyKernel,
                        config->sampleRate,
                        config->bitRate,
                        elInfo->nChannelsInEl,
                        tnsMask,
                        hAacEnc->bandwidth90dB);
  }
  
  
  /* allocate the Q&C Out structure */
  
  if (!error) {
      
    error = QCOutNew(hAacEnc->aac_ram, &hAacEnc->qcOut,
                     elInfo->nChannelsInEl);
  }
  
  /* allocate the Q&C kernel */
  
  if (!error) {
      
    error = QCNew(&hAacEnc->qcKernel);
  }
  
  
  if (!error) {
    struct QC_INIT qcInit;
    
      
    qcInit.elInfo = &hAacEnc->elInfo;
    
      
    qcInit.maxBits = MAX_CHANNEL_BITS * elInfo->nChannelsInEl;
    
    
    qcInit.bitRes = qcInit.maxBits;
    
      
    qcInit.averageBits = (config->bitRate * FRAME_LEN_LONG) / config->sampleRate;
    
    
    qcInit.padding.paddingRest = config->sampleRate;
    
      
    qcInit.meanPe = 10.0f * FRAME_LEN_LONG * hAacEnc->bandwidth90dB/(config->sampleRate/2.0f);
    
       
    qcInit.maxBitFac = (float)((MAX_CHANNEL_BITS-744)*elInfo->nChannelsInEl) /
      (float)(qcInit.averageBits?qcInit.averageBits:1);
    
    
    qcInit.bitrate = config->bitRate;
    
      
    error = QCInit(hAacEnc->aac_ram, &hAacEnc->qcKernel, &qcInit);
  }
  
  /* init bitstream encoder */
  
  if (!error) {

     
      hAacEnc->bseInit.nChannels   = elInfo->nChannelsInEl;
      hAacEnc->bseInit.bitrate     = config->bitRate;
      hAacEnc->bseInit.sampleRate  = config->sampleRate;
      hAacEnc->bseInit.profile     = profile;
  }
  
  
  /* common things */
  
  if (!error) {
    
       
    hAacEnc->downmix = (config->nChannelsIn==2 && config->nChannelsOut==1);
    
      
    hAacEnc->downmixFac = (hAacEnc->downmix) ? config->nChannelsIn : 1;
  }
  
  
  
  
  if(!error) {

    /*
      decide if stereo preprocessing should be activated 
    */
    
        
    if ( elInfo->elType == ID_CPE &&
         (config->sampleRate <= 24000 && (config->bitRate/elInfo->nChannelsInEl*2) < 60000) ) {

      float scfUsedRatio = (float) hAacEnc->psyKernel.psyConfLong.sfbActive / hAacEnc->psyKernel.psyConfLong.sfbCnt ;
      
       /* counting previous operation */
      
      
      error = InitStereoPreProcessing(&(hAacEnc->stereoPrePro),
                                      elInfo->nChannelsInEl, 
                                      config->bitRate,
                                      config->sampleRate,
                                      scfUsedRatio);
    }
  }
  
  
  if (error) {
    
    
    AacEncClose(hAacEnc);
    
    
  }
  
  
  
  
  
  return error;
}


int AacEncEncode(struct AAC_ENCODER *aacEnc,          /*!< an encoder handle */
                 float  *timeSignal,                  /*!< BLOCKSIZE*nChannels audio samples */
                 unsigned int timeInStride,
                 const unsigned char *ancBytes,       /*!< pointer to ancillary data bytes */
                 unsigned int        *numAncBytes,    /*!< number of ancillary Data Bytes */
                 unsigned int        *outBytes,       /*!< pointer to output buffer            */
                 int                 *numOutBytes     /*!< number of bytes in output buffer */
                 )
{
  ELEMENT_INFO* elInfo = &aacEnc->elInfo;
  int globUsedBits;
  int ancDataBytes, ancDataBytesLeft;
  
  

    
  aacEnc->hBitStream = CreateBitBuffer(&aacEnc->bitStream, 
                                       (unsigned char*) outBytes, 
                                       (MAX_CHANNEL_BITS/8)*MAX_CHANNELS);

    /* counting previous operation */

  
  ancDataBytes = ancDataBytesLeft =  *numAncBytes;


  /* advance psychoacoustic */

   
  if (elInfo->elType == ID_CPE) {

      
    ApplyStereoPreProcess(&aacEnc->stereoPrePro,
                          timeInStride,
                          elInfo,
                          timeSignal,
                          FRAME_LEN_LONG);
  }

    
  psyMain(aacEnc->fftctx,
          timeInStride,
          elInfo,
          timeSignal,
          &aacEnc->psyKernel.psyData[elInfo->ChannelIndex[0]],
          &aacEnc->psyKernel.tnsData[elInfo->ChannelIndex[0]],
          &aacEnc->psyKernel.psyConfLong,
          &aacEnc->psyKernel.psyConfShort,
          &aacEnc->psyOut.psyOutChannel[elInfo->ChannelIndex[0]],
          &aacEnc->psyOut.psyOutElement,
          aacEnc->psyKernel.pScratchTns);
  

    
  AdjustBitrate(&aacEnc->qcKernel,
                aacEnc->config->bitRate,
                aacEnc->config->sampleRate);

   /* 
                  aacEnc->qcKernel.elementBits
                  aacEnc->qcKernel.adjThr.adjThrStateElem,
                  aacEnc->psyOut.psyOutElement,
                  aacEnc->qcOut.qcElement
               */

     /* min() */   
  QCMain( aacEnc->aac_ram,
          &aacEnc->qcKernel,
          elInfo->nChannelsInEl,
          &aacEnc->qcKernel.elementBits,
          &aacEnc->qcKernel.adjThr.adjThrStateElem,
          &aacEnc->psyOut.psyOutChannel[elInfo->ChannelIndex[0]],
          &aacEnc->psyOut.psyOutElement,
          &aacEnc->qcOut.qcChannel[elInfo->ChannelIndex[0]],
          &aacEnc->qcOut.qcElement,
          min(ancDataBytesLeft,ancDataBytes));

    
   
  if ( elInfo->elType == ID_CPE ) {

      
    UpdateStereoPreProcess(&aacEnc->psyOut.psyOutChannel[elInfo->ChannelIndex[0]],
                           &aacEnc->qcOut.qcElement,
                           &aacEnc->stereoPrePro,
                           aacEnc->psyOut.psyOutElement.weightMsLrPeRatio);
  }

  
  ancDataBytesLeft-=ancDataBytes;
   
  
    
  FinalizeBitConsumption( &aacEnc->qcKernel,
                          &aacEnc->qcOut);

    
  WriteBitstreamData( aacEnc->hBitStream,
                  *elInfo,
                  &aacEnc->qcOut,
                  &aacEnc->psyOut,
                  &globUsedBits,
                  ancBytes);

    
  UpdateBitres(&aacEnc->qcKernel,
               &aacEnc->qcOut);

  /* write out the bitstream */
     
  *numOutBytes = GetBitsAvail(aacEnc->hBitStream)/8;
  
  /* assert this frame is not too large */
  assert(*numOutBytes*8 <= MAX_CHANNEL_BITS * elInfo->nChannelsInEl);
  
  
  
  return 0;
}


/*---------------------------------------------------------------------------

    functionname:AacEncClose
    description: deallocate an encoder instance

  ---------------------------------------------------------------------------*/

void AacEncClose (struct AAC_ENCODER* hAacEnc)
{
  int error=0;

  

   /* counting previous operation */

  
  if (hAacEnc) {
    
      
    QCDelete(&hAacEnc->qcKernel);
    
      
    QCOutDelete(&hAacEnc->qcOut);
    
      
    error = PsyDelete(&hAacEnc->psyKernel);
    
      
    error = PsyOutDelete(&hAacEnc->psyOut);
    
    if(hAacEnc->hBitStream)
        DeleteBitBuffer(&hAacEnc->hBitStream);
    
    
  }
  
  
}
