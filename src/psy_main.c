/*
  Psychoaccoustic main module
*/
#include <stdlib.h>
#include <string.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

/*                                    long       start       short       stop */
static int blockType2windowShape[] = {KBD_WINDOW,SINE_WINDOW,SINE_WINDOW,KBD_WINDOW};

/*
   forward definitions
*/
static int advancePsychLong(PSY_DATA* psyData,
                            TNS_DATA* tnsData,
                            PSY_CONFIGURATION_LONG *psyConfLong,
                            PSY_OUT_CHANNEL* psyOutChannel,
                            float *pScratchTns,
                            const TNS_DATA *tnsData2,
                            const int ch);

static int advancePsychLongMS (PSY_DATA  psyData[MAX_CHANNELS],
                               PSY_CONFIGURATION_LONG *psyConfLong);

static int advancePsychShort(PSY_DATA* psyData,
                             TNS_DATA* tnsData,
                             PSY_CONFIGURATION_SHORT *psyConfShort,
                             PSY_OUT_CHANNEL* psyOutChannel,
                             float *pScratchTns,
                             const TNS_DATA *tnsData2,
                             const int ch);

static int advancePsychShortMS (PSY_DATA  psyData[MAX_CHANNELS],
                                PSY_CONFIGURATION_SHORT *psyConfShort);


/*****************************************************************************

    functionname: PsyNew
    description:  allocates memory for psychoacoustic
    returns:      an error code
    input:        pointer to a psych handle

*****************************************************************************/
int PsyNew(AACRam_t *aacram, SBRRam_t *sbrram, PSY_KERNEL  *hPsy, int nChan)
{
    int i;

    

     /* hPsy->psyData[]
                    mdctDelayBuffer[]
                    mdctSpectrum[]
                 */
    
    for (i=0; i<nChan; i++){
      /*
        reserve memory for mdct delay buffer
      */
      
      hPsy->psyData[i].mdctDelayBuffer = &aacram->mdctDelayBuffer[i*BLOCK_SWITCHING_OFFSET];

      
      /*
        reserve memory for mdct spectum
      */
      
      
      hPsy->psyData[i].mdctSpectrum = &sbrram->sbr_envRBuffer[i*FRAME_LEN_LONG];

    }
  
     
    hPsy->pScratchTns = sbrram->sbr_envIBuffer;

    
 
    return 0;
}

/*****************************************************************************

    functionname: PsyDelete
    description:  allocates memory for psychoacoustic
    returns:      an error code

*****************************************************************************/
int PsyDelete(PSY_KERNEL  *hPsy)
{
  

  /*
    nothing to do
  */
  
  hPsy=NULL;

  
  return 0;
}


/*****************************************************************************

    functionname: PsyOutNew
    description:  allocates memory for psyOut struc
    returns:      an error code
    input:        pointer to a psych handle

*****************************************************************************/
int PsyOutNew(PSY_OUT *hPsyOut)
{
  
  

      
  memset(hPsyOut,0,sizeof(PSY_OUT));
 
  /*
    alloc some more stuff, tbd
  */

  
  return 0;
}

/*****************************************************************************

    functionname: PsyOutDelete
    description:  allocates memory for psychoacoustic
    returns:      an error code

*****************************************************************************/
int PsyOutDelete(PSY_OUT *hPsyOut)
{
  

  
  hPsyOut=NULL;

  
  return 0;
}


/*****************************************************************************

    functionname: psyMainInit
    description:  initializes psychoacoustic
    returns:      an error code

*****************************************************************************/

int psyMainInit(PSY_KERNEL *hPsy,
                int sampleRate,
                int bitRate,
                int channels,
                int tnsMask,
                int bandwidth)
{
  int ch, err;

  

    
  err = InitPsyConfiguration(bitRate/channels, sampleRate, bandwidth, &(hPsy->psyConfLong));

  
  if (!err)
  {
      
    err = InitTnsConfiguration(bitRate, sampleRate, channels,
                               &hPsy->psyConfLong.tnsConf, hPsy->psyConfLong, (int)(tnsMask&2));
  }

  
  if (!err)
  {
     
    err = InitPsyConfigurationShort(bitRate/channels, sampleRate, bandwidth, &hPsy->psyConfShort);
  }

  
  if (!err)
  {
      
    err =
      InitTnsConfigurationShort(bitRate, sampleRate, channels,
                                &hPsy->psyConfShort.tnsConf, hPsy->psyConfShort,(int)(tnsMask&1));
  }

  
  if (!err)
  {
     /* hPsy->psyData[] */
    
    for(ch=0;ch < channels;ch++){
  
       
      InitBlockSwitching(&hPsy->psyData[ch].blockSwitchingControl,
                         bitRate, channels);

       
      InitPreEchoControl(   hPsy->psyData[ch].sfbThresholdnm1,
                            hPsy->psyConfLong.sfbCnt,
                            hPsy->psyConfLong.sfbThresholdQuiet);
  }
  }
 
  

  return(err);
}







/*****************************************************************************

    functionname: psyMain
    description:  psychoacoustic
    returns:      an error code

        This function assumes that enough input data is in the modulo buffer.

*****************************************************************************/

int psyMain(FFTWFContext_t *fftctx,
            int   timeInStride,
            ELEMENT_INFO *elemInfo,
            float *timeSignal, 
            PSY_DATA  psyData[MAX_CHANNELS],
            TNS_DATA  tnsData[MAX_CHANNELS],
            PSY_CONFIGURATION_LONG  *psyConfLong,
            PSY_CONFIGURATION_SHORT *psyConfShort,
            PSY_OUT_CHANNEL  psyOutChannel[MAX_CHANNELS],
            PSY_OUT_ELEMENT* psyOutElement,
            float* pScratchTns)
{
  int maxSfbPerGroup[MAX_CHANNELS];
  int groupedSfbOffset[MAX_CHANNELS][MAX_GROUPED_SFB+1];  /* plus one for last dummy offset ! */
  float groupedSfbMinSnr[MAX_CHANNELS][MAX_GROUPED_SFB];
  
  int ch;   /* counts through channels          */
  int sfb;  /* counts through scalefactor bands */
  int line; /* counts through lines             */
  int channels = elemInfo->nChannelsInEl;
  
  
  
  
  /* block switching */
  
   /* pointers for psyData[],
                  elemInfo->ChannelIndex[]
               */
  
  for(ch = 0; ch < channels; ch++) {
    
     
    BlockSwitching (&psyData[ch].blockSwitchingControl,
                    timeSignal+elemInfo->ChannelIndex[ch],
                    timeInStride);
  }
  
  /* synch left and right block type */
   
  SyncBlockSwitching( &psyData[0].blockSwitchingControl,
                      &psyData[1].blockSwitchingControl,
                      channels);
  
  
  /* transform */
   /* pointers for psyData[],
                  elemInfo->ChannelIndex[]
               */
  
  for(ch = 0; ch < channels; ch++) {

     
    if(psyData[ch].blockSwitchingControl.windowSequence != SHORT_WINDOW){
      
       
      Transform_Real( fftctx,
                      psyData[ch].mdctDelayBuffer,
                      timeSignal+elemInfo->ChannelIndex[ch],
                      timeInStride,
                      psyData[ch].mdctSpectrum,
                      psyData[ch].blockSwitchingControl.windowSequence);
    }
    else {
      
       
      Transform_Real( fftctx,
                      psyData[ch].mdctDelayBuffer,
                      timeSignal+elemInfo->ChannelIndex[ch],
                      timeInStride,
                      psyData[ch].mdctSpectrum,
                      SHORT_WINDOW);
    }
  }
  
   /* pointers for psyData[],
                  tnsData[],
                  psyOutChannel[],
                  maxSfbPerGroup[]
               */
  
  for(ch = 0; ch < channels; ch++) {
    
     
    if(psyData[ch].blockSwitchingControl.windowSequence != SHORT_WINDOW) {

       
      advancePsychLong(   &psyData[ch],
                          &tnsData[ch],
                          psyConfLong,
                          &psyOutChannel[ch],
                          pScratchTns,
                          &tnsData[1-ch],
                          ch);
      
      /* determine maxSfb */
       /* pointers for psyConfLong->sfbOffset[] */
      
      for (sfb = psyConfLong->sfbCnt-1; sfb >= 0; sfb--) {
        
         /* pointers for psyData[ch].mdctSpectrum[] */
        
        for (line = psyConfLong->sfbOffset[sfb+1]-1; line >= psyConfLong->sfbOffset[sfb]; line--) {
          
          
          if (psyData[ch].mdctSpectrum[line] != 0) break;
        }
         
        if (line >= psyConfLong->sfbOffset[sfb]) break;
      }
      
       
      maxSfbPerGroup[ch] = sfb + 1;
      
      /* Calc bandwise energies for mid and side channel
         Do it only if 2 channels exist */
       
      if (ch==1) {

        
        advancePsychLongMS(psyData, psyConfLong);
      }
      
    }
    else {

       
      advancePsychShort(  &psyData[ch],
                          &tnsData[ch],
                          psyConfShort,
                          &psyOutChannel[ch],
                          pScratchTns,
                          &tnsData[1-ch],
                          ch);
      
      /* Calc bandwise energies for mid and side channel
         Do it only if 2 channels exist */
       
      if (ch==1) {

        
        advancePsychShortMS (psyData, psyConfShort);
      }
      
    }
  }
  
  /* group short data (maxSfb for short blocks is determined here) */
   /* pointers for psyData[],
                  tnsData[],
                  psyOutChannel[],
                  maxSfbPerGroup[],
                  groupedSfbOffset[]
               */
  
  for(ch=0;ch<channels;ch++) {

     
    if(psyData[ch].blockSwitchingControl.windowSequence == SHORT_WINDOW) {
        
        
      groupShortData( psyData[ch].mdctSpectrum,
                      pScratchTns,
                      &psyData[ch].sfbThreshold,
                      &psyData[ch].sfbEnergy,
                      &psyData[ch].sfbEnergyMS,
                      &psyData[ch].sfbSpreadedEnergy,
                      psyConfShort->sfbCnt,
                      psyConfShort->sfbOffset,
                      psyConfShort->sfbMinSnr,
                      groupedSfbOffset[ch],
                      &maxSfbPerGroup[ch],
                      groupedSfbMinSnr[ch],
                      psyData[ch].blockSwitchingControl.noOfGroups,
                      psyData[ch].blockSwitchingControl.groupLen);
    }
    
  }
  
#if (MAX_CHANNELS>1)
  
  /*
    stereo Processing
  */
   
  if(channels == 2) {

      
    psyOutElement->toolsInfo.msDigest = MS_NONE;
    
      
    maxSfbPerGroup[0] = maxSfbPerGroup[1] = max(maxSfbPerGroup[0], maxSfbPerGroup[1]);

     
    if(psyData[0].blockSwitchingControl.windowSequence != SHORT_WINDOW) {
      
        
      MsStereoProcessing( psyData[0].sfbEnergy.Long,
                          psyData[1].sfbEnergy.Long,
                          psyData[0].sfbEnergyMS.Long,
                          psyData[1].sfbEnergyMS.Long,
                          psyData[0].mdctSpectrum,
                          psyData[1].mdctSpectrum,
                          psyData[0].sfbThreshold.Long,
                          psyData[1].sfbThreshold.Long,
                          psyData[0].sfbSpreadedEnergy.Long,
                          psyData[1].sfbSpreadedEnergy.Long,
                          &psyOutElement->toolsInfo.msDigest,
                          psyOutElement->toolsInfo.msMask,
                          psyConfLong->sfbCnt,
                          psyConfLong->sfbCnt,
                          maxSfbPerGroup[0],
                          psyConfLong->sfbOffset,
                          &psyOutElement->weightMsLrPeRatio);
    }
    else {
      
         
      MsStereoProcessing( psyData[0].sfbEnergy.Long,
                          psyData[1].sfbEnergy.Long,
                          psyData[0].sfbEnergyMS.Long,
                          psyData[1].sfbEnergyMS.Long,
                          psyData[0].mdctSpectrum,
                          psyData[1].mdctSpectrum,
                          psyData[0].sfbThreshold.Long,
                          psyData[1].sfbThreshold.Long,
                          psyData[0].sfbSpreadedEnergy.Long,
                          psyData[1].sfbSpreadedEnergy.Long,
                          &psyOutElement->toolsInfo.msDigest,
                          psyOutElement->toolsInfo.msMask,
                          psyData[0].blockSwitchingControl.noOfGroups*psyConfShort->sfbCnt,
                          psyConfShort->sfbCnt,
                          maxSfbPerGroup[0],
                          groupedSfbOffset[0],
                          &psyOutElement->weightMsLrPeRatio);
    }
  }
  
#endif /* (MAX_CHANNELS>1) */
  
  /*
    build output
  */
   /* pointers for psyData[],
                  psyOutChannel[],
                  maxSfbPerGroup[],
                  groupedSfbOffset[],
                  groupedSfbMinSnr[]
               */
  
  for(ch=0;ch<channels;ch++) {

     
    if(psyData[ch].blockSwitchingControl.windowSequence != SHORT_WINDOW) {

         
      BuildInterface( psyData[ch].mdctSpectrum,
                      &psyData[ch].sfbThreshold,
                      &psyData[ch].sfbEnergy,
                      &psyData[ch].sfbSpreadedEnergy,
                      psyData[ch].sfbEnergySum,
                      psyData[ch].sfbEnergySumMS,
                      psyData[ch].blockSwitchingControl.windowSequence,
                      blockType2windowShape[psyData[ch].blockSwitchingControl.windowSequence],
                      psyConfLong->sfbCnt,
                      psyConfLong->sfbOffset,
                      maxSfbPerGroup[ch],
                      psyConfLong->sfbMinSnr,
                      psyData[ch].blockSwitchingControl.noOfGroups,
                      psyData[ch].blockSwitchingControl.groupLen,
                      &psyOutChannel[ch]);
    }
    else {

         
      BuildInterface( psyData[ch].mdctSpectrum,
                      &psyData[ch].sfbThreshold,
                      &psyData[ch].sfbEnergy,
                      &psyData[ch].sfbSpreadedEnergy,
                      psyData[ch].sfbEnergySum,
                      psyData[ch].sfbEnergySumMS,
                      SHORT_WINDOW,
                      SINE_WINDOW,
                      psyData[0].blockSwitchingControl.noOfGroups*psyConfShort->sfbCnt,
                      groupedSfbOffset[ch],
                      maxSfbPerGroup[ch],
                      groupedSfbMinSnr[ch],
                      psyData[ch].blockSwitchingControl.noOfGroups,
                      psyData[ch].blockSwitchingControl.groupLen,
                      &psyOutChannel[ch]);
    }
  }
  
  
  
  return 0; /* no error */
}


/*****************************************************************************

    functionname: advancePsychLong
    description:  psychoacoustic for long blocks

*****************************************************************************/

static int advancePsychLong(PSY_DATA* psyData,
                            TNS_DATA* tnsData,
                            PSY_CONFIGURATION_LONG *psyConfLong,
                            PSY_OUT_CHANNEL* psyOutChannel,
                            float *pScratchTns,
                            const TNS_DATA* tnsData2,
                            const int ch)
{
    int i;
  
    /* the energy Shift is actually needed since we do not multiply the spectrum by two after the cmdct */
    float energyShift = 0.25f;

    float clipEnergy = psyConfLong->clipEnergy * energyShift;

    

       /* counting previous operations */
 
    /* low pass */
     /* pointers for psyData->mdctSpectrum[] */
     
    for(i=psyConfLong->lowpassLine; i<FRAME_LEN_LONG; i++){
        
        psyData->mdctSpectrum[i] = 0.0f;
    }

    /* Calc sfb-bandwise mdct-energies for left and right channel */
      
    CalcBandEnergy( psyData->mdctSpectrum,
                    psyConfLong->sfbOffset,
                    psyConfLong->sfbActive,
                    psyData->sfbEnergy.Long,
                    &psyData->sfbEnergySum.Long);


    /*
        TNS
    */
     
    TnsDetect(  tnsData,
                psyConfLong->tnsConf,
                pScratchTns,
                psyConfLong->sfbOffset,
                psyData->mdctSpectrum,
                0,
                (int)psyData->blockSwitchingControl.windowSequence,
                psyData->sfbEnergy.Long);

    /*  TnsSync */
     
    if (ch==1) {
        
       TnsSync(tnsData,
               tnsData2,
               psyConfLong->tnsConf,
               0,
               (int)psyData->blockSwitchingControl.windowSequence);
    }

      
    TnsEncodeData(  &psyOutChannel->tnsInfo,
                tnsData,
                psyConfLong->sfbCnt,
                psyConfLong->tnsConf,
                psyConfLong->lowpassLine,
                psyData->mdctSpectrum,
                0,
                (int)psyData->blockSwitchingControl.windowSequence);

    /* first part of threshold calculation */
     /* pointers for psyData->sfbThreshold.Long[],
                                 psyData->sfbEnergy.Long[]
                 */
     
    for (i=0; i<psyConfLong->sfbCnt; i++) {

        /* multiply sfbEnergy by ratio */
          
        psyData->sfbThreshold.Long[i] = psyData->sfbEnergy.Long[i] * psyConfLong->ratio;

        /* limit threshold to avoid clipping */
          
        psyData->sfbThreshold.Long[i] = min(psyData->sfbThreshold.Long[i], clipEnergy);
    }

    /* Calc sfb-bandwise mdct-energies for left and right channel again,
       if tns has modified the spectrum */
      
    if (psyOutChannel->tnsInfo.tnsActive[0]==1) {

          
        CalcBandEnergy( psyData->mdctSpectrum,
                        psyConfLong->sfbOffset,
                        psyConfLong->sfbActive,
                        psyData->sfbEnergy.Long,
                        &psyData->sfbEnergySum.Long);
    }


    /* spreading */
     
    SpreadingMax(   psyConfLong->sfbCnt,
                    psyConfLong->sfbMaskLowFactor,
                    psyConfLong->sfbMaskHighFactor,
                    psyData->sfbThreshold.Long);

    /* threshold in quiet */
     /* pointers for psyData->sfbThreshold.Long[],
                                 psyConfLong.sfbThresholdQuiet[]
                 */
     
    for (i=0; i<psyConfLong->sfbCnt;i++)
    {
           
        psyData->sfbThreshold.Long[i] = max(psyData->sfbThreshold.Long[i],
                                            (psyConfLong->sfbThresholdQuiet[i] * energyShift));
    }


    /* preecho control */
      
    if(psyData->blockSwitchingControl.windowSequence == STOP_WINDOW)
    {
        /* prevent PreEchoControl from comparing stop
           thresholds with short thresholds */
         /* pointer for psyData->sfbThresholdnm1[] */
         
        for (i=0; i<psyConfLong->sfbCnt;i++)
        {
            
            psyData->sfbThresholdnm1[i] = 1.0e20f;
        }
    }

     
    PreEchoControl( psyData->sfbThresholdnm1,
                    psyConfLong->sfbCnt,
                    psyConfLong->maxAllowedIncreaseFactor,
                    psyConfLong->minRemainingThresholdFactor,
                    psyData->sfbThreshold.Long);

      
    if(psyData->blockSwitchingControl.windowSequence == START_WINDOW)
    {
        /* prevent PreEchoControl in next frame to compare start
           thresholds with short thresholds */
         /* pointer for psyData->sfbThresholdnm1[] */
         
        for (i=0; i<psyConfLong->sfbCnt;i++)
        {
            
            psyData->sfbThresholdnm1[i] = 1.0e20f;
        }
    }

    /* apply tns mult table on cb thresholds */
    
    if (psyOutChannel->tnsInfo.tnsActive[0]) {
       
      ApplyTnsMultTableToRatios(  psyConfLong->tnsConf.tnsRatioPatchLowestCb,
                                  psyConfLong->tnsConf.tnsStartBand,
                                  psyData->sfbThreshold.Long);
    }

    /* spreaded energy for avoid hole detection */
     /* pointers for psyData->sfbSpreadedEnergy.Long[],
                                 psyData->sfbEnergy.Long[]
                 */
     
    for (i=0; i<psyConfLong->sfbCnt; i++)
    {
      
      psyData->sfbSpreadedEnergy.Long[i] = psyData->sfbEnergy.Long[i];
    }

     
    SpreadingMax(psyConfLong->sfbCnt,
                 psyConfLong->sfbMaskLowFactorSprEn, 
                 psyConfLong->sfbMaskHighFactorSprEn,
                 psyData->sfbSpreadedEnergy.Long);

    

    return 0;
}


static int advancePsychLongMS (PSY_DATA psyData[MAX_CHANNELS],
                               PSY_CONFIGURATION_LONG *psyConfLong)
{
    

      
    CalcBandEnergyMS(psyData[0].mdctSpectrum,
                     psyData[1].mdctSpectrum,
                     psyConfLong->sfbOffset,
                     psyConfLong->sfbActive,
                     psyData[0].sfbEnergyMS.Long,
                     &psyData[0].sfbEnergySumMS.Long,
                     psyData[1].sfbEnergyMS.Long,
                     &psyData[1].sfbEnergySumMS.Long);

    

    return 0;
}


/*****************************************************************************

    functionname: advancePsychShort
    description:  psychoacoustic for short blocks

*****************************************************************************/

static int advancePsychShort(PSY_DATA* psyData,
                             TNS_DATA* tnsData,
                             PSY_CONFIGURATION_SHORT *psyConfShort,
                             PSY_OUT_CHANNEL* psyOutChannel,
                             float *pScratchTns,
                             const TNS_DATA *tnsData2,
                             const int ch)
{
    int w;

    /* the energy Shift is actually needed since we do not multiply the spectrum by two after the cmdct */
    float energyShift = 0.25f;
    float clipEnergy = psyConfShort->clipEnergy *energyShift;

    

       /* counting previous operations */

     /* pointers for psyData->sfbEnergy.Short[w],
                                 psyData->sfbEnergySum.Short[w],
                                 psyData->sfbThreshold.Short[w],
                                 psyData->sfbSpreadedEnergy.Short[w],
                                 tnsData->dataRaw.Short.subBlockInfo[w],
                                 tnsData->dataRaw.Short.ratioMultTable[w]
                 */
    
    for(w = 0; w < TRANS_FAC; w++)
    {
        int wOffset=w*FRAME_LEN_SHORT;
        int i;

         /* counting previous operation */

        /* low pass */
         /* pointers for psyData->mdctSpectrum[] */
         
        for(i=psyConfShort->lowpassLine; i<FRAME_LEN_SHORT; i++){
            
            psyData->mdctSpectrum[i+wOffset] = 0.0f;
        }

        /* Calc sfb-bandwise mdct-energies for left and right channel */
         
        CalcBandEnergy( psyData->mdctSpectrum+wOffset,
                        psyConfShort->sfbOffset,
                        psyConfShort->sfbActive,
                        psyData->sfbEnergy.Short[w],
                        &psyData->sfbEnergySum.Short[w]);

        /*
        TNS
        */
         
        TnsDetect(  tnsData,
                    psyConfShort->tnsConf,
                    pScratchTns,
                    psyConfShort->sfbOffset,
                    psyData->mdctSpectrum+wOffset,
                    w,
                    (int)psyData->blockSwitchingControl.windowSequence,
                    psyData->sfbEnergy.Short[w]);

        /*  TnsSync */
         
        if (ch==1) {
             
            TnsSync(tnsData,
                    tnsData2,
                    psyConfShort->tnsConf,
                    w,
                    (int)psyData->blockSwitchingControl.windowSequence);
        }

          
        TnsEncodeData(  &psyOutChannel->tnsInfo,
                    tnsData,
                    psyConfShort->sfbCnt,
                    psyConfShort->tnsConf,
                    psyConfShort->lowpassLine,
                    psyData->mdctSpectrum+wOffset,
                    w,
                    (int)psyData->blockSwitchingControl.windowSequence);

        /* first part of threshold calculation */
         /* pointers for psyData->sfbThreshold.Short[][],
                                     psyData->sfbEnergy.Short[][]
                     */
         
        for (i=0; i<psyConfShort->sfbCnt; i++) {

            /* multiply sfbEnergy by ratio */
              
            psyData->sfbThreshold.Short[w][i] = psyData->sfbEnergy.Short[w][i] * psyConfShort->ratio;

            /* limit threshold to avoid clipping */
              
            psyData->sfbThreshold.Short[w][i] = min(psyData->sfbThreshold.Short[w][i], clipEnergy);
        }

        /* Calc sfb-bandwise mdct-energies for left and right channel again,
           if tns has modified the spectrum */

         
        if (psyOutChannel->tnsInfo.tnsActive[w]) {
           
          CalcBandEnergy( psyData->mdctSpectrum+wOffset,
                          psyConfShort->sfbOffset,
                          psyConfShort->sfbActive,
                          psyData->sfbEnergy.Short[w],
                          &psyData->sfbEnergySum.Short[w]);

        }

        /* spreading */
         
        SpreadingMax(   psyConfShort->sfbCnt,
                        psyConfShort->sfbMaskLowFactor,
                        psyConfShort->sfbMaskHighFactor,
                        psyData->sfbThreshold.Short[w]);


        /* threshold in quiet */
         /* pointers for psyData->sfbThreshold.Short[w][i] 
                                     psyConfShort.sfbThresholdQuiet[i]
                     */
         
        for(i=0;i<psyConfShort->sfbCnt;i++)
        {
               
            psyData->sfbThreshold.Short[w][i] = max(psyData->sfbThreshold.Short[w][i],
                                                    (psyConfShort->sfbThresholdQuiet[i] * energyShift));
        }


        /* preecho */
         
        PreEchoControl( psyData->sfbThresholdnm1,
                        psyConfShort->sfbCnt,
                        psyConfShort->maxAllowedIncreaseFactor,
                        psyConfShort->minRemainingThresholdFactor,
                        psyData->sfbThreshold.Short[w]);

        /* apply tns mult table on cb thresholds */
        
        if (psyOutChannel->tnsInfo.tnsActive[w]) {
           
          ApplyTnsMultTableToRatios( psyConfShort->tnsConf.tnsRatioPatchLowestCb,
                                     psyConfShort->tnsConf.tnsStartBand,
                                     psyData->sfbThreshold.Short[w]);
        }
        
        /* spreaded energy for avoid hole detection */
         /* pointers for psyData->sfbSpreadedEnergy.Short[][],
                                     psyData->sfbEnergy.Short[][]
                     */
         
        for (i=0; i<psyConfShort->sfbCnt; i++)
        {
          
          psyData->sfbSpreadedEnergy.Short[w][i] = psyData->sfbEnergy.Short[w][i];
        }

         
        SpreadingMax(psyConfShort->sfbCnt,
                     psyConfShort->sfbMaskLowFactorSprEn, 
                     psyConfShort->sfbMaskHighFactorSprEn,
                     psyData->sfbSpreadedEnergy.Short[w]);

    } /* for TRANS_FAC */

    

    return 0;
}


static int advancePsychShortMS (PSY_DATA psyData[MAX_CHANNELS],
                                PSY_CONFIGURATION_SHORT *psyConfShort)
{
    int w;

    

     /* pointers for psyData[0].sfbEnergyMS.Short[w],
                                 psyData[0].sfbEnergySumMS.Short[w],
                                 psyData[1].sfbEnergyMS.Short[w],
                                 psyData[1].sfbEnergySumMS.Short[w]
                 */
    
    for(w = 0; w < TRANS_FAC; w++) {
        int   wOffset=w*FRAME_LEN_SHORT;

         /* counting previous operation */

         
        CalcBandEnergyMS(   psyData[0].mdctSpectrum+wOffset,
                            psyData[1].mdctSpectrum+wOffset,
                            psyConfShort->sfbOffset,
                            psyConfShort->sfbActive,
                            psyData[0].sfbEnergyMS.Short[w],
                            &psyData[0].sfbEnergySumMS.Short[w],
                            psyData[1].sfbEnergyMS.Short[w],
                            &psyData[1].sfbEnergySumMS.Short[w]);
    }

    

    return 0;
}
