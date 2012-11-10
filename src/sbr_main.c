/*
  SBR encoder top level processing
*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#define INVALID_TABLE_IDX -1

/*
   tuningTable
*/
static const struct
{
  unsigned int    bitrateFrom ;
  unsigned int    bitrateTo ;

  unsigned int    sampleRate ;
  unsigned int    numChannels ;

  unsigned int    startFreq ;
  unsigned int    stopFreq ;

  int numNoiseBands ;
  int noiseFloorOffset ;
  int noiseMaxLevel ;
  SBR_STEREO_MODE stereoMode ;
  int freqScale ;

} tuningTable[] = {


  /*** mono ***/
  { 10000, 12000,  16000, 1,  1,  3, 1, 0, 6, SBR_MONO, 3 }, /* nominal: 10 kbit/s */
  { 12000, 16000,  16000, 1,  2,  0, 1, 0, 6, SBR_MONO, 3 }, /* nominal: 12 kbit/s */
  { 16000, 20000,  16000, 1,  2,  3, 1, 0, 6, SBR_MONO, 3 }, /* nominal: 18 kbit/s */

  { 12000, 18000,  22050, 1,  1,  4, 1, 0, 6, SBR_MONO, 3 }, /* nominal: 14 kbit/s */
  { 18000, 22000,  22050, 1,  3,  5, 2, 0, 6, SBR_MONO, 2 }, /* nominal: 20 kbit/s */
  { 22000, 28000,  22050, 1,  7,  8, 2, 0, 6, SBR_MONO, 2 }, /* nominal: 24 kbit/s */
  { 28000, 36000,  22050, 1, 11,  9, 2, 0, 3, SBR_MONO, 2 }, /* nominal: 32 kbit/s */
  { 36000, 44001,  22050, 1, 11,  9, 2, 0, 3, SBR_MONO, 1 }, /* nominal: 40 kbit/s */

  { 12000, 18000,  24000, 1,  1,  4, 1, 0, 6, SBR_MONO, 3 }, /* nominal: 14 kbit/s */
  { 18000, 22000,  24000, 1,  3,  5, 2, 0, 6, SBR_MONO, 2 }, /* nominal: 20 kbit/s */
  { 22000, 28000,  24000, 1,  7,  8, 2, 0, 6, SBR_MONO, 2 }, /* nominal: 24 kbit/s */
  { 28000, 36000,  24000, 1, 10,  9, 2, 0, 3, SBR_MONO, 2 }, /* nominal: 32 kbit/s */
  { 36000, 44001,  24000, 1, 11,  9, 2, 0, 3, SBR_MONO, 1 }, /* nominal: 40 kbit/s */

  /*** stereo ***/
  { 18000, 24000,  16000, 2,  4,  2, 1, 0, -3, SBR_SWITCH_LRC, 3 }, /* nominal: 18 kbit/s */

  { 24000, 28000,  22050, 2,  5,  6, 1, 0, -3, SBR_SWITCH_LRC, 3 }, /* nominal: 24 kbit/s */
  { 28000, 36000,  22050, 2,  7,  8, 2, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 32 kbit/s */
  { 36000, 44000,  22050, 2, 10,  9, 2, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 40 kbit/s */
  { 44000, 52000,  22050, 2, 12,  9, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 48 kbit/s */
  { 52000, 60000,  22050, 2, 12,  9, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 56 kbit/s */
  { 60000, 68000,  22050, 2, 14, 10, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 64 kbit/s */
  { 68000, 72001,  22050, 2, 14, 10, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 72 kbit/s */

  { 24000, 28000,  24000, 2,  5,  6, 1, 0, -3, SBR_SWITCH_LRC, 3 }, /* nominal: 24 kbit/s */
  { 28000, 36000,  24000, 2,  7,  8, 2, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 32 kbit/s */
  { 36000, 44000,  24000, 2, 10,  9, 2, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 40 kbit/s */
  { 44000, 52000,  24000, 2, 12,  9, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 48 kbit/s */
  { 52000, 60000,  24000, 2, 12,  9, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 56 kbit/s */
  { 60000, 68000,  24000, 2, 14, 10, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 64 kbit/s */
  { 68000, 72001,  24000, 2, 14, 10, 3, 0, -3, SBR_SWITCH_LRC, 2 }, /* nominal: 72 kbit/s */

};


/***************************************************************************/
/*!

  \brief  Selects the SBR tuning.
  \return Index

****************************************************************************/
static int
getSbrTuningTableIndex(unsigned int bitrate,
                       unsigned int numChannels,
                       unsigned int sampleRate
                       )
{
  int i, paramSets = sizeof (tuningTable) / sizeof (tuningTable [0]);

  

   /* counting previous operation */

   /* tuningTable[] */
  
  for (i = 0 ; i < paramSets ; i++)  {
     
    if (numChannels == tuningTable [i].numChannels) {

        
      if ((sampleRate == tuningTable [i].sampleRate) &&
          (bitrate >= tuningTable [i].bitrateFrom) &&
          (bitrate < tuningTable [i].bitrateTo)) {
        
        return i ;
      }
    }
  }

  

  return INVALID_TABLE_IDX;
}

/***************************************************************************/
/*!

  \brief  Is SBR setting available.
  \return a flag indicating success: yes (1) or no (0)

****************************************************************************/
unsigned int
IsSbrSettingAvail (unsigned int bitrate,
                   unsigned int numOutputChannels,
                   unsigned int sampleRateInput,
                   unsigned int *sampleRateCore
                   )
{
  int idx = INVALID_TABLE_IDX;

  

   /* counting previous operation */

   
  if (sampleRateInput < 32000)
  {
    
    return 0;
  }

   
  *sampleRateCore = sampleRateInput/2;

  
  idx = getSbrTuningTableIndex(bitrate,numOutputChannels, *sampleRateCore);

    /* counting post-operation */

  

  return (idx == INVALID_TABLE_IDX) ? 0 : 1;
}


/***************************************************************************/
/*!

  \brief  Adjusts the SBR settings.
  \return success:

****************************************************************************/
unsigned int
AdjustSbrSettings (const sbrConfigurationPtr config,
                   unsigned int bitRate,
                   unsigned int numChannels,
                   unsigned int fsCore,
                   unsigned int transFac,
                   unsigned int standardBitrate
                   )
{
  int idx = INVALID_TABLE_IDX;
  unsigned int sampleRate;

  

   /* counting previous operation */

  /* set the codec settings */
   
  config->codecSettings.bitRate         = bitRate;
  config->codecSettings.nChannels       = numChannels;
  config->codecSettings.sampleFreq      = fsCore;
  config->codecSettings.transFac        = transFac;
  config->codecSettings.standardBitrate = standardBitrate;

  
  sampleRate  = fsCore * 2;

  
  idx = getSbrTuningTableIndex(bitRate,numChannels,fsCore);

   
  if (idx != INVALID_TABLE_IDX) {

     
    config->startFreq       = tuningTable [idx].startFreq ;
    config->stopFreq        = tuningTable [idx].stopFreq ;
    config->sbr_noise_bands = tuningTable [idx].numNoiseBands ;

     
    config->noiseFloorOffset= tuningTable[idx].noiseFloorOffset;

     
    config->ana_max_level   = tuningTable [idx].noiseMaxLevel ;
    config->stereoMode      = tuningTable[idx].stereoMode ;
    config->freqScale       = tuningTable[idx].freqScale ;

     
    if (bitRate <= 20000) {

       
      config->parametricCoding = 0;
      config->useSpeechConfig  = 1;
    }

#ifndef MONO_ONLY
     
    if(config->usePs)
    {
        
      config->psMode = GetPsMode(bitRate);
    }
#endif

    

    return 1 ;
  }

  

  return 0 ;
}


/*****************************************************************************

 description:  initializes the SBR confifuration
 returns:      error status

*****************************************************************************/
unsigned int
InitializeSbrDefaults (sbrConfigurationPtr config

                       )
{
  






   
  config->SendHeaderDataTime     = 500;
  config->crcSbr                 = 0;
  config->tran_thr               = 13000;
  config->detectMissingHarmonics = 1;
  config->parametricCoding       = 1;
  config->useSpeechConfig        = 0;



   
  config->sbr_data_extra         = 0;
  config->amp_res                = SBR_AMP_RES_3_0 ;
  config->tran_fc                = 0 ;
  config->tran_det_mode          = 1 ;
  config->spread                 = 1 ;
  config->stat                   = 0 ;
  config->e                      = 1 ;
  config->deltaTAcrossFrames     = 1 ;
  config->dF_edge_1stEnv         = 0.3f ;
  config->dF_edge_incr           = 0.3f ;

   
  config->sbr_invf_mode   = INVF_SWITCHED;
  config->sbr_xpos_mode   = XPOS_LC;
  config->sbr_xpos_ctrl   = SBR_XPOS_CTRL_DEFAULT;
  config->sbr_xpos_level  = 0;

   
  config->usePs           = 0;
  config->psMode          = -1;

   
  config->stereoMode             = SBR_SWITCH_LRC;
  config->ana_max_level          = 6;
  config->noiseFloorOffset       = 0;
  config->startFreq              = 5;
  config->stopFreq               = 9;

   
  config->freqScale       = SBR_FREQ_SCALE_DEFAULT;
  config->alterScale      = SBR_ALTER_SCALE_DEFAULT;
  config->sbr_noise_bands = SBR_NOISE_BANDS_DEFAULT;

   
  config->sbr_limiter_bands    = SBR_LIMITER_BANDS_DEFAULT;
  config->sbr_limiter_gains    = SBR_LIMITER_GAINS_DEFAULT;
  config->sbr_interpol_freq    = SBR_INTERPOL_FREQ_DEFAULT;
  config->sbr_smoothing_length = SBR_SMOOTHING_LENGTH_DEFAULT;

  

  return 1;
}


/*****************************************************************************

 description:  frees memory of one SBR channel
 returns:      void

*****************************************************************************/
static void
deleteEnvChannel (HANDLE_ENV_CHANNEL hEnvCut)
{

  

  
  if (hEnvCut) {

      
    deleteFrameInfoGenerator (&hEnvCut->SbrEnvFrame);

      
    deleteQmfBank (&hEnvCut->sbrQmf);

      
    deleteSbrCodeEnvelope (&hEnvCut->sbrCodeEnvelope);

      
    deleteSbrCodeEnvelope (&hEnvCut->sbrCodeNoiseFloor);


      
    deleteSbrTransientDetector (&hEnvCut->sbrTransientDetector);

      
    deleteExtractSbrEnvelope (&hEnvCut->sbrExtractEnvelope);


      
    DeleteTonCorrParamExtr(&hEnvCut->TonCorr);

  }

  
}


/*****************************************************************************

 description:  close the envelope coding handle
 returns:      void

*****************************************************************************/
void
EnvClose (HANDLE_SBR_ENCODER hEnvEnc)
{
  int i;

  

  
  if (hEnvEnc != NULL) {

     /* hEnvEnc->hEnvChannel[] */
    
    for (i = 0; i < MAX_CHANNELS; i++) {

      
      if (hEnvEnc->hEnvChannel[i] != NULL) {

        
        deleteEnvChannel (hEnvEnc->hEnvChannel[i]);
        free(hEnvEnc->hEnvChannel[i]);
        
        hEnvEnc->hEnvChannel[i] = NULL;
      }
    }

    
#ifndef MONO_ONLY
     
    if (hEnvEnc->hSynthesisQmfBank)
    {
       
      DeleteSynthesisQmfBank ((HANDLE_SBR_QMF_FILTER_BANK*)&hEnvEnc->hSynthesisQmfBank);
      free(hEnvEnc->hSynthesisQmfBank);
    }

     
    if (hEnvEnc->hPsEnc)
    {
       
      DeletePsEnc(&hEnvEnc->hPsEnc);
      free(hEnvEnc->hPsEnc);
    }
#endif /* #ifndef MONO_ONLY */

  }

  
}

/*****************************************************************************

 description:  updates vk_master
 returns:      void

*****************************************************************************/
static int updateFreqBandTable(HANDLE_SBR_CONFIG_DATA  sbrConfigData,
                               HANDLE_SBR_HEADER_DATA  sbrHeaderData,
                               int noQmfChannels)
{


  int k0, k2;

  

     
  if(FindStartAndStopBand(sbrConfigData->sampleFreq,
                          noQmfChannels,
                          sbrHeaderData->sbr_start_frequency,
                          sbrHeaderData->sbr_stop_frequency,
                          sbrHeaderData->sampleRateMode,
                          &k0, &k2))
  {
    
    return(1);
  }



     
  if(UpdateFreqScale(sbrConfigData->v_k_master, &sbrConfigData->num_Master,
                     k0, k2, sbrHeaderData->freqScale,
                     sbrHeaderData->alterScale))
  {
    
    return(1);
  }


   
  sbrHeaderData->sbr_xover_band=0;


     
  if(UpdateHiRes(sbrConfigData->freqBandTable[HI],
                 &sbrConfigData->nSfb[HI],
                 sbrConfigData->v_k_master,
                 sbrConfigData->num_Master ,
                 &sbrHeaderData->sbr_xover_band,
                 sbrHeaderData->sampleRateMode,
                 noQmfChannels))
  {
    
    return(1);
  }


    
  UpdateLoRes(sbrConfigData->freqBandTable[LO],
              &sbrConfigData->nSfb[LO],
              sbrConfigData->freqBandTable[HI],
              sbrConfigData->nSfb[HI]);

       
  sbrConfigData->xOverFreq = (sbrConfigData->freqBandTable[LOW_RES][0] * sbrConfigData->sampleFreq / noQmfChannels+1)>>1;


  

  return (0);
}

/*****************************************************************************

 description: performs the sbr envelope calculation
 returns:

*****************************************************************************/
int
EnvEncodeFrame (HANDLE_SBR_ENCODER hEnvEncoder,
                float *samples,
                float *pCoreBuffer,
                unsigned int timeInStride,
                unsigned int  *numAncBytes,
                unsigned char *ancData)
{
  

  
  if (hEnvEncoder != NULL) {
    /* header bitstream handling */

    HANDLE_SBR_BITSTREAM_DATA sbrBitstreamData = &hEnvEncoder->sbrBitstreamData;

     /* counting previous operation */

     
    sbrBitstreamData->HeaderActive = 0;

     
    if (sbrBitstreamData->CountSendHeaderData == 0)
    {
       
      sbrBitstreamData->HeaderActive = 1;
    }

     
    if (sbrBitstreamData->NrSendHeaderData == 0) {
       
      sbrBitstreamData->CountSendHeaderData = 1;
    }
    else {
        
      sbrBitstreamData->CountSendHeaderData++;

        
      sbrBitstreamData->CountSendHeaderData %= sbrBitstreamData->NrSendHeaderData;
    }



       
     InitSbrBitstream(&hEnvEncoder->CmonData,
                      (unsigned char*) hEnvEncoder->sbrPayload,
                      sizeof(hEnvEncoder->sbrPayload),
                      hEnvEncoder->sbrBitstreamData.CRCActive);


       
     extractSbrEnvelope(hEnvEncoder->fftctx,
                        samples,
                        pCoreBuffer,
                        timeInStride,
                        &hEnvEncoder->sbrConfigData,
                        &hEnvEncoder->sbrHeaderData,
                        &hEnvEncoder->sbrBitstreamData,
                        hEnvEncoder->hEnvChannel,
                        hEnvEncoder->hPsEnc,
                        hEnvEncoder->hSynthesisQmfBank,
                        &hEnvEncoder->CmonData);


      
    AssembleSbrBitstream(&hEnvEncoder->CmonData);


    assert(GetBitsAvail(&hEnvEncoder->CmonData.sbrBitbuf) % 8 == 0);


        
    hEnvEncoder->sbrPayloadSize = GetBitsAvail(&hEnvEncoder->CmonData.sbrBitbuf) / 8;

      
    if(hEnvEncoder->sbrPayloadSize > MAX_PAYLOAD_SIZE)
    {
      
      hEnvEncoder->sbrPayloadSize=0;
    }

   
   if(ancData){

       
      *numAncBytes = hEnvEncoder->sbrPayloadSize;

          
      memcpy(ancData,hEnvEncoder->sbrPayload,hEnvEncoder->sbrPayloadSize);
    }
  }

 
 return (0);
}

/*****************************************************************************

 description:  initializes parameters and allocates memory
 returns:      error status

*****************************************************************************/

static int
createEnvChannel (SBRRam_t *sbrram,
                  int chan,
                  HANDLE_SBR_CONFIG_DATA sbrConfigData,
                  HANDLE_SBR_HEADER_DATA sbrHeaderData,
                  HANDLE_ENV_CHANNEL     hEnv,
                  sbrConfigurationPtr params)
{

  int e;
  int tran_fc;

  int startIndex;
  int noiseBands[2] = { 3, 3 };

  

   /* counting previous operations */

   
  e = 1 << params->e;

   
  hEnv->encEnvData.freq_res_fixfix = FREQ_RES_HIGH;

   
  hEnv->encEnvData.sbr_xpos_mode = (XPOS_MODE)params->sbr_xpos_mode;
  hEnv->encEnvData.sbr_xpos_ctrl = params->sbr_xpos_ctrl;


     
  if(createQmfBank (sbrram, chan, &hEnv->sbrQmf)){

    
    return (1); /* initialisation failed */
  }

   
  startIndex = 576;


     
  if(CreateTonCorrParamExtr (sbrram,
                             chan,
                             &hEnv->TonCorr,



                             sbrConfigData->sampleFreq,
                             sbrConfigData->freqBandTable[LOW_RES][sbrConfigData->nSfb[LO]],
                             64,
                             params->sbr_xpos_ctrl,
                             sbrConfigData->freqBandTable[LOW_RES][0],
                             0,
                             sbrConfigData->v_k_master,
                             sbrConfigData->num_Master,
                             params->ana_max_level,
                             sbrConfigData->freqBandTable,
                             sbrConfigData->nSfb,
                             sbrHeaderData->sbr_noise_bands,
                             params->noiseFloorOffset,
                             params->useSpeechConfig))
  {
    
    return(1);
  }

   
  hEnv->encEnvData.noOfnoisebands = hEnv->TonCorr.sbrNoiseFloorEstimate.noNoiseBands;

   
  noiseBands[0] = hEnv->encEnvData.noOfnoisebands;
  noiseBands[1] = hEnv->encEnvData.noOfnoisebands;

   
  hEnv->encEnvData.sbr_invf_mode = (INVF_MODE)params->sbr_invf_mode;

   
  if (hEnv->encEnvData.sbr_invf_mode == INVF_SWITCHED) {

     
    hEnv->encEnvData.sbr_invf_mode = INVF_MID_LEVEL;
    hEnv->TonCorr.switchInverseFilt = TRUE;
  }
  else {

     
    hEnv->TonCorr.switchInverseFilt = FALSE;
  }


   
  tran_fc  = params->tran_fc;

  
  if (tran_fc == 0)
  {
        
    tran_fc = min (5000, getSbrStartFreqRAW (sbrHeaderData->sbr_start_frequency,64,sbrConfigData->sampleFreq));
  }


      
  tran_fc = (tran_fc*4*64/sbrConfigData->sampleFreq + 1)>>1;


     
  if(CreateExtractSbrEnvelope (sbrram,
                               chan,
                               &hEnv->sbrExtractEnvelope,

                               startIndex

                               ))
  {
    
    return(1);
  }

     
  if(CreateSbrCodeEnvelope (&hEnv->sbrCodeEnvelope,
                            sbrConfigData->nSfb,
                            params->deltaTAcrossFrames,
                            params->dF_edge_1stEnv,
                            params->dF_edge_incr))
  {
    
    return(1);
  }

     
  if(CreateSbrCodeEnvelope (&hEnv->sbrCodeNoiseFloor,
                            noiseBands,
                            params->deltaTAcrossFrames,
                            0, 0))
  {
    
    return(1);
  }

     
  if(InitSbrHuffmanTables (&hEnv->encEnvData,
                           &hEnv->sbrCodeEnvelope,
                           &hEnv->sbrCodeNoiseFloor,
                           sbrHeaderData->sbr_amp_res))
  {
   
   return(1);
  }

    
  CreateFrameInfoGenerator (&hEnv->SbrEnvFrame,
                            params->spread,
                            e,
                            params->stat,

                            hEnv->encEnvData.freq_res_fixfix);

     
  if(CreateSbrTransientDetector (sbrram,
                                 chan,
                                 &hEnv->sbrTransientDetector,

                                 sbrConfigData->sampleFreq,
                                 params->codecSettings.standardBitrate *
                                 params->codecSettings.nChannels,
                                 params->codecSettings.bitRate,
                                 params->tran_thr,
                                 params->tran_det_mode,
                                 tran_fc


                                 ))
  {
    
    return(1);
  }

   
  sbrConfigData->xposCtrlSwitch = params->sbr_xpos_ctrl;


     
    hEnv->encEnvData.noHarmonics = sbrConfigData->nSfb[HI];
    hEnv->encEnvData.syntheticCoding = sbrConfigData->detectMissingHarmonics;
    hEnv->encEnvData.addHarmonicFlag = 0;

  

  return (0);
}


/*****************************************************************************

 description:  initializes parameters and allocates memory
 returns:      error status

*****************************************************************************/
int
EnvOpen (SBRRam_t *sbrram,
         HANDLE_SBR_ENCODER   hEnvEnc,
         float *pCoreBuffer,
         sbrConfigurationPtr params,
         int      *coreBandWith)
{
  int ch;

  

  

  

   /* hEnvEnc->hEnvChannel[]
                  EnvChannel[]
               */
  
  for (ch=0; ch<MAX_CHANNELS; ch++)
  {
    
    hEnvEnc->hEnvChannel[ch] = calloc(1, sizeof(struct ENV_CHANNEL));
  }

     
  if ((params->codecSettings.nChannels < 1) || (params->codecSettings.nChannels > MAX_CHANNELS)){

    
    EnvClose (hEnvEnc);

    
    return(1);
  }


   
  hEnvEnc->sbrConfigData.freqBandTable[LO] =  sbrram->sbr_freqBandTableLO;

      
  memset(hEnvEnc->sbrConfigData.freqBandTable[LO],0,sizeof(unsigned char)*MAX_FREQ_COEFFS/2+1);

   
  hEnvEnc->sbrConfigData.freqBandTable[HI] =  sbrram->sbr_freqBandTableHI;

      
  memset(hEnvEnc->sbrConfigData.freqBandTable[HI],0,sizeof(unsigned char)*MAX_FREQ_COEFFS+1);

   
  hEnvEnc->sbrConfigData.v_k_master =  sbrram->sbr_v_k_master;

      
  memset(hEnvEnc->sbrConfigData.v_k_master,0,sizeof(unsigned char)*MAX_FREQ_COEFFS+1);

   
  if (hEnvEnc->CmonData.sbrBitbuf.isValid == 0) {

      
    CreateBitBuffer(&hEnvEnc->CmonData.sbrBitbuf,
                    (unsigned char*) hEnvEnc->sbrPayload,
                    sizeof(hEnvEnc->sbrPayload));
  }

   
  if (hEnvEnc->CmonData.sbrBitbufPrev.isValid == 0) {

      
    CreateBitBuffer(&hEnvEnc->CmonData.sbrBitbufPrev,
                    (unsigned char*) hEnvEnc->sbrPayloadPrevious,
                    sizeof(hEnvEnc->sbrPayload));
  }

   
  hEnvEnc->sbrConfigData.nChannels = params->codecSettings.nChannels;

    
  if(params->codecSettings.nChannels == 2)
  {
      
     hEnvEnc->sbrConfigData.stereoMode = params->stereoMode;
  }
  else
  {
      
     hEnvEnc->sbrConfigData.stereoMode = SBR_MONO;
  }






    
  if (params->codecSettings.sampleFreq <= 24000) {

     
    hEnvEnc->sbrHeaderData.sampleRateMode = DUAL_RATE;

      
    hEnvEnc->sbrConfigData.sampleFreq = 2 * params->codecSettings.sampleFreq;
  }
  else {

     
    hEnvEnc->sbrHeaderData.sampleRateMode = SINGLE_RATE;
    hEnvEnc->sbrConfigData.sampleFreq = params->codecSettings.sampleFreq;
  }


   
  hEnvEnc->sbrBitstreamData.CountSendHeaderData = 0;

   
  if (params->SendHeaderDataTime > 0 ) {

       
    hEnvEnc->sbrBitstreamData.NrSendHeaderData = (int)(params->SendHeaderDataTime * 0.001*
                                             hEnvEnc->sbrConfigData.sampleFreq / 2048
                                             );

      
    hEnvEnc->sbrBitstreamData.NrSendHeaderData = max(hEnvEnc->sbrBitstreamData.NrSendHeaderData,1);
  }
  else {

    
   hEnvEnc->sbrBitstreamData.NrSendHeaderData = 0;
  }

   
  hEnvEnc->sbrHeaderData.sbr_data_extra = params->sbr_data_extra;
  hEnvEnc->sbrBitstreamData.CRCActive = params->crcSbr;
  hEnvEnc->sbrBitstreamData.HeaderActive = 0;
  hEnvEnc->sbrHeaderData.sbr_start_frequency = params->startFreq;
  hEnvEnc->sbrHeaderData.sbr_stop_frequency  = params->stopFreq;
  hEnvEnc->sbrHeaderData.sbr_xover_band = 0;

      
    if (params->sbr_xpos_ctrl!= SBR_XPOS_CTRL_DEFAULT)
    {
        
       hEnvEnc->sbrHeaderData.sbr_data_extra = 1;
    }

    
   hEnvEnc->sbrHeaderData.protocol_version = SI_SBR_PROTOCOL_VERSION_ID;

    
   hEnvEnc->sbrHeaderData.sbr_amp_res = (AMP_RES)params->amp_res;



   
  hEnvEnc->sbrHeaderData.freqScale  = params->freqScale;
  hEnvEnc->sbrHeaderData.alterScale = params->alterScale;
  hEnvEnc->sbrHeaderData.sbr_noise_bands = params->sbr_noise_bands;
  hEnvEnc->sbrHeaderData.header_extra_1 = 0;

    
  if ((params->freqScale != SBR_FREQ_SCALE_DEFAULT) ||
      (params->alterScale != SBR_ALTER_SCALE_DEFAULT) ||
      (params->sbr_noise_bands != SBR_NOISE_BANDS_DEFAULT)) {

   
   hEnvEnc->sbrHeaderData.header_extra_1 = 1;
  }

   
  hEnvEnc->sbrHeaderData.sbr_limiter_bands = params->sbr_limiter_bands;
  hEnvEnc->sbrHeaderData.sbr_limiter_gains = params->sbr_limiter_gains;
  hEnvEnc->sbrHeaderData.sbr_interpol_freq = params->sbr_interpol_freq;
  hEnvEnc->sbrHeaderData.sbr_smoothing_length = params->sbr_smoothing_length;
  hEnvEnc->sbrHeaderData.header_extra_2 = 0;

    
  if ((params->sbr_limiter_bands != SBR_LIMITER_BANDS_DEFAULT) ||
      (params->sbr_limiter_gains != SBR_LIMITER_GAINS_DEFAULT) ||
      (params->sbr_interpol_freq != SBR_INTERPOL_FREQ_DEFAULT) ||
      (params->sbr_smoothing_length != SBR_SMOOTHING_LENGTH_DEFAULT)) {

     
     hEnvEnc->sbrHeaderData.header_extra_2 = 1;
  }

   
  hEnvEnc->sbrConfigData.detectMissingHarmonics    = params->detectMissingHarmonics;
  hEnvEnc->sbrConfigData.useParametricCoding       = params->parametricCoding;


     
  if(updateFreqBandTable(&hEnvEnc->sbrConfigData,
                         &hEnvEnc->sbrHeaderData,
                         QMF_CHANNELS)){
    
    EnvClose (hEnvEnc);

    
    return(1);
  }


   /* hEnvEnc->hEnvChannel[] */
    
  for ( ch = 0; ch < hEnvEnc->sbrConfigData.nChannels; ch++ ) {

     
    if(createEnvChannel(sbrram,
                        ch,
                        &hEnvEnc->sbrConfigData,
                        &hEnvEnc->sbrHeaderData,
                        hEnvEnc->hEnvChannel[ch],
                        params)){

      
      EnvClose (hEnvEnc);

      
      return(1);
    }
  }


   
  hEnvEnc->hPsEnc = NULL;


#ifndef MONO_ONLY

   
  if (params->usePs)
  {
       
    if(createQmfBank (sbrram, 1, &hEnvEnc->hEnvChannel[1]->sbrQmf)){
      
      return (1);
    }
       
    if(CreateExtractSbrEnvelope (sbrram,
                                 1,
                                 &hEnvEnc->hEnvChannel[1]->sbrExtractEnvelope,
                                 576
                                 )) {
      
      return(1);
    }

     
    hEnvEnc->hSynthesisQmfBank = calloc(1, sizeof(SBR_QMF_FILTER_BANK));

       
    if(CreateSynthesisQmfBank (sbrram, hEnvEnc->hSynthesisQmfBank)){

      
      DeleteSynthesisQmfBank ((HANDLE_SBR_QMF_FILTER_BANK*)&hEnvEnc->hSynthesisQmfBank);

      
      return 1;
    }

     
    hEnvEnc->hPsEnc  = calloc(1, sizeof(struct PS_ENC));

       
    if(CreatePsEnc (sbrram,
                    hEnvEnc->hPsEnc,
                    params->psMode)){

      
      DeletePsEnc(&hEnvEnc->hPsEnc);

      
      return 1;
    }
  }
#endif /* #ifndef MONO_ONLY */

   
  hEnvEnc->CmonData.sbrNumChannels  = hEnvEnc->sbrConfigData.nChannels;

   
  hEnvEnc->sbrPayloadSize = 0;

   
  *coreBandWith = hEnvEnc->sbrConfigData.xOverFreq;

  

  return 0;
}

