/*
  fast aac coder interface library functions
*/

#ifndef _aacenc_h_
#define _aacenc_h_

struct AAC_ENCODER {
//#ifdef _FFTW3
  FFTWFContext_t *fftctx;
//#endif

  aacplusEncConfiguration *config;
  AACRam_t *aac_ram;
  SBRRam_t *sbr_ram;

  ELEMENT_INFO elInfo;

  QC_STATE qcKernel;
  QC_OUT   qcOut;

  PSY_OUT    psyOut;
  PSY_KERNEL psyKernel;

  struct BITSTREAMENCODER_INIT bseInit;

  struct STEREO_PREPRO stereoPrePro;

  struct BIT_BUF  bitStream;
  HANDLE_BIT_BUF  hBitStream;

  /* lifetime vars */
  int downmix;
  int downmixFac;
  int dualMono;
  int bandwidth90dB;
};


/*-----------------------------------------------------------------------------

     functionname: AacInitDefaultConfig
     description:  gives reasonable default configuration
     returns:      ---

 ------------------------------------------------------------------------------*/
void AacInitDefaultConfig(aacplusEncConfiguration *config);

/*---------------------------------------------------------------------------

    functionname:AacEncOpen
    description: allocate and initialize a new encoder instance
    returns:     AACENC_OK if success

  ---------------------------------------------------------------------------*/

int  AacEncOpen
(  struct AAC_ENCODER*      hAacEnc,       /* pointer to an encoder handle, initialized on return */
   aacplusEncConfiguration    *config          /* pre-initialized config struct */
);

int AacEncEncode(struct AAC_ENCODER  *hAacEnc,
                 float               *timeSignal,
                 unsigned int        timeInStride,
                 const unsigned char *ancBytes,      /*!< pointer to ancillary data bytes */
                 unsigned int        *numAncBytes,   /*!< number of ancillary Data Bytes, send as fill element  */
                 unsigned int        *outBytes,      /*!< pointer to output buffer            */
                 int                 *numOutBytes    /*!< number of bytes in output buffer */
                 );


/*---------------------------------------------------------------------------

    functionname:AacEncClose
    description: deallocate an encoder instance

  ---------------------------------------------------------------------------*/

void AacEncClose (struct AAC_ENCODER* hAacEnc); /* an encoder handle */

#endif /* _aacenc_h_ */
