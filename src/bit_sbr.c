/*
  SBR bit writing routines
*/
#include <stdlib.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */



#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif


typedef enum {
  SBR_ID_SCE = 1,
  SBR_ID_CPE
} SBR_ELEMENT_TYPE;


static int encodeSbrData (HANDLE_SBR_ENV_DATA   sbrEnvDataLeft,
                          HANDLE_SBR_ENV_DATA   sbrEnvDataRight,
                          HANDLE_COMMON_DATA    cmonData,
                          SBR_ELEMENT_TYPE      sbrElem,
                          int                   sampleRateMode,
                          int                   data_extra,
                          struct PS_ENC*     h_ps_e,
                          int                   bHeaderActive,
                          int                   coupling);

static int encodeSbrHeader (HANDLE_SBR_HEADER_DATA     sbrHeaderData,
                            HANDLE_SBR_BITSTREAM_DATA  sbrBitstreamData,
                            HANDLE_COMMON_DATA         cmonData,
                            SBR_ELEMENT_TYPE           sbrElem);


static int encodeSbrHeaderData (HANDLE_SBR_HEADER_DATA sbrHeaderData,
                                HANDLE_BIT_BUF         hBitStream,
                                SBR_ELEMENT_TYPE       sbrElem);

static int encodeSbrSingleChannelElement (HANDLE_SBR_ENV_DATA    sbrEnvData,
                                          HANDLE_BIT_BUF         hBitStream,
                                          int                    data_extra);


static int encodeSbrChannelPairElement (HANDLE_SBR_ENV_DATA  sbrEnvDataLeft,
                                        HANDLE_SBR_ENV_DATA  sbrEnvDataRight,
                                        HANDLE_BIT_BUF        hBitStream,
                                        int                   data_extra,
                                        int                   coupling);


static int encodeSbrGrid (HANDLE_SBR_ENV_DATA   sbrEnvData,
                          HANDLE_BIT_BUF        hBitStream);

static int encodeSbrDtdf (HANDLE_SBR_ENV_DATA   sbrEnvData,
                          HANDLE_BIT_BUF        hBitStream);

static int writeNoiseLevelData (HANDLE_SBR_ENV_DATA   sbrEnvData,
                                HANDLE_BIT_BUF        hBitStream,
                                int                   coupling);

static int writeEnvelopeData (HANDLE_SBR_ENV_DATA    sbrEnvData,
                              HANDLE_BIT_BUF         hBitStream,
                              int                    coupling);

static int writeSyntheticCodingData (HANDLE_SBR_ENV_DATA  sbrEnvData,
                                     HANDLE_BIT_BUF       hBitStream);


static void encodeExtendedData (HANDLE_SBR_ENV_DATA  sbrEnvDataLeft,
                                HANDLE_SBR_ENV_DATA  sbrEnvDataRight,
                                struct PS_ENC*    h_ps_e,
                                int                  bHeaderActive,
                                HANDLE_BIT_BUF       hBitStreamPrev,
                                int*                 sbrHdrBits,
                                HANDLE_BIT_BUF       hBitStream,
                                int* payloadBitsReturn);



static int
getSbrExtendedDataSize (HANDLE_SBR_ENV_DATA sbrEnvDataLeft,
                        HANDLE_SBR_ENV_DATA sbrEnvDataRight,
                        struct PS_ENC*      h_ps_e,
                        int                 bHeaderActive
);


/*****************************************************************************

    description:  writes SBR single channel data element
    returns:      number of bits written

*****************************************************************************/
int
WriteEnvSingleChannelElement(HANDLE_SBR_HEADER_DATA     sbrHeaderData,
                             HANDLE_SBR_BITSTREAM_DATA  sbrBitstreamData,
                             HANDLE_SBR_ENV_DATA        sbrEnvData,
                             struct PS_ENC*             h_ps_e,
                             HANDLE_COMMON_DATA         cmonData)

{
  int payloadBits = 0;

  

   /* counting previous operation */

   
  cmonData->sbrHdrBits  = 0;
  cmonData->sbrDataBits = 0;
  cmonData->sbrCrcLen   = 0;

  /* write pure sbr data */
  
  if (sbrEnvData != NULL) {

    /* write header */
      
    payloadBits += encodeSbrHeader (sbrHeaderData,
                                    sbrBitstreamData,
                                    cmonData,
                                    SBR_ID_SCE);


    /* write data */
      
    payloadBits += encodeSbrData (sbrEnvData,
                                  NULL,
                                  cmonData,
                                  SBR_ID_SCE,
                                  sbrHeaderData->sampleRateMode,
                                  sbrHeaderData->sbr_data_extra,
                                  h_ps_e,
                                  sbrBitstreamData->HeaderActive,
                                  0);


  }

  

  return payloadBits;
}

/*****************************************************************************

    description:  writes SBR channel pair data element
    returns:      number of bits written

*****************************************************************************/
int
WriteEnvChannelPairElement (HANDLE_SBR_HEADER_DATA     sbrHeaderData,
                            HANDLE_SBR_BITSTREAM_DATA  sbrBitstreamData,
                            HANDLE_SBR_ENV_DATA        sbrEnvDataLeft,
                            HANDLE_SBR_ENV_DATA        sbrEnvDataRight,
                            HANDLE_COMMON_DATA         cmonData)

{
  int payloadBits = 0;

  

   /* counting previous operation */

   
  cmonData->sbrHdrBits  = 0;
  cmonData->sbrDataBits = 0;
  cmonData->sbrCrcLen   = 0;

  /* write pure sbr data */
   
  if ((sbrEnvDataLeft != NULL) && (sbrEnvDataRight != NULL)) {

    /* write header */
     
    payloadBits += encodeSbrHeader (sbrHeaderData,
                                    sbrBitstreamData,
                                    cmonData,
                                    SBR_ID_CPE);


    /* write data */
      
    payloadBits += encodeSbrData (sbrEnvDataLeft,
                                  sbrEnvDataRight,
                                  cmonData,
                                  SBR_ID_CPE,
                                  sbrHeaderData->sampleRateMode,
                                  sbrHeaderData->sbr_data_extra,
                                  NULL,
                                  0,
                                  sbrHeaderData->coupling);


  }
  

  return payloadBits;
}

int
CountSbrChannelPairElement (HANDLE_SBR_HEADER_DATA     sbrHeaderData,
                            HANDLE_SBR_BITSTREAM_DATA  sbrBitstreamData,
                            HANDLE_SBR_ENV_DATA        sbrEnvDataLeft,
                            HANDLE_SBR_ENV_DATA        sbrEnvDataRight,
                            HANDLE_COMMON_DATA         cmonData)

{
  int payloadBits;
  struct BIT_BUF bitBufTmp = cmonData->sbrBitbuf;

  

    /* counting previous operation */

  
  payloadBits = WriteEnvChannelPairElement(sbrHeaderData,
                                           sbrBitstreamData,
                                           sbrEnvDataLeft,
                                           sbrEnvDataRight,
                                           cmonData);

   
  cmonData->sbrBitbuf = bitBufTmp;

  

  return payloadBits;
}

/*****************************************************************************

    description:  encodes SBR Header information
    returns:      number of bits written

*****************************************************************************/
static int
encodeSbrHeader (HANDLE_SBR_HEADER_DATA     sbrHeaderData,
                 HANDLE_SBR_BITSTREAM_DATA  sbrBitstreamData,
                 HANDLE_COMMON_DATA         cmonData,
                 SBR_ELEMENT_TYPE           sbrElem)

{
  int payloadBits = 0;

  

   /* counting previous operation */

   
  if (sbrBitstreamData->CRCActive) {

     
    cmonData->sbrCrcLen = 1;
  }
  else {

     
    cmonData->sbrCrcLen = 0;
  }


  
  if (sbrBitstreamData->HeaderActive) {

       
    payloadBits += WriteBits (&cmonData->sbrBitbuf, 1, 1);

       
    payloadBits += encodeSbrHeaderData (sbrHeaderData,
                                        &cmonData->sbrBitbuf,
                                        sbrElem);
  }
  else {

       
    payloadBits += WriteBits (&cmonData->sbrBitbuf, 0, 1);
  }


   
  cmonData->sbrHdrBits = payloadBits;

  

  return payloadBits;
}



/*****************************************************************************

    functionname: encodeSbrHeaderData
    description:  writes sbr_header()
    returns:      number of bits written

*****************************************************************************/
static int
encodeSbrHeaderData (HANDLE_SBR_HEADER_DATA sbrHeaderData,
                     HANDLE_BIT_BUF hBitStream,
                     SBR_ELEMENT_TYPE sbrElem)

{
  int payloadBits = 0;

  

   /* counting previous operation */

  
  if (sbrHeaderData != NULL) {

      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_amp_res,
                              SI_SBR_AMP_RES_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_start_frequency,
                              SI_SBR_START_FREQ_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_stop_frequency,
                              SI_SBR_STOP_FREQ_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_xover_band,
                              SI_SBR_XOVER_BAND_BITS);

     
    payloadBits += WriteBits (hBitStream, 0,
                              SI_SBR_RESERVED_BITS);



      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->header_extra_1,
                              SI_SBR_HEADER_EXTRA_1_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrHeaderData->header_extra_2,
                              SI_SBR_HEADER_EXTRA_2_BITS);


     
    if (sbrHeaderData->header_extra_1) {

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->freqScale,
                                SI_SBR_FREQ_SCALE_BITS);

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->alterScale,
                                SI_SBR_ALTER_SCALE_BITS);

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_noise_bands,
                                SI_SBR_NOISE_BANDS_BITS);
    } /* sbrHeaderData->header_extra_1 */

     
    if (sbrHeaderData->header_extra_2) {

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_limiter_bands,
                                SI_SBR_LIMITER_BANDS_BITS);

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_limiter_gains,
                                SI_SBR_LIMITER_GAINS_BITS);

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_interpol_freq,
                                SI_SBR_INTERPOL_FREQ_BITS);

        
      payloadBits += WriteBits (hBitStream, sbrHeaderData->sbr_smoothing_length,
                                SI_SBR_SMOOTHING_LENGTH_BITS);

    } /* sbrHeaderData->header_extra_2 */
  } /* sbrHeaderData != NULL */

  

  return payloadBits;
}


/*****************************************************************************

    description:  encodes sbr Data information
    returns:      number of bits written

*****************************************************************************/
static int
encodeSbrData (HANDLE_SBR_ENV_DATA   sbrEnvDataLeft,
               HANDLE_SBR_ENV_DATA   sbrEnvDataRight,
               HANDLE_COMMON_DATA    cmonData,
               SBR_ELEMENT_TYPE      sbrElem,
               int                   sampleRateMode,
               int                   data_extra,
               struct PS_ENC*        h_ps_e,
               int                   bHeaderActive,
               int                   coupling)
{
  int payloadBits = 0;

  

   /* counting previous operation */

  
  switch (sbrElem) {
  case SBR_ID_SCE:

       
    payloadBits += encodeSbrSingleChannelElement (sbrEnvDataLeft,
                                                  &cmonData->sbrBitbuf,
                                                  data_extra);

      
    encodeExtendedData(sbrEnvDataLeft,
                       NULL,
                       h_ps_e,
                       bHeaderActive,
                       &cmonData->sbrBitbufPrev,
                       &cmonData->sbrHdrBits,
                       &cmonData->sbrBitbuf,
                       &payloadBits);


    break;
  case SBR_ID_CPE:

       
    payloadBits += encodeSbrChannelPairElement (sbrEnvDataLeft, sbrEnvDataRight, &cmonData->sbrBitbuf,
                                                data_extra, coupling);


      
    encodeExtendedData(sbrEnvDataLeft,
                       sbrEnvDataRight,
                       NULL,
                       0,
                       NULL,
                       0,
                       &cmonData->sbrBitbuf,
                       &payloadBits);

    break;
  default:
    /* we never should apply SBR to any other element type */
    assert (0);
  }


   
  cmonData->sbrDataBits = payloadBits;

  

  return payloadBits;
}



/*****************************************************************************

    description:  encodes sbr SCE information
    returns:      number of bits written

*****************************************************************************/
static int
encodeSbrSingleChannelElement (HANDLE_SBR_ENV_DATA   sbrEnvData,
                               HANDLE_BIT_BUF        hBitStream,
                               int                   data_extra)

{
  int payloadBits = 0;

  

   /* counting previous operation */

     
    payloadBits += WriteBits (hBitStream, 0, 1); /* no reserved bits */



   
  payloadBits += encodeSbrGrid (sbrEnvData, hBitStream);

   
  payloadBits += encodeSbrDtdf (sbrEnvData, hBitStream);

  {
    int i;

     /* sbrEnvData->sbr_invf_mode_vec[] */
     
    for (i = 0; i < sbrEnvData->noOfnoisebands; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvData->sbr_invf_mode_vec[i], SI_SBR_INVF_MODE_BITS);
    }
  }


   
  payloadBits += writeEnvelopeData (sbrEnvData, hBitStream, 0);

   
  payloadBits += writeNoiseLevelData (sbrEnvData, hBitStream, 0);

   
  payloadBits += writeSyntheticCodingData (sbrEnvData,hBitStream);


  

  return payloadBits;
}


/*****************************************************************************

    description:  encodes sbr CPE information

*****************************************************************************/
static int
encodeSbrChannelPairElement (HANDLE_SBR_ENV_DATA   sbrEnvDataLeft,
                             HANDLE_SBR_ENV_DATA   sbrEnvDataRight,
                             HANDLE_BIT_BUF        hBitStream,
                             int                   data_extra,
                             int                   coupling)
{
  int payloadBits = 0;
  int i = 0;

  

   /* counting previous operation */

     
    payloadBits += WriteBits (hBitStream, 0, 1); /* no reserved bits */

   
  payloadBits += WriteBits (hBitStream, coupling, SI_SBR_COUPLING_BITS);

  
  if (coupling) {

     
    payloadBits += encodeSbrGrid (sbrEnvDataLeft, hBitStream);

     
    payloadBits += encodeSbrDtdf (sbrEnvDataLeft, hBitStream);

     
    payloadBits += encodeSbrDtdf (sbrEnvDataRight, hBitStream);

     /* sbrEnvDataLeft->sbr_invf_mode_vec[] */
     
    for (i = 0; i < sbrEnvDataLeft->noOfnoisebands; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvDataLeft->sbr_invf_mode_vec[i], SI_SBR_INVF_MODE_BITS);
    }




     
    payloadBits += writeEnvelopeData  (sbrEnvDataLeft,  hBitStream,1);

     
    payloadBits += writeNoiseLevelData (sbrEnvDataLeft,  hBitStream,1);

     
    payloadBits += writeEnvelopeData  (sbrEnvDataRight, hBitStream,1);

     
    payloadBits += writeNoiseLevelData (sbrEnvDataRight, hBitStream,1);




     
    payloadBits += writeSyntheticCodingData (sbrEnvDataLeft,hBitStream);

     
    payloadBits += writeSyntheticCodingData (sbrEnvDataRight,hBitStream);


  } else { /* no coupling */

     
    payloadBits += encodeSbrGrid (sbrEnvDataLeft,  hBitStream);

     
    payloadBits += encodeSbrGrid (sbrEnvDataRight, hBitStream);

     
    payloadBits += encodeSbrDtdf (sbrEnvDataLeft,  hBitStream);

     
    payloadBits += encodeSbrDtdf (sbrEnvDataRight, hBitStream);

     /* sbrEnvDataLeft->sbr_invf_mode_vec[] */
     
    for (i = 0; i < sbrEnvDataLeft->noOfnoisebands; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvDataLeft->sbr_invf_mode_vec[i],
                                SI_SBR_INVF_MODE_BITS);
    }

     /* sbrEnvDataRight->sbr_invf_mode_vec[] */
     
    for (i = 0; i < sbrEnvDataRight->noOfnoisebands; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvDataRight->sbr_invf_mode_vec[i],
                                SI_SBR_INVF_MODE_BITS);
    }



     
    payloadBits += writeEnvelopeData  (sbrEnvDataLeft,  hBitStream,0);

     
    payloadBits += writeEnvelopeData  (sbrEnvDataRight, hBitStream,0);

     
    payloadBits += writeNoiseLevelData (sbrEnvDataLeft,  hBitStream,0);

     
    payloadBits += writeNoiseLevelData (sbrEnvDataRight, hBitStream,0);


     
    payloadBits += writeSyntheticCodingData (sbrEnvDataLeft,hBitStream);

     
    payloadBits += writeSyntheticCodingData (sbrEnvDataRight,hBitStream);

  } /* coupling */



  

  return payloadBits;
}

static int ceil_ln2(int x)
{
  int tmp=-1;

  

   /* counting previous operation */

  
  while((1<<++tmp) < x);
  

  

  return(tmp);
}


/*****************************************************************************

    description:  Encode SBR grid information

*****************************************************************************/
static int
encodeSbrGrid (HANDLE_SBR_ENV_DATA sbrEnvData, HANDLE_BIT_BUF hBitStream)
{
  int payloadBits = 0;
  int i, temp;



  

   /* counting previous operations */

    
  payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->frameClass,
                            SBR_CLA_BITS);

   
  switch (sbrEnvData->hSbrBSGrid->frameClass) {
  case FIXFIX:

     
    temp = ceil_ln2(sbrEnvData->hSbrBSGrid->bs_num_env);

     
    payloadBits += WriteBits (hBitStream, temp, SBR_ENV_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->freq_res_fixfix, SBR_RES_BITS);
    break;

  case FIXVAR:
  case VARFIX:

      
    if (sbrEnvData->hSbrBSGrid->frameClass == FIXVAR)
    {
       
      temp = sbrEnvData->hSbrBSGrid->bs_abs_bord - 16;
    }
    else
    {
      
      temp = sbrEnvData->hSbrBSGrid->bs_abs_bord;
    }

     
    payloadBits += WriteBits (hBitStream, temp, SBR_ABS_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->n, SBR_NUM_BITS);

     /* sbrEnvData->hSbrBSGrid->bs_rel_bord[] */
     
    for (i = 0; i < sbrEnvData->hSbrBSGrid->n; i++) {

       
      temp = (sbrEnvData->hSbrBSGrid->bs_rel_bord[i] - 2) >> 1;

       
      payloadBits += WriteBits (hBitStream, temp, SBR_REL_BITS);
    }

      
    temp = ceil_ln2(sbrEnvData->hSbrBSGrid->n + 2);

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->p, temp);

     /* sbrEnvData->hSbrBSGrid->v_f[] */
      
    for (i = 0; i < sbrEnvData->hSbrBSGrid->n + 1; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->v_f[i],
                                SBR_RES_BITS);
    }
    break;

  case VARVAR:

    
    temp = sbrEnvData->hSbrBSGrid->bs_abs_bord_0;

     
    payloadBits += WriteBits (hBitStream, temp, SBR_ABS_BITS);

     
    temp = sbrEnvData->hSbrBSGrid->bs_abs_bord_1 - 16;

     
    payloadBits += WriteBits (hBitStream, temp, SBR_ABS_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->bs_num_rel_0, SBR_NUM_BITS);

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->bs_num_rel_1, SBR_NUM_BITS);

     /* sbrEnvData->hSbrBSGrid->bs_rel_bord_0[] */
     
    for (i = 0; i < sbrEnvData->hSbrBSGrid->bs_num_rel_0; i++) {

       
      temp = (sbrEnvData->hSbrBSGrid->bs_rel_bord_0[i] - 2) >> 1;

       
      payloadBits += WriteBits (hBitStream, temp, SBR_REL_BITS);
    }

     /* sbrEnvData->hSbrBSGrid->bs_rel_bord_1[] */
     
    for (i = 0; i < sbrEnvData->hSbrBSGrid->bs_num_rel_1; i++) {

       
      temp = (sbrEnvData->hSbrBSGrid->bs_rel_bord_1[i] - 2) >> 1;

       
      payloadBits += WriteBits (hBitStream, temp, SBR_REL_BITS);
    }

      
    temp = ceil_ln2(sbrEnvData->hSbrBSGrid->bs_num_rel_0 +
                             sbrEnvData->hSbrBSGrid->bs_num_rel_1 + 2);

      
    payloadBits +=  WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->p, temp);

     
    temp = sbrEnvData->hSbrBSGrid->bs_num_rel_0 +
           sbrEnvData->hSbrBSGrid->bs_num_rel_1 + 1;

     /* sbrEnvData->hSbrBSGrid->v_fLR[] */
     
    for (i = 0; i < temp; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvData->hSbrBSGrid->v_fLR[i],
                                SBR_RES_BITS);
    }
    break;
  }

  

  return payloadBits;
}


/*****************************************************************************

    description:  writes bits that describes the direction of the envelopes of a frame
    returns:      number of bits written

*****************************************************************************/
static int
encodeSbrDtdf (HANDLE_SBR_ENV_DATA sbrEnvData, HANDLE_BIT_BUF hBitStream)
{
  int i, payloadBits = 0, noOfNoiseEnvelopes;

  

   /* counting previous operations */

     
  noOfNoiseEnvelopes = (sbrEnvData->noOfEnvelopes > 1) ? 2 : 1;

   /* sbrEnvData->domain_vec[] */
   
  for (i = 0; i < sbrEnvData->noOfEnvelopes; ++i) {

     
    payloadBits += WriteBits (hBitStream, sbrEnvData->domain_vec[i], SBR_DIR_BITS);
  }

   /* sbrEnvData->domain_vec_noise[] */
  
  for (i = 0; i < noOfNoiseEnvelopes; ++i) {

     
    payloadBits +=  WriteBits (hBitStream, sbrEnvData->domain_vec_noise[i], SBR_DIR_BITS);
  }

  

  return payloadBits;
}


/*****************************************************************************

    description:  writes bits corresponding to the noise-floor-level
    returns:      number of bits written

*****************************************************************************/
static int
writeNoiseLevelData (HANDLE_SBR_ENV_DATA sbrEnvData, HANDLE_BIT_BUF hBitStream, int coupling)
{
  int j, i, payloadBits = 0;
  int nNoiseEnvelopes = ((sbrEnvData->noOfEnvelopes > 1) ? 2 : 1);

  

      /* .. > .. ? */ 

   /* sbrEnvData->domain_vec_noise[]
                  sbrEnvData->sbr_noise_levels[i * sbrEnvData->noOfnoisebands]
               */
  
  for (i = 0; i < nNoiseEnvelopes; i++) {

    
    switch (sbrEnvData->domain_vec_noise[i]) {
    case FREQ:

        
      if (coupling && sbrEnvData->balance) {

          
        payloadBits += WriteBits (hBitStream,
                                  sbrEnvData->sbr_noise_levels[i * sbrEnvData->noOfnoisebands],
                                  sbrEnvData->si_sbr_start_noise_bits_balance);
      } else {

          
        payloadBits += WriteBits (hBitStream,
                                  sbrEnvData->sbr_noise_levels[i * sbrEnvData->noOfnoisebands],
                                  sbrEnvData->si_sbr_start_noise_bits);
      }

       /* sbrEnvData->sbr_noise_levels[] */
         
      for (j = 1 + i * sbrEnvData->noOfnoisebands; j < (sbrEnvData->noOfnoisebands * (1 + i)); j++) {
        
        if (coupling) {

           
          if (sbrEnvData->balance) {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableNoiseBalanceFreqC[sbrEnvData->sbr_noise_levels[j] +
                                                                            CODE_BOOK_SCF_LAV_BALANCE11],
                                      sbrEnvData->hufftableNoiseBalanceFreqL[sbrEnvData->sbr_noise_levels[j] +
                                                                            CODE_BOOK_SCF_LAV_BALANCE11]);
          } else {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableNoiseLevelFreqC[sbrEnvData->sbr_noise_levels[j] +
                                                                          CODE_BOOK_SCF_LAV11],
                                      sbrEnvData->hufftableNoiseLevelFreqL[sbrEnvData->sbr_noise_levels[j] +
                                                                          CODE_BOOK_SCF_LAV11]);
          }
        } else {

            
          payloadBits += WriteBits (hBitStream,
                                    sbrEnvData->hufftableNoiseFreqC[sbrEnvData->sbr_noise_levels[j] +
                                                                   CODE_BOOK_SCF_LAV11],
                                    sbrEnvData->hufftableNoiseFreqL[sbrEnvData->sbr_noise_levels[j] +
                                                                   CODE_BOOK_SCF_LAV11]);
        }
      }
      break;

    case TIME:

       /* sbrEnvData->sbr_noise_levels[] */
         
      for (j = i * sbrEnvData->noOfnoisebands; j < (sbrEnvData->noOfnoisebands * (1 + i)); j++) {

        
        if (coupling) {

           
          if (sbrEnvData->balance) {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableNoiseBalanceTimeC[sbrEnvData->sbr_noise_levels[j] +
                                                                            CODE_BOOK_SCF_LAV_BALANCE11],
                                      sbrEnvData->hufftableNoiseBalanceTimeL[sbrEnvData->sbr_noise_levels[j] +
                                                                            CODE_BOOK_SCF_LAV_BALANCE11]);
          } else {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableNoiseLevelTimeC[sbrEnvData->sbr_noise_levels[j] +
                                                                          CODE_BOOK_SCF_LAV11],
                                      sbrEnvData->hufftableNoiseLevelTimeL[sbrEnvData->sbr_noise_levels[j] +
                                                                          CODE_BOOK_SCF_LAV11]);
          }
        } else {

            
          payloadBits += WriteBits (hBitStream,
                                    sbrEnvData->hufftableNoiseLevelTimeC[sbrEnvData->sbr_noise_levels[j] +
                                                                        CODE_BOOK_SCF_LAV11],
                                    sbrEnvData->hufftableNoiseLevelTimeL[sbrEnvData->sbr_noise_levels[j] +
                                                                        CODE_BOOK_SCF_LAV11]);
        }
      }
      break;
    }
  }

  

  return payloadBits;
}


/*****************************************************************************

    description:  writes bits corresponding to the envelope data
    returns:      number of bits written

*****************************************************************************/
static int
writeEnvelopeData (HANDLE_SBR_ENV_DATA sbrEnvData, HANDLE_BIT_BUF hBitStream, int coupling)
{
  int payloadBits = 0, j, i, delta;

  

   /* sbrEnvData->domain_vec[]
                  sbrEnvData->ienvelope[][]
                  sbrEnvData->noScfBands[]
                */
   
  for (j = 0; j < sbrEnvData->noOfEnvelopes; j++) { /* loop over all envelopes */

     
    if (sbrEnvData->domain_vec[j] == FREQ) {

        
      if (coupling && sbrEnvData->balance) {

          
        payloadBits += WriteBits (hBitStream, sbrEnvData->ienvelope[j][0], sbrEnvData->si_sbr_start_env_bits_balance);
      } else {

          
        payloadBits += WriteBits (hBitStream, sbrEnvData->ienvelope[j][0], sbrEnvData->si_sbr_start_env_bits);
      }
    }

    
    for (i = 1 - sbrEnvData->domain_vec[j]; i < sbrEnvData->noScfBands[j]; i++) {

      
      delta = sbrEnvData->ienvelope[j][i];

      if (coupling && sbrEnvData->balance) {
        assert (abs (delta) <= sbrEnvData->codeBookScfLavBalance);
      } else {
        assert (abs (delta) <= sbrEnvData->codeBookScfLav);
      }

      
      if (coupling) {

         
        if (sbrEnvData->balance) {

          
          if (sbrEnvData->domain_vec[j]) {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableBalanceTimeC[delta + sbrEnvData->codeBookScfLavBalance],
                                      sbrEnvData->hufftableBalanceTimeL[delta + sbrEnvData->codeBookScfLavBalance]);
          } else {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableBalanceFreqC[delta + sbrEnvData->codeBookScfLavBalance],
                                      sbrEnvData->hufftableBalanceFreqL[delta + sbrEnvData->codeBookScfLavBalance]);
          }
        } else {

          
          if (sbrEnvData->domain_vec[j]) {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableLevelTimeC[delta + sbrEnvData->codeBookScfLav],
                                      sbrEnvData->hufftableLevelTimeL[delta + sbrEnvData->codeBookScfLav]);
          } else {

              
            payloadBits += WriteBits (hBitStream,
                                      sbrEnvData->hufftableLevelFreqC[delta + sbrEnvData->codeBookScfLav],
                                      sbrEnvData->hufftableLevelFreqL[delta + sbrEnvData->codeBookScfLav]);
          }
        }
      } else {

        
        if (sbrEnvData->domain_vec[j]) {

            
          payloadBits += WriteBits (hBitStream,
                                    sbrEnvData->hufftableTimeC[delta + sbrEnvData->codeBookScfLav],
                                    sbrEnvData->hufftableTimeL[delta + sbrEnvData->codeBookScfLav]);
        } else {

            
          payloadBits += WriteBits (hBitStream,
                                    sbrEnvData->hufftableFreqC[delta + sbrEnvData->codeBookScfLav],
                                    sbrEnvData->hufftableFreqL[delta + sbrEnvData->codeBookScfLav]);
        }
      }
    }
  }

  

  return payloadBits;
}


/*****************************************************************************

    description:  writes bits corresponding to the extended data
    returns:      number of bits written

*****************************************************************************/
static void
encodeExtendedData (HANDLE_SBR_ENV_DATA  sbrEnvDataLeft,
                    HANDLE_SBR_ENV_DATA  sbrEnvDataRight,
                    struct PS_ENC*    h_ps_e,
                    int                  bHeaderActive,
                    HANDLE_BIT_BUF       hBitStreamPrev,
                    int                  *sbrHdrBits,
                    HANDLE_BIT_BUF       hBitStream,
                    int* payloadBitsReturn)
{
  int extDataSize;
  int payloadBitsIn = *payloadBitsReturn;
  int payloadBits = 0;

  

   /* counting previous operations */

  
  extDataSize = getSbrExtendedDataSize(sbrEnvDataLeft,
                                       sbrEnvDataRight,
                                       h_ps_e,
                                       bHeaderActive
                                       );



  
  if (extDataSize != 0) {

#ifndef MONO_ONLY

      
    if (h_ps_e && AppendPsBS (h_ps_e, NULL, NULL, 0)) {

       
      *payloadBitsReturn = AppendPsBS (h_ps_e, hBitStream, hBitStreamPrev, sbrHdrBits);
    }
    else {

#endif /*#ifndef MONO_ONLY */

      int maxExtSize = (1<<SI_SBR_EXTENSION_SIZE_BITS) - 1;

        /* counting previous operation */

       
      payloadBits += WriteBits (hBitStream, 1, SI_SBR_EXTENDED_DATA_BITS);

      assert(extDataSize <= SBR_EXTENDED_DATA_MAX_CNT);

       
      if (extDataSize < maxExtSize) {

         
        payloadBits += WriteBits (hBitStream, extDataSize, SI_SBR_EXTENSION_SIZE_BITS);
      } else {

         
        payloadBits += WriteBits (hBitStream, maxExtSize, SI_SBR_EXTENSION_SIZE_BITS);

         
        payloadBits += WriteBits (hBitStream, extDataSize - maxExtSize, SI_SBR_EXTENSION_ESC_COUNT_BITS);
      }

       
      *payloadBitsReturn = payloadBits + payloadBitsIn;

#ifndef MONO_ONLY
    }
#endif /*#ifndef MONO_ONLY */

  }
  else {

     
    payloadBits += WriteBits (hBitStream, 0, SI_SBR_EXTENDED_DATA_BITS);

     
    *payloadBitsReturn = payloadBits + payloadBitsIn;
  }

  
}


/*****************************************************************************

    description:  writes bits corresponding to the "synthetic-coding"-extension
    returns:      number of bits written

*****************************************************************************/
static int writeSyntheticCodingData (HANDLE_SBR_ENV_DATA  sbrEnvData,
                                     HANDLE_BIT_BUF       hBitStream)

{
  int i;
  int payloadBits = 0;

  

   /* counting previous operation */

      
    payloadBits += WriteBits (hBitStream, sbrEnvData->addHarmonicFlag, 1);

  
  if (sbrEnvData->addHarmonicFlag) {

     /* sbrEnvData->addHarmonic[] */
    
    for (i = 0; i < sbrEnvData->noHarmonics; i++) {

       
      payloadBits += WriteBits (hBitStream, sbrEnvData->addHarmonic[i], 1);
    }
  }

  

  return payloadBits;
}


/*****************************************************************************

    description:  counts the number of bits needed

    returns:      number of bits needed for the extended data

*****************************************************************************/
static int
getSbrExtendedDataSize (HANDLE_SBR_ENV_DATA sbrEnvDataLeft,
                        HANDLE_SBR_ENV_DATA sbrEnvDataRight,
                        struct PS_ENC*      h_ps_e,
                        int                 bHeaderActive)

{
  int extDataBits = 0;

  

   /* counting previous operation */

#ifndef MONO_ONLY

  
  if (h_ps_e) {

     
    extDataBits += WritePsData(h_ps_e, bHeaderActive);
  }

#endif /*#ifndef MONO_ONLY */

  /*
    no extended data
  */

  
  if (extDataBits != 0)
  {
    
    extDataBits += SI_SBR_EXTENSION_ID_BITS;
  }

    /* counting post-operation */

  

  return (extDataBits+7) >> 3;
}


