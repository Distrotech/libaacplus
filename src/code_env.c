/*
  DPCM envelope coding
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif


/*****************************************************************************

 functionname: InitSbrHuffmanTables
 description:  initializes Huffman Tables dependent on chosen amp_res
 returns:      error handle

*****************************************************************************/
int
InitSbrHuffmanTables (HANDLE_SBR_ENV_DATA       sbrEnvData,
                      HANDLE_SBR_CODE_ENVELOPE  henv,
                      HANDLE_SBR_CODE_ENVELOPE  hnoise,
                      AMP_RES                   amp_res)
{
  

   
  if ( (!henv)  ||  (!hnoise)  || (!sbrEnvData) ) {
    

    return (1);
  }

   
  sbrEnvData->init_sbr_amp_res = amp_res;

  
  switch (amp_res) {
  case  SBR_AMP_RES_3_0:



     
    sbrEnvData->hufftableLevelTimeC   = v_Huff_envelopeLevelC11T;
    sbrEnvData->hufftableLevelTimeL   = v_Huff_envelopeLevelL11T;
    sbrEnvData->hufftableBalanceTimeC = bookSbrEnvBalanceC11T;
    sbrEnvData->hufftableBalanceTimeL = bookSbrEnvBalanceL11T;

     
    sbrEnvData->hufftableLevelFreqC   = v_Huff_envelopeLevelC11F;
    sbrEnvData->hufftableLevelFreqL   = v_Huff_envelopeLevelL11F;
    sbrEnvData->hufftableBalanceFreqC = bookSbrEnvBalanceC11F;
    sbrEnvData->hufftableBalanceFreqL = bookSbrEnvBalanceL11F;


     
    sbrEnvData->hufftableTimeC        = v_Huff_envelopeLevelC11T;
    sbrEnvData->hufftableTimeL        = v_Huff_envelopeLevelL11T;
    sbrEnvData->hufftableFreqC        = v_Huff_envelopeLevelC11F;
    sbrEnvData->hufftableFreqL        = v_Huff_envelopeLevelL11F;

     
    sbrEnvData->codeBookScfLavBalance  = CODE_BOOK_SCF_LAV_BALANCE11;
    sbrEnvData->codeBookScfLav         = CODE_BOOK_SCF_LAV11;

     
    sbrEnvData->si_sbr_start_env_bits           = SI_SBR_START_ENV_BITS_AMP_RES_3_0;
    sbrEnvData->si_sbr_start_env_bits_balance   = SI_SBR_START_ENV_BITS_BALANCE_AMP_RES_3_0;
    break;

  case SBR_AMP_RES_1_5:


     
    sbrEnvData->hufftableLevelTimeC   = v_Huff_envelopeLevelC10T;
    sbrEnvData->hufftableLevelTimeL   = v_Huff_envelopeLevelL10T;
    sbrEnvData->hufftableBalanceTimeC = bookSbrEnvBalanceC10T;
    sbrEnvData->hufftableBalanceTimeL = bookSbrEnvBalanceL10T;

     
    sbrEnvData->hufftableLevelFreqC   = v_Huff_envelopeLevelC10F;
    sbrEnvData->hufftableLevelFreqL   = v_Huff_envelopeLevelL10F;
    sbrEnvData->hufftableBalanceFreqC = bookSbrEnvBalanceC10F;
    sbrEnvData->hufftableBalanceFreqL = bookSbrEnvBalanceL10F;


     
    sbrEnvData->hufftableTimeC        = v_Huff_envelopeLevelC10T;
    sbrEnvData->hufftableTimeL        = v_Huff_envelopeLevelL10T;
    sbrEnvData->hufftableFreqC        = v_Huff_envelopeLevelC10F;
    sbrEnvData->hufftableFreqL        = v_Huff_envelopeLevelL10F;

     
    sbrEnvData->codeBookScfLavBalance = CODE_BOOK_SCF_LAV_BALANCE10;
    sbrEnvData->codeBookScfLav = CODE_BOOK_SCF_LAV10;

     
    sbrEnvData->si_sbr_start_env_bits           = SI_SBR_START_ENV_BITS_AMP_RES_1_5;
    sbrEnvData->si_sbr_start_env_bits_balance   = SI_SBR_START_ENV_BITS_BALANCE_AMP_RES_1_5;
    break;

  default:
    
    return (1);
    break;
  }


   
  sbrEnvData->hufftableNoiseLevelTimeC   = v_Huff_NoiseLevelC11T;
  sbrEnvData->hufftableNoiseLevelTimeL   = v_Huff_NoiseLevelL11T;
  sbrEnvData->hufftableNoiseBalanceTimeC = bookSbrNoiseBalanceC11T;
  sbrEnvData->hufftableNoiseBalanceTimeL = bookSbrNoiseBalanceL11T;

   
  sbrEnvData->hufftableNoiseLevelFreqC   = v_Huff_envelopeLevelC11F;
  sbrEnvData->hufftableNoiseLevelFreqL   = v_Huff_envelopeLevelL11F;
  sbrEnvData->hufftableNoiseBalanceFreqC = bookSbrEnvBalanceC11F;
  sbrEnvData->hufftableNoiseBalanceFreqL = bookSbrEnvBalanceL11F;


   
  sbrEnvData->hufftableNoiseTimeC        = v_Huff_NoiseLevelC11T;
  sbrEnvData->hufftableNoiseTimeL        = v_Huff_NoiseLevelL11T;
  sbrEnvData->hufftableNoiseFreqC        = v_Huff_envelopeLevelC11F;
  sbrEnvData->hufftableNoiseFreqL        = v_Huff_envelopeLevelL11F;

   
  sbrEnvData->si_sbr_start_noise_bits         = SI_SBR_START_NOISE_BITS_AMP_RES_3_0;
  sbrEnvData->si_sbr_start_noise_bits_balance = SI_SBR_START_NOISE_BITS_BALANCE_AMP_RES_3_0;


   
  henv->codeBookScfLavBalanceTime = sbrEnvData->codeBookScfLavBalance;
  henv->codeBookScfLavBalanceFreq = sbrEnvData->codeBookScfLavBalance;
  henv->codeBookScfLavLevelTime = sbrEnvData->codeBookScfLav;
  henv->codeBookScfLavLevelFreq = sbrEnvData->codeBookScfLav;
  henv->codeBookScfLavTime = sbrEnvData->codeBookScfLav;
  henv->codeBookScfLavFreq = sbrEnvData->codeBookScfLav;

   
  henv->hufftableLevelTimeL = sbrEnvData->hufftableLevelTimeL;
  henv->hufftableBalanceTimeL = sbrEnvData->hufftableBalanceTimeL;
  henv->hufftableTimeL = sbrEnvData->hufftableTimeL;
  henv->hufftableLevelFreqL = sbrEnvData->hufftableLevelFreqL;
  henv->hufftableBalanceFreqL = sbrEnvData->hufftableBalanceFreqL;
  henv->hufftableFreqL = sbrEnvData->hufftableFreqL;

   
  henv->codeBookScfLavFreq = sbrEnvData->codeBookScfLav;
  henv->codeBookScfLavTime = sbrEnvData->codeBookScfLav;

   
  henv->start_bits = sbrEnvData->si_sbr_start_env_bits;
  henv->start_bits_balance = sbrEnvData->si_sbr_start_env_bits_balance;

   
  hnoise->codeBookScfLavBalanceTime = CODE_BOOK_SCF_LAV_BALANCE11;
  hnoise->codeBookScfLavBalanceFreq = CODE_BOOK_SCF_LAV_BALANCE11;
  hnoise->codeBookScfLavLevelTime = CODE_BOOK_SCF_LAV11;
  hnoise->codeBookScfLavLevelFreq = CODE_BOOK_SCF_LAV11;
  hnoise->codeBookScfLavTime = CODE_BOOK_SCF_LAV11;
  hnoise->codeBookScfLavFreq = CODE_BOOK_SCF_LAV11;

   
  hnoise->hufftableLevelTimeL = sbrEnvData->hufftableNoiseLevelTimeL;
  hnoise->hufftableBalanceTimeL = sbrEnvData->hufftableNoiseBalanceTimeL;
  hnoise->hufftableTimeL = sbrEnvData->hufftableNoiseTimeL;
  hnoise->hufftableLevelFreqL = sbrEnvData->hufftableNoiseLevelFreqL;
  hnoise->hufftableBalanceFreqL = sbrEnvData->hufftableNoiseBalanceFreqL;
  hnoise->hufftableFreqL = sbrEnvData->hufftableNoiseFreqL;


   
  hnoise->start_bits = sbrEnvData->si_sbr_start_noise_bits;
  hnoise->start_bits_balance = sbrEnvData->si_sbr_start_noise_bits_balance;

   
  henv->upDate = 0;
  hnoise->upDate = 0;

  

  return  (0);
}

/*******************************************************************************
 Functionname:  indexLow2High
 *******************************************************************************

 Description:   patch-functions in order to cope with non-factor-2
                ratios between high-res and low-res

 Arguments:     int offset, int index, FREQ_RES res

 Return:        int

*******************************************************************************/
static int indexLow2High(int offset, int index, FREQ_RES res)
{
  

   
  if(res == FREQ_RES_LOW)
  {
    
    if (offset >= 0)
    {
         
        if (index < offset) {
          
          return(index);
        }
        else {
            /* counting post-operation */
          
          return(2*index - offset);
        }
    }
    else
    {
        
        offset = -offset;

         
        if (index < offset) {
            /* counting post-operation */
          
          return(2*index+index);
        }
        else {
            /* counting post-operation */
          
          return(2*index + offset);
        }
    }
  }
  else {
    
    return(index);
  }
}



/*******************************************************************************
 Functionname:  mapLowResEnergyVal
 *******************************************************************************


 Arguments:     int currVal,int* prevData, int offset, int index, FREQ_RES res

 Return:        none

*******************************************************************************/
static void mapLowResEnergyVal(int currVal,int* prevData, int offset, int index, FREQ_RES res)
{

  

   
  if(res == FREQ_RES_LOW)
  {
    
    if (offset >= 0)
    {
         
        if(index < offset) {
             

            prevData[index] = currVal;
        }
        else
        {
               

            prevData[2*index - offset] = currVal;
            prevData[2*index+1 - offset] = currVal;
        }
    }
    else
    {
        
        offset = -offset;

         
        if (index < offset)
        {
              

            prevData[3*index] = currVal;
            prevData[3*index+1] = currVal;
            prevData[3*index+2] = currVal;
        }
        else
        {
              

            prevData[2*index + offset] = currVal;
            prevData[2*index + 1 + offset] = currVal;
        }
    }
  }
  else {
     

    prevData[index] = currVal;
  }

  
}



/*******************************************************************************
 Functionname:  computeBits
 *******************************************************************************

 Description:

 Arguments:     int delta,
                int codeBookScfLavLevel,
                int codeBookScfLavBalance,
                const unsigned char * hufftableLevel,
                const unsigned char * hufftableBalance, int coupling, int channel)

 Return:        int

*******************************************************************************/
static int
computeBits (int delta,
             int codeBookScfLavLevel,
             int codeBookScfLavBalance,
             const unsigned char * hufftableLevel,
             const unsigned char * hufftableBalance, int coupling, int channel)
{
  int index;
  int delta_bits = 0;

  

   /* counting previous operations */

  
  if (coupling) {

     
    if (channel == 1)
      {
           
        index =
          (delta < 0) ? max (delta, -codeBookScfLavBalance) : min (delta,
                                                                   codeBookScfLavBalance);

         
        if (index != delta) {
          assert(0);

          
          return (10000);
        }

         
        delta_bits = hufftableBalance[index + codeBookScfLavBalance];
      }
    else {

         
      index =
        (delta < 0) ? max (delta, -codeBookScfLavLevel) : min (delta,
                                                               codeBookScfLavLevel);

       
      if (index != delta) {
        assert(0);

        
        return (10000);
      }

       
      delta_bits = hufftableLevel[index + codeBookScfLavLevel];
    }
  }
  else {

       
    index =
      (delta < 0) ? max (delta, -codeBookScfLavLevel) : min (delta,
                                                             codeBookScfLavLevel);

     
    if (index != delta) {
      assert(0);

      
      return (10000);
    }

     
    delta_bits = hufftableLevel[index + codeBookScfLavLevel];
  }

  

  return (delta_bits);
}




/*******************************************************************************
 Functionname:  codeEnvelope
 *******************************************************************************

 Arguments:     int *sfb_nrg,
                const FREQ_RES *freq_res,
                SBR_CODE_ENVELOPE * h_sbrCodeEnvelope,
                int *directionVec, int scalable, int nEnvelopes, int channel,
                int headerActive)

*******************************************************************************/
void
codeEnvelope (int               *sfb_nrg,
              const FREQ_RES    *freq_res,
              SBR_CODE_ENVELOPE * h_sbrCodeEnvelope,
              int *directionVec,
              int coupling,
              int nEnvelopes,
              int channel,
              int headerActive)
{
  int i, no_of_bands, band, last_nrg, curr_nrg;
  int *ptr_nrg;

  int codeBookScfLavLevelTime;
  int codeBookScfLavLevelFreq;
  int codeBookScfLavBalanceTime;
  int codeBookScfLavBalanceFreq;
  const unsigned char *hufftableLevelTimeL;
  const unsigned char *hufftableBalanceTimeL;
  const unsigned char *hufftableLevelFreqL;
  const unsigned char *hufftableBalanceFreqL;

  int offset = h_sbrCodeEnvelope->offset;
  int envDataTableCompFactor;

  int delta_F_bits = 0, delta_T_bits = 0;
  int use_dT;

  int delta_F[MAX_FREQ_COEFFS];
  int delta_T[MAX_FREQ_COEFFS];

  float dF_edge_1stEnv = h_sbrCodeEnvelope->dF_edge_1stEnv +
    h_sbrCodeEnvelope->dF_edge_incr * h_sbrCodeEnvelope->dF_edge_incr_fac;

  

      /* counting previous operations */

  
  if (coupling) {

     
    codeBookScfLavLevelTime = h_sbrCodeEnvelope->codeBookScfLavLevelTime;
    codeBookScfLavLevelFreq = h_sbrCodeEnvelope->codeBookScfLavLevelFreq;
    codeBookScfLavBalanceTime = h_sbrCodeEnvelope->codeBookScfLavBalanceTime;
    codeBookScfLavBalanceFreq = h_sbrCodeEnvelope->codeBookScfLavBalanceFreq;
    hufftableLevelTimeL = h_sbrCodeEnvelope->hufftableLevelTimeL;
    hufftableBalanceTimeL = h_sbrCodeEnvelope->hufftableBalanceTimeL;
    hufftableLevelFreqL = h_sbrCodeEnvelope->hufftableLevelFreqL;
    hufftableBalanceFreqL = h_sbrCodeEnvelope->hufftableBalanceFreqL;
  }
  else {

     
    codeBookScfLavLevelTime = h_sbrCodeEnvelope->codeBookScfLavTime;
    codeBookScfLavLevelFreq = h_sbrCodeEnvelope->codeBookScfLavFreq;
    codeBookScfLavBalanceTime = h_sbrCodeEnvelope->codeBookScfLavTime;
    codeBookScfLavBalanceFreq = h_sbrCodeEnvelope->codeBookScfLavFreq;
    hufftableLevelTimeL = h_sbrCodeEnvelope->hufftableTimeL;
    hufftableBalanceTimeL = h_sbrCodeEnvelope->hufftableTimeL;
    hufftableLevelFreqL = h_sbrCodeEnvelope->hufftableFreqL;
    hufftableBalanceFreqL = h_sbrCodeEnvelope->hufftableFreqL;
  }

     
  if(coupling == 1 && channel == 1)
    envDataTableCompFactor = 1;
  else
    envDataTableCompFactor = 0;


   
  if (h_sbrCodeEnvelope->deltaTAcrossFrames == 0) {

     
    h_sbrCodeEnvelope->upDate = 0;
  }

  
  if (headerActive) {

     
    h_sbrCodeEnvelope->upDate = 0;
  }


   /* pointers for freq_res[],
                               directionVec[]
               */
  
  for (i = 0; i < nEnvelopes; i++)
  {
       
    if (freq_res[i] == FREQ_RES_HIGH)
      no_of_bands = h_sbrCodeEnvelope->nSfb[FREQ_RES_HIGH];
    else
      no_of_bands = h_sbrCodeEnvelope->nSfb[FREQ_RES_LOW];


    
    ptr_nrg = sfb_nrg;
    curr_nrg = *ptr_nrg;

     
    delta_F[0] = curr_nrg >> envDataTableCompFactor;

        
    if (coupling && channel == 1)
      delta_F_bits = h_sbrCodeEnvelope->start_bits_balance;
    else
      delta_F_bits = h_sbrCodeEnvelope->start_bits;


     
    if(h_sbrCodeEnvelope->upDate != 0)
    {
         
      delta_T[0] = (curr_nrg - h_sbrCodeEnvelope->sfb_nrg_prev[0]) >> envDataTableCompFactor;

      
      delta_T_bits = computeBits (delta_T[0],
                                  codeBookScfLavLevelTime,
                                  codeBookScfLavBalanceTime,
                                  hufftableLevelTimeL,
                                  hufftableBalanceTimeL, coupling, channel);
    }


     
    mapLowResEnergyVal(curr_nrg, h_sbrCodeEnvelope->sfb_nrg_prev, offset, 0, freq_res[i]);



    /* ensure that nrg difference is not higher than codeBookScfLavXXXFreq */
      
    if ( coupling && channel == 1 ) {
       
      for (band = no_of_bands - 1; band > 0; band--) {
          
        if ( ptr_nrg[band] - ptr_nrg[band-1] > codeBookScfLavBalanceFreq ) {
           
          ptr_nrg[band-1] = ptr_nrg[band] - codeBookScfLavBalanceFreq;
        }
      }
       
      for (band = 1; band < no_of_bands; band++) {
          
        if ( ptr_nrg[band-1] - ptr_nrg[band] > codeBookScfLavBalanceFreq ) {
           
          ptr_nrg[band] = ptr_nrg[band-1] - codeBookScfLavBalanceFreq;
        }
      }
    }
    else {
       
      for (band = no_of_bands - 1; band > 0; band--) {
          
        if ( ptr_nrg[band] - ptr_nrg[band-1] > codeBookScfLavLevelFreq ) {
           
          ptr_nrg[band-1] = ptr_nrg[band] - codeBookScfLavLevelFreq;
        }
      }
       
      for (band = 1; band < no_of_bands; band++) {
          
        if ( ptr_nrg[band-1] - ptr_nrg[band] > codeBookScfLavLevelFreq ) {
           
          ptr_nrg[band] = ptr_nrg[band-1] - codeBookScfLavLevelFreq;
        }
      }
    }



     /* pointers for delta_F[band],
                                 delta_T[band]
                 */
    
    for (band = 1; band < no_of_bands; band++)
    {
      
      last_nrg = (*ptr_nrg);

      
      ptr_nrg++;

      
      curr_nrg = (*ptr_nrg);


        
      delta_F[band] = (curr_nrg - last_nrg) >> envDataTableCompFactor;

       
      delta_F_bits += computeBits (delta_F[band],
                                   codeBookScfLavLevelFreq,
                                   codeBookScfLavBalanceFreq,
                                   hufftableLevelFreqL,
                                   hufftableBalanceFreqL, coupling, channel);

       
      if(h_sbrCodeEnvelope->upDate != 0)
      {
           
        delta_T[band] = curr_nrg - h_sbrCodeEnvelope->sfb_nrg_prev[indexLow2High(offset, band, freq_res[i])];

         
        delta_T[band] = delta_T[band] >> envDataTableCompFactor;
      }

       
      mapLowResEnergyVal(curr_nrg, h_sbrCodeEnvelope->sfb_nrg_prev, offset, band, freq_res[i]);

       
      if(h_sbrCodeEnvelope->upDate != 0)
      {

         
        delta_T_bits += computeBits (delta_T[band],
                                     codeBookScfLavLevelTime,
                                     codeBookScfLavBalanceTime,
                                     hufftableLevelTimeL,
                                     hufftableBalanceTimeL, coupling, channel);
      }
    }

    
    if (i == 0) {

         
      use_dT = (h_sbrCodeEnvelope->upDate != 0 && (delta_F_bits > delta_T_bits * (1 + dF_edge_1stEnv)));
    }
    else {

      
      use_dT = (delta_F_bits > delta_T_bits);
    }

    
    if (use_dT) {

      
      directionVec[i] = TIME;

          
      memcpy (sfb_nrg, delta_T, no_of_bands * sizeof (int));
    }
    else {

      
      directionVec[i] = FREQ;

          
      memcpy (sfb_nrg, delta_F, no_of_bands * sizeof (int));
    }

    
    sfb_nrg += no_of_bands;

     
    h_sbrCodeEnvelope->upDate = 1;
  }

  
}




/*******************************************************************************
 Functionname:  CreateSbrCodeEnvelope
 *******************************************************************************

 Description:

 Arguments:

 Return:

*******************************************************************************/
int
CreateSbrCodeEnvelope (HANDLE_SBR_CODE_ENVELOPE  h_sbrCodeEnvelope,
                       int *nSfb,
                       int deltaTAcrossFrames,
                       float dF_edge_1stEnv,
                       float dF_edge_incr)
{
  

      
  memset(h_sbrCodeEnvelope,0,sizeof(SBR_CODE_ENVELOPE));

   
  h_sbrCodeEnvelope->deltaTAcrossFrames = deltaTAcrossFrames;
  h_sbrCodeEnvelope->dF_edge_1stEnv = dF_edge_1stEnv;
  h_sbrCodeEnvelope->dF_edge_incr = dF_edge_incr;
  h_sbrCodeEnvelope->dF_edge_incr_fac = 0;
  h_sbrCodeEnvelope->upDate = 0;
  h_sbrCodeEnvelope->nSfb[FREQ_RES_LOW] = nSfb[FREQ_RES_LOW];
  h_sbrCodeEnvelope->nSfb[FREQ_RES_HIGH] = nSfb[FREQ_RES_HIGH];

     
  h_sbrCodeEnvelope->offset = 2*h_sbrCodeEnvelope->nSfb[FREQ_RES_LOW] - h_sbrCodeEnvelope->nSfb[FREQ_RES_HIGH];

  

  return (0);
}





/*******************************************************************************
 Functionname:  deleteSbrCodeEnvelope
 *******************************************************************************

 Description:

 Arguments:

 Return:

*******************************************************************************/
void
deleteSbrCodeEnvelope (HANDLE_SBR_CODE_ENVELOPE h_sbrCodeEnvelope)
{

  
  /*
    nothing to do
  */

  
}
