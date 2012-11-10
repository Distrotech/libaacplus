/*
  Framing generator
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static SBR_FRAME_INFO
frameInfo1_2048 = {1,
                   {0, 16},
                   {FREQ_RES_HIGH},
                   0,
                   1,
                   {0, 16}};

static SBR_FRAME_INFO
frameInfo2_2048 = {2,
                   {0, 8, 16},
                   {FREQ_RES_HIGH, FREQ_RES_HIGH},
                   0,
                   2,
                   {0, 8, 16}};

static SBR_FRAME_INFO
frameInfo4_2048 = {4,
                   {0, 4, 8, 12, 16},
                   {FREQ_RES_HIGH,
                    FREQ_RES_HIGH,
                    FREQ_RES_HIGH,
                    FREQ_RES_HIGH},
                   0,
                   2,
                   {0, 8, 16}};



static void
addFreqLeft (FREQ_RES *vector, int *length_vector, FREQ_RES value)
{
  int i;

  

   /* vector[] */
  
  for (i = *length_vector; i > 0; i--)
  {
    
    vector[i] = vector[i - 1];
  }

  
  vector[0] = value;

  
  (*length_vector)++;

  
}


static void
addFreqVecLeft (FREQ_RES *dst, int *length_dst, FREQ_RES *src, int length_src)
{
  int i;

  

   /* src[] */
  
  for (i = length_src - 1; i >= 0; i--)
  {
    
    addFreqLeft (dst, length_dst, src[i]);
  }

  
}



static void
addFreqRight (FREQ_RES *vector, int *length_vector, FREQ_RES value)
{
  

   
  vector[*length_vector] = value;

  
  (*length_vector)++;

  
}


/***************************************************************************/
/*!
  \brief     Add mandatory borders

  \return    void

****************************************************************************/
static void
fillFrameTran (int *v_bord, int *length_v_bord,
               FREQ_RES *v_freq, int *length_v_freq,
               int *bmin, int *bmax,
               int tran,
               int *v_tuningSegm, FREQ_RES *v_tuningFreq)
{
  int bord, i;

  

  
  *length_v_bord = 0;
  *length_v_freq = 0;

  
  if (v_tuningSegm[0]) {

     
    AddRight (v_bord, length_v_bord, (tran - v_tuningSegm[0]));

    
    addFreqRight (v_freq, length_v_freq, v_tuningFreq[0]);
  }

  
  bord = tran;

  
  AddRight (v_bord, length_v_bord, tran);

  
  if (v_tuningSegm[1]) {

    
    bord += v_tuningSegm[1];

    
    AddRight (v_bord, length_v_bord, bord);

    
    addFreqRight (v_freq, length_v_freq, v_tuningFreq[1]);
  }


  
  if (v_tuningSegm[2] != 0) {

    
    bord += v_tuningSegm[2];


    
    AddRight (v_bord, length_v_bord, bord);

    
    addFreqRight (v_freq, length_v_freq, v_tuningFreq[2]);
  }

  
  addFreqRight (v_freq, length_v_freq, FREQ_RES_HIGH);


  
  *bmin = v_bord[0];

   /* v_bord[] */
  
  for (i = 0; i < *length_v_bord; i++)
  {
     
    if (v_bord[i] < *bmin)
    {
      
      *bmin = v_bord[i];
    }
  }

  
  *bmax = v_bord[0];

   /* v_bord[] */
  
  for (i = 0; i < *length_v_bord; i++)
  {
     
    if (v_bord[i] > *bmax)
    {
      
      *bmax = v_bord[i];
    }
  }

  
}


/***************************************************************************/
/*!
  \brief     Add borders before mandatory borders, if needed

  \return    void

****************************************************************************/
static void
fillFramePre (int dmax,
              int *v_bord, int *length_v_bord,
              FREQ_RES *v_freq, int *length_v_freq,
              int bmin, int rest)
{

  int parts, d, j, S, s = 0, segm, bord;

  

   /* counting previous operation */


  
  parts = 1;
  d = rest;

  
  while (d > dmax) {
    
    parts++;

    
    segm = rest / parts;

      
    S = (int) floor ((segm - 2) * 0.5);

       
    s = min (8, 2 * S + 2);

     
    d = rest - (parts - 1) * s;
  }


  
  bord = bmin;

  
  for (j = 0; j <= parts - 2; j++) {

    
    bord = bord - s;

    
    AddLeft (v_bord, length_v_bord, bord);

    
    addFreqLeft (v_freq, length_v_freq, FREQ_RES_HIGH);
  }

  
}

/***************************************************************************/
/*!
  \brief Add borders after mandatory borders

  \return    void

****************************************************************************/
static void
fillFramePost (int *parts, int *d, int dmax, int *v_bord, int *length_v_bord,
               FREQ_RES *v_freq, int *length_v_freq, int bmax,
               int fmax)
{
  int j, rest, segm, S, s = 0, bord;

  

   /* counting previous operation */

  
  rest = 32 - bmax;

  
  *d = rest;

  
  if (*d > 0) {

    
    *parts = 1;           /* start with one envelope */

    
    while (*d > dmax) {

       
      *parts = *parts + 1;

      
      segm = rest / (*parts);

        
      S = (int) floor ((segm - 2) * 0.5);

         
      s = min (fmax, 2 * S + 2);

        
      *d = rest - (*parts - 1) * s;
    }

    
    bord = bmax;

     
    for (j = 0; j <= *parts - 2; j++) {

      
      bord += s;

      
      AddRight (v_bord, length_v_bord, bord);

      
      addFreqRight (v_freq, length_v_freq, FREQ_RES_HIGH);
    }
  }
  else {
    
    *parts = 1;

     
    *length_v_bord = *length_v_bord - 1;
    *length_v_freq = *length_v_freq - 1;

  }

  
}


/***************************************************************************/
/*!
  \brief    Fills in the envelopes between transients

  \return   void

****************************************************************************/
static void
fillFrameInter (int *nL,
                int *v_tuningSegm,
                int *v_bord,
                int *length_v_bord,
                int bmin,
                FREQ_RES *v_freq,
                int *length_v_freq,
                int *v_bordFollow,
                int *length_v_bordFollow,
                FREQ_RES *v_freqFollow,
                int *length_v_freqFollow,
                int i_fillFollow,
                int dmin,
                int dmax)
{
  int middle, b_new, numBordFollow, bordMaxFollow, i;

  


     
    if (i_fillFollow >= 1) {

      
      *length_v_bordFollow = i_fillFollow;
      *length_v_freqFollow = i_fillFollow;
    }

    
    numBordFollow = *length_v_bordFollow;

     
    bordMaxFollow = v_bordFollow[numBordFollow - 1];

    
    middle = bmin - bordMaxFollow;

     /* v_bordFollow[] */
    
    while (middle < 0) {

      
      numBordFollow--;

      
      bordMaxFollow = v_bordFollow[numBordFollow - 1];

      
      middle = bmin - bordMaxFollow;
    }

    
    *length_v_bordFollow = numBordFollow;
    *length_v_freqFollow = numBordFollow;

     
    *nL = numBordFollow - 1;

    
    b_new = *length_v_bord;


     
    if (middle <= dmax) {

       
      if (middle >= dmin) {

        
        AddVecLeft (v_bord, length_v_bord, v_bordFollow, *length_v_bordFollow);

        
        addFreqVecLeft (v_freq, length_v_freq, v_freqFollow, *length_v_freqFollow);
      }

      else {

        
        if (v_tuningSegm[0] != 0) {

           
          *length_v_bord = b_new - 1;

          
          AddVecLeft (v_bord, length_v_bord, v_bordFollow,
              *length_v_bordFollow);

           
          *length_v_freq = b_new - 1;

           
          addFreqVecLeft (v_freq + 1, length_v_freq, v_freqFollow,
                          *length_v_freqFollow);

        }
        else {

           
          if (*length_v_bordFollow > 1) {

             
            AddVecLeft (v_bord, length_v_bord, v_bordFollow,
                        *length_v_bordFollow - 1);

             
            addFreqVecLeft (v_freq, length_v_freq, v_freqFollow,
                            *length_v_bordFollow - 1);

             
            *nL = *nL - 1;
          }
          else {


             /* v_bord[] */
            
            for (i = 0; i < *length_v_bord - 1; i++)
            {
              
              v_bord[i] = v_bord[i + 1];
            }

             /* v_freq[] */
            
            for (i = 0; i < *length_v_freq - 1; i++)
            {
              
              v_freq[i] = v_freq[i + 1];
            }

             
            *length_v_bord = b_new - 1;
            *length_v_freq = b_new - 1;

            
            AddVecLeft (v_bord, length_v_bord, v_bordFollow,
                        *length_v_bordFollow);

            
            addFreqVecLeft (v_freq, length_v_freq, v_freqFollow,
                            *length_v_freqFollow);
          }
        }
      }
    }
    else {

      
      fillFramePre (dmax, v_bord, length_v_bord, v_freq, length_v_freq, bmin,
                    middle);

      
      AddVecLeft (v_bord, length_v_bord, v_bordFollow, *length_v_bordFollow);

      
      addFreqVecLeft (v_freq, length_v_freq, v_freqFollow, *length_v_freqFollow);
    }

  
}



/***************************************************************************/
/*!
  \brief    Calculates frame class for current frame

  \return   void

****************************************************************************/
static void
calcFrameClass (FRAME_CLASS *frameClass,
                FRAME_CLASS *frameClassOld,
                int tranFlag,
                int *pSpreadFlag
                )
{
  

  
  switch (*frameClassOld) {

  case FIXFIX:

     
    if (tranFlag)
      *frameClass = FIXVAR;
    else
      *frameClass = FIXFIX;
    break;

  case FIXVAR:

    
    if (tranFlag) {

      
      *frameClass = VARVAR;
      *pSpreadFlag = 0;
    }
    else {

       
      if (*pSpreadFlag)
        *frameClass = VARVAR;
      else
        *frameClass = VARFIX;
    }
    break;

  case VARFIX:

     
    if (tranFlag)
      *frameClass = FIXVAR;
    else
      *frameClass = FIXFIX;
    break;

  case VARVAR:

    
    if (tranFlag) {

      
      *frameClass = VARVAR;
      *pSpreadFlag = 0;
    }
    else {

       
      if (*pSpreadFlag)
        *frameClass = VARVAR;
      else
        *frameClass = VARFIX;
    }
    break;

  default:
    assert(0);
  }

  
  *frameClassOld = *frameClass;

  
}


/***************************************************************************/
/*!
  \brief    Handles special case for the framing generator

  \return   void

****************************************************************************/
static void
specialCase (int *pSpreadFlag,
             int allowSpread,
             int *aBorders,
             int *pBorderVecLen,
             FREQ_RES *aFreqRes,
             int *pFreqResVecLen,
             int *pNumParts,
             int d
             )
{
  int L;

  

  
  L = *pBorderVecLen;

  
  if (allowSpread) {

    
    *pSpreadFlag = 1;

      
    AddRight (aBorders, pBorderVecLen, aBorders[L - 1] + 8);

    
    addFreqRight (aFreqRes, pFreqResVecLen, FREQ_RES_HIGH);

     
    (*pNumParts)++;
  }
  else {

     
    if (d == 1) {

       
      *pBorderVecLen = L - 1;
      *pFreqResVecLen = L - 1;
    }
    else {

        
      if ((aBorders[L - 1] - aBorders[L - 2]) > 2) {

         
        aBorders[L - 1] = aBorders[L - 1] - 2;

         
        aFreqRes[*pFreqResVecLen - 1] = FREQ_RES_LOW;
      }
    }
  }

  
}



/***************************************************************************/
/*!
  \brief    Calculates a common border

  \return   void

****************************************************************************/
static void
calcCmonBorder (int *pCmonBorderIdx,  /*!< */
                int *pTranIdx,        /*!< */
                int *aBorders,        /*!< */
                int *pBorderVecLen,   /*!< */
                int tran             /*!< */


                )
{
  int i;

  

   /* aBorders[] */
  
  for (i = 0; i < *pBorderVecLen; i++)
  {
     
    if (aBorders[i] >= 16) {

      
      *pCmonBorderIdx = i;
      break;
    }
  }

   /* aBorders[] */
  
  for (i = 0; i < *pBorderVecLen; i++)
  {
     
    if (aBorders[i] >= tran) {

      
      *pTranIdx = i;
      break;
    }
    else
    {
      
      *pTranIdx = EMPTY;
    }
  }

  
}



/***************************************************************************/
/*!
  \brief    Stores values for following frame

  \return   void

****************************************************************************/
static void
keepForFollowUp (int *aBordersFollow,
                 int *pBorderVecLenFollow,
                 FREQ_RES *aFreqResFollow,
                 int *pFreqResVecLenFollow,
                 int *pTranIdxFollow,
                 int *pFillIdxFollow,
                 int *aBorders,
                 int *pBorderVecLen,
                 FREQ_RES *aFreqRes,
                 int cmonBorderIdx,
                 int tranIdx,
                 int numParts

                 )
{
  int L, i, j;

  

  
  L = *pBorderVecLen;

  
  (*pBorderVecLenFollow) = 0;
  (*pFreqResVecLenFollow) = 0;

   /* aBordersFollow[]
                  aBorders[]
                  aFreqResFollow[]
                  aFreqRes[];
               */
  
  for (j = 0, i = cmonBorderIdx; i < L; i++, j++) {

     
    aBordersFollow[j] = aBorders[i] - 16;

    
    aFreqResFollow[j] = aFreqRes[i];

     
    (*pBorderVecLenFollow)++;
    (*pFreqResVecLenFollow)++;
  }

   
  if (tranIdx != EMPTY)
  {
     
    *pTranIdxFollow = tranIdx - cmonBorderIdx;
  }
  else
  {
    
    *pTranIdxFollow = EMPTY;
  }

   
  *pFillIdxFollow = L - (numParts - 1) - cmonBorderIdx;

  
}


/***************************************************************************/
/*!
  \brief    Calculates the control signal

  \return   void

****************************************************************************/
static void
calcCtrlSignal (HANDLE_SBR_GRID hSbrGrid,
                FRAME_CLASS frameClass, int *v_bord, int length_v_bord, FREQ_RES *v_freq,
                int length_v_freq, int i_cmon, int i_tran, int spreadFlag,
                int nL)
{


  int i, r, a, n, p, b, aL, aR, ntot, nmax, nR;

  FREQ_RES *v_f = hSbrGrid->v_f;
  FREQ_RES *v_fLR = hSbrGrid->v_fLR;
  int *v_r = hSbrGrid->bs_rel_bord;
  int *v_rL = hSbrGrid->bs_rel_bord_0;
  int *v_rR = hSbrGrid->bs_rel_bord_1;

  int length_v_r = 0;
  int length_v_rR = 0;
  int length_v_rL = 0;

  

     /* counting previous operations */

  
  switch (frameClass) {
  case FIXVAR:

    
    a = v_bord[i_cmon];

    
    length_v_r = 0;
    i = i_cmon;

     /* v_bord[i] */
    
    while (i >= 1) {

      
      r = v_bord[i] - v_bord[i - 1];

       
      AddRight (v_r, &length_v_r, r);

      i--;
    }


    
    n = length_v_r;

     /* v_f[i]
                    v_freq[i_cmon - 1 - i]
                 */
    
    for (i = 0; i < i_cmon; i++)
    {
      
      v_f[i] = v_freq[i_cmon - 1 - i];
    }

     
    v_f[i_cmon] = FREQ_RES_HIGH;

    /* pointer: */
      
    if (i_cmon >= i_tran && i_tran != EMPTY)
    {
      
      p = i_cmon - i_tran + 1;
    }
    else
    {
      
      p = 0;
    }

     
    hSbrGrid->frameClass = frameClass;
    hSbrGrid->bs_abs_bord = a;
    hSbrGrid->n = n;
    hSbrGrid->p = p;

    break;
  case VARFIX:

    
    a = v_bord[0];

    
    length_v_r = 0;

     /* v_bord[] */
    
    for (i = 1; i < length_v_bord; i++) {

      
      r = v_bord[i] - v_bord[i - 1];

       
      AddRight (v_r, &length_v_r, r);
    }


    
    n = length_v_r;

        
    memcpy (v_f, v_freq, length_v_freq * sizeof (int));

    /* pointer: */
      
    if (i_tran >= 0 && i_tran != EMPTY)
    {
      
      p = i_tran + 1;
    }
    else
    {
      
      p = 0;
    }

     
    hSbrGrid->frameClass = frameClass;
    hSbrGrid->bs_abs_bord = a;
    hSbrGrid->n = n;
    hSbrGrid->p = p;

    break;
  case VARVAR:

    
    if (spreadFlag) {

      
      b = length_v_bord;

      
      aL = v_bord[0];

       
      aR = v_bord[b - 1];


      
      ntot = b - 2;

      
      nmax = 2;

       
      if (ntot > nmax) {

        
        nL = nmax;
        nR = ntot - nmax;
      }
      else {

        
        nL = ntot;
        nR = 0;
      }

      
      length_v_rL = 0;

       /* v_bord[] */
      
      for (i = 1; i <= nL; i++) {

        
        r = v_bord[i] - v_bord[i - 1];

         
        AddRight (v_rL, &length_v_rL, r);
      }

      
      length_v_rR = 0;

      
      i = b - 1;

       /* v_bord[] */
       
      while (i >= b - nR) {

        
        r = v_bord[i] - v_bord[i - 1];

         
        AddRight (v_rR, &length_v_rR, r);

        i--;
      }


        
      if (i_tran > 0 && i_tran != EMPTY)
      {
        
        p = b - i_tran;
      }
      else
      {
        
        p = 0;
      }


       /* v_fLR[i]
                      v_freq[i]
                   */
      
      for (i = 0; i < b - 1; i++)
      {
        
        v_fLR[i] = v_freq[i];
      }
    }
    else {

      
      length_v_bord = i_cmon + 1;
      length_v_freq = i_cmon + 1;


      
      b = length_v_bord;

      
      aL = v_bord[0];

       
      aR = v_bord[b - 1];

      
      ntot = b - 2;
      nR = ntot - nL;


      
      length_v_rL = 0;

       /* v_bord[i] */
      
      for (i = 1; i <= nL; i++) {

        
        r = v_bord[i] - v_bord[i - 1];

         
        AddRight (v_rL, &length_v_rL, r);
      }



      
      length_v_rR = 0;

      
      i = b - 1;

       /* v_bord[i] */
       
      while (i >= b - nR) {

        
        r = v_bord[i] - v_bord[i - 1];

         
        AddRight (v_rR, &length_v_rR, r);

        i--;
      }



        
      if (i_cmon >= i_tran && i_tran != EMPTY)
      {
        
        p = i_cmon - i_tran + 1;
      }
      else
      {
        
        p = 0;
      }




       /* v_fLR[i]
                      v_freq[i]
                   */
      
      for (i = 0; i < b - 1; i++)
      {
        
        v_fLR[i] = v_freq[i];
      }

    }


     
    hSbrGrid->frameClass = frameClass;
    hSbrGrid->bs_abs_bord_0 = aL;
    hSbrGrid->bs_abs_bord_1 = aR;
    hSbrGrid->bs_num_rel_0 = nL;
    hSbrGrid->bs_num_rel_1 = nR;
    hSbrGrid->p = p;

    break;

  default:
    /* do nothing */
    break;
  }

  
}


/***************************************************************************/
/*!
  \brief    Creates a frame info default

  \return   void

****************************************************************************/
static void
createDefFrameInfo(HANDLE_SBR_FRAME_INFO hSbrFrameInfo,
                   int nEnv

                   )
{
  

  
  switch (nEnv) {
  case 1:



          
      memcpy (hSbrFrameInfo, &frameInfo1_2048, sizeof (SBR_FRAME_INFO));


    break;

  case 2:



          
      memcpy (hSbrFrameInfo, &frameInfo2_2048, sizeof (SBR_FRAME_INFO));


    break;

  case 4:



          
      memcpy (hSbrFrameInfo, &frameInfo4_2048, sizeof (SBR_FRAME_INFO));


    break;

  default:
    assert(0);
  }

  
}

/***************************************************************************/
/*!
  \brief    Translates frame_info struct to control signal

  \return   void

****************************************************************************/
static void
ctrlSignal2FrameInfo (HANDLE_SBR_GRID hSbrGrid,
                      HANDLE_SBR_FRAME_INFO hSbrFrameInfo,
                      FREQ_RES freq_res_fixfix
                      )
{
  int nEnv = 0, border = 0, i, k, p;
  int *v_r = hSbrGrid->bs_rel_bord;
  FREQ_RES *v_f = hSbrGrid->v_f;

  FRAME_CLASS frameClass = hSbrGrid->frameClass;



  

     /* counting previous operations */

  
  switch (frameClass) {
  case FIXFIX:
     
    createDefFrameInfo(hSbrFrameInfo, hSbrGrid->bs_num_env);

     
    if (freq_res_fixfix == FREQ_RES_LOW) {

       /* hSbrFrameInfo->freqRes[] */
       
      for (i = 0; i < hSbrFrameInfo->nEnvelopes; i++) {
        
        hSbrFrameInfo->freqRes[i] = FREQ_RES_LOW;
      }
    }
    break;

  case FIXVAR:
  case VARFIX:

     
    nEnv = hSbrGrid->n + 1;

    assert(nEnv <= MAX_ENVELOPES_FIXVAR_VARFIX);

     
    hSbrFrameInfo->nEnvelopes = nEnv;

     
    border = hSbrGrid->bs_abs_bord;

       
    if (nEnv == 1)
      hSbrFrameInfo->nNoiseEnvelopes = 1;
    else
      hSbrFrameInfo->nNoiseEnvelopes = 2;

    break;

  default:
    /* do nothing */
    break;
  }

  
  switch (frameClass) {
  case FIXVAR:

     
    hSbrFrameInfo->borders[0] = 0;

     
    hSbrFrameInfo->borders[nEnv] = border;

     /* v_r[k];
                    hSbrFrameInfo->borders[i]
                 */
     
    for (k = 0, i = nEnv - 1; k < nEnv - 1; k++, i--) {
      
      border -= v_r[k];

      
      hSbrFrameInfo->borders[i] = border;
    }

     
    p = hSbrGrid->p;

    
    if (p == 0) {

       
      hSbrFrameInfo->shortEnv = 0;
    } else {

        
      hSbrFrameInfo->shortEnv = nEnv + 1 - p;
    }

     /* hSbrFrameInfo->freqRes[i]
                    v_f[k]
                 */
     
    for (k = 0, i = nEnv - 1; k < nEnv; k++, i--) {
      
      hSbrFrameInfo->freqRes[i] = v_f[k];
    }

      
    if (p == 0 || p == 1) {
       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
    } else {
       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
    }

    break;

  case VARFIX:
     
    hSbrFrameInfo->borders[0] = border;

     /* v_r[]
                    hSbrFrameInfo->borders[]
                 */
     
    for (k = 0; k < nEnv - 1; k++) {
      
      border += v_r[k];

      
      hSbrFrameInfo->borders[k + 1] = border;
    }

     
    hSbrFrameInfo->borders[nEnv] = 16;

     
    p = hSbrGrid->p;

      
    if (p == 0 || p == 1) {
       
      hSbrFrameInfo->shortEnv = 0;
    } else {
        
      hSbrFrameInfo->shortEnv = p - 1;
    }

     /* hSbrFrameInfo->freqRes[k]
                    v_f[k]
                 */
    
    for (k = 0; k < nEnv; k++) {
      
      hSbrFrameInfo->freqRes[k] = v_f[k];
    }

    
    switch (p) {
    case 0:
       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[1];
      break;
    case 1:
       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
      break;
    default:
       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
      break;
    }
    break;

  case VARVAR:
     
    nEnv = hSbrGrid->bs_num_rel_0 + hSbrGrid->bs_num_rel_1 + 1;

    assert(nEnv <= MAX_ENVELOPES_VARVAR); /* just to be sure */

     
    hSbrFrameInfo->nEnvelopes = nEnv;

      
    hSbrFrameInfo->borders[0] = border = hSbrGrid->bs_abs_bord_0;

     /* hSbrGrid->bs_rel_bord_0[k]
                    hSbrFrameInfo->borders[i]
                 */
     
    for (k = 0, i = 1; k < hSbrGrid->bs_num_rel_0; k++, i++) {
      
      border += hSbrGrid->bs_rel_bord_0[k];

      
      hSbrFrameInfo->borders[i] = border;
    }

     
    border = hSbrGrid->bs_abs_bord_1;

     
    hSbrFrameInfo->borders[nEnv] = border;

     /* hSbrGrid->bs_rel_bord_1[k]
                    hSbrFrameInfo->borders[i]
                 */
      
    for (k = 0, i = nEnv - 1; k < hSbrGrid->bs_num_rel_1; k++, i--) {
      
      border -= hSbrGrid->bs_rel_bord_1[k];

      
      hSbrFrameInfo->borders[i] = border;
    }

     
    p = hSbrGrid->p;

    
    if (p == 0) {

       
      hSbrFrameInfo->shortEnv = 0;
    } else {

        
      hSbrFrameInfo->shortEnv = nEnv + 1 - p;
    }

     /* hSbrFrameInfo->freqRes[k]
                    hSbrGrid->v_fLR[k]
                 */
    
    for (k = 0; k < nEnv; k++) {
      
      hSbrFrameInfo->freqRes[k] = hSbrGrid->v_fLR[k];
    }

     
    if (nEnv == 1) {

       
      hSbrFrameInfo->nNoiseEnvelopes = 1;
      hSbrFrameInfo->bordersNoise[0] = hSbrGrid->bs_abs_bord_0;
      hSbrFrameInfo->bordersNoise[1] = hSbrGrid->bs_abs_bord_1;
    } else {

       
      hSbrFrameInfo->nNoiseEnvelopes = 2;
      hSbrFrameInfo->bordersNoise[0] = hSbrGrid->bs_abs_bord_0;

        
      if (p == 0 || p == 1) {

         
        hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
      } else {

         
        hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
      }

       
      hSbrFrameInfo->bordersNoise[2] = hSbrGrid->bs_abs_bord_1;
    }
    break;

  default:
    /* do nothing */
    break;
  }

    
  if (frameClass == VARFIX || frameClass == FIXVAR) {

     
    hSbrFrameInfo->bordersNoise[0] = hSbrFrameInfo->borders[0];

     
    if (nEnv == 1) {

       
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv];
    } else {

       
      hSbrFrameInfo->bordersNoise[2] = hSbrFrameInfo->borders[nEnv];
    }
  }

  
}


/***************************************************************************/
/*!
  \brief    Produces the FRAME_INFO struct for the current frame

  \return   The frame info handle for the current frame

****************************************************************************/
HANDLE_SBR_FRAME_INFO
frameInfoGenerator (HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame,
                    int *v_pre_transient_info,
                    int *v_transient_info,
                    int *v_tuning)
{
  int numEnv, tran=0, bmin=0, bmax=0, parts, d, i_cmon, i_tran, nL;
  int fmax = 0;

  int *v_bord = hSbrEnvFrame->v_bord;
  FREQ_RES *v_freq = hSbrEnvFrame->v_freq;
  int *v_bordFollow = hSbrEnvFrame->v_bordFollow;
  FREQ_RES *v_freqFollow = hSbrEnvFrame->v_freqFollow;


  int *length_v_bordFollow = &hSbrEnvFrame->length_v_bordFollow;
  int *length_v_freqFollow = &hSbrEnvFrame->length_v_freqFollow;
  int *length_v_bord = &hSbrEnvFrame->length_v_bord;
  int *length_v_freq = &hSbrEnvFrame->length_v_freq;
  int *spreadFlag = &hSbrEnvFrame->spreadFlag;
  int *i_tranFollow = &hSbrEnvFrame->i_tranFollow;
  int *i_fillFollow = &hSbrEnvFrame->i_fillFollow;
  FRAME_CLASS *frameClassOld = &hSbrEnvFrame->frameClassOld;
  FRAME_CLASS frameClass;


  int allowSpread = hSbrEnvFrame->allowSpread;
  int numEnvStatic = hSbrEnvFrame->numEnvStatic;
  int staticFraming = hSbrEnvFrame->staticFraming;
  int dmin = hSbrEnvFrame->dmin;
  int dmax = hSbrEnvFrame->dmax;





  int tranPos = v_transient_info[0];
  int tranFlag = v_transient_info[1];

  int *v_tuningSegm = v_tuning;
  FREQ_RES *v_tuningFreq = (FREQ_RES*) (v_tuning + 3);

  FREQ_RES freq_res_fixfix = hSbrEnvFrame->freq_res_fixfix;

  

     /* counting previous operation */


  
  if (staticFraming) {

     
    frameClass = FIXFIX;
    numEnv = numEnvStatic;
    *frameClassOld = FIXFIX;
    hSbrEnvFrame->SbrGrid.bs_num_env = numEnv;
    hSbrEnvFrame->SbrGrid.frameClass = frameClass;
  }
  else {
     
    calcFrameClass (&frameClass, frameClassOld, tranFlag, spreadFlag);

    
    if (tranFlag) {

       
      if (tranPos < 4)
      {
        
        fmax = 6;
      }
      else {
         
      if (tranPos == 4 || tranPos == 5)
        fmax = 4;
      else
        fmax = 8;
      }

      
      tran = tranPos + 4;

       
      fillFrameTran (v_bord, length_v_bord, v_freq, length_v_freq, &bmin,
                     &bmax, tran, v_tuningSegm, v_tuningFreq);
    }

    
    switch (frameClass) {
    case FIXVAR:

      
      fillFramePre (dmax, v_bord, length_v_bord, v_freq, length_v_freq,
                    bmin, bmin);


       
      fillFramePost (&parts, &d, dmax, v_bord, length_v_bord, v_freq,
                     length_v_freq, bmax, fmax);

        
      if (parts == 1 && d < dmin)
      {
         
        specialCase (spreadFlag, allowSpread, v_bord, length_v_bord,
                     v_freq, length_v_freq, &parts, d);
      }

       
      calcCmonBorder (&i_cmon, &i_tran, v_bord, length_v_bord, tran
                      );
      
      keepForFollowUp (v_bordFollow, length_v_bordFollow, v_freqFollow,
                       length_v_freqFollow, i_tranFollow, i_fillFollow,
                       v_bord, length_v_bord, v_freq, i_cmon, i_tran, parts);  /* FH 00-06-26 */

        
      calcCtrlSignal (&hSbrEnvFrame->SbrGrid, frameClass,
                      v_bord, *length_v_bord, v_freq, *length_v_freq,
                      i_cmon, i_tran, *spreadFlag, DC);
      break;
    case VARFIX:
        
      calcCtrlSignal (&hSbrEnvFrame->SbrGrid, frameClass,
                      v_bordFollow, *length_v_bordFollow, v_freqFollow,
                      *length_v_freqFollow, DC, *i_tranFollow,
                      *spreadFlag, DC);
      break;
    case VARVAR:

      
      if (*spreadFlag) {
          
        calcCtrlSignal (&hSbrEnvFrame->SbrGrid,
                        frameClass, v_bordFollow, *length_v_bordFollow,
                        v_freqFollow, *length_v_freqFollow, DC,
                        *i_tranFollow, *spreadFlag, DC);

        
        *spreadFlag = 0;

          
        v_bordFollow[0] = hSbrEnvFrame->SbrGrid.bs_abs_bord_1 - 16;

        
        v_freqFollow[0] = FREQ_RES_HIGH;
        *length_v_bordFollow = 1;
        *length_v_freqFollow = 1;

         
        *i_tranFollow = -DC;
        *i_fillFollow = -DC;
      }
      else {
         
        fillFrameInter (&nL, v_tuningSegm, v_bord, length_v_bord, bmin,
                        v_freq, length_v_freq, v_bordFollow,
                        length_v_bordFollow, v_freqFollow,
                        length_v_freqFollow, *i_fillFollow, dmin, dmax
                        );

         
        fillFramePost (&parts, &d, dmax, v_bord, length_v_bord, v_freq,
                       length_v_freq, bmax, fmax);

          
        if (parts == 1 && d < dmin)
        {
           
          specialCase (spreadFlag, allowSpread, v_bord, length_v_bord,
                       v_freq, length_v_freq, &parts, d);
        }

         
        calcCmonBorder (&i_cmon, &i_tran, v_bord, length_v_bord, tran
                        );

        
        keepForFollowUp (v_bordFollow, length_v_bordFollow,
                         v_freqFollow, length_v_freqFollow,
                         i_tranFollow, i_fillFollow, v_bord,
                         length_v_bord, v_freq, i_cmon, i_tran, parts);

          
        calcCtrlSignal (&hSbrEnvFrame->SbrGrid,
                        frameClass, v_bord, *length_v_bord, v_freq,
                        *length_v_freq, i_cmon, i_tran, 0, nL);
      }
      break;
    case FIXFIX:
       
      if (tranPos == 0)
        numEnv = 1;
      else
        numEnv = 2;

       
      hSbrEnvFrame->SbrGrid.bs_num_env = numEnv;
      hSbrEnvFrame->SbrGrid.frameClass = frameClass;

      break;
    default:
      assert(0);
    }
  }


    
  ctrlSignal2FrameInfo (&hSbrEnvFrame->SbrGrid,
                        &hSbrEnvFrame->SbrFrameInfo,
                        freq_res_fixfix);

   

  

  return &hSbrEnvFrame->SbrFrameInfo;
}

/***************************************************************************/
/*!
  \brief    Creates an instance of the SBR framing handle

  \return  void

****************************************************************************/
void
CreateFrameInfoGenerator (HANDLE_SBR_ENVELOPE_FRAME  hSbrEnvFrame,
                          int allowSpread,
                          int numEnvStatic,
                          int staticFraming,

                          FREQ_RES freq_res_fixfix)

{

  

   
  hSbrEnvFrame->allowSpread = allowSpread;
  hSbrEnvFrame->numEnvStatic = numEnvStatic;
  hSbrEnvFrame->staticFraming = staticFraming;
  hSbrEnvFrame->freq_res_fixfix = freq_res_fixfix;

   
  hSbrEnvFrame->dmin = 4;
  hSbrEnvFrame->dmax = 12;

  
}

/***************************************************************************/
/*!
  \brief    Deletes an instance of the SBR framing handle

  \return   void

****************************************************************************/
void
deleteFrameInfoGenerator (HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame)
{
  

  /*
    nothing to do
  */

  
}


