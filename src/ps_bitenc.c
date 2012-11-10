/*
  Parametric stereo bitstream encoder
*/

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static int
appendBitstream( HANDLE_BIT_BUF hBitBufWrite,
                 HANDLE_BIT_BUF hBitBufRead,
                 int nBits )
{
  int i;
  unsigned long b;

  

  
  for (i=0; i< nBits; i++)
  {
    
    b = ReadBits (hBitBufRead, 1);

    
    WriteBits (hBitBufWrite, b, 1);
  }

  

  return nBits;
}

int
WritePsData (HANDLE_PS_ENC h_ps_e,
             int           bHeaderActive)
{
  int temp, gr;
  const int *aaHuffBookIidC;
  const short *aaHuffBookIccC;
  const char *aaHuffBookIidL;
  const char *aaHuffBookIccL;
  int *aaDeltaIid;
  int *aaDeltaIcc;

  int index, lastIndex;
  int noBitsF = 0;
  int noBitsT = 0;

  int aaDeltaIidT[NO_IID_BINS];
  int aaDeltaIccT[NO_ICC_BINS];

  int aaDeltaIidF[NO_IID_BINS];
  int aaDeltaIccF[NO_ICC_BINS];

  int abDtFlagIid;
  int abDtFlagIcc;

  int bSendHeader;

  unsigned int bZeroIid = 1;
  unsigned int bZeroIcc = 1;
  unsigned int bKeepParams = 1;

  HANDLE_BIT_BUF bb = &h_ps_e->psBitBuf;

  

     /* counting previous operations */

  
  temp = GetBitsAvail(bb);

  /* bit buffer shall be empty */
  
  if (temp != 0)
  {
    
    return -1;
  }

  
  if (bHeaderActive) {

    
    bKeepParams = 0;
  }

  
  lastIndex = 0;

   /* h_ps_e->aaaIIDDataBuffer[][]
                  h_ps_e->aaaICCDataBuffer[][]
                  h_ps_e->aLastIidIndex[]
               */
   
  for (gr = 0; gr < h_ps_e->iidIccBins; gr++) {
    float panValue = h_ps_e->aaaIIDDataBuffer[gr][SYSTEMLOOKAHEAD];

     /* counting previous operation */

     
    if (panValue >= -panClass[0] &&
        panValue <=  panClass[0]) {

      
      index = 0;
    }
    else  {

      
      if (panValue < 0) {

         /* panClass[] */
        
        for (index = NO_IID_STEPS-1; panValue > -panClass[index]; index--) {
           /* for() condition */
        }

        
        index = -index - 1;
      }
      else  {
         /* panClass[] */
        
        for (index = NO_IID_STEPS-1; panValue <  panClass[index]; index--) {
           /* for() condition */
        }

        
        index++;
      }

      
      bZeroIid = 0;
    }

    
    if (gr == 0) {

      
      aaDeltaIidF[gr] = index;
      noBitsT = 0;

       
      noBitsF = aBookPsIidFreqLength[index + CODE_BOOK_LAV_IID];
    }
    else {

       
      aaDeltaIidF[gr] = index - lastIndex;

       
      noBitsF += aBookPsIidFreqLength[aaDeltaIidF[gr] + CODE_BOOK_LAV_IID];
    }

    
    lastIndex = index;

     
    aaDeltaIidT[gr] = index - h_ps_e->aLastIidIndex[gr];

    
    h_ps_e->aLastIidIndex[gr] = index;

     
    noBitsT += aBookPsIidTimeLength[aaDeltaIidT[gr] + CODE_BOOK_LAV_IID];

    
    if (aaDeltaIidT[gr] != 0) {

      
      bKeepParams = 0;
    }
  } /* gr */

    
  if (noBitsT < noBitsF && !bHeaderActive) {

    
    abDtFlagIid    = 1;
    aaDeltaIid     = aaDeltaIidT;
    aaHuffBookIidC = aBookPsIidTimeCode;
    aaHuffBookIidL = aBookPsIidTimeLength;
  }
  else  {

    
    abDtFlagIid    = 0;
    aaDeltaIid     = aaDeltaIidF;
    aaHuffBookIidC = aBookPsIidFreqCode;
    aaHuffBookIidL = aBookPsIidFreqLength;
  }

  
  lastIndex = 0;

   /* h_ps_e->aaaICCDataBuffer[][]
                  h_ps_e->aLastIccIndex[]
                  aaDeltaIccF[]
                  aaDeltaIccT[]
               */
   
  for (gr = 0; gr < h_ps_e->iidIccBins; gr++) {

    float saValue = h_ps_e->aaaICCDataBuffer[gr][SYSTEMLOOKAHEAD];

     /* counting previous operation */

     
    if (saValue <= saClass[0]) {

      
      index = 0;
    }
    else  {
       /* saClass[] */
      
      for (index = NO_ICC_STEPS-2;  saValue < saClass[index]; index--) {
         /* for() condition */
      }

      
      index++;

      
      bZeroIcc = 0;
    }

    
    if (gr == 0) {

      
      aaDeltaIccF[gr] = index;

       
      noBitsF = aBookPsIccFreqLength[index + CODE_BOOK_LAV_ICC];

      
      noBitsT = 0;
    }
    else  {
       
      aaDeltaIccF[gr] = index - lastIndex;

       
      noBitsF += aBookPsIccFreqLength[aaDeltaIccF[gr] + CODE_BOOK_LAV_ICC];
    }

    
    lastIndex = index;

     
    aaDeltaIccT[gr] = index - h_ps_e->aLastIccIndex[gr];

    
    h_ps_e->aLastIccIndex[gr] = index;

     
    noBitsT += aBookPsIccTimeLength[aaDeltaIccT[gr] + CODE_BOOK_LAV_ICC];

    
    if (aaDeltaIccT[gr] != 0) {

      
      bKeepParams = 0;
    }
  } /* gr */

    
  if (noBitsT < noBitsF && !bHeaderActive) {

    
    abDtFlagIcc    = 1;
    aaDeltaIcc     = aaDeltaIccT;
    aaHuffBookIccC = aBookPsIccTimeCode;
    aaHuffBookIccL = aBookPsIccTimeLength;
  }
  else {

    
    abDtFlagIcc    = 0;
    aaDeltaIcc     = aaDeltaIccF;
    aaHuffBookIccC = aBookPsIccFreqCode;
    aaHuffBookIccL = aBookPsIccFreqLength;
  }

  {
    static int initheader = 0;

     /* counting previous operation */

     
    if (!initheader || bHeaderActive) {

       
      initheader = 1;
      h_ps_e->bEnableHeader = 1;
    }
    else {

       
      h_ps_e->bEnableHeader = 0;
    }
  }

    
  bSendHeader = h_ps_e->bEnableHeader            ||
                h_ps_e->bPrevZeroIid != bZeroIid ||
                h_ps_e->bPrevZeroIcc != bZeroIcc;

  
  WriteBits (bb, bSendHeader, 1);

  
  if (bSendHeader) {

     
    WriteBits (bb, !bZeroIid, 1);

    
    if (!bZeroIid)
    {
       
      WriteBits (bb, (h_ps_e->bHiFreqResIidIcc)?1:0, 3);
    }

     
    WriteBits (bb, !bZeroIcc, 1);

    
    if (!bZeroIcc)
    {
       
      WriteBits (bb, (h_ps_e->bHiFreqResIidIcc)?1:0, 3);
    }

    
    WriteBits (bb, 0, 1);
  }

  
  WriteBits (bb, 0, 1);

   
  WriteBits (bb, 1-bKeepParams, 2);

   /* h_ps_e->iidIccBins
                */
  
  if (!bKeepParams) {

     
    if (!bZeroIid) {

      
      WriteBits (bb, abDtFlagIid, 1);

      
      for (gr = 0; gr < h_ps_e->iidIccBins; gr++) {

         
        WriteBits (bb,
                   aaHuffBookIidC[aaDeltaIid[gr] + CODE_BOOK_LAV_IID],
                   aaHuffBookIidL[aaDeltaIid[gr] + CODE_BOOK_LAV_IID]);
      } /* gr */
    }  /* if (!bZeroIid) */
  }

   /* h_ps_e->iidIccBins
               */
  
  if (!bKeepParams) {

     
    if (!bZeroIcc) {

      
      WriteBits (bb, abDtFlagIcc, 1);

      
      for (gr = 0; gr < h_ps_e->iidIccBins; gr++) {

         
        WriteBits (bb,
                   aaHuffBookIccC[aaDeltaIcc[gr] + CODE_BOOK_LAV_ICC],
                   aaHuffBookIccL[aaDeltaIcc[gr] + CODE_BOOK_LAV_ICC]);
      } /* gr */
    }  /* if (!bZeroIcc) */
  }

   
  h_ps_e->bPrevZeroIid = bZeroIid;
  h_ps_e->bPrevZeroIcc = bZeroIcc;

   /* counting post-operation */

  

  return  GetBitsAvail(bb);

} /* writePsData */


/***************************************************************************/
/*!

  \brief  appends the parametric stereo bitstream portion to the output
          bitstream

  \return Number of bits needed for parametric stereo coding
          of 0 if no EXTENSION_ID_PS_CODING element should be transmitted

****************************************************************************/
int
AppendPsBS (HANDLE_PS_ENC h_ps_e,
            HANDLE_BIT_BUF   hBitStream,
            HANDLE_BIT_BUF   hBitStreamPrev,
            int*             sbrHdrBits)
{
  struct BIT_BUF bitbufTmp;
  unsigned char tmp[MAX_PAYLOAD_SIZE];

  

  
  if (!h_ps_e)
  {
    
    return 0;
  }

  
  if (!hBitStream) {

       /* counting post-operation */
    

    return GetBitsAvail (&h_ps_e->psBitBuf);
  }
  else {
    int writtenNoBits = 0;
    int maxExtSize = (1<<SI_SBR_EXTENSION_SIZE_BITS) - 1;
    int numBits = GetBitsAvail (&h_ps_e->psBitBuf);
    int extDataSize = (numBits+SI_SBR_EXTENSION_ID_BITS+7)>>3;

          /* counting previous operations */

     
    if ( GetBitsAvail(hBitStreamPrev) == 0) {

       
      h_ps_e->hdrBitsPrevFrame = *sbrHdrBits;

      
      CopyBitBuf(hBitStream, hBitStreamPrev);
    }
    else {

      int tmpBits;

       
      CreateBitBuffer (&bitbufTmp, tmp, sizeof(tmp));

      
      tmpBits = *sbrHdrBits;

       
      *sbrHdrBits = h_ps_e->hdrBitsPrevFrame;
      h_ps_e->hdrBitsPrevFrame = tmpBits;

      
      CopyBitBuf (hBitStreamPrev, &bitbufTmp);

      
      CopyBitBuf (hBitStream,  hBitStreamPrev);

      
      CopyBitBuf (&bitbufTmp,  hBitStream);
    }

    
    WriteBits (hBitStream, 1, SI_SBR_EXTENDED_DATA_BITS);

     
    if (extDataSize < maxExtSize) {

      
      WriteBits (hBitStream, extDataSize, SI_SBR_EXTENSION_SIZE_BITS);
    } else {
      
      WriteBits (hBitStream, maxExtSize, SI_SBR_EXTENSION_SIZE_BITS);

       
      WriteBits (hBitStream, extDataSize - maxExtSize, SI_SBR_EXTENSION_ESC_COUNT_BITS);
    }

     
    writtenNoBits += WriteBits (hBitStream, EXTENSION_ID_PS_CODING, SI_SBR_EXTENSION_ID_BITS);

       
    writtenNoBits += appendBitstream( hBitStream, &h_ps_e->psBitBuf, numBits );

    
    writtenNoBits = writtenNoBits%8;

    
    if(writtenNoBits)
    {
       
      WriteBits(hBitStream, 0, (8 - writtenNoBits));
    }

      /* counting post-operation */
    

    return GetBitsAvail(hBitStream)-(*sbrHdrBits)-SI_FILL_EXTENTION_BITS;
  }

  
}
