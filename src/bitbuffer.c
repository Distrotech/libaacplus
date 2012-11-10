/*
  Bit Buffer Management
*/

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

/*
  The pointer will be incremented if the parameter cnt is positive. Otherwise,
  the read pointer will be decremented.

  return  none
*/
static void updateBitBufWordPtr(HANDLE_BIT_BUF hBitBuf,       /* handle to bit buffer structure */
                                unsigned char **pBitBufWord,  /* pointer to pointer to bitsteam buffer */
                                long cnt)                     /* number of words */
{
  

   
  *pBitBufWord += cnt;
   
   /* hBitBuf->pBitBufEnd
                  hBitBuf->pBitBufBase
               */

   
  if(*pBitBufWord > hBitBuf->pBitBufEnd) {

     
    *pBitBufWord -= (hBitBuf->pBitBufEnd - hBitBuf->pBitBufBase + 1);
  }

   
  if(*pBitBufWord < hBitBuf->pBitBufBase) {

     
    *pBitBufWord += (hBitBuf->pBitBufEnd - hBitBuf->pBitBufBase + 1);
  }

  
}


/*
  Creates and initializes a bit buffer instance;

  returns pointer to bit buffer instance
*/
HANDLE_BIT_BUF CreateBitBuffer(HANDLE_BIT_BUF hBitBuf,      /* handle to bit buffer structure */
                               unsigned char *pBitBufBase,  /* pointer to bitstream buffer */
                               unsigned int bitBufSize)     /* size of bitstream buffer in words */
{
  

   
  hBitBuf->pBitBufBase = pBitBufBase;

    
  hBitBuf->pBitBufEnd = pBitBufBase+bitBufSize-1;

   
  hBitBuf->pReadNext = pBitBufBase;
  hBitBuf->pWriteNext = pBitBufBase;

   
  hBitBuf->wBitPos = 7;
  hBitBuf->rBitPos = 7;
  
   
  hBitBuf->cntBits = 0;

    
  hBitBuf->size    = bitBufSize * 8;

   
  hBitBuf->isValid = 1;

  

  return hBitBuf;
}


/*
  Deletes a bit buffer instance
*/
void DeleteBitBuffer(HANDLE_BIT_BUF *hBitBuf)  /* handle to pointer to bit buffer structure */
{
  

   
  (*hBitBuf)->isValid = 0;

  
  *hBitBuf = NULL;

  
}

/*
  Resets elements of a bit buffer instance
*/
void ResetBitBuf(HANDLE_BIT_BUF hBitBuf,      /* handle to bit buffer structure */
                 unsigned char *pBitBufBase,  /* pointer to bitstream buffer */
                 unsigned int bitBufSize)     /* size of bitstream buffer in words */
{
  

   
  hBitBuf->pBitBufBase = pBitBufBase;

    
  hBitBuf->pBitBufEnd = pBitBufBase+bitBufSize-1;

   
  hBitBuf->pReadNext = pBitBufBase;
  hBitBuf->pWriteNext = pBitBufBase;

   
  hBitBuf->rBitPos = 7;
  hBitBuf->wBitPos = 7;
    
   
  hBitBuf->cntBits = 0;

  
}


/*
  Copy source bit buffer instance to destination bit buffer instance
*/
void CopyBitBuf(HANDLE_BIT_BUF hBitBufSrc,  /* handle to source bit buffer structure */
                HANDLE_BIT_BUF hBitBufDst)  /* handle to destination bit buffer structure */
{
  int i;
  int bytesToGoSrc = (hBitBufSrc->pBitBufEnd - hBitBufSrc->pBitBufBase);
  int bytesToGoDst = (hBitBufDst->pBitBufEnd - hBitBufDst->pBitBufBase);

  

    /* counting previous operations */

  assert (bytesToGoSrc==bytesToGoDst);

   /* hBitBufDst->pBitBufBase[]
                  hBitBufSrc->pBitBufBase[]
               */
  
  for (i=0; i<bytesToGoSrc;i++)
  {
    
    hBitBufDst->pBitBufBase[i] = hBitBufSrc->pBitBufBase[i];
  }

   
  hBitBufDst->pReadNext  = hBitBufSrc->pReadNext;
  hBitBufDst->pWriteNext = hBitBufSrc->pWriteNext;

   
  hBitBufDst->rBitPos    = hBitBufSrc->rBitPos;
  hBitBufDst->wBitPos    = hBitBufSrc->wBitPos;
  
   
  hBitBufDst->cntBits    = hBitBufSrc->cntBits;

   
  hBitBufDst->isValid    = hBitBufSrc->isValid;

  
}


int GetBitsAvail(HANDLE_BIT_BUF hBitBuf)
{
  return hBitBuf->cntBits;
}
/*
  Read a certain number of bits from the bitstream buffer. The read direction is from left to right.
  The element cntBits will be decremented with the number of read bits. The element rBitPos will
  be set to the new "bit position"

  returns number of bits read
*/
unsigned long ReadBits(HANDLE_BIT_BUF hBitBuf,      /* handle to bit buffer structure */
                       unsigned char noBitsToRead)  /* number of bits to read */
{
  unsigned long returnValue;

  

  /* return value is of type unsigned int, it can hold up to 32 bits OETELAAR TODO wat is dit een 64 platform is??
     this optimized code can read upto 25 Bits a time*/
   
  if (noBitsToRead >= 25) {
    
    return 0;
  }

 
    
  hBitBuf->cntBits -= noBitsToRead;
  hBitBuf->rBitPos -= noBitsToRead;

   /* hBitBuf->rBitPos
                  hBitBuf->pReadNext
                  hBitBuf->pBitBufEnd
                  hBitBuf->pBitBufBase
               */

  
  returnValue = (unsigned long)*hBitBuf->pReadNext;
  
  
  while (hBitBuf->rBitPos < 0)
  {
     
    hBitBuf->rBitPos += 8;
    hBitBuf->pReadNext++;

     
    if(hBitBuf->pReadNext > hBitBuf->pBitBufEnd) {

      
      hBitBuf->pReadNext = hBitBuf->pBitBufBase;
    }

    
    returnValue <<= 8;

    
    returnValue |= (unsigned long)*hBitBuf->pReadNext;
  } 
    
   
  returnValue = returnValue << ((LongSize-1) - noBitsToRead - hBitBuf->rBitPos) >> (LongSize - noBitsToRead);

  
  return (returnValue);
}


/*
  Write a certain number of bits to the bitstream buffer. The write direction is from left to right.
  The element cntBits will be incremented with the number of written bits. It is actually irrelevant
  if the bits are really written to the bitstream buffer or only to the wCache.

  returns number of bits to write
*/
unsigned char WriteBits(HANDLE_BIT_BUF hBitBuf,       /* handle to bit buffer structure */
                        unsigned long writeValue,     /* write bits in word right bounded */
                        unsigned char noBitsToWrite)  /* number of bits to write */
{
  int bitsToWrite;
  unsigned char bitsWritten = noBitsToWrite;

  assert(noBitsToWrite <= LongSize);

  

   /* counting previous operation */

    
  hBitBuf->cntBits += noBitsToWrite;

  /* Bit Buffer Management: do not write more bits to input buffer than possible */
  assert ( hBitBuf->cntBits <= (hBitBuf->pBitBufEnd - hBitBuf->pBitBufBase + 1) * 8);
  
   /* hBitBuf->wBitPos
                  hBitBuf->pWriteNext
                  hBitBuf->pBitBufEnd
                  hBitBuf->pBitBufBase
               */
  
  while (noBitsToWrite) {
    unsigned char tmp,msk;

      
    bitsToWrite = min(hBitBuf->wBitPos + 1, noBitsToWrite);
	
     
    tmp = (unsigned char) ( writeValue << (LongSize - noBitsToWrite) >> (LongSize - bitsToWrite) << (hBitBuf->wBitPos + 1 - bitsToWrite) );

       /* (hBitBuf->wBitPos + 1 - bitsToWrite) --> already calculated */
    msk = ~(((1 << bitsToWrite) - 1) << (hBitBuf->wBitPos + 1 - bitsToWrite));

     
    *hBitBuf->pWriteNext &= msk;
    *hBitBuf->pWriteNext |= tmp;

     
    hBitBuf->wBitPos -= bitsToWrite;

    
    noBitsToWrite    -= bitsToWrite;

    
    if (hBitBuf->wBitPos < 0) {

       
      hBitBuf->wBitPos += 8;
      hBitBuf->pWriteNext++;

       
      if (hBitBuf->pWriteNext > hBitBuf->pBitBufEnd) {

        
        hBitBuf->pWriteNext= hBitBuf->pBitBufBase;
      }
    }

  }

  

  return(bitsWritten);
}


/*
  The read pointer will be winded a certain number of bits in forward or backward direction. The forward direction
  is chosen if the offset is positive, and the backward direction is chosen of the offset is negative. The pointer
  pReadNext and the elements rCache and rBitsLeft will be updated. The element cntBits will be decremented if the
  read pointer is updated in forward direction, and incremented if the read pointer is updated in backward direction.
*/
void WindBitBufferBidirectional(HANDLE_BIT_BUF hBitBuf,  /* handle to bit buffer structure */
                                long offset)             /* positive number => wind offset in forward direction
                                                            negative number => wind offset in backward direction */
{
  

  
  if (offset != 0)
  {
    int bOffset;
    
      
    hBitBuf->rBitPos -= offset;

    
    bOffset           = hBitBuf->rBitPos >> 3;

      
    hBitBuf->rBitPos -= bOffset << 3;

    
    if (bOffset) {

         
      updateBitBufWordPtr(hBitBuf, &hBitBuf->pReadNext, -bOffset);
    }

      
    hBitBuf->cntBits -= offset;
  }

  
}
