/*
  Remaining SBR Bit Writing Routines
*/

#include <stdlib.h>
#include <string.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

#ifndef min
#define min(a,b) ( a < b ? a:b)
#endif

#ifndef max
#define max(a,b) ( a > b ? a:b)
#endif

/* ***************************** crcAdvance **********************************/
/**
 * @brief    updates crc data register
 *
 */
static void crcAdvance(unsigned short crcPoly,
                       unsigned short crcMask,
                       unsigned short *crc,
                       unsigned long bValue,
                       int           bBits
                )
{
  int i;

  

   
  for (i=bBits-1; i>=0; i--)
  {
    unsigned short flag = (*crc) & crcMask ? 1:0;

       /* counting previous operation */

       /* .. ? */ 
    flag ^= (bValue & (1<<i) ? 1 : 0);

     
    (*crc)<<=1;

    
    if(flag)
    {
       
      (*crc) ^= crcPoly;
    }
  }

  
}


/* ***************************** InitSbrBitstream **********************************/
/**
 * @brief    Inittialisation of sbr bitstream, write of dummy header and CRC
 * @return   none
 *
 */

void InitSbrBitstream(HANDLE_COMMON_DATA  hCmonData,
                      unsigned char *memoryBase,
                      int            memorySize,
                      int CRCActive)
{
  

    
  ResetBitBuf(&hCmonData->sbrBitbuf,memoryBase,memorySize);

   
  hCmonData->tmpWriteBitbuf = hCmonData->sbrBitbuf;

    
  WriteBits (&hCmonData->sbrBitbuf, 0, SI_FILL_EXTENTION_BITS);

  
  if(CRCActive)
  {
      
    WriteBits (&hCmonData->sbrBitbuf, 0,SI_SBR_CRC_BITS);
  }



  
}


/* ************************** AssembleSbrBitstream *******************************/
/**
 * @brief    Formats the SBR payload
 * @return   nothing
 *
 */

void
AssembleSbrBitstream( HANDLE_COMMON_DATA  hCmonData)
{

  unsigned short crcReg =  SBR_CRCINIT;
  int numCrcBits,i;
  int sbrLoad=0;

  

   /* counting previous operations */

  
  if ( hCmonData==NULL )
  {
    
    return;
  }


   
  sbrLoad = hCmonData->sbrHdrBits + hCmonData->sbrDataBits;

  
  sbrLoad += SI_FILL_EXTENTION_BITS;

   
  if ( hCmonData->sbrCrcLen )
  {
    
    sbrLoad += SI_SBR_CRC_BITS;
  }

     
  hCmonData->sbrFillBits = (8 - (sbrLoad) % 8) % 8;

  
  sbrLoad += hCmonData->sbrFillBits;

   
  WriteBits(&hCmonData->sbrBitbuf, 0,  hCmonData->sbrFillBits );


   
  if ( hCmonData->sbrCrcLen ){
    struct BIT_BUF  tmpCRCBuf = hCmonData->sbrBitbuf;

     /* counting previous operation */

   
  ReadBits (&tmpCRCBuf, SI_FILL_EXTENTION_BITS);

   
  ReadBits (&tmpCRCBuf, SI_SBR_CRC_BITS);


     
    numCrcBits = hCmonData->sbrHdrBits + hCmonData->sbrDataBits + hCmonData->sbrFillBits;

    
    for(i=0;i<numCrcBits;i++){
      int bit;

       
      bit=ReadBits(&tmpCRCBuf,1);

       
      crcAdvance(SBR_CRC_POLY,SBR_CRC_MASK,&crcReg,bit,1);
    }

    
    crcReg &= (SBR_CRC_RANGE);
  }


     
    if ( hCmonData->sbrCrcLen ) {

        
      WriteBits (&hCmonData->tmpWriteBitbuf, AAC_SI_FIL_SBR_CRC, SI_FILL_EXTENTION_BITS);

      
      WriteBits (&hCmonData->tmpWriteBitbuf, crcReg,SI_SBR_CRC_BITS);

    }
    else {

        
      WriteBits (&hCmonData->tmpWriteBitbuf, AAC_SI_FIL_SBR, SI_FILL_EXTENTION_BITS);
    }




    
 }

