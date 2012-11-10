/*
 Huffman Bitcounter & coder
*/
#include <stdlib.h>

#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


#define HI_LTAB(a) (a>>8)
#define LO_LTAB(a) (a & 0xff)

/*
  needed for fast bit counter, since
  huffman cw length tables are packed into 16 
  bit now
*/
#define HI_EXPLTAB(a) (a>>16)
#define LO_EXPLTAB(a) (a & 0xffff)
#define EXPAND(a)  ((((int)(a&0xff00))<<8)|(int)(a&0xff)) 


/*****************************************************************************


    functionname: count1_2_3_4_5_6_7_8_9_10_11
    description:  counts tables 1-11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 1-11

*****************************************************************************/

static void count1_2_3_4_5_6_7_8_9_10_11(const short *values,
                                         const int  width,
                                         int       *bitCount)
{

  int i;
  int bc1_2,bc3_4,bc5_6,bc7_8,bc9_10,bc11,sc;
  int t0,t1,t2,t3;

  

  
  bc1_2=0;
  bc3_4=0;
  bc5_6=0;
  bc7_8=0;
  bc9_10=0;
  bc11=0;
  sc=0;

   /* pointer for values[] */
  
  for(i=0;i<width;i+=4){
    
    
    t0= values[i+0];
    t1= values[i+1];
    t2= values[i+2];
    t3= values[i+3];
  
    /* 1,2 */

      /* EXPAND() */  
    bc1_2+=EXPAND(huff_ltab1_2[t0+1][t1+1][t2+1][t3+1]);

    /* 5,6 */
      /* EXPAND() */  
    bc5_6+=EXPAND(huff_ltab5_6[t0+4][t1+4]);

      /* EXPAND() */  
    bc5_6+=EXPAND(huff_ltab5_6[t2+4][t3+4]);

    
    t0=abs(t0);
    t1=abs(t1);
    t2=abs(t2);
    t3=abs(t3);

    
      /* EXPAND() */  
    bc3_4+= EXPAND(huff_ltab3_4[t0][t1][t2][t3]);
    
      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t0][t1]);

      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t2][t3]);
    
      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t0][t1]);

      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t2][t3]);
    
      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t2][t3]);
   
     
    sc+=(t0>0)+(t1>0)+(t2>0)+(t3>0);
  }
  
   /* HI_EXPLTAB() */ 
  bitCount[1]=HI_EXPLTAB(bc1_2);

   /* LO_EXPLTAB() */ 
  bitCount[2]=LO_EXPLTAB(bc1_2);

   /* HI_EXPLTAB() */  
  bitCount[3]=HI_EXPLTAB(bc3_4)+sc;

   /* LO_EXPLTAB() */  
  bitCount[4]=LO_EXPLTAB(bc3_4)+sc;

   /* HI_EXPLTAB() */ 
  bitCount[5]=HI_EXPLTAB(bc5_6);

   /* LO_EXPLTAB() */ 
  bitCount[6]=LO_EXPLTAB(bc5_6);

   /* HI_EXPLTAB() */  
  bitCount[7]=HI_EXPLTAB(bc7_8)+sc;

   /* LO_EXPLTAB() */  
  bitCount[8]=LO_EXPLTAB(bc7_8)+sc;

   /* HI_EXPLTAB() */  
  bitCount[9]=HI_EXPLTAB(bc9_10)+sc;

   /* LO_EXPLTAB() */  
  bitCount[10]=LO_EXPLTAB(bc9_10)+sc;

   
  bitCount[11]=bc11+sc;
  
  
}


/*****************************************************************************

    functionname: count3_4_5_6_7_8_9_10_11
    description:  counts tables 3-11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 3-11

*****************************************************************************/

static void count3_4_5_6_7_8_9_10_11(const short *values,
                                     const int  width,
                                     int       *bitCount)
{

  int i;
  int bc3_4,bc5_6,bc7_8,bc9_10,bc11,sc;
  int t0,t1,t2,t3;
  
  

  
  bc3_4=0;
  bc5_6=0;
  bc7_8=0;
  bc9_10=0;
  bc11=0;
  sc=0;

   /* pointer for values[] */
  
  for(i=0;i<width;i+=4){

    
    t0= values[i+0];
    t1= values[i+1];
    t2= values[i+2];
    t3= values[i+3];
    
    /*
      5,6
    */
      /* EXPAND() */  
    bc5_6+=EXPAND(huff_ltab5_6[t0+4][t1+4]);

      /* EXPAND() */  
    bc5_6+=EXPAND(huff_ltab5_6[t2+4][t3+4]);

    
    t0=abs(t0);
    t1=abs(t1);
    t2=abs(t2);
    t3=abs(t3);


      /* EXPAND() */  
    bc3_4+= EXPAND(huff_ltab3_4[t0][t1][t2][t3]);
    
      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t0][t1]);

      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t2][t3]);
    
      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t0][t1]);

      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t2][t3]);
    
      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t2][t3]);

     
    sc+=(t0>0)+(t1>0)+(t2>0)+(t3>0);
      
   
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;

   /* HI_EXPLTAB() */  
  bitCount[3]=HI_EXPLTAB(bc3_4)+sc;

   /* LO_EXPLTAB() */  
  bitCount[4]=LO_EXPLTAB(bc3_4)+sc;

   /* HI_EXPLTAB() */ 
  bitCount[5]=HI_EXPLTAB(bc5_6);

   /* LO_EXPLTAB() */ 
  bitCount[6]=LO_EXPLTAB(bc5_6);

   /* HI_EXPLTAB() */  
  bitCount[7]=HI_EXPLTAB(bc7_8)+sc;

   /* LO_EXPLTAB() */  
  bitCount[8]=LO_EXPLTAB(bc7_8)+sc;

   /* HI_EXPLTAB() */  
  bitCount[9]=HI_EXPLTAB(bc9_10)+sc;

   /* LO_EXPLTAB() */  
  bitCount[10]=LO_EXPLTAB(bc9_10)+sc;

   
  bitCount[11]=bc11+sc;
  
  
}



/*****************************************************************************

    functionname: count5_6_7_8_9_10_11
    description:  counts tables 5-11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 5-11

*****************************************************************************/


static void count5_6_7_8_9_10_11(const short *values,
                                 const int  width,
                                 int       *bitCount)
{

  int i;
  int bc5_6,bc7_8,bc9_10,bc11,sc;
  int t0,t1;

  

  
  bc5_6=0;
  bc7_8=0;
  bc9_10=0;
  bc11=0;
  sc=0;

   /* pointer for values[] */
  
  for(i=0;i<width;i+=2){

    
    t0 = values[i+0];
    t1 = values[i+1];

      /* EXPAND() */  
    bc5_6+=EXPAND(huff_ltab5_6[t0+4][t1+4]);

    
    t0=abs(t0);
    t1=abs(t1);
     
      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t0][t1]);

      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t0][t1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);
    
     
    sc+=(t0>0)+(t1>0);
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;
  bitCount[3]=INVALID_BITCOUNT;
  bitCount[4]=INVALID_BITCOUNT;

   /* HI_EXPLTAB() */ 
  bitCount[5]=HI_EXPLTAB(bc5_6);

   /* LO_EXPLTAB() */ 
  bitCount[6]=LO_EXPLTAB(bc5_6);

   /* HI_EXPLTAB() */  
  bitCount[7]=HI_EXPLTAB(bc7_8)+sc;

   /* LO_EXPLTAB() */  
  bitCount[8]=LO_EXPLTAB(bc7_8)+sc;

   /* HI_EXPLTAB() */  
  bitCount[9]=HI_EXPLTAB(bc9_10)+sc;

   /* LO_EXPLTAB() */  
  bitCount[10]=LO_EXPLTAB(bc9_10)+sc;

   
  bitCount[11]=bc11+sc;
  
  
}


/*****************************************************************************

    functionname: count7_8_9_10_11
    description:  counts tables 7-11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 7-11

*****************************************************************************/

static void count7_8_9_10_11(const short *values,
                             const int  width,
                             int       *bitCount)
{

  int i;
  int bc7_8,bc9_10,bc11,sc;
  int t0,t1;

  
  
  
  bc7_8=0;
  bc9_10=0;
  bc11=0;
  sc=0;

   /* pointer for values[] */
  
  for(i=0;i<width;i+=2){

    
    
    t0=abs(values[i+0]);
    t1=abs(values[i+1]);

      /* EXPAND() */  
    bc7_8+=EXPAND(huff_ltab7_8[t0][t1]);

      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t0][t1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);

     
    sc+=(t0>0)+(t1>0);
   
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;
  bitCount[3]=INVALID_BITCOUNT;
  bitCount[4]=INVALID_BITCOUNT;
  bitCount[5]=INVALID_BITCOUNT;
  bitCount[6]=INVALID_BITCOUNT;

   /* HI_EXPLTAB() */  
  bitCount[7]=HI_EXPLTAB(bc7_8)+sc;

   /* LO_EXPLTAB() */  
  bitCount[8]=LO_EXPLTAB(bc7_8)+sc;

   /* HI_EXPLTAB() */  
  bitCount[9]=HI_EXPLTAB(bc9_10)+sc;

   /* LO_EXPLTAB() */  
  bitCount[10]=LO_EXPLTAB(bc9_10)+sc;

   
  bitCount[11]=bc11+sc;
  
  
}

/*****************************************************************************

    functionname: count9_10_11
    description:  counts tables 9-11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 9-11

*****************************************************************************/



static void count9_10_11(const short *values,
                         const int  width,
                         int       *bitCount)
{

  int i;
  int bc9_10,bc11,sc;
  int t0,t1;

  

  
  bc9_10=0;
  bc11=0;
  sc=0;

   /* pointer for values[] */
  
  for(i=0;i<width;i+=2){


    
    t0=abs(values[i+0]);
    t1=abs(values[i+1]);
    

      /* EXPAND() */  
    bc9_10+=EXPAND(huff_ltab9_10[t0][t1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);

     
    sc+=(t0>0)+(t1>0);
   
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;
  bitCount[3]=INVALID_BITCOUNT;
  bitCount[4]=INVALID_BITCOUNT;
  bitCount[5]=INVALID_BITCOUNT;
  bitCount[6]=INVALID_BITCOUNT;
  bitCount[7]=INVALID_BITCOUNT;
  bitCount[8]=INVALID_BITCOUNT;

   /* HI_EXPLTAB() */  
  bitCount[9]=HI_EXPLTAB(bc9_10)+sc;

   /* LO_EXPLTAB() */  
  bitCount[10]=LO_EXPLTAB(bc9_10)+sc;

   
  bitCount[11]=bc11+sc;
  
  
}
 
/*****************************************************************************

    functionname: count11
    description:  counts table 11 
    returns:      
    input:        quantized spectrum
    output:       bitCount for table 11

*****************************************************************************/
 
static void count11(const short *values,
                    const int  width,
                    int        *bitCount)
{

  int i;
  int bc11,sc;
  int t0,t1;

  

  
  bc11=0;
  sc=0;

   /* pointer for values[] */
   
  for(i=0;i<width;i+=2){

    
    t0=abs(values[i+0]);
    t1=abs(values[i+1]);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t0][t1]);

     
    sc+=(t0>0)+(t1>0);
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;
  bitCount[3]=INVALID_BITCOUNT;
  bitCount[4]=INVALID_BITCOUNT;
  bitCount[5]=INVALID_BITCOUNT;
  bitCount[6]=INVALID_BITCOUNT;
  bitCount[7]=INVALID_BITCOUNT;
  bitCount[8]=INVALID_BITCOUNT;
  bitCount[9]=INVALID_BITCOUNT;
  bitCount[10]=INVALID_BITCOUNT;

   
  bitCount[11]=bc11+sc;

  
}

/*****************************************************************************

    functionname: countEsc
    description:  counts table 11 (with Esc) 
    returns:      
    input:        quantized spectrum
    output:       bitCount for tables 11 (with Esc)

*****************************************************************************/

static void countEsc(const short *values,
                     const int  width,
                     int       *bitCount)
{

  int i;
  int bc11,ec,sc;
  int t0,t1,t00,t01;

  

  
  bc11=0;
  sc=0;
  ec=0;

   /* pointer for values[] */
   
  for(i=0;i<width;i+=2){

    
    t0=abs(values[i+0]);
    t1=abs(values[i+1]);
    
     
    sc+=(t0>0)+(t1>0);

      
    t00 = min(t0,16);
    t01 = min(t1,16);

      /* EXPAND() */  
    bc11+=EXPAND(huff_ltab11[t00][t01]);
    
     
    if(t0>=16){

      
      ec+=5;

         /* while() condition */
      while((t0>>=1) >= 16)
      {
        
        ec+=2;
      }
    }
    
     
    if(t1>=16){

      
      ec+=5;

         /* while() condition */
      while((t1>>=1) >= 16)
      {

        
        ec+=2;
      }
    }
  }

  
  bitCount[1]=INVALID_BITCOUNT;
  bitCount[2]=INVALID_BITCOUNT;
  bitCount[3]=INVALID_BITCOUNT;
  bitCount[4]=INVALID_BITCOUNT;
  bitCount[5]=INVALID_BITCOUNT;
  bitCount[6]=INVALID_BITCOUNT;
  bitCount[7]=INVALID_BITCOUNT;
  bitCount[8]=INVALID_BITCOUNT;
  bitCount[9]=INVALID_BITCOUNT;
  bitCount[10]=INVALID_BITCOUNT;

   
  bitCount[11]=bc11+sc+ec;

  
}


typedef void (*COUNT_FUNCTION)(const short *values,
                               const int  width,
                               int       *bitCount);

static COUNT_FUNCTION countFuncTable[CODE_BOOK_ESC_LAV+1] =
{

 count1_2_3_4_5_6_7_8_9_10_11,  /* 0  */
 count1_2_3_4_5_6_7_8_9_10_11,  /* 1  */
 count3_4_5_6_7_8_9_10_11,      /* 2  */
 count5_6_7_8_9_10_11,          /* 3  */
 count5_6_7_8_9_10_11,          /* 4  */
 count7_8_9_10_11,              /* 5  */
 count7_8_9_10_11,              /* 6  */
 count7_8_9_10_11,              /* 7  */
 count9_10_11,                  /* 8  */
 count9_10_11,                  /* 9  */
 count9_10_11,                  /* 10 */
 count9_10_11,                  /* 11 */
 count9_10_11,                  /* 12 */
 count11,                       /* 13 */
 count11,                       /* 14 */
 count11,                       /* 15 */
 countEsc                       /* 16 */
};



int    bitCount(const short      *values,
                const int         width,
                      int         maxVal,
                int              *bitCount)
{

  /*
    check if we can use codebook 0
  */

  

   
  if(maxVal == 0)
    bitCount[0] = 0;
  else
    bitCount[0] = INVALID_BITCOUNT;

    
  maxVal = min(maxVal,CODE_BOOK_ESC_LAV);

   
  countFuncTable[maxVal](values,width,bitCount);

  
  return(0);
}







int codeValues(short *values, int width, int codeBook,  HANDLE_BIT_BUF hBitstream)
{
  int i,t0,t1,t2,t3,t00,t01;
  int codeWord,codeLength;
  int sign,signLength;

  


  
  switch(codeBook){

  case CODE_BOOK_ZERO_NO:
    break;

  case CODE_BOOK_1_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=4) {

      
      t0         = values[i+0];
      t1         = values[i+1];
      t2         = values[i+2];
      t3         = values[i+3];

      
      codeWord   = huff_ctab1[t0+1][t1+1][t2+1][t3+1];

       /* HI_LTAB() */ 
      codeLength = HI_LTAB(huff_ltab1_2[t0+1][t1+1][t2+1][t3+1]);

      
      WriteBits(hBitstream,codeWord,codeLength);        
    }
    break;

  case CODE_BOOK_2_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=4) {

      
      t0         = values[i+0];
      t1         = values[i+1];
      t2         = values[i+2];
      t3         = values[i+3];

      
      codeWord   = huff_ctab2[t0+1][t1+1][t2+1][t3+1];

       /* LO_TAB() */ 
      codeLength = LO_LTAB(huff_ltab1_2[t0+1][t1+1][t2+1][t3+1]);

      
      WriteBits(hBitstream,codeWord,codeLength);
      
    }
    break;

  case CODE_BOOK_3_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=4) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }

      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      t2 = values[i+2];

      
      if(t2 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t2 < 0){

          
          sign|=1;

          
          t2=abs(t2);
        }
      }

      
      t3 = values[i+3];

      
      if(t3 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t3 < 0){

          
          sign|=1;

          
          t3=abs(t3);
        }
      }
      
      
      codeWord   = huff_ctab3[t0][t1][t2][t3];

       /* HI_LTAB() */ 
      codeLength = HI_LTAB(huff_ltab3_4[t0][t1][t2][t3]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_4_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=4) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }

      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      t2 = values[i+2];

      
      if(t2 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t2 < 0){

          
          sign|=1;

          
          t2=abs(t2);
        }
      }

      
      t3 = values[i+3];

      
      if(t3 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t3 < 0){

          
          sign|=1;

          
          t3=abs(t3);
        }
      }

      
      codeWord   = huff_ctab4[t0][t1][t2][t3];

       /* LO_TAB() */ 
      codeLength = LO_LTAB(huff_ltab3_4[t0][t1][t2][t3]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_5_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      t0         = values[i+0];
      t1         = values[i+1];

      
      codeWord   = huff_ctab5[t0+4][t1+4];

       /* HI_LTAB() */ 
      codeLength = HI_LTAB(huff_ltab5_6[t0+4][t1+4]);

      
      WriteBits(hBitstream,codeWord,codeLength);
    

    }
    break;

  case CODE_BOOK_6_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      t0         = values[i+0];
      t1         = values[i+1];

      
      codeWord   = huff_ctab6[t0+4][t1+4];

       /* LO_TAB() */ 
      codeLength = LO_LTAB(huff_ltab5_6[t0+4][t1+4]);

      
      WriteBits(hBitstream,codeWord,codeLength);
      
    

    }
    break;

  case CODE_BOOK_7_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2){

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }
      
      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      codeWord   = huff_ctab7[t0][t1];

       /* HI_LTAB() */ 
      codeLength = HI_LTAB(huff_ltab7_8[t0][t1]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_8_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }
      
      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      codeWord   = huff_ctab8[t0][t1];

       /* LO_TAB() */ 
      codeLength = LO_LTAB(huff_ltab7_8[t0][t1]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_9_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }

      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      codeWord   = huff_ctab9[t0][t1];
    
       /* HI_LTAB() */ 
      codeLength = HI_LTAB(huff_ltab9_10[t0][t1]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_10_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }
      
      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

      
      codeWord   = huff_ctab10[t0][t1];

       /* LO_TAB() */ 
      codeLength = LO_LTAB(huff_ltab9_10[t0][t1]);

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);
    }
    break;

  case CODE_BOOK_ESC_NO:

     /* pointer for values[] */
    
    for(i=0; i<width; i+=2) {

      
      sign=0;
      signLength=0;
      t0 = values[i+0];

      
      if(t0 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t0 < 0){

          
          sign|=1;

          
          t0=abs(t0);
        }
      }

      
      t1 = values[i+1];

      
      if(t1 != 0){

        
        signLength++;

        
        sign<<=1;

        
        if(t1 < 0){

          
          sign|=1;

          
          t1=abs(t1);
        }
      }

        
      t00 = min(t0,16);
      t01 = min(t1,16);
      
       
      codeWord   = huff_ctab11[t00][t01];
      codeLength = huff_ltab11[t00][t01];

      
      WriteBits(hBitstream,codeWord,codeLength);

      
      WriteBits(hBitstream,sign,signLength);

       
      if(t0 >=16){
        int n,p;

        
        n=0;
        p=t0;

          
        while((p>>=1) >=16){

          
          WriteBits(hBitstream,1,1);

          
          n++;
        }

        
        WriteBits(hBitstream,0,1);

          
        WriteBits(hBitstream,t0-(1<<(n+4)),n+4);
      }
      if(t1 >=16){
        int n,p;

        
        n=0;
        p=t1;

          
        while((p>>=1) >=16){

          
          WriteBits(hBitstream,1,1);

          
          n++;
        }

        
        WriteBits(hBitstream,0,1);

          
        WriteBits(hBitstream,t1-(1<<(n+4)),n+4);
      }
    }
    break;
  
  default:
    break;
  }

  

  return(0);
}

int bitCountScalefactorDelta(signed int delta)
{
   
   
   

  return(huff_ltabscf[delta+CODE_BOOK_SCF_LAV]);
}

int codeScalefactorDelta(signed int delta, HANDLE_BIT_BUF hBitstream)
{

  int codeWord,codeLength;

  
  
    
  if(labs(delta) >CODE_BOOK_SCF_LAV)
  {
    
    return(1);
  }
  
   
  codeWord   = huff_ctabscf[delta+CODE_BOOK_SCF_LAV];
  codeLength = huff_ltabscf[delta+CODE_BOOK_SCF_LAV];

  
  WriteBits(hBitstream,codeWord,codeLength);

  

  return(0);
}



