/*
  QMF analysis
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static void fct3_4(float *x)
{
  float tmp00, tmp01, tmp10, tmp11, wc, ws;

  

   
  x[1] *= 0.7071067812f;

  
  tmp00 = x[0] + x[1];
  tmp01 = x[0] - x[1];

  wc = 0.923879532511287f;
  ws = 0.382683432365090f;

    
  tmp10 = x[2] * wc + x[3] * ws;
  tmp11 = x[2] * ws - x[3] * wc;

   
  x[0] = tmp00 + tmp10;
  x[3] = tmp00 - tmp10;
  x[1] = tmp01 + tmp11;
  x[2] = tmp01 - tmp11;

  
}

static void fst3_4r(float *x)
{
  float tmp00, tmp01, tmp10, tmp11, wc, ws;

  

   
  x[2] *= 0.7071067812f;

  
  tmp00 = x[3] + x[2];
  tmp01 = x[3] - x[2];

  wc = 0.923879532511287f;
  ws = 0.382683432365090f;

    
  tmp10 = x[1] * wc + x[0] * ws;
  tmp11 = x[0] * wc - x[1] * ws;

   
  x[3] = tmp00 + tmp10;
  x[0] = tmp10 - tmp00;
  x[2] = tmp11 - tmp01;
  x[1] = tmp11 + tmp01;

  
}


static void fct4_4r(float *x)
{
  float tmp00, tmp01, tmp10, tmp11, wc, ws;

  

  wc = 0.7071067812f;

   
  x[1] *= wc;

  
  tmp00 = x[0] + x[1];
  tmp01 = x[0] - x[1];

   
  x[2] *= wc;

  
  tmp11 = x[3] - x[2];
  tmp10 = x[3] + x[2];

  wc = 0.980785280403230f;
  ws = 0.195090322016128f;

     
  x[3] = tmp00 * wc + tmp10 * ws;
  x[0] = tmp00 * ws - tmp10 * wc;

  wc = 0.831469612302545f;
  ws = 0.555570233019602f;

     
  x[2] = tmp01 * wc - tmp11 * ws;
  x[1] = tmp01 * ws + tmp11 * wc;

  
}

static void fst4_4(float *x)
{
  float tmp00, tmp01, tmp10, tmp11, wc, ws;

  

  wc = 0.7071067812f;

   
  x[1] *= wc;

  
  tmp10 = x[0] + x[1];
  tmp11 = x[0] - x[1];

   
  x[2] *= wc;

  
  tmp01 = x[3] - x[2];
  tmp00 = x[3] + x[2];

  wc = 0.980785280403230f;
  ws = 0.195090322016128f;

     
  x[0] = tmp00 * wc + tmp10 * ws;
  x[3] = tmp10 * wc - tmp00 * ws;

  wc = 0.831469612302545f;
  ws = 0.555570233019602f;

     
  x[1] = tmp11 * ws - tmp01 * wc;
  x[2] = tmp01 * ws + tmp11 * wc;

  
}

static void fct3_64(float *a)
{
  int k;
  const float *tPtr;
  float wc, ws;
  float xp;

  

  /* bit reversal */
  
  xp =  a[1];  a[1] = a[32]; a[32] = xp;
  xp =  a[2];  a[2] = a[16]; a[16] = xp;
  xp =  a[3];  a[3] = a[48]; a[48] = xp;
  xp =  a[4];  a[4] =  a[8];  a[8] = xp;
  xp =  a[5];  a[5] = a[40]; a[40] = xp;
  xp =  a[6];  a[6] = a[24]; a[24] = xp;
  xp =  a[7];  a[7] = a[56]; a[56] = xp;
  xp =  a[9];  a[9] = a[36]; a[36] = xp;
  xp = a[10]; a[10] = a[20]; a[20] = xp;
  xp = a[11]; a[11] = a[52]; a[52] = xp;
  xp = a[13]; a[13] = a[44]; a[44] = xp;
  xp = a[14]; a[14] = a[28]; a[28] = xp;
  xp = a[15]; a[15] = a[60]; a[60] = xp;
  xp = a[17]; a[17] = a[34]; a[34] = xp;
  xp = a[19]; a[19] = a[50]; a[50] = xp;
  xp = a[21]; a[21] = a[42]; a[42] = xp;
  xp = a[22]; a[22] = a[26]; a[26] = xp;
  xp = a[23]; a[23] = a[58]; a[58] = xp;
  xp = a[25]; a[25] = a[38]; a[38] = xp;
  xp = a[27]; a[27] = a[54]; a[54] = xp;
  xp = a[29]; a[29] = a[46]; a[46] = xp;
  xp = a[31]; a[31] = a[62]; a[62] = xp;
  xp = a[35]; a[35] = a[49]; a[49] = xp;
  xp = a[37]; a[37] = a[41]; a[41] = xp;
  xp = a[39]; a[39] = a[57]; a[57] = xp;
  xp = a[43]; a[43] = a[53]; a[53] = xp;
  xp = a[47]; a[47] = a[61]; a[61] = xp;
  xp = a[55]; a[55] = a[59]; a[59] = xp;


    

  xp = a[33] + a[62]; a[62] -= a[33]; a[33] = xp;

  xp = a[34] + a[60]; a[60] -= a[34]; a[34] = xp;
  xp = a[35] + a[61]; a[61] -= a[35]; a[35] = xp;

  xp = a[36] + a[56]; a[56] -= a[36]; a[36] = xp;
  xp = a[37] + a[57]; a[57] -= a[37]; a[37] = xp;
  xp = a[38] + a[58]; a[58] -= a[38]; a[38] = xp;
  xp = a[39] + a[59]; a[59] -= a[39]; a[39] = xp;

  xp = a[40] + a[48]; a[48] -= a[40]; a[40] = xp;
  xp = a[41] + a[49]; a[49] -= a[41]; a[41] = xp;
  xp = a[42] + a[50]; a[50] -= a[42]; a[42] = xp;
  xp = a[43] + a[51]; a[51] -= a[43]; a[43] = xp;
  xp = a[44] + a[52]; a[52] -= a[44]; a[44] = xp;
  xp = a[45] + a[53]; a[53] -= a[45]; a[45] = xp;
  xp = a[46] + a[54]; a[54] -= a[46]; a[46] = xp;
  xp = a[47] + a[55]; a[55] -= a[47]; a[47] = xp;


    

  xp = a[17] + a[30]; a[30] -= a[17]; a[17] = xp;

  xp = a[18] + a[28]; a[28] -= a[18]; a[18] = xp;
  xp = a[19] + a[29]; a[29] -= a[19]; a[19] = xp;

  xp = a[20] + a[24]; a[24] -= a[20]; a[20] = xp;
  xp = a[21] + a[25]; a[25] -= a[21]; a[21] = xp;
  xp = a[22] + a[26]; a[26] -= a[22]; a[22] = xp;
  xp = a[23] + a[27]; a[27] -= a[23]; a[23] = xp;


    


  xp =  a[9] + a[14]; a[14] -=  a[9];  a[9] = xp;
  xp = a[10] + a[12]; a[12] -= a[10]; a[10] = xp;
  xp = a[11] + a[13]; a[13] -= a[11]; a[11] = xp;

  xp = a[41] + a[46]; a[46] -= a[41]; a[41] = xp;
  xp = a[42] + a[44]; a[44] -= a[42]; a[42] = xp;
  xp = a[43] + a[45]; a[45] -= a[43]; a[43] = xp;

  xp = a[49] + a[54]; a[49] -= a[54]; a[54] = xp;
  xp = a[50] + a[52]; a[50] -= a[52]; a[52] = xp;
  xp = a[51] + a[53]; a[51] -= a[53]; a[53] = xp;


    


  xp =  a[5] +  a[6]; a[6]  -=  a[5];  a[5] = xp;
  xp = a[21] + a[22]; a[22] -= a[21]; a[21] = xp;
  xp = a[25] + a[26]; a[25] -= a[26]; a[26] = xp;
  xp = a[37] + a[38]; a[38] -= a[37]; a[37] = xp;
  xp = a[57] + a[58]; a[57] -= a[58]; a[58] = xp;

  
  fct3_4(a);
   
  fct4_4r(a+4);
   
  fct3_4(a+8);
   
  fst3_4r(a+12);

   
  fct3_4(a+16);
   
  fct4_4r(a+20);
   
  fst4_4(a+24);
   
  fst3_4r(a+28);

   
  fct3_4(a+32);
   
  fct4_4r(a+36);
   
  fct3_4(a+40);
   
  fst3_4r(a+44);

   
  fct3_4(a+48);
   
  fst3_4r(a+52);
   
  fst4_4(a+56);
   
  fst3_4r(a+60);


  
  for (k=0; k<4; k++) {

      

    xp      =    a[k] +  a[7-k];
    a[7-k]  =    a[k] -  a[7-k];
    a[k]    = xp;
    xp      = a[16+k] + a[23-k];
    a[23-k] = a[16+k] - a[23-k];
    a[16+k] = xp;
    xp      = a[24+k] + a[31-k];
    a[24+k] = a[24+k] - a[31-k];
    a[31-k] = xp;
    xp      = a[32+k] + a[39-k];
    a[39-k] = a[32+k] - a[39-k];
    a[32+k] = xp;
    xp      = a[56+k] + a[63-k];
    a[56+k] = a[56+k] - a[63-k];
    a[63-k] = xp;
  }

  tPtr = trigData_fct4_8;

  
  for (k=0; k<4; k++) {
    wc = *tPtr++;
    ws = *tPtr++;

        

    xp      =  a[8+k]*wc + a[15-k]*ws;
    a[8+k]  =  a[8+k]*ws - a[15-k]*wc;
    a[15-k] = xp;
    xp      = a[40+k]*wc + a[47-k]*ws;
    a[40+k] = a[40+k]*ws - a[47-k]*wc;
    a[47-k] = xp;
    xp      = a[48+k]*ws + a[55-k]*wc;
    a[55-k] = a[48+k]*wc - a[55-k]*ws;
    a[48+k] = xp;
  }

  
  for (k=0; k<8; k++) {

      

    xp      =    a[k] + a[15-k];
    a[15-k] =    a[k] - a[15-k];
    a[k]    = xp;
    xp      = a[32+k] + a[47-k];
    a[47-k] = a[32+k] - a[47-k];
    a[32+k] = xp;
    xp      = a[48+k] + a[63-k];
    a[48+k] = a[48+k] - a[63-k];
    a[63-k] = xp;
  }

  tPtr = trigData_fct4_16;

  
  for (k=0; k<8; k++) {
    wc = *tPtr++;
    ws = *tPtr++;

        

    xp      = a[16+k]*wc + a[31-k]*ws;
    a[16+k] = a[16+k]*ws - a[31-k]*wc;
    a[31-k] = xp;
  }


  
  for (k=0; k<16; k++) {

      

    xp      = a[k] + a[31-k];
    a[31-k] = a[k] - a[31-k];
    a[k]    = xp;
  }

  tPtr = trigData_fct4_32;

  
  for (k=0; k<16; k++) {
    wc = *tPtr++;
    ws = *tPtr++;

        

    xp = a[32+k]*wc + a[63-k]*ws;
    a[32+k] = a[32+k]*ws - a[63-k]*wc;
    a[63-k] = xp;
  }


  
  for (k=0; k<32; k++) {

      

    xp      = a[k] + a[63-k];
    a[63-k] = a[k] - a[63-k];
    a[k]    = xp;
  }

  
}

static void fst3_64(float *a)
{
  int k;
  float xp;

  

  /* reorder input */
  
  for (k=0; k<32; k++) {

    
    xp = a[k];
    a[k] = a[63-k];
    a[63-k] = xp;
  }

  
  fct3_64(a);

  
  for (k=1; k<64; k+=2) {

     
    a[k] = -a[k];
  }

  
}

#define INV_SQRT2    7.071067811865475e-1f
#define COS_PI_DIV8  9.238795325112867e-1f
#define COS_3PI_DIV8 3.826834323650898e-1f
#define SQRT2PLUS1   2.414213562373095f
#define SQRT2MINUS1  4.142135623730952e-1f

static void fft16(float *vec)
{
  float temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17,
    temp18, temp19, temp110, temp111, temp112, temp113, temp114, temp115;
  float temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27,
    temp28, temp29, temp210, temp211, temp212, temp213, temp214, temp215;
  float vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7,
    vec8, vec9, vec10, vec11, vec12, vec13, vec14, vec15;

  

  
  vec0 = vec[0] + vec[16];
  vec1 = vec[1] + vec[17];
  vec2 = vec[2] + vec[18];
  vec3 = vec[3] + vec[19];
  vec4 = vec[4] + vec[20];
  vec5 = vec[5] + vec[21];
  vec6 = vec[6] + vec[22];
  vec7 = vec[7] + vec[23];
  vec8 = vec[8] + vec[24];
  vec9 = vec[9] + vec[25];
  vec10 = vec[10] + vec[26];
  vec11 = vec[11] + vec[27];
  vec12 = vec[12] + vec[28];
  vec13 = vec[13] + vec[29];
  vec14 = vec[14] + vec[30];
  vec15 = vec[15] + vec[31];

  
  temp10  = vec0 + vec8;
  temp12  = vec0 - vec8;
  temp11  = vec1 + vec9;
  temp13  = vec1 - vec9;
  temp14  = vec2 + vec10;
  temp16  = vec2 - vec10;
  temp15  = vec3 + vec11;
  temp17  = vec3 - vec11;
  temp18  = vec4 + vec12;
  temp110 = vec4 - vec12;
  temp19  = vec5 + vec13;
  temp111 = vec5 - vec13;
  temp112 = vec6 + vec14;
  temp114 = vec6 - vec14;
  temp113 = vec7 + vec15;
  temp115 = vec7 - vec15;


  
  temp20  =  temp10 + temp18;
  temp24  =  temp10 - temp18;
  temp21  =  temp11 + temp19;
  temp25  =  temp11 - temp19;
  temp28  =  temp12 - temp111;
  temp210 =  temp12 + temp111;
  temp29  =  temp13 + temp110;
  temp211 =  temp13 - temp110;
  temp22  =  temp14 + temp112;
  temp27  =  temp14 - temp112;
  temp23  =  temp15 + temp113;
  temp26  =  temp113- temp15;

  
  temp11  =  temp16 + temp114;
  temp12  =  temp16 - temp114;
  temp10  =  temp17 + temp115;
  temp13  =  temp17 - temp115;

   
  temp212 = (temp10 + temp12) *  INV_SQRT2;
  temp214 = (temp10 - temp12) *  INV_SQRT2;
  temp213 = (temp13 - temp11) *  INV_SQRT2;
  temp215 = (temp11 + temp13) * -INV_SQRT2;


  
  vec0 = vec[0] - vec[16];
  vec1 = vec[1] - vec[17];
  vec2 = vec[2] - vec[18];
  vec3 = vec[3] - vec[19];
  vec4 = vec[4] - vec[20];
  vec5 = vec[5] - vec[21];
  vec6 = vec[6] - vec[22];
  vec7 = vec[7] - vec[23];
  vec8 = vec[8] - vec[24];
  vec9 = vec[9] - vec[25];
  vec10 = vec[10] - vec[26];
  vec11 = vec[11] - vec[27];
  vec12 = vec[12] - vec[28];
  vec13 = vec[13] - vec[29];
  vec14 = vec[14] - vec[30];
  vec15 = vec[15] - vec[31];

   
  temp19  = (vec2 + vec14) * -COS_3PI_DIV8;
  temp110 = (vec2 - vec14) * COS_PI_DIV8;
  temp18  = (vec3 + vec15) * COS_3PI_DIV8;
  temp111 = (vec3 - vec15) *  COS_PI_DIV8;
  temp15  = (vec4 + vec12) * -INV_SQRT2;
  temp16  = (vec4 - vec12) * INV_SQRT2;
  temp14  = (vec5 + vec13) * INV_SQRT2;
  temp17  = (vec5 - vec13) *  INV_SQRT2;
  temp113 = (vec6 + vec10) * -COS_PI_DIV8;
  temp114 = (vec6 - vec10) * COS_3PI_DIV8;
  temp112 = (vec7 + vec11) * COS_PI_DIV8;
  temp115 = (vec7 - vec11) *  COS_3PI_DIV8;

   
  vec2 = temp18  * SQRT2PLUS1  - temp112 * SQRT2MINUS1;
  vec3 = temp19  * SQRT2PLUS1  - temp113 * SQRT2MINUS1;
  vec4 = temp110 * SQRT2MINUS1 - temp114 * SQRT2PLUS1;
  vec5 = temp111 * SQRT2MINUS1 - temp115 * SQRT2PLUS1;


  
  temp18 += temp112;
  temp19 += temp113;
  temp110+= temp114;
  temp111+= temp115;

  
  vec6   = vec0  + temp14;
  vec10  = vec0  - temp14;
  vec7   = vec1  + temp15;
  vec11  = vec1  - temp15;

  
  vec12  = temp16 - vec9;
  vec14  = temp16 + vec9;
  vec13  = vec8  + temp17;
  vec15  = vec8  - temp17;

  
  temp10  = vec6  - vec14;
  temp12  = vec6  + vec14;
  temp11  = vec7  + vec15;
  temp13  = vec7  - vec15;
  temp14  = vec10 + vec12;
  temp16  = vec10 - vec12;
  temp15  = vec11 + vec13;
  temp17  = vec11 - vec13;

  
  vec10  = temp18 + temp110;
  temp110 = temp18 - temp110;
  vec11  = temp19 + temp111;
  temp111 = temp19 - temp111;

  
  temp112 = vec2  + vec4;
  temp114 = vec2  - vec4;
  temp113 = vec3  + vec5;
  temp115 = vec3  - vec5;


   
  *vec++ = temp20 + temp22;
  *vec++ = temp21 + temp23;
  *vec++ = temp12 + vec10;
  *vec++ = temp13 + vec11;
  *vec++ = temp210+ temp212;
  *vec++ = temp211+ temp213;
  *vec++ = temp10 + temp112;
  *vec++ = temp11 + temp113;
  *vec++ = temp24 - temp26;
  *vec++ = temp25 - temp27;
  *vec++ = temp16 + temp114;
  *vec++ = temp17 + temp115;
  *vec++ = temp28 + temp214;
  *vec++ = temp29 + temp215;
  *vec++ = temp14 + temp110;
  *vec++ = temp15 + temp111;
  *vec++ = temp20 - temp22;
  *vec++ = temp21 - temp23;
  *vec++ = temp12 - vec10;
  *vec++ = temp13 - vec11;
  *vec++ = temp210- temp212;
  *vec++ = temp211- temp213;
  *vec++ = temp10 - temp112;
  *vec++ = temp11 - temp113;
  *vec++ = temp24 + temp26;
  *vec++ = temp25 + temp27;
  *vec++ = temp16 - temp114;
  *vec++ = temp17 - temp115;
  *vec++ = temp28 - temp214;
  *vec++ = temp29 - temp215;
  *vec++ = temp14 - temp110;
  *vec++ = temp15 - temp111;

  
}

/*******************************************************************************
 Functionname:  cosMod
 *******************************************************************************

 Description: Performs cosine modulation.
 Return:      none

*******************************************************************************/
static void
cosMod (float *subband, HANDLE_SBR_QMF_FILTER_BANK qmfBank)
{
  int i;
  float wim, wre;
  float re1, im1, re2, im2;


  



    /* subband[2 * i]
                   subband[32 - 2 * i]
                   qmfBank->sin_twiddle[i]
                   qmfBank->cos_twiddle[i]
                   qmfBank->sin_twiddle[15 - i]
                   qmfBank->cos_twiddle[15 - i]
                */
  
  for (i = 0; i < 8; i++) {

    
    re1 = subband[2 * i];
    im2 = subband[2 * i + 1];
    re2 = subband[30 - 2 * i];
    im1 = subband[31 - 2 * i];

    
    wim = qmfBank->sin_twiddle[i];
    wre = qmfBank->cos_twiddle[i];

     
    subband[2 * i]     = im1 * wim + re1 * wre;
     
    subband[2 * i + 1] = im1 * wre - re1 * wim;

    
    wim = qmfBank->sin_twiddle[15 - i];
    wre = qmfBank->cos_twiddle[15 - i];

     
    subband[30 - 2 * i] = im2 * wim + re2 * wre;
     
    subband[31 - 2 * i] = im2 * wre - re2 * wim;
  }

  
  fft16(subband);

    /* subband[2 * i]
                   subband[32 - 2 * i]
                   qmfBank->alt_sin_twiddle[i]
                   qmfBank->sin_twiddle[15 - i]
                */

  
  wim = qmfBank->alt_sin_twiddle[0];
  wre = qmfBank->alt_sin_twiddle[16];

  
  for (i = 0; i < 8; i++) {

    
    re1 = subband[2 * i];
    im1 = subband[2 * i + 1];
    re2 = subband[30 - 2 * i];
    im2 = subband[31 - 2 * i];

     
    subband[2 * i]      = re1 * wre + im1 * wim;
     
    subband[31 - 2 * i] = re1 * wim - im1 * wre;

    
    wim = qmfBank->alt_sin_twiddle[i + 1];
    wre = qmfBank->alt_sin_twiddle[15 - i];

     
    subband[30 - 2 * i] = re2 * wim + im2 * wre;
     
    subband[2 * i + 1]  = re2 * wre - im2 * wim;
  }

  
}

/*******************************************************************************
 Functionname:  sinMod
 *******************************************************************************

 Description: Performs sine modulation.
 Return:      none

*******************************************************************************/
static void
sinMod (float *subband, HANDLE_SBR_QMF_FILTER_BANK qmfBank)
{
  int i;
  float wre, wim;
  float re1, im1, re2, im2;



  



    /* subband[2 * i]
                   subband[32 - 2 * i]
                   qmfBank->sin_twiddle[i]
                   qmfBank->cos_twiddle[i]
                   qmfBank->sin_twiddle[15 - i]
                   qmfBank->cos_twiddle[15 - i]
                */
  
  for (i = 0; i < 8; i++) {

    
    re1 = subband[2 * i];
    im2 = subband[2 * i + 1];
    re2 = subband[30 - 2 * i];
    im1 = subband[31 - 2 * i];

    
    wre = qmfBank->sin_twiddle[i];
    wim = qmfBank->cos_twiddle[i];

     
    subband[2 * i + 1] = im1 * wim + re1 * wre;
     
    subband[2 * i]     = im1 * wre - re1 * wim;

    
    wre = qmfBank->sin_twiddle[15 - i];
    wim = qmfBank->cos_twiddle[15 - i];

     
    subband[31 - 2 * i] = im2 * wim + re2 * wre;
     
    subband[30 - 2 * i] = im2 * wre - re2 * wim;
  }

  
  fft16(subband);

    /* subband[2 * i],
                   subband[32 - 2 * i],
                   qmfBank->alt_sin_twiddle[i],
                   qmfBank->sin_twiddle[15 - i]
                */

  
  wim = qmfBank->alt_sin_twiddle[0];
  wre = qmfBank->alt_sin_twiddle[16];

  
  for (i = 0; i < 8; i++) {

    
    re1 = subband[2 * i];
    im1 = subband[2 * i + 1];
    re2 = subband[30 - 2 * i];
    im2 = subband[31 - 2 * i];


     
    subband[31 - 2 * i] = -(re1 * wre + im1 * wim);
     
    subband[2 * i]      = -(re1 * wim - im1 * wre);

    
    wim = qmfBank->alt_sin_twiddle[i + 1];
    wre = qmfBank->alt_sin_twiddle[15 - i];

     
    subband[2 * i + 1]  = -(re2 * wim + im2 * wre);
     
    subband[30 - 2 * i] = -(re2 * wre - im2 * wim);

  }

  
}


/*******************************************************************************
 Functionname:  inverseModulation
 *******************************************************************************

 Description: Performs inverse modulation.
 Return:      none

*******************************************************************************/
static void
inverseModulation (const float *sbrReal,
                   const float *sbrImag,
                   float *timeOut,
                   HANDLE_SBR_QMF_FILTER_BANK qmfBank)
{
  int i;


  float gain = 0.015625f;
  float r1, i1, r2, i2;

  

   /* counting previous operations */

   /* timeOut[i]
                  timeOut[32 + i]
                  sbrReal[i]
                  sbrImag[i]
               */
  
  for (i = 0; i < 32; i++) {
     
    timeOut[i]      = gain * sbrReal[i];
    timeOut[32 + i] = gain * sbrImag[i];
  }

  
  cosMod (timeOut, qmfBank);
  
  sinMod (timeOut + 32, qmfBank);

    /* timeOut[i]
                   timeOut[63 - i]
                   timeOut[31 - i]
                   timeOut[32 + i]
                */
  
  for (i = 0; i < 16; i++) {

    
    r1 = timeOut[i];
    i2 = timeOut[63 - i];
    r2 = timeOut[31 - i];
    i1 = timeOut[32 + i];

     
    timeOut[i] = r1 - i1;

      
    timeOut[63 - i] = -(r1 + i1);

     
    timeOut[31 - i] = r2 - i2;

      
    timeOut[32 + i] = -(r2 + i2);

  }

  
}


/*!

  \brief  Performs forward modulation.
  \return none

*/

static void
forwardModulation (const float *timeIn,
                   float *rSubband,
                   float *iSubband,
                   HANDLE_SBR_QMF_FILTER_BANK qmfBank
                   )
{
  int i;


  



  
  rSubband[0] = timeIn[0];

   /* pointer for rSubband[],
                              iSubband[],
                              timeIn[128 - i]
                */
  
  for (i = 1; i < 64; i++) {

     
    rSubband[i]   = timeIn[i] - timeIn[128 - i];
    iSubband[i-1] = timeIn[i] + timeIn[128 - i];
  }

  
  iSubband[63] = timeIn[64];

  
  fct3_64(rSubband);

  
  fst3_64(iSubband);

  
}


/*!

  \brief      Performs complex-valued subband filtering.
  \return     none

*/
void
sbrAnalysisFiltering (const float *timeIn,
                      int   timeInStride,
                      float **rAnalysis,
                      float **iAnalysis,
                      HANDLE_SBR_QMF_FILTER_BANK qmfBank
                      )
{
  int i, p, k;



  float syn_buffer[2*QMF_CHANNELS];
  const float *ptr_pf;

  



   /* pointer for timeIn[(i * 64 + k)*timeInStride] */
  
  for (i = 0; i < 32; i++) {
    float accu;

     
    ptr_pf = qmfBank->p_filter;

        
    memmove(qmfBank->qmf_states_buffer, qmfBank->qmf_states_buffer + 64, 576 * sizeof(float));

     /* pointer for qmfBank->qmf_states_buffer[] */
    
    for(k = 0; k < 64; k++) {
       
      qmfBank->qmf_states_buffer[576 + k] = timeIn ? timeIn[(i * 64 + k)*timeInStride] : 0;
    }


     /* pointer for syn_buffer[] */
    
    for (k = 0; k < 64; k++){

      
      accu = 0.0;

       /* pointer for qmfBank->qmf_states_buffer[] */
      
      for (p = 0; p < 5; p++) {

        
        accu += *ptr_pf++ * qmfBank->qmf_states_buffer[p * 128 + k];
      }

      
      syn_buffer[127 - k] = accu;
    }

    
    accu = 0.0;

     /* pointer for qmfBank->qmf_states_buffer[] */
    
    for (p = 0; p < 5; p++) {

      
      accu += *ptr_pf++ * qmfBank->qmf_states_buffer[p * 128 + 127];
    }

    
    syn_buffer[0] = accu;

    
    ptr_pf -= 10;

     /* pointer for syn_buffer[] */
    
    for (k = 0; k < 63; k++){

      
      accu = 0.0;

       /* pointer for qmfBank->qmf_states_buffer[] */
      
      for (p = 0; p < 5; p++) {

        
        accu += *--ptr_pf * qmfBank->qmf_states_buffer[p * 128 + k + 64];
      }

      
      syn_buffer[63 - k] = accu;
    }

    
    forwardModulation (syn_buffer,
                       *(rAnalysis + i),
                       *(iAnalysis + i),
                       qmfBank);
  }

  
}


/*!

  \brief      Calculates energy form real and imaginary part
              of the QMF subsamples
  \return     none

*/
void
getEnergyFromCplxQmfData(float **energyValues,
                         float **realValues,
                         float **imagValues


                         )
{
  int j, k;

  

  
  for (k = 0; k < 16; k++) {

    
    for (j = 0; j < 64; j++) {

        
      energyValues[k][j] = (  (realValues[2*k][j] * realValues[2*k][j] +
                               imagValues[2*k][j] * imagValues[2*k][j]) +
                              (realValues[2*k+1][j] * realValues[2*k+1][j] +
                               imagValues[2*k+1][j] * imagValues[2*k+1][j])  ) * 0.5f;

    }
  }

  
}


/*!

  \brief      Creates QMF filter bank instance
  \return     1 in case of an error, otherwise 0

*/

int
createQmfBank (SBRRam_t *sbrram,
               int chan,
               HANDLE_SBR_QMF_FILTER_BANK   h_sbrQmf


               )
{


  

      
  memset (h_sbrQmf, 0, sizeof (SBR_QMF_FILTER_BANK));

   

  h_sbrQmf->p_filter = sbr_qmf_64_640;

     
  h_sbrQmf->qmf_states_buffer = sbrram->sbr_QmfStatesAnalysis + chan * 640;

  

  return (0);
}

void
deleteQmfBank (
               HANDLE_SBR_QMF_FILTER_BANK  h_sbrQmf /*!<pointer to QMF filter bank instance */
               )
{

  
  /*
    nothing to do
  */
  
}

/*******************************************************************************
 Functionname:  cplxSynthesisQmfFiltering
 *******************************************************************************

 Description: Performs inverse complex-valued subband filtering.
 Return:      void

*******************************************************************************/
void
SynthesisQmfFiltering (float **sbrReal,
                       float **sbrImag,
                       float *timeOut,
                       HANDLE_SBR_QMF_FILTER_BANK qmfBank)
{
  int k, j, p;



  const float *p_filter = qmfBank->p_filter + 1;

  

    /* counting previous operations */

   /* qmfBank->workBuffer

               */
  
  for (k = 0; k < 32; k++) {

     
    inverseModulation ( *(sbrReal + k), *(sbrImag + k), qmfBank->workBuffer, qmfBank);

     /* qmfBank->timeBuffer[p * 64 + j]
                    qmfBank->workBuffer[63 - j]
                    p_filter[2 * (p * 64 + j)]
                 */
    
    for (p = 0; p < 4; p++) {
      
      for (j = 0; j < 64; j++)
      {
         
        qmfBank->timeBuffer[p * 64 + j] +=
          p_filter[2 * (p * 64 + j)] * qmfBank->workBuffer[63 - j];
      }
    }

     /* qmfBank->timeBuffer[p * L2 + j]
                    qmfBank->workBuffer[L2 - 1 - j]
                    p_filter[SYNTH_STRIDE * (p * L2 + j)]
                 */
    
    for (j = 0; j < 32; j++)
    {
       
      qmfBank->timeBuffer[p * 64 + j] +=
        p_filter[2 * (p * 64 + j)] * qmfBank->workBuffer[63 - j];
    }

     /* timeOut[31 - j]
                    qmfBank->timeBuffer[p * 64 + 32 + j]
                    qmfBank->workBuffer[31 - j]
                    p_filter[2 * (p * 64 + 32 + j)]
                 */
    
    for (j = 0; j < 32; j++)
    {
        
      timeOut[31 - j] = qmfBank->timeBuffer[p * 64 + 32 + j] +
        p_filter[2 * (p * 64 + 32 + j)] * qmfBank->workBuffer[31 - j];
    }

    timeOut += 32;

        
    memmove (qmfBank->timeBuffer + 32, qmfBank->timeBuffer,
             288 * sizeof (float));

        
    memset (qmfBank->timeBuffer, 0, 32 * sizeof (float));

  }

  
}

/*******************************************************************************
 Functionname:  createCplxSynthesisQmfBank
 *******************************************************************************

 Description: Creates QMF filter bank instance
 Return:      errorCode.

*******************************************************************************/
int
CreateSynthesisQmfBank (SBRRam_t *sbrram,HANDLE_SBR_QMF_FILTER_BANK hs)
{
  float *ptr;

  

  
  ptr = &sbrram->PsBuf5[0];

      
  memset (hs, 0, sizeof (SBR_QMF_FILTER_BANK));

   

  hs->p_filter = p_64_640_qmf;

   
  hs->cos_twiddle     = sbr_cos_twiddle;
  hs->sin_twiddle     = sbr_sin_twiddle;
  hs->alt_sin_twiddle = sbr_alt_sin_twiddle;

   
  hs->workBuffer = ptr;ptr+=64;
  hs->timeBuffer = ptr;ptr+=320;

  

  return 0;

}

/*******************************************************************************
 Functionname:  DeleteSynthesisQmfBank
 *******************************************************************************

 Description: Destroys QMF filter bank instance
 Return:      none

 Written:     Per Ekstrand, CTS AB
*******************************************************************************/
void
DeleteSynthesisQmfBank (HANDLE_SBR_QMF_FILTER_BANK *h_sbrQmf)
{
  
  

  return;
}

