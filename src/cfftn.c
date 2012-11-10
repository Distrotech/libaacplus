/*
  Complex FFT core for transforms
*/
#include "cfftn.h"
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _FFTW3
#include <fftw3.h>
#else
#include <math.h>
#include <string.h>


 /* the 3GPP instrumenting tool */

#define MAX_FACTORS  23
#define MAX_PERM    209
#define NFACTOR	     11
#ifndef M_PI
#define M_PI	3.14159265358979323846264338327950288
#endif
#define SIN60	0.86602540378443865
#define COS72	0.30901699437494742
#define SIN72	0.95105651629515357
#endif

#ifdef _FFTW3
void init_plans(FFTWFContext_t *ctx)
{
	fftwf_complex fft_data;

	ctx->plan4 = fftwf_plan_dft_1d(4, &fft_data, &fft_data, FFTW_BACKWARD, FFTW_ESTIMATE);
	ctx->plan8 = fftwf_plan_dft_1d(8, &fft_data, &fft_data, FFTW_BACKWARD, FFTW_ESTIMATE);
	ctx->plan64 = fftwf_plan_dft_1d(64, &fft_data, &fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
	ctx->plan512 = fftwf_plan_dft_1d(512, &fft_data, &fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
}

void destroy_plans(FFTWFContext_t *ctx)
{
	fftwf_destroy_plan(ctx->plan512);
	fftwf_destroy_plan(ctx->plan64);
	fftwf_destroy_plan(ctx->plan8);
	fftwf_destroy_plan(ctx->plan4);
}
#endif

#ifndef _FFTW3
int cfftn(float Re[],
          float Im[],
          int  nTotal,
          int  nPass,
          int  nSpan,
          int  iSign)
{

  int ii, mfactor, kspan, ispan, inc;
  int j, jc, jf, jj, k, k1, k2, k3=0, k4, kk, kt, nn, ns, nt;
  
  double radf;
  double c1, c2=0.0, c3=0.0, cd;
  double s1, s2=0.0, s3=0.0, sd;
  float  ak,bk,akp,bkp,ajp,bjp,ajm,bjm,akm,bkm,aj,bj,aa,bb;

  float Rtmp[MAX_FACTORS];
  float Itmp[MAX_FACTORS];
  double Cos[MAX_FACTORS];
  double Sin[MAX_FACTORS];
  
  int Perm[MAX_PERM];
  int factor [NFACTOR];
  
  double s60 = SIN60;
  double c72 = COS72;
  double s72 = SIN72;
  double pi2 = M_PI;
  
  

   /* pointer for Rtmp[],
                              Itmp[],
                              Cos[],
                              Sin[],
                              Perm[],
                              factor[] */

  
  Re--;
  Im--;
  
   
  if (nPass < 2) {
    
    return 0;
  }
  
  
  inc = iSign;
  
  if( iSign < 0 ) {
    
    s72 = -s72;
    s60 = -s60;
    pi2 = -pi2;
    inc = -inc;
  }
  
  
  nt = inc * nTotal;
  ns = inc * nSpan;
  
  kspan = ns;
  
  
  nn = nt - inc;

  
  jc = ns / nPass;

  
  radf = pi2 * (double) jc;
  pi2 *= 2.0;
  
  
  ii = 0;
  jf = 0;

  
  mfactor = 0;
  k = nPass;

   /* pointer for factor[] */
  
  while (k % 16 == 0) {
    
    mfactor++;

    
    factor [mfactor - 1] = 4;

    
    k /= 16;
  }

  
  j = 3;
  jj = 9;

   /* pointer for factor[] */
  
  do {
    
    while (k % jj == 0) {
      
      mfactor++;

      
      factor [mfactor - 1] = j;

      
      k /= jj;
    }

    
    j += 2;

    
    jj = j * j;
  } while (jj <= k);

   
  if (k <= 4) {

    
    kt = mfactor;
    factor [mfactor] = k;

     
    if (k != 1) {
      
      mfactor++;
    }
  } else {

      
    if (k - (k / 4 << 2) == 0) {
      
      mfactor++;
      
      factor [mfactor - 1] = 2;
      
      k /= 4;
    }

    
    kt = mfactor;
    j = 2;

    
    do {

      
      if (k % j == 0) {
        
        mfactor++;

        
        factor [mfactor - 1] = j;

        
        k /= j;
      }
       
      j = ((j + 1) / 2 << 1) + 1;
    } while (j <= k);
  }

  
  if (kt) {

    
    j = kt;

     /* pointer for factor[] */
    
    do {

      
      mfactor++;

      
      factor [mfactor - 1] = factor [j - 1];

      
      j--;
    } while (j);
  }
  
   
  if (mfactor > NFACTOR) {
    
    return(0);
  }
  
  
  for (;;) {
    
    sd = radf / (double) kspan;

    
    cd = sin(sd);

    
    cd = 2.0 * cd * cd;

    
    sd = sin(sd + sd);

    
    kk = 1;

    
    ii++;
    
    
    switch (factor [ii - 1]) {
    case 2:
      
      kspan /= 2;

      
      k1 = kspan + 2;

      
      do {

         /* pointer for Re[k2], Re[kk],
                                    Im[k2], Im[kk]  */
        
        do {
          
          k2 = kk + kspan;

          
          ak = Re [k2];
          bk = Im [k2];

           
          Re [k2] = Re [kk] - ak;
          Im [k2] = Im [kk] - bk;
          Re [kk] += ak;
          Im [kk] += bk;

          
          kk = k2 + kspan;

        } while (kk <= nn);

        
        kk -= nn;
      } while (kk <= jc);

       
      if (kk > kspan)
        goto Permute_Results_Label;


       /* pointers for Re[k2], Re[kk],
                                   Im[k2], Im[kk]  */
      
      do {

        
        c1 = 1.0 - cd;

        
        s1 = sd;

        
        do {

          
          do {

            
            do {
              
              k2 = kk + kspan;
              ak = Re [kk] - Re [k2];
              bk = Im [kk] - Im [k2];

               
              Re [kk] += Re [k2];
              Im [kk] += Im [k2];

                
              Re [k2] = (float)(c1 * ak - s1 * bk);

                
              Im [k2] = (float)(s1 * ak + c1 * bk);

              
              kk = k2 + kspan;
            } while (kk < nt);

            
            k2 = kk - nt;

            
            c1 = -c1;

            
            kk = k1 - k2;

          } while (kk > k2);

            
          ak = (float)(c1 - (cd * c1 + sd * s1));

           
          s1 = sd * c1 - cd * s1 + s1;

            
          c1 = 2.0 - (ak * ak + s1 * s1);

          
          s1 *= c1;
          c1 *= ak;

          
          kk += jc;
        } while (kk < k2);

        
        k1 += inc + inc;

         
        kk = (k1 - kspan) / 2 + jc;
      } while (kk <= jc + jc);
      break;
      
    case 4:
      
      ispan = kspan;

      
      kspan /= 4;
      

       /* pointer for Re[k2], Re[kk],
                                  Im[k2], Im[kk]  */
      
      do {

        
        c1 = 1.0;
        s1 = 0.0;

        
        do {
          
          do {
            
            k1 = kk + kspan;
            k2 = k1 + kspan;
            k3 = k2 + kspan;

            
            akp = Re [kk] + Re [k2];
            akm = Re [kk] - Re [k2];
            ajp = Re [k1] + Re [k3];
            ajm = Re [k1] - Re [k3];
            bkp = Im [kk] + Im [k2];
            bkm = Im [kk] - Im [k2];
            bjp = Im [k1] + Im [k3];
            bjm = Im [k1] - Im [k3];

             
            Re [kk] = akp + ajp;
            Im [kk] = bkp + bjp;

            
            ajp = akp - ajp;
            bjp = bkp - bjp;

            
            if (iSign < 0) {

              
              akp = akm + bjm;
              bkp = bkm - ajm;
              akm -= bjm;
              bkm += ajm;

            } else {

              
              akp = akm - bjm;
              bkp = bkm + ajm;
              akm += bjm;
              bkm -= ajm;
            }

            
            if (s1 != 0.0) {

                
              Re [k1] = (float)(akp * c1 - bkp * s1);
              Re [k2] = (float)(ajp * c2 - bjp * s2);
              Re [k3] = (float)(akm * c3 - bkm * s3);

                
              Im [k1] = (float)(akp * s1 + bkp * c1);
              Im [k2] = (float)(ajp * s2 + bjp * c2);
              Im [k3] = (float)(akm * s3 + bkm * c3);
            } else {

              
              Re [k1] = akp;
              Re [k2] = ajp;
              Re [k3] = akm;
              Im [k1] = bkp;
              Im [k2] = bjp;
              Im [k3] = bkm;
            }
            
            kk = k3 + kspan;
          } while (kk <= nt);
          
            
          c2 = c1 - (cd * c1 + sd * s1);

           
          s1 = sd * c1 - cd * s1 + s1;

            
          c1 = 2.0 - (c2 * c2 + s1 * s1);

          
          s1 *= c1;
          c1 *= c2;

           
          c2 = c1 * c1 - s1 * s1;

          
          s2 = 2.0 * c1 * s1;

           
          c3 = c2 * c1 - s2 * s1;

           
          s3 = c2 * s1 + s2 * c1;

          
          kk = kk - nt + jc;
        } while (kk <= kspan);

        
        kk = kk - kspan + inc;

      } while (kk <= jc);

       
      if (kspan == jc)
        goto Permute_Results_Label;
      break;
      
    default:
      
      k = factor [ii - 1];
      
      ispan = kspan;
      
      kspan /= k;
      
      
      switch (k) {
      case 3:


         /* pointer for Re[k1], Re[k2],
                                    Im[k1], Im[k2]  */
        
        do {
          
          do {

            
            k1 = kk + kspan;
            k2 = k1 + kspan;

            
            ak = Re [kk];
            bk = Im [kk];

            
            aj = Re [k1] + Re [k2];
            bj = Im [k1] + Im [k2];

             
            Re [kk] = ak + aj;
            Im [kk] = bk + bj;

             
            ak -= 0.5f * aj;
            bk -= 0.5f * bj;

             
            aj = (float)((Re [k1] - Re [k2]) * s60);
            bj = (float)((Im [k1] - Im [k2]) * s60);

             
            Re [k1] = ak - bj;
            Re [k2] = ak + bj;
            Im [k1] = bk + aj;
            Im [k2] = bk - aj;

            
            kk = k2 + kspan;
          } while (kk < nn);

          
          kk -= nn;
        } while (kk <= kspan);
        break;
        
      case 5:

         
        c2 = c72 * c72 - s72 * s72;

        
        s2 = 2.0 * c72 * s72;


         /* pointer for Re[k1], Re[k2], Re[k3], Re[k4], Re[k4],
                                     Im[k1], Im[k2], Im[k3], Im[k4], Im[kk] */
        
        do {

          
          do {

            
            k1 = kk + kspan;
            k2 = k1 + kspan;
            k3 = k2 + kspan;
            k4 = k3 + kspan;

            
            akp = Re [k1] + Re [k4];
            akm = Re [k1] - Re [k4];
            bkp = Im [k1] + Im [k4];
            bkm = Im [k1] - Im [k4];
            ajp = Re [k2] + Re [k3];
            ajm = Re [k2] - Re [k3];
            bjp = Im [k2] + Im [k3];
            bjm = Im [k2] - Im [k3];

            
            aa = Re [kk];
            bb = Im [kk];

             
            Re [kk] = aa + akp + ajp;
            Im [kk] = bb + bkp + bjp;

              
            ak = (float)(akp * c72 + ajp * c2 + aa);
            bk = (float)(bkp * c72 + bjp * c2 + bb);

             
            aj = (float)(akm * s72 + ajm * s2);
            bj = (float)(bkm * s72 + bjm * s2);

             
            Re [k1] = ak - bj;
            Re [k4] = ak + bj;
            Im [k1] = bk + aj;
            Im [k4] = bk - aj;

              
            ak = (float)(akp * c2 + ajp * c72 + aa);
            bk = (float)(bkp * c2 + bjp * c72 + bb);

             
            aj = (float)(akm * s2 - ajm * s72);
            bj = (float)(bkm * s2 - bjm * s72);

             
            Re [k2] = ak - bj;
            Re [k3] = ak + bj;
            Im [k2] = bk + aj;
            Im [k3] = bk - aj;

            
            kk = k4 + kspan;
          } while (kk < nn);

          
          kk -= nn;

        } while (kk <= kspan);
        break;
        
      default:

         
        if (k != jf) {

          
          jf = k;

          
          s1 = pi2 / (double) k;

          
          c1 = cos(s1);
          s1 = sin(s1);

           
          if (jf > MAX_FACTORS ) {
            
            return(0);
          }
          
          
          Cos [jf - 1] = 1.0;
          Sin [jf - 1] = 0.0;

          
          j = 1;

           /* pointers for Cos[j], Cos[k],
                                       Sin[j], Sin[k] */
          
          do {
              
            Cos [j - 1] = Cos [k - 1] * c1 + Sin [k - 1] * s1;

              
            Sin [j - 1] = Cos [k - 1] * s1 - Sin [k - 1] * c1;

            
            k--;

            
            Cos [k - 1] = Cos [j - 1];

             
            Sin [k - 1] = -Sin [j - 1];

            
            j++;
          } while (j < k);
        }

         /* pointers for Re[k1], Re[k2], Re[kk], Rtmp[], Cos[]
                                      Im[k1], Im[k2], Im[kk], Itmp[], Sin[]  */
        
        do {
          
          do {

            
            k1 = kk;

            
            k2 = kk + ispan;

            
            ak = aa = Re [kk];
            bk = bb = Im [kk];

            
            j = 1;

            
            k1 += kspan;

            
            do {

              
              k2 -= kspan;

              
              j++;

               
              Rtmp [j - 1] = Re [k1] + Re [k2];

              
              ak += Rtmp [j - 1];

               
              Itmp [j - 1] = Im [k1] + Im [k2];

              
              bk += Itmp [j - 1];

              
              j++;

               
              Rtmp [j - 1] = Re [k1] - Re [k2];
              Itmp [j - 1] = Im [k1] - Im [k2];

              
              k1 += kspan;

            } while (k1 < k2);

            
            Re [kk] = ak;
            Im [kk] = bk;

            
            k1 = kk;

            
            k2 = kk + ispan;

            
            j = 1;

            
            do {

              
              k1 += kspan;
              k2 -= kspan;

              
              jj = j;
              ak = aa;
              bk = bb;
              aj = 0.0;
              bj = 0.0;
              k = 1;

              
              do {

                
                k++;

                
                ak += (float)(Rtmp [k - 1] * Cos [jj - 1]);
                bk += (float)(Itmp [k - 1] * Cos [jj - 1]);

                
                k++;

                
                aj += (float)(Rtmp [k - 1] * Sin [jj - 1]);
                bj += (float)(Itmp [k - 1] * Sin [jj - 1]);

                
                jj += j;

                 
                if (jj > jf) {
                  
                  jj -= jf;
                }
              } while (k < jf);

              
              k = jf - j;

               
              Re [k1] = ak - bj;
              Im [k1] = bk + aj;
              Re [k2] = ak + bj;
              Im [k2] = bk - aj;

              
              j++;
            } while (j < k);

            
            kk += ispan;

          } while (kk <= nn);

          
          kk -= nn;

        } while (kk <= kspan);
        break;
      }

       
      if (ii == mfactor)
        goto Permute_Results_Label;

      
      kk = jc + 1;

       /* pointers for Re[kk],
                                   Im[kk]  */
      
      do {

        
        c2 = 1.0 - cd;

        
        s1 = sd;

        
        do {

          
          c1 = c2;
          s2 = s1;

          
          kk += kspan;

          
          do {

            
            do {

              
              ak = Re [kk];

                
              Re [kk] = (float)(c2 * ak - s2 * Im [kk]);

                
              Im [kk] = (float)(s2 * ak + c2 * Im [kk]);

              
              kk += ispan;

            } while (kk <= nt);

            
            ak = (float)(s1 * s2);

             
            s2 = s1 * c2 + c1 * s2;

             
            c2 = c1 * c2 - ak;

            
            kk = kk - nt + kspan;

          } while (kk <= ispan);

            
          c2 = c1 - (cd * c1 + sd * s1);

           
          s1 += sd * c1 - cd * s1;

            
          c1 = 2.0 - (c2 * c2 + s1 * s1);

          
          s1 *= c1;
          c2 *= c1;

          
          kk = kk - ispan + jc;

        } while (kk <= kspan);

        
        kk = kk - kspan + jc + inc;

      } while (kk <= jc + jc);
      break;
    }
  }
  
Permute_Results_Label:

  
  Perm [0] = ns;

  
  if (kt) {

    
    k = kt + kt + 1;

     
    if (mfactor < k) {
      
      k--;
    }

    
    j = 1;

    
    Perm [k] = jc;

     /* pointers for Perm[j], Perm[k]  */
    
    do {
       
      Perm [j] = Perm [j - 1] / factor [j - 1];

       
      Perm [k - 1] = Perm [k] * factor [j - 1];

      
      j++;
	     k--;
    } while (j < k);

    
    k3 = Perm [k];
    kspan = Perm [1];

    
    kk = jc + 1;
    k2 = kspan + 1;

    
    j = 1;

     
    if (nPass != nTotal) {
Permute_Multi_Label:

     /* pointers for Re[kk], Re[k2],
                                 Im[kk], Im[k2],
                                 Perm[]          */
    
    do {

      
      do {

        
        k = kk + jc;

        
        do {
            
		        ak = Re [kk]; Re [kk] = Re [k2]; Re [k2] = ak;
            bk = Im [kk]; Im [kk] = Im [k2]; Im [k2] = bk;

            
            kk += inc;
            k2 += inc;

        } while (kk < k);

        
        kk += ns - jc;
        k2 += ns - jc;

      } while (kk < nt);

      
      k2 = k2 - nt + kspan;
      kk = kk - nt + jc;

    } while (k2 < ns);

    
    do {

      
      do {

        
        k2 -= Perm [j - 1];

        
        j++;

        
        k2 = Perm [j] + k2;


      } while (k2 > Perm [j - 1]);

      
      j = 1;

      
      do {

         
        if (kk < k2)
		        goto Permute_Multi_Label;

        
        kk += jc;
        k2 += kspan;

      } while (k2 < ns);
    } while (kk < ns);
    } else {
Permute_Single_Label:

     /* pointers for Re[kk], Re[k2],
                                 Im[kk], Im[k2],
                                 Perm[]          */
    
    do {
      
      ak = Re [kk]; Re [kk] = Re [k2]; Re [k2] = ak;
      bk = Im [kk]; Im [kk] = Im [k2]; Im [k2] = bk;

      
      kk += inc;
      k2 += kspan;

    } while (k2 < ns);

    
    do {

      
      do {

        
        k2 -= Perm [j - 1];

        
        j++;

        
        k2 = Perm [j] + k2;

      } while (k2 > Perm [j - 1]);

      
      j = 1;

      
      do {

         
        if (kk < k2)
		        goto Permute_Single_Label;

        
        kk += inc;
        k2 += kspan;

      } while (k2 < ns);
    } while (kk < ns);
    }

    
    jc = k3;
  }
  
    
  if ((kt << 1) + 1 >= mfactor) {
    
    return 1;
  }

  
  ispan = Perm [kt];

  
  j = mfactor - kt;

  
  factor [j] = 1;

   /* pointers for factor[] */
  
  do {

     
    factor [j - 1] *= factor [j];

    
    j--;
  } while (j != kt);

  
  kt++;

  
  nn = factor [kt - 1] - 1;

   
  if (nn > MAX_PERM) {
    
    return(0);
  }
  
  
  j = jj = 0;

   /* pointers for factor[kt], factor[k],
                               Perm[]                  */
  
  for (;;) {

    
    k = kt + 1;

    
    k2 = factor [kt - 1];
    kk = factor [k - 1];

    
    j++;

     
    if (j > nn)
      break;

    
    jj += kk;

    
    while (jj >= k2) {

      
      jj -= k2;

      
      k2 = kk;

      
      k++;

      
      kk = factor [k - 1];

      
      jj += kk;
    }

    
    Perm [j - 1] = jj;
  }

  
  j = 0;

   /* pointers for Perm[j], Perm[k] */
  
  for (;;) {

    
    do {

      
      j++;

      
      kk = Perm [j - 1];

    } while (kk < 0);

     
    if (kk != j) {

      
      do {

        
        k = kk;
        kk = Perm [k - 1];

         
        Perm [k - 1] = -kk;
      } while (kk != j);

      
      k3 = kk;

    } else {

       
      Perm [j - 1] = -j;

       
      if (j == nn)
        break;
    }
  }
  
   /* pointers for Re[k1], Re[k2], Rtmp[],
                               Im[k1], Im[k2], Itmp[],
                               Perm[]                  */
  
  for (;;){

    
    j = k3 + 1;

    
    nt -= ispan;

    
    ii = nt - inc + 1;

    
    if (nt < 0)
      break;

    
    do {

      
      do {

        
        j--;
      } while (Perm [j - 1] < 0);

      
      jj = jc;

      
      do {

        
        kspan = jj;

          
        if (jj > MAX_FACTORS * inc) {

          
          kspan = MAX_FACTORS * inc;
        }

        
        jj -= kspan;

        
        k = Perm [j - 1];

         
        kk = jc * k + ii + jj;

        
        k1 = kk + kspan;

        
        k2 = 0;

        
        do {

          
          k2++;

          
          Rtmp [k2 - 1] = Re [k1];
          Itmp [k2 - 1] = Im [k1];

          
          k1 -= inc;

        } while (k1 != kk);

        
        do {

          
          k1 = kk + kspan;

           
          k2 = k1 - jc * (k + Perm [k - 1]);

          
          k = -Perm [k - 1];

          
          do {

            
            Re [k1] = Re [k2];
            Im [k1] = Im [k2];

            
            k1 -= inc;
            k2 -= inc;

          } while (k1 != kk);

          
          kk = k2;

        } while (k != j);

        
        k1 = kk + kspan;

        
        k2 = 0;

        
        do {

          
          k2++;

          
          Re [k1] = Rtmp [k2 - 1];
          Im [k1] = Itmp [k2 - 1];

          
          k1 -= inc;

        } while (k1 != kk);
      } while (jj);
    } while (j != 1);
  }

  
  return 1;
}
#endif

/*
   computes complex fourier transform of length len
                  
   returns status
*/
int CFFTN(FFTWFContext_t *ctx, float *afftData,int len, int isign)
{
#ifdef _FFTW3
	switch(len) {
	case 4:
		fftwf_execute_dft(ctx->plan4, (fftwf_complex*)afftData, (fftwf_complex*)afftData);
		break;
	case 8:
		fftwf_execute_dft(ctx->plan8, (fftwf_complex*)afftData, (fftwf_complex*)afftData);
		break;
	case 64:
		fftwf_execute_dft(ctx->plan64, (fftwf_complex*)afftData, (fftwf_complex*)afftData);
		break;
	case 512:
		fftwf_execute_dft(ctx->plan512, (fftwf_complex*)afftData, (fftwf_complex*)afftData);
		break;
	default:
		printf("non standard len for FFT: %d\nWill now die", len);
		exit(1);
	}

	return 1;
#else
  return(cfftn(afftData,afftData+1,len,len,len,2*isign));
#endif
}
