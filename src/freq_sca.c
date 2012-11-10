/*
  frequency scale
*/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static int getStartFreq(int fs, int start_freq);
static int getStopFreq(int fs, int stop_freq);

static int  numberOfBands(int b_p_o, int start, int stop, float warp_factor);
static void CalcBands(int * diff, int start , int stop , int num_bands);
static int  modifyBands(int max_band, int * diff, int length);
static void cumSum(int start_value, int* diff, int length, unsigned char  *start_adress);



/*******************************************************************************
 Functionname:  getSbrStartFreqRAW
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/

int
getSbrStartFreqRAW (int startFreq, int QMFbands, int fs)
{
  int result;

  

    
  if ( startFreq < 0 || startFreq > 15)
  {
    
    return -1;
  }

  
  result = getStartFreq(fs, startFreq);

     
  result =   (result*fs/QMFbands+1)>>1;

  

  return (result);

} /* End getSbrStartFreqRAW */


/*******************************************************************************
 Functionname:  getSbrStopFreq
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
int getSbrStopFreqRAW  (int stopFreq, int QMFbands, int fs)
{
  int result;

  

    
  if ( stopFreq < 0 || stopFreq > 13)
  {
    
    return -1;
  }


  
  result = getStopFreq( fs, stopFreq);

     
  result =   (result*fs/QMFbands+1)>>1;

  

  return (result);
} /* End getSbrStopFreq */


/*******************************************************************************
 Functionname:  getStartFreq
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
static int
getStartFreq(int fs, int start_freq)
{
  int k0_min;

  

   
  switch(fs){
  case 16000: k0_min = 24;
    break;
  case 22050: k0_min = 17;
    break;
  case 24000: k0_min = 16;
    break;
  case 32000: k0_min = 16;
    break;
  case 44100: k0_min = 12;
    break;
  case 48000: k0_min = 11;
    break;
  default:
    k0_min=11; /* illegal fs */

  }

     /* counting post-operations */
  

  switch (fs) {

  case 16000:
    {
      int v_offset[]= {-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
      return (k0_min + v_offset[start_freq]);
    }
    break;

  case 22050:
    {
      int v_offset[]= {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13};
      return (k0_min + v_offset[start_freq]);
    }
    break;

  case 24000:
    {
      int v_offset[]= {-5, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16};
      return (k0_min + v_offset[start_freq]);
    }
    break;

  case 32000:
    {
      int v_offset[]= {-6, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16};
      return (k0_min + v_offset[start_freq]);
    }
    break;

  case 44100:
  case 48000:
    {
      int v_offset[]= {-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20};
      return (k0_min + v_offset[start_freq]);
    }
    break;

  default:
    {
      int v_offset[]= {0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20, 24, 28, 33};
      return (k0_min + v_offset[start_freq]);
    }

  }

}/* End getStartFreq */


/*******************************************************************************
 Functionname:  getStopFreq
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
static int
getStopFreq(int fs, int stop_freq)
{
  int result,i;
  int *v_stop_freq = NULL;
  int k1_min;
  int v_dstop[13];

  int v_stop_freq_32[14] =  {32,34,36,38,40,42,44,46,49,52,55,58,61,64};
  int v_stop_freq_44[14] =  {23,25,27,29,32,34,37,40,43,47,51,55,59,64};
  int v_stop_freq_48[14] =  {21,23,25,27,30,32,35,38,42,45,49,54,59,64};

  

   /* counting previous operations */

  
  switch(fs){
  case 32000:
    
    k1_min = 32;

    
    v_stop_freq =v_stop_freq_32;
    break;

  case 44100:
    
    k1_min = 23;

    
    v_stop_freq =v_stop_freq_44;
    break;

  case 48000:
    
    k1_min = 21;

    
    v_stop_freq =v_stop_freq_48;
    break;

  default:
    
    k1_min = 21; /* illegal fs  */
  }



   /* v_dstop[]
                  v_stop_freq[]
               */
  
  for(i = 0; i <= 12; i++) {
     
    v_dstop[i] = v_stop_freq[i+1] - v_stop_freq[i];
  }

  
  Shellsort_int(v_dstop, 13);

  
  result = k1_min;

   /* v_dstop[] */
  
  for(i = 0; i < stop_freq; i++) {
    
    result = result + v_dstop[i];
  }

  

  return result;

}/* End getStopFreq */


/*******************************************************************************
 Functionname:  FindStartAndStopBand
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
int
FindStartAndStopBand(const int samplingFreq,
                     const int noChannels,
                     const int startFreq,
                     const int stopFreq,
                     const SR_MODE sampleRateMode,
                     int *k0,
                     int *k2)
{

  


   
  *k0 = getStartFreq(samplingFreq, startFreq);


     
  if( ( sampleRateMode == 1 ) &&
      ( samplingFreq*noChannels  <
        2**k0 * samplingFreq) ) {
    

    return (1);

  }


   
  if ( stopFreq < 14 ) {
     
    *k2 = getStopFreq(samplingFreq, stopFreq);
  } else {
     
    if( stopFreq == 14 ) {
     
    *k2 = 2 * *k0;
  } else {
     
    *k2 = 3 * *k0;
  }
  }


   
  if (*k2 > noChannels) {
    
    *k2 = noChannels;
  }


    
  if (*k2 - *k0 > noChannels / 2 - 4) {

    
    return (1);
  }

   
  if (*k2 > noChannels - 2) {

    
    return (1);
  }


   
  if ((*k2 - *k0) > MAX_FREQ_COEFFS)
  {
    
    return (1);
  }

   
  if ((*k2 - *k0) < 0)
  {
    
    return (1);
  }

  

  return(0);
}

/*******************************************************************************
 Functionname:  UpdateFreqScale
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
int
UpdateFreqScale(unsigned char  *v_k_master, int *h_num_bands,
                const int k0, const int k2,
                const int freqScale,
                const int alterScale)

{

  int     b_p_o = 0;
  float   warp = 0;
  int     dk = 0;

  /* Internal variables */
  int     two_regions = 0;
  int     k1 = 0, i;
  int     num_bands0;
  int     num_bands1;
  int     diff_tot[MAX_OCTAVE + MAX_SECOND_REGION];
  int     *diff0 = diff_tot;
  int     *diff1 = diff_tot+MAX_OCTAVE;
  int     k2_achived;
  int     k2_diff;
  int     incr = 0;

  

    /* counting previous operations */


   
  if(freqScale==1)
  {
    
    b_p_o=12;
  }

   
  if(freqScale==2)
  {
    
    b_p_o=10;
  }

   
  if(freqScale==3)
  {
    
    b_p_o=8;
  }


  
  if(freqScale > 0)
    {
      
      if(alterScale==0)
      {
        
        warp=1.0;
      }
      else
      {
        
        warp=1.3f;
      }


        
      if(4*k2 >= 9*k0)
        {
          
          two_regions=1;

          
          k1=2*k0;

          
          num_bands0=numberOfBands(b_p_o, k0, k1, 1.0);

          
          num_bands1=numberOfBands(b_p_o, k1, k2, warp);

          
          CalcBands(diff0, k0, k1, num_bands0);

          
          Shellsort_int( diff0, num_bands0);

          
          if (diff0[0] == 0)
          {
            

            return (1);

          }

          
          cumSum(k0, diff0, num_bands0, v_k_master);

          
          CalcBands(diff1, k1, k2, num_bands1);

          
          Shellsort_int( diff1, num_bands1);

            
          if(diff0[num_bands0-1] > diff1[0])
            {
                
              if(modifyBands(diff0[num_bands0-1],diff1, num_bands1))
              {
                
                return(1);
              }

            }


            
          cumSum(k1, diff1, num_bands1, &v_k_master[num_bands0]);
           
          *h_num_bands=num_bands0+num_bands1;

        }
      else
        {
          
          two_regions=0;
          k1=k2;

          
          num_bands0=numberOfBands(b_p_o, k0, k1, 1.0);

          
          CalcBands(diff0, k0, k1, num_bands0);

          
          Shellsort_int( diff0, num_bands0);

          
          if (diff0[0] == 0)
          {
            
            return (1);

          }

          
          cumSum(k0, diff0, num_bands0, v_k_master);

          
          *h_num_bands=num_bands0;

        }
    }
  else
    {
      
      if (alterScale==0) {

        
        dk = 1;

          
        num_bands0 = 2 * ((k2 - k0)/2);
      } else {

        
        dk = 2;

          
        num_bands0 = 2 * (((k2 - k0)/dk +1)/2);
      }

       
      k2_achived = k0 + num_bands0*dk;

      
      k2_diff = k2 - k2_achived;

      
      
      for(i=0;i<num_bands0;i++)
      {
        
        diff_tot[i] = dk;
      }



      
      if (k2_diff < 0) {

          
          incr = 1;
          i = 0;
      }



      
      if (k2_diff > 0) {

          
          incr = -1;

          
          i = num_bands0-1;
      }


       /* diff_tot[] */
      
      while (k2_diff != 0) {

         
        diff_tot[i] = diff_tot[i] - incr;

        /* ADD(1): pointer increment */
        i = i + incr;

        
        k2_diff = k2_diff + incr;
      }

      
      cumSum(k0, diff_tot, num_bands0, v_k_master);

      
      *h_num_bands=num_bands0;

    }

   
  if (*h_num_bands < 1)
  {
    
    return(1);
  }

  

  return (0);
}/* End UpdateFreqScale */


static int
numberOfBands(int b_p_o, int start, int stop, float warp_factor)
{
  int result=0;

  

   /* counting previous operation */

     
  result = 2* (int) ( b_p_o * log( (float) (stop)/start) / (2.0*log(2.0)*warp_factor) +0.5);

  

  return(result);
}


static void
CalcBands(int * diff, int start , int stop , int num_bands)
{
  int i;
  int previous;
  int current;

  

  
  previous=start;

   /* diff[] */
  
  for(i=1; i<= num_bands; i++)
    {
         
      current=(int) ( (start * pow( (float)stop/start, (float)i/num_bands)) + 0.5f);

       
      diff[i-1]=current-previous;

      
      previous=current;
    }

  

}/* End CalcBands */


static void
cumSum(int start_value, int* diff, int length,  unsigned char *start_adress)
{
  int i;

  

  
  start_adress[0]=start_value;

   /* start_adress[]
                  diff[]
               */
  
  for(i=1;i<=length;i++)
  {
     
    start_adress[i]=start_adress[i-1]+diff[i-1];
  }

  

} /* End cumSum */


static int
modifyBands(int max_band_previous, int * diff, int length)
{
  int change=max_band_previous-diff[0];

  

   /* counting previous operation */


     
  if ( change > (diff[length-1] - diff[0]) / 2 )
  {
    
    change = (diff[length-1] - diff[0]) / 2;
  }

   
  diff[0] += change;

   
  diff[length-1] -= change;

  
  Shellsort_int(diff, length);

  

  return(0);
}/* End modifyBands */


/*******************************************************************************
 Functionname:  UpdateHiRes
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
int
UpdateHiRes(unsigned char *h_hires, int *num_hires,unsigned char * v_k_master,
            int num_master , int *xover_band, SR_MODE drOrSr,
            int noQMFChannels)
{
  int i;
  int divider;
  int max1,max2;

  


    
  divider = (drOrSr == DUAL_RATE) ? 2 : 1;

      
  if( (v_k_master[*xover_band] > (noQMFChannels/divider) ) ||
      ( *xover_band > num_master ) )  {



    
    max1=0;
    max2=num_master;

     /* v_k_master[] */
     
    while( (v_k_master[max1+1] < (noQMFChannels/divider)) &&
           ( (max1+1) < max2) )
      {
          /* while()-condition */

        
        max1++;
      }

    
    *xover_band=max1;
  }

   
  *num_hires = num_master - *xover_band;

   /* h_hires[]
                  v_k_master[]
               */
  
  for(i = *xover_band; i <= num_master; i++)
    {
      
      h_hires[i - *xover_band] = v_k_master[i];
    }

  

  return (0);
}/* End UpdateHiRes */


/*******************************************************************************
 Functionname:  UpdateLoRes
 *******************************************************************************
 Description:

 Arguments:

 Return:
 *******************************************************************************/
void
UpdateLoRes(unsigned char * h_lores, int *num_lores, unsigned char * h_hires, int num_hires)
{
  int i;

  

   
  if(num_hires%2 == 0)
    {
       
      *num_lores=num_hires/2;


       /* h_lores[]
                      h_hires[]
                   */
      
      for(i=0;i<=*num_lores;i++)
      {
        
        h_lores[i]=h_hires[i*2];
      }

    }
  else
    {
        
      *num_lores=(num_hires+1)/2;


      
      h_lores[0]=h_hires[0];

       /* h_lores[]
                      h_hires[]
                   */
      
      for(i=1;i<=*num_lores;i++)
        {
          
          h_lores[i]=h_hires[i*2-1];
        }
    }

  

}/* End UpdateLoRes */
