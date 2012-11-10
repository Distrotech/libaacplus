/*
  Quantization submodule
*/
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */





/*****************************************************************************

    functionname:quantizeLines 
    description: quantizes spectrum lines  
                 quaSpectrum = mdctSpectrum^3/4*2^(-(3/16)*gain)    
    input: global gain, number of lines to process, spectral data         
    output: quantized spectrum

*****************************************************************************/
static void quantizeLines(const int gain,
                          const int noOfLines,
                          const float *mdctSpectrum,
                          signed short *quaSpectrum)
{
  float  quantizer;
  float k = (-0.0946f+ 0.5f)/8192.0f;
  int   line;

  

   /* counting previous operation */

   /* mdctSpectrum[line]
                  quaSpectrum[line]
               */
  
  for (line = 0; line < noOfLines; line++)
  {
    float tmp  = mdctSpectrum[line];

     /* counting previous operation */

    
    k = -0.0946f + 0.5f;

        
    quantizer=quantTableE[(gain>>4)+8] * quantTableQ[gain & 15];

    
    if (tmp < 0.0f) {

       
      tmp = (float)sqrt(-tmp);

       
      tmp *= (float)sqrt(tmp); /* x^(3/4) */
      
        
      quaSpectrum[line] = -(int)(k + quantizer * tmp);
    }
    else {
      
      tmp = (float)sqrt(tmp);

       
      tmp *= (float)sqrt(tmp); /* x^(3/4) */

        
      quaSpectrum[line] = (int)(k + quantizer * tmp);
    }
  }

  
}


/*****************************************************************************

    functionname: QuantizeSpectrum
    description: quantizes the entire spectrum
    returns:
    input: number of scalefactor bands to be quantized, ...
    output: quantized spectrum

*****************************************************************************/
void QuantizeSpectrum(int sfbCnt,
                      int maxSfbPerGroup,
                      int sfbPerGroup,
                      int *sfbOffset,
                      float *mdctSpectrum,
                      int globalGain,
                      short *scalefactors,
                      short *quantizedSpectrum)
{
  int sfbOffs,sfb;

  

  
  for(sfbOffs=0;sfbOffs<sfbCnt;sfbOffs+=sfbPerGroup)
  {
   /* scalefactors[]
                  sfbOffset[]
               */
  
  for (sfb = 0; sfb < maxSfbPerGroup; sfb++)
  {
    int scalefactor = scalefactors[sfbOffs+sfb];

     /* counting previous operation */

      
    quantizeLines(globalGain - scalefactor,
                  sfbOffset[sfbOffs+sfb+1] - sfbOffset[sfbOffs+sfb],
                  mdctSpectrum + sfbOffset[sfbOffs+sfb],
                  quantizedSpectrum + sfbOffset[sfbOffs+sfb]);
  }
  }

  
}


void calcExpSpec(float *expSpec, float *mdctSpectrum, int noOfLines)
{
   int i;
   float tmp;

   

    /* pointers for mdctSpectrum[],
                                expSpec[]
                */
   
   for (i=0; i<noOfLines; i++) {
      /* x^(3/4) */

      
      tmp = (float)fabs(mdctSpectrum[i]);

      
      tmp = (float)sqrt(tmp);

        
      expSpec[i] = (float)(tmp * sqrt(tmp)); 
   }

   
}

/*****************************************************************************

    functionname:calcSfbDist 
    description: quantizes and requantizes lines to calculate distortion
    input:  number of lines to be quantized, ...
    output: distortion

*****************************************************************************/
float calcSfbDist(const float *spec,
                  const float *expSpec,
                  short *quantSpec,
                   int  sfbWidth,
                   int  gain)
{
  int i;
  float dist = 0;
  float k =  -0.0946f + 0.5f;
  float quantizer = quantTableE[(gain>>4)+8]*quantTableQ[gain & 15]; 
  float invQuantizer = invQuantTableE[(gain>>4)+8]*invQuantTableQ[gain & 15];
  
  

      /* counting previous operations */

   /* pointer for quantSpec[],
                              expSpec[],
                              spec[]
               */
  
  for(i=0;i<sfbWidth;i++){
    float iqval;
    float diff;

     
    quantSpec[i] = (int)(k + quantizer * expSpec[i]);

     
    if (quantSpec[i] < MAX_POW4_3_TABLE) {
       
      iqval = pow4_3_tab[quantSpec[i]] * invQuantizer;
    }
    else {
        
      iqval = (float) ( (pow((float)abs(quantSpec[i]), 4.0f/3.0f)) * invQuantizer );
    }

      
    diff = (float)fabs(spec[i])-iqval;

    
    dist+= diff * diff;

    }


  

  return dist;
}
