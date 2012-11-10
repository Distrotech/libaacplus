/*
  calculation of band energies
*/
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

void CalcBandEnergy(const float *mdctSpectrum,
                    const int   *bandOffset,
                    const int    numBands,
                    float      *bandEnergy,
                    float       *bandEnergySum) {
  int i, j;

  

  
  j = 0;
  *bandEnergySum = 0.0f;
 
   /* pointers for bandEnergy[],
                               bandOffset[],
                               mdctSpectrum[]
                */
  
  for(i=0; i<numBands; i++) {

    
    bandEnergy[i] = 0.0f;

    
    while (j < bandOffset[i+1]) {

       
      bandEnergy[i] += mdctSpectrum[j] * mdctSpectrum[j];

      j++;
    }
     
    *bandEnergySum += bandEnergy[i];
  }

  
}



void CalcBandEnergyMS(const float *mdctSpectrumLeft,
                      const float *mdctSpectrumRight,
                      const int   *bandOffset,
                      const int    numBands,
                      float      *bandEnergyMid,
                      float       *bandEnergyMidSum,
                      float      *bandEnergySide,
                      float       *bandEnergySideSum) {

  int i, j;

  

  
  j = 0;
	*bandEnergyMidSum = 0.0f;
	*bandEnergySideSum = 0.0f;

   /* pointers for bandEnergyMid[],
                               bandEnergySide[],
                               bandOffset[],
                               mdctSpectrumLeft[],
                               mdctSpectrumRight[]
               */
  
  for(i=0; i<numBands; i++) {

    
    bandEnergyMid[i] = 0.0f;
    bandEnergySide[i] = 0.0f;

    
    while (j < bandOffset[i+1]) {
        float specm, specs;

         
        specm = 0.5f * (mdctSpectrumLeft[j] + mdctSpectrumRight[j]);
        specs = 0.5f * (mdctSpectrumLeft[j] - mdctSpectrumRight[j]);

         
        bandEnergyMid[i]  += specm * specm;
        bandEnergySide[i] += specs * specs;

        j++;
    }
     
    *bandEnergyMidSum += bandEnergyMid[i];
		*bandEnergySideSum += bandEnergySide[i];

  }

  
}

