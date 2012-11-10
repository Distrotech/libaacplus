/*
  Transient detector
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */



static const float ABS_THRES=1.28e5;
static const float STD_FAC=1.0f;


/*******************************************************************************
 Functionname:  spectralChange
 *******************************************************************************
 \author  Klaus Peichl
 \brief   Calculates a measure for the spectral change within the frame

 \return  calculated value
*******************************************************************************/
static float spectralChange(float Energies[16][MAX_FREQ_COEFFS],
                            float EnergyTotal,
                            int nSfb,
                            int start,
                            int border,
                            int stop)
{
  int i,j;
  int len1 = border-start;
  int len2 = stop-border;
  float lenRatio = (float)len1 / (float)len2;
  float delta;
  float delta_sum= 0.0f;
  float nrg1[MAX_FREQ_COEFFS];
  float nrg2[MAX_FREQ_COEFFS];


  float pos_weight = (0.5f - (float)len1 / (float)(len1+len2));

  

     /* counting previous operations */

   
  pos_weight = 1.0f - 4.0f * pos_weight * pos_weight;


   /* pointer for nrg1[j],
                              nrg2[j]
               */
  
  for ( j=0; j<nSfb; j++ ) {


     
    nrg1[j] = 1.0e6f * len1;
    nrg2[j] = 1.0e6f * len2;

     /* pointer for Energies[i][j] */
    
    for (i=start; i<border; i++) {

       
      nrg1[j] += Energies[i][j];
    }

     /* pointer for Energies[i][j] */
    
    for (i=border; i<stop; i++) {

       
      nrg2[j] += Energies[i][j];
    }
  }

   /* pointer for nrg1[j],
                              nrg2[j]
               */
  
  for ( j=0; j<nSfb; j++ ) {

       
    delta = (float) fabs( log( nrg2[j] / nrg1[j] * lenRatio) );

       
    delta_sum += (float) (sqrt( (nrg1[j] + nrg2[j]) / EnergyTotal) * delta);
  }

   /* counting post-operations */

  

  return delta_sum * pos_weight;
}


/*******************************************************************************
 Functionname:  addLowbandEnergies
 *******************************************************************************
 \author  Klaus Peichl
 \brief   Calculates total lowband energy
 \return  total energy in the lowband
*******************************************************************************/
static float addLowbandEnergies(float **Energies,
                                unsigned char * freqBandTable,
                                int nSfb,
                                int slots)
{
  int i,k;
  float nrgTotal=1.0f;

  

   /* counting previous operations */

   
  for (k = 0; k < freqBandTable[0]; k++) {

     /* pointer for Energies[][] */
    
    for (i=0; i<slots; i++) {

        
        nrgTotal += Energies[(i+slots/2)/2][k];
    }
  }

  

  return(nrgTotal);
}


/*******************************************************************************
 Functionname:  addHighbandEnergies
 *******************************************************************************
 \author  Klaus Peichl
 \brief   Add highband energies
 \return  total energy in the highband
*******************************************************************************/
static float addHighbandEnergies(float **Energies, /*!< input */
                                 float EnergiesM[16][MAX_FREQ_COEFFS], /*!< Combined output */
                                 unsigned char * freqBandTable,
                                 int nSfb,
                                 int sbrSlots,
                                 int timeStep)
{
  int i,j,k,slotIn,slotOut;
  int li,ui;
  float nrgTotal=1.0;

  

   /* counting previous operations */

  
  for (slotOut=0; slotOut<sbrSlots; slotOut++) {

    
    slotIn = 2*slotOut;

     /* pointers for Energies[][],
                                 freqBandTable[]
                 */
    
    for ( j=0; j<nSfb; j++ ) {

      
      EnergiesM[slotOut][j] = 0;

      
      li = freqBandTable[j];
      ui = freqBandTable[j + 1];

      
      for (k = li; k < ui; k++) {

         /* pointers for Energies[][] */
        
        for (i=0; i<timeStep; i++) {

          
         EnergiesM[slotOut][j] += Energies[(slotIn+i)/2][k];
        }
      }
    }
  }

  
  for (slotOut=0; slotOut<sbrSlots; slotOut++) {

     /* pointers for Energies[][] */
    
    for (j=0; j<nSfb; j++) {

      
      nrgTotal += EnergiesM[slotOut][j];
    }
  }

  

  return(nrgTotal);
}


/*******************************************************************************
 Functionname:  frameSplitter
 *******************************************************************************
 \author  Klaus Peichl
 \brief   Decides if a frame shall be splitted into 2 envelopes
*******************************************************************************/
void
frameSplitter(float **Energies,
              HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector,
              unsigned char * freqBandTable,
              int nSfb,
              int timeStep,
              int no_cols,
              int *tran_vector)
{
  

  
  if (tran_vector[1]==0) { /* no transient was detected */
    float delta;
    float EnergiesM[16][MAX_FREQ_COEFFS];
    float EnergyTotal, newLowbandEnergy;
    int border;
    int sbrSlots = no_cols / timeStep;

    assert( sbrSlots * timeStep == no_cols );

     /* counting previous operations */

    
    newLowbandEnergy = addLowbandEnergies(Energies,
                                          freqBandTable,
                                          nSfb,
                                          no_cols);

      
    EnergyTotal = 0.5f * (newLowbandEnergy + h_sbrTransientDetector->prevLowBandEnergy);

      
    h_sbrTransientDetector->totalHighBandEnergy = addHighbandEnergies( Energies,
                                                  EnergiesM,
                                                  freqBandTable,
                                                  nSfb,
                                                  sbrSlots,
                                                  timeStep);

     
    EnergyTotal += h_sbrTransientDetector->totalHighBandEnergy;

       
    h_sbrTransientDetector->totalHighBandEnergy /= (sbrSlots*nSfb);

     
    border = (sbrSlots+1) / 2;

    
    delta = spectralChange( EnergiesM,
                            EnergyTotal,
                            nSfb,
                            0,
                            border,
                            sbrSlots);

       
    if (delta > h_sbrTransientDetector->split_thr)
      tran_vector[0] = 1; /* Set flag for splitting */
    else
      tran_vector[0] = 0;

     
    h_sbrTransientDetector->prevLowBandEnergy = newLowbandEnergy;
  }
  
}



static void
calculateThresholds (float **Energies,
                     int     noCols,
                     int     noRows,
                     float  *thresholds)
{
  int i, j;
  float mean_val, std_val, temp;
  float i_noCols  = (float)1.0f / (float)(noCols+noCols/2);
  float i_noCols1 = (float)1.0f / (float)(noCols+noCols/2 - 1);

  

     /* counting previous operations */

   /* pointer for Energies[],
                              thresholds[]
               */
  
  for (i = 0; i < noRows; i++)
  {

    
    mean_val = std_val = 0;

    
    for (j = noCols/2; j < 2*noCols; j++) {

      
      mean_val += (float)Energies[j/2][i];
    }

    
    mean_val = mean_val * i_noCols;       /* average */

    
    for (j = noCols/2; j < 2*noCols;j++) {

      
      temp = mean_val - (float)Energies[j/2][i];

      
      std_val += temp * temp;
    }

     
    std_val = (float) sqrt (std_val * i_noCols1);

       /* max() */  
    thresholds[i] =
      max (ABS_THRES, (float)0.66f*thresholds[i] + (float)0.34f*STD_FAC * std_val);

  }

  

}


static void
extractTransientCandidates (float **Energies,
                            int     noCols,
                            int     start_band,
                            int     stop_band,
                            float  *thresholds,
                            float  *transients,
                            int     bufferLength)

{
  int i, j;
  float delta_1, delta_2, delta_3, i_thres;
  int bufferMove = bufferLength / 2;

  

      
  memmove(transients, transients + noCols, bufferMove * sizeof (float));

      
  memset (transients + bufferMove, 0, (bufferLength-bufferMove) * sizeof (float));


   /* pointer for thresholds[]
               */
    
  for (i = start_band; i < stop_band; i++) {

    
    i_thres = (float)1.0f / thresholds[i];

     /* pointer for Energies[],
                                transients[]
                  */
    
    for (j = 0; j < noCols+noCols/2 - 3; j++) {

      
      delta_1 = Energies[(j + noCols/2 + 1)/2][i] - Energies[(j + noCols/2 - 1)/2][i];

       
      if(delta_1 > thresholds[i]) {
         
        transients[j + bufferMove] += delta_1 * i_thres - (float)1.0f;
      }

      
      delta_2 = Energies[(j + noCols/2 + 2)/2][i] - Energies[(j + noCols/2 - 2)/2][i] + delta_1;

       
      if(delta_2 > thresholds[i]) {

         
        transients[j + bufferMove] += delta_2 * i_thres - (float)1.0f;
      }

      
      delta_3 = Energies[(j + noCols/2 + 3)/2][i] - Energies[(j + noCols/2 - 3)/2][i] + delta_2;

       
      if(delta_3 > thresholds[i]) {

         
        transients[j + bufferMove] += delta_3 * i_thres - (float)1.0f;
      }


    }
  }

  
}


void
transientDetect (float **Energies,
                 HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTran,
                 int *tran_vector,
                 int timeStep
                 )
{

  int i;
  int cond;
  int no_cols = h_sbrTran->no_cols;
  int qmfStartSample = no_cols + timeStep * 4;
  float int_thres = h_sbrTran->tran_thr / (float) h_sbrTran->no_rows;

  

      /* counting previous operations */

   
  calculateThresholds (Energies,
                       h_sbrTran->no_cols,
                       h_sbrTran->no_rows,
                       h_sbrTran->thresholds);


   
  extractTransientCandidates (Energies,
                              h_sbrTran->no_cols,
                              0,
                              h_sbrTran->no_rows,
                              h_sbrTran->thresholds,
                              h_sbrTran->transients,
                              h_sbrTran->buffer_length);


  
  tran_vector[0] = 0;
  tran_vector[1] = 0;


  
  for (i = qmfStartSample; i < qmfStartSample + no_cols; i++){

      
      cond = (h_sbrTran->transients[i] < (float)0.9f * h_sbrTran->transients[i - 1]) &&
           (h_sbrTran->transients[i - 1] > int_thres);


      
      if (cond) {

           
	      tran_vector[0] =
	        (int) floor ( (i - qmfStartSample) / timeStep );

        
	      tran_vector[1] = 1;
        break;
      }
  }

  
}





int
CreateSbrTransientDetector (SBRRam_t *sbrram,
                            int chan,
                            HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector,

                            int   sampleFreq,
                            int   totalBitrate,
                            int   codecBitrate,
                            int   tran_thr,
                            int   mode,
                            int   tran_fc


                            )
{

  float bitrateFactor;
  float frame_dur = (float) 2048 / (float) sampleFreq;
  float temp = frame_dur - (float)0.010f;

  

    /* counting previous operations */

      
  memset(h_sbrTransientDetector,0,sizeof(SBR_TRANSIENT_DETECTOR));

  
  if(codecBitrate)
  {
    
    bitrateFactor = (float)totalBitrate / (float)codecBitrate;
  }
  else
  {
    
    bitrateFactor = (float)1.0f;
  }

   
  if (temp < (float)0.0001f)
  {
    
    temp = (float)0.0001f;
  }

   
  temp = (float)0.000075f / (temp * temp);





   
  h_sbrTransientDetector->no_cols = 32;
  h_sbrTransientDetector->tran_thr = (float) tran_thr;
  h_sbrTransientDetector->tran_fc = tran_fc;

    
  h_sbrTransientDetector->split_thr = temp * bitrateFactor;
  h_sbrTransientDetector->buffer_length = 96;

   
  h_sbrTransientDetector->no_rows = 64;
  h_sbrTransientDetector->mode = mode;
  h_sbrTransientDetector->prevLowBandEnergy = 0;

   
  h_sbrTransientDetector->thresholds = &sbrram->sbr_thresholds[chan*QMF_CHANNELS];

      
  memset(h_sbrTransientDetector->thresholds,0,sizeof(float)*QMF_CHANNELS);

    
  h_sbrTransientDetector->transients =  &sbrram->sbr_transients[chan*h_sbrTransientDetector->buffer_length];

      
  memset(h_sbrTransientDetector->transients,0,sizeof(float)*h_sbrTransientDetector->buffer_length);

  

  return 0;
}



void
deleteSbrTransientDetector (HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector)
{
  

  
  if(h_sbrTransientDetector){
    /*
      nothing to do
    */

  }

  

}
