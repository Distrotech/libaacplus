/*
  Sbr miscellaneous helper functions
*/
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */


/* Sorting routine */
void Shellsort_int (int *in, int n)
{

  int i, j, v;
  int inc = 1;

  

   /* counting previous operation */

  
  do
  {
     
    inc = 3 * inc + 1;
  }
  while (inc <= n);

  
  do {

    
    inc = inc / 3;

     /* pointers for in[] */
     
    for (i = inc + 1; i <= n; i++) {

      
      v = in[i-1];
      j = i;

       /* pointers for in[j-inc-1],
                                   in[j-1]
                   */
      
      while (in[j-inc-1] > v) {

        
        in[j-1] = in[j-inc-1];

        
        j -= inc;

         
        if (j <= inc)
          break;
      }

      
      in[j-1] = v;
    }
  } while (inc > 1);

  
}



/*******************************************************************************
 Functionname:  AddVecLeft
 *******************************************************************************

 Arguments:   int* dst, int* length_dst, int* src, int length_src
 Return:      none
*******************************************************************************/
void
AddVecLeft (int *dst, int *length_dst, int *src, int length_src)
{
  int i;

  

   /* src[i] */
   
  for (i = length_src - 1; i >= 0; i--)
  {
    
    AddLeft (dst, length_dst, src[i]);
  }

  
}


/*******************************************************************************
 Functionname:  AddLeft
 *******************************************************************************

 Arguments:   int* vector, int* length_vector, int value
 Return:      none
*******************************************************************************/
void
AddLeft (int *vector, int *length_vector, int value)
{
  int i;

  

   /* vector[] */
  
  for (i = *length_vector; i > 0; i--)
  {
    
    vector[i] = vector[i - 1];
  }

  
  vector[0] = value;

  
  (*length_vector)++;

  
}


/*******************************************************************************
 Functionname:  AddRight
 *******************************************************************************

 Arguments:   int* vector, int* length_vector, int value
 Return:      none
*******************************************************************************/
void
AddRight (int *vector, int *length_vector, int value)
{
  

   
  vector[*length_vector] = value;

  
  (*length_vector)++;

  
}



/*******************************************************************************
 Functionname:  AddVecRight
 *******************************************************************************

 Arguments:   int* dst, int* length_dst, int* src, int length_src)
 Return:      none
*******************************************************************************/
void
AddVecRight (int *dst, int *length_dst, int *src, int length_src)
{
  int i;

  

   /* src[] */
  
  for (i = 0; i < length_src; i++)
  {
    
    AddRight (dst, length_dst, src[i]);
  }

  
}

