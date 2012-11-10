/*
  Implementation of math functions
*/
   
#include <math.h>
#include <assert.h>
#include "aacplusenc.h"

 /* the 3GPP instrumenting tool */

static float logDualisTable[LOG_DUALIS_TABLE_SIZE];

/*
  Creates lookup tables for some arithmetic functions
*/
void FloatFR_Init(void)
{
  int i;

  

   
  logDualisTable[0] = -1.0f; /* actually, ld 0 is not defined */

   /* logDualisTable[] */
  
  for (i=1; i<LOG_DUALIS_TABLE_SIZE; i++) {
      /* xxx * (1 / 0.30103) */ 
    logDualisTable[i] = (float) ( log(i)/log(2.0f) );
  }
  
  
}


/*
  The table must have been created before by FloatFR_Init().
  The valid range for a is 1 to LOG_DUALIS_TABLE_SIZE.
  For a=0, the result will be -1 (should be -inf).

  returns ld(a)
*/
float FloatFR_logDualis(int a)  /* Index for logarithm table */
{
  assert( a>=0 && a<LOG_DUALIS_TABLE_SIZE );

  
  
  

  return logDualisTable[a];
}


/*
  The function FloatFR_Init() must have been called before use.
  The valid range for a and b is 1 to LOG_DUALIS_TABLE_SIZE.

  returns ld(a/b)
*/
float FloatFR_getNumOctaves(int a, /* lower band */
                            int b) /* upper band */
{
  
    
  

  return (FloatFR_logDualis(b) - FloatFR_logDualis(a));
}

