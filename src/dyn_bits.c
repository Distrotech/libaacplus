/*
  Noiseless coder module
*/
#include <stdlib.h>
#include <limits.h>
#include "aacplusenc.h"
 /* the 3GPP instrumenting tool */


static int
calcSideInfoBits(int sfbCnt,
                 int blockType)
{
  int seg_len_bits = (blockType == SHORT_WINDOW ? 3 : 5);
  int escape_val = (blockType == SHORT_WINDOW ? 7 : 31);
  int sideInfoBits, tmp;

  

     /* counting previous operations */

  
  sideInfoBits = CODE_BOOK_BITS;
  tmp = sfbCnt;

  
  while (tmp >= 0)
  {
    
    sideInfoBits += seg_len_bits;
    tmp -= escape_val;
  }

  

  return (sideInfoBits);
}



/*
   count bits using all possible tables
 */

static void
buildBitLookUp(const short *quantSpectrum,
               const int maxSfb,
               const int *sfbOffset,
               const unsigned short *sfbMax,
               int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
               SECTION_INFO * section)
{
  int i;

  

   /* pointers for section[],
                               sfbOffset[],
                               sfbMax[],
                               bitLookUp[]
               */
  
  for (i = 0; i < maxSfb; i++)
  {
    int sfbWidth, maxVal;

    
    section[i].sfbCnt = 1;
    section[i].sfbStart = i;
    section[i].sectionBits = INVALID_BITCOUNT;
    section[i].codeBook = -1;

    
    sfbWidth = sfbOffset[i + 1] - sfbOffset[i];

    
    maxVal = sfbMax[i];

     
    bitCount(quantSpectrum + sfbOffset[i], sfbWidth, maxVal, bitLookUp[i]);
  }

  
}

/*
   essential helper functions
 */

static int
findBestBook(const int *bc, int *book)
{
  int minBits = INVALID_BITCOUNT, j;

  

   /* counting previous operations */

   /* pointer for bc[] */
  
  for (j = 0; j <= CODE_BOOK_ESC_NDX; j++)
  {
     
    if (bc[j] < minBits)
    {

      
      minBits = bc[j];
      *book = j;
    }
  }

  

  return (minBits);
}

static int
findMinMergeBits(const int *bc1, const int *bc2)
{
  int minBits = INVALID_BITCOUNT, j;

  

   /* counting previous operations */

   /* pointers for bc1[],
                               bc2[]
               */
  
  for (j = 0; j <= CODE_BOOK_ESC_NDX; j++)
  {
     
    if (bc1[j] + bc2[j] < minBits)
    {
      
      minBits = bc1[j] + bc2[j];
    }
  }

  

  return (minBits);
}

static void
mergeBitLookUp(int *bc1, const int *bc2)
{
  int j;

  

   /* pointers for bc1[],
                               bc2[]
               */
  
  for (j = 0; j <= CODE_BOOK_ESC_NDX; j++)
  {
      
    bc1[j] = min(bc1[j] + bc2[j], INVALID_BITCOUNT);
  }

  
}

static int
findMaxMerge(const int mergeGainLookUp[MAX_SFB_LONG],
             const SECTION_INFO *section,
             const int maxSfb, int *maxNdx)
{
  int i, maxMergeGain = 0;

  

   /* counting previous operations */

   /* pointers for section[],
                               mergeGainLookUp[]
               */
  
  for (i = 0; i + section[i].sfbCnt < maxSfb; i += section[i].sfbCnt)
  {
     
    if (mergeGainLookUp[i] > maxMergeGain)
    {
      
      maxMergeGain = mergeGainLookUp[i];
      *maxNdx = i;
    }
  }

  

  return (maxMergeGain);
}



static int
CalcMergeGain(const SECTION_INFO * section,
              int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
              const int *sideInfoTab,
              const int ndx1,
              const int ndx2)
{
  int SplitBits;
  int MergeBits;
  int MergeGain;

  

   /* pointers for section[ndx1],
                               section[ndx2],
                               bitLookUp[ndx1],
                               bitLookUp[ndx2]
               */

  /*
     Bit amount for splitted sections
   */
  
  SplitBits = section[ndx1].sectionBits + section[ndx2].sectionBits;

    
  MergeBits = sideInfoTab[section[ndx1].sfbCnt + section[ndx2].sfbCnt] + findMinMergeBits(bitLookUp[ndx1], bitLookUp[ndx2]);

  
  MergeGain = SplitBits - MergeBits;

  

  return (MergeGain);
}


static void
gmStage0(SECTION_INFO * section,
         int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
         const int maxSfb)
{
  int i;

  

   /* pointers for section[],
                               bitLookUp[]
               */
  
  for (i = 0; i < maxSfb; i++)
  {
     
    if (section[i].sectionBits == INVALID_BITCOUNT)
    {
        
      section[i].sectionBits = findBestBook(bitLookUp[i], &(section[i].codeBook));

    }
  }

  
}


static void
gmStage1(SECTION_INFO * section,
         int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
         const int maxSfb,
         const int *sideInfoTab)
{
  int mergeStart = 0, mergeEnd;

  

   /* counting previous operations */

  
  do
  {

     /* pointers for section[mergeStart]
                                 section[mergeEnd]
                                 bitLookUp[mergeStart]
                                 bitLookUp[mergeEnd]
                 */
     
    for (mergeEnd = mergeStart + 1; mergeEnd < maxSfb; mergeEnd++)
    {
       
      if (section[mergeStart].codeBook != section[mergeEnd].codeBook)
        break;

       
      section[mergeStart].sfbCnt++;

       
      section[mergeStart].sectionBits += section[mergeEnd].sectionBits;

      
      mergeBitLookUp(bitLookUp[mergeStart], bitLookUp[mergeEnd]);
    }

      
    section[mergeStart].sectionBits += sideInfoTab[section[mergeStart].sfbCnt];

    
    section[mergeEnd - 1].sfbStart = section[mergeStart].sfbStart;

    
    mergeStart = mergeEnd;

  } while (mergeStart < maxSfb);

  
}


static void
gmStage2(SECTION_INFO * section,
         int mergeGainLookUp[MAX_SFB_LONG],
         int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
         const int maxSfb,
         const int *sideInfoTab)
{
  int i;

  

   /* pointers for section[],
                               mergeGainLookUp[]
               */
  
  for (i = 0; i + section[i].sfbCnt < maxSfb; i += section[i].sfbCnt)
  {
      
    mergeGainLookUp[i] = CalcMergeGain(section,
                                       bitLookUp,
                                       sideInfoTab,
                                       i,
                                       i + section[i].sfbCnt);

  }

  
  while (TRUE)
  {
    int maxMergeGain, maxNdx, maxNdxNext, maxNdxLast;

     
    maxMergeGain = findMaxMerge(mergeGainLookUp, section, maxSfb, &maxNdx);

    /*
       exit while loop if no more gain is possible
     */
    
    if (maxMergeGain <= 0)
      break;


     /* pointers for section[maxNdx],
                                 bitLookUp[maxNdx],
                                 mergeGainLookUp[maxNdx]
                 */
    
    maxNdxNext = maxNdx + section[maxNdx].sfbCnt;

     /* pointers for section[maxNdxNext],
                                 bitLookUp[maxNdxNext]
                 */

    
    section[maxNdx].sfbCnt += section[maxNdxNext].sfbCnt;

    
    section[maxNdx].sectionBits += section[maxNdxNext].sectionBits - maxMergeGain;

    
    mergeBitLookUp(bitLookUp[maxNdx], bitLookUp[maxNdxNext]);

    
    if (maxNdx != 0)
    {
      
      maxNdxLast = section[maxNdx - 1].sfbStart;

        
      mergeGainLookUp[maxNdxLast] = CalcMergeGain(section,
                                                  bitLookUp,
                                                  sideInfoTab,
                                                  maxNdxLast,
                                                  maxNdx);
    }

    
    maxNdxNext = maxNdx + section[maxNdx].sfbCnt;

     /* pointers for section[maxNdxNext] */

    
    section[maxNdxNext - 1].sfbStart = section[maxNdx].sfbStart;

     
    if (maxNdxNext < maxSfb)
    {
       
      mergeGainLookUp[maxNdx] = CalcMergeGain(section,
                                              bitLookUp,
                                              sideInfoTab,
                                              maxNdx,
                                              maxNdxNext);
    }
  }

  
}

/*
   count bits used by the noiseless coder
 */
static void
noiselessCounter(AACRam_t *aacram, SECTION_DATA * sectionData,
                 int mergeGainLookUp[MAX_SFB_LONG],
                 int bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
                 const short *quantSpectrum,
                 const unsigned short *maxValueInSfb,
                 const int *sfbOffset,
                 const int blockType)
{
  int grpNdx, i;
  int *sideInfoTab = NULL;
  SECTION_INFO *section;

  

   /* counting previous operations */

  
  switch (blockType)
  {
  case LONG_WINDOW:
  case START_WINDOW:
  case STOP_WINDOW:

    
    sideInfoTab = aacram->sideInfoTabLong;
    break;
  case SHORT_WINDOW:

    
    sideInfoTab = aacram->sideInfoTabShort;
    break;
  }


   
  sectionData->noOfSections = 0;
  sectionData->huffmanBits = 0;
  sectionData->sideInfoBits = 0;



   
  if (sectionData->maxSfbPerGroup == 0) {

    
    return;
  }

   
  for (grpNdx = 0; grpNdx < sectionData->sfbCnt; grpNdx += sectionData->sfbPerGroup)
  {

     
    section = sectionData->section + sectionData->noOfSections;

      
    buildBitLookUp(quantSpectrum,
                   sectionData->maxSfbPerGroup,
                   sfbOffset + grpNdx,
                   maxValueInSfb + grpNdx,
                   bitLookUp,
                   section);

     
    gmStage0(section, bitLookUp, sectionData->maxSfbPerGroup);

    
    gmStage1(section, bitLookUp, sectionData->maxSfbPerGroup, sideInfoTab);

     
    gmStage2(section,
             mergeGainLookUp,
             bitLookUp,
             sectionData->maxSfbPerGroup,
             sideInfoTab);


     /* pointers for section[],
                                 bitLookUp[]
                 */
     
    for (i = 0; i < sectionData->maxSfbPerGroup; i += section[i].sfbCnt)
    {
       
      findBestBook(bitLookUp[i], &(section[i].codeBook));

       
      section[i].sfbStart += grpNdx;

        
      sectionData->huffmanBits += section[i].sectionBits - sideInfoTab[section[i].sfbCnt];

        
      sectionData->sideInfoBits += sideInfoTab[section[i].sfbCnt];

        
      sectionData->section[sectionData->noOfSections++] = section[i];
    }
  }
  
}


static void scfCount(const short *scalefacGain,
                     const unsigned short *maxValueInSfb,
                     SECTION_DATA * sectionData)

{
  /* counter */
  int i = 0; /* section counter */
  int j = 0; /* sfb counter */
  int k = 0; /* current section auxiliary counter */
  int m = 0; /* other section auxiliary counter */
  int n = 0; /* other sfb auxiliary counter */

  /* further variables */
  int lastValScf     = 0;
  int deltaScf       = 0;
  int found          = 0;
  int scfSkipCounter = 0;

  

   /* counting previous operations */

   
  sectionData->scalefacBits = 0;
  
  
  if (scalefacGain == NULL) {
    
    return;
  }

   
  lastValScf = 0;
  sectionData->firstScf = 0;
  
   /* sectionData->section[] */
   
  for (i = 0; i < sectionData->noOfSections; i++) {

     
    if (sectionData->section[i].codeBook != CODE_BOOK_ZERO_NO) {

       
      sectionData->firstScf = sectionData->section[i].sfbStart;
      
       
      lastValScf = scalefacGain[sectionData->firstScf];
      break;
    }
  }
  
   /* sectionData->section[] */
   
  for (i = 0; i < sectionData->noOfSections; i++) {

      
    if ((sectionData->section[i].codeBook != CODE_BOOK_ZERO_NO) &&
        (sectionData->section[i].codeBook != CODE_BOOK_PNS_NO)) {

       /* maxValueInSfb[]
                      scalefacGain[]
                   */
       
      for (j = sectionData->section[i].sfbStart;
           j < sectionData->section[i].sfbStart + sectionData->section[i].sfbCnt;
           j++) {
        
        
        if (maxValueInSfb[j] == 0) {

          
          found = 0;
          
          
          if (scfSkipCounter == 0) {
            
             
            if (j == (sectionData->section[i].sfbStart + sectionData->section[i].sfbCnt - 1) ) {
              
              
              found = 0;
            }
            else {
              
               /* maxValueInSfb[]
                              scalefacGain[]
                           */
              
              for (k = (j+1); k < sectionData->section[i].sfbStart + sectionData->section[i].sfbCnt; k++) {
                
                if (maxValueInSfb[k] != 0) {
                  
                  found = 1;
                  
                    
                  if ( (abs(scalefacGain[k] - lastValScf)) < CODE_BOOK_SCF_LAV) {
                    
                    deltaScf = 0;
                  }
                  else {
                     
                    deltaScf = -(scalefacGain[j] - lastValScf);
                    
                    
                    lastValScf = scalefacGain[j];
                    scfSkipCounter = 0;
                  }
                  break;
                }
                /* count scalefactor skip */
                
                scfSkipCounter = scfSkipCounter + 1;
              }
            }
            
            /* search for the next maxValueInSfb[] != 0 in all other sections */
             /* sectionData->section[] */
             
            for (m = (i+1); (m < sectionData->noOfSections) && (found == 0); m++) {

                
              if ((sectionData->section[m].codeBook != CODE_BOOK_ZERO_NO) &&
                  (sectionData->section[m].codeBook != CODE_BOOK_PNS_NO)) {
                 /* maxValueInSfb[]
                                scalefacGain[]
                             */
                
                for (n = sectionData->section[m].sfbStart;
                     n < sectionData->section[m].sfbStart + sectionData->section[m].sfbCnt;
                     n++) {
                  
                  if (maxValueInSfb[n] != 0) {
                    
                    found = 1;
                    
                      
                    if ( (abs(scalefacGain[n] - lastValScf)) < CODE_BOOK_SCF_LAV) {
                      
                      deltaScf = 0;
                    }
                    else {
                      
                       
                      deltaScf = -(scalefacGain[j] - lastValScf);
                      
                      
                      lastValScf = scalefacGain[j];
                      scfSkipCounter = 0;
                    }
                    break;
                  }
                  
                  
                  scfSkipCounter = scfSkipCounter + 1;
                }
              }
            }
            
            
            if (found == 0)   {
              
              deltaScf = 0;
              scfSkipCounter = 0;
            }
          }
          else {

            
            deltaScf = 0;
            
            
            scfSkipCounter = scfSkipCounter - 1;
          }
        }
        else {
           
          deltaScf = -(scalefacGain[j] - lastValScf);
          
          
          lastValScf = scalefacGain[j];
        }
        
           
        sectionData->scalefacBits += bitCountScalefactorDelta(deltaScf);
      }
    }
  }
  
  
}


typedef int (*lookUpTable)[CODE_BOOK_ESC_NDX + 1];

int
dynBitCount(AACRam_t *aacram,
            const short          *quantSpectrum,
            const unsigned short   *maxValueInSfb,
            const signed short     *scalefac,
            const int             blockType,
            const int             sfbCnt,
            const int             maxSfbPerGroup,
            const int             sfbPerGroup,
            const int            *sfbOffset,
            SECTION_DATA         *sectionData)
{
  int bitLookUp[MAX_SFB_LONG*(CODE_BOOK_ESC_NDX+1)];
  int mergeGainLookUp[MAX_SFB_LONG];

  

   
  sectionData->blockType      = blockType;
  sectionData->sfbCnt         = sfbCnt;
  sectionData->sfbPerGroup    = sfbPerGroup;

    
  sectionData->noOfGroups     = sfbCnt / sfbPerGroup;

  
  sectionData->maxSfbPerGroup = maxSfbPerGroup;

  
  noiselessCounter(aacram, sectionData,
                   mergeGainLookUp,
                   (lookUpTable)bitLookUp,
                   quantSpectrum,
                   maxValueInSfb,
                   sfbOffset,
                   blockType);

  
  scfCount(scalefac,
           maxValueInSfb,
           sectionData);

    /* counting post-operations */

  

  return (sectionData->huffmanBits +
          sectionData->sideInfoBits +
          sectionData->scalefacBits);
}


int
BCInit(AACRam_t *aacram)
{
  int i;

  

   /* sideInfoTabLong[] */
  
  for (i = 0; i <= MAX_SFB_LONG; i++)
  {
     
      aacram->sideInfoTabLong[i] = calcSideInfoBits(i, LONG_WINDOW);
  }

   /* sideInfoTabShort[] */
  for (i = 0; i <= MAX_SFB_SHORT; i++)
  {
     
      aacram->sideInfoTabShort[i] = calcSideInfoBits(i, SHORT_WINDOW);
  }

  

  return 0;
}
