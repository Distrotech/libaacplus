/*
  Bitstream encoder
*/
#ifndef _BITENC_H
#define _BITENC_H

#include "aacplusenc.h"

struct BITSTREAMENCODER_INIT
{
  int nChannels;
  int bitrate;
  int sampleRate;
  int profile;
};



int WriteBitstreamData (HANDLE_BIT_BUF hBitstream,
                    ELEMENT_INFO elInfo,
                    QC_OUT* qcOut,
                    PSY_OUT* psyOut,
                    int* globUsedBits,
                    const unsigned char *ancBytes
                    );

#endif /* _BITENC_H */
