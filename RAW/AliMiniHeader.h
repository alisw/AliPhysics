#ifndef ALIMINIHEADER_H
#define ALIMINIHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


struct AliMiniHeader {
  AliMiniHeader() {
    fMagicWord[0] = (kMagicWord >> 0) & 0xFF;
    fMagicWord[1] = (kMagicWord >> 8) & 0xFF;
    fMagicWord[2] = (kMagicWord >> 16) & 0xFF;
    fVersion = 1;
  }

  enum {kMagicWord = 0x123456};

  UInt_t    fSize;              // size of the raw data in bytes
  UChar_t   fDetectorID;        // unique detector number
  UChar_t   fMagicWord[3];      // hexadecimal word 123456 (used to detect byte swapping)
  UChar_t   fVersion;           // mini header version
  UChar_t   fCompressionFlag;   // compressed (=1) or uncompressed (=0)
  UShort_t  fDDLID;             // unique DDL number
};

#endif
