// @(#) $Id$

#ifndef _DIGITDATA_H_
#define _DIGITDATA_H_

#include "AliL3RootTypes.h" 

struct AliL3DigitData
{
#ifdef do_mc
  Int_t fTrackID[3];
#endif
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliL3DigitData AliL3DigitData;

struct AliL3DigitRowData
{
  UInt_t fNDigit;
  UInt_t fRow;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliL3DigitData fDigitData[1];
#else
  AliL3DigitData fDigitData[0];
#endif
};
typedef struct AliL3DigitRowData AliL3DigitRowData;

struct AliL3RandomDigitData{
  UChar_t fRow;
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliL3RandomDigitData AliL3RandomDigitData;
#endif /* _DIGITDATA_H_ */
