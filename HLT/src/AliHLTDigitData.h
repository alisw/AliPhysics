// @(#) $Id$

#ifndef _DIGITDATA_H_
#define _DIGITDATA_H_

#include "AliHLTRootTypes.h" 

struct AliHLTDigitData
{
#ifdef do_mc
  Int_t fTrackID[3];
#endif
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
#ifdef IA64
  UChar_t dummy1;
  UChar_t dummy2;
#endif
};
typedef struct AliHLTDigitData AliHLTDigitData;
typedef AliHLTDigitData AliL3DigitData;

struct AliHLTDigitRowData
{
  UInt_t fNDigit;
  UInt_t fRow;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTDigitData fDigitData[1];
#else
  AliHLTDigitData fDigitData[0];
#endif
};
typedef struct AliHLTDigitRowData AliHLTDigitRowData;
typedef AliHLTDigitRowData AliL3DigitRowData;

struct AliHLTRandomDigitData{
  UChar_t fRow;
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliHLTRandomDigitData AliHLTRandomDigitData;
typedef AliHLTRandomDigitData AliL3RandomDigitData;
#endif /* _DIGITDATA_H_ */
