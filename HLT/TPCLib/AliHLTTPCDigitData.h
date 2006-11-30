// @(#) $Id$
// Original: AliHLTDigitData.h,v 1.5 2004/05/12 11:51:27 loizides 

#ifndef _ALIHLTTPCDIGITDATA_H_
#define _ALIHLTTPCDIGITDATA_H_

#include "AliHLTTPCRootTypes.h" 

struct AliHLTTPCDigitData
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
typedef struct AliHLTTPCDigitData AliHLTTPCDigitData;

struct AliHLTTPCDigitRowData
{
  UInt_t fNDigit;
  UInt_t fRow;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCDigitData fDigitData[1];
#else
  AliHLTTPCDigitData fDigitData[0];
#endif
};
typedef struct AliHLTTPCDigitRowData AliHLTTPCDigitRowData;

struct AliHLTTPCRandomDigitData{
  UChar_t fRow;
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliHLTTPCRandomDigitData AliHLTTPCRandomDigitData;
#endif /* _ALIHLTTPCDIGITDATA_H_ */
