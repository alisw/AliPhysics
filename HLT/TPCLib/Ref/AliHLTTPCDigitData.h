// @(#) $Id$

#ifndef _DIGITDATA_H_
#define _DIGITDATA_H_

#include "AliHLTTPCRootTypes.h" 

struct AliHLTTPCDigitData
{
#ifdef do_mc
  Int_t fTrackID[3];
#endif
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliHLTTPCDigitData AliHLTTPCDigitData;

struct AliHLTTPCDigitRowData
{
  UInt_t fNDigit;
  UInt_t fRow;
  AliHLTTPCDigitData fDigitData[0];
};
typedef struct AliHLTTPCDigitRowData AliHLTTPCDigitRowData;

struct AliHLTTPCRandomDigitData{
  UChar_t fRow;
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliHLTTPCRandomDigitData AliHLTTPCRandomDigitData;
#endif /* _DIGITDATA_H_ */
