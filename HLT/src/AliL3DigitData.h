#ifndef _DIGITDATA_H_
#define _DIGITDATA_H_

#include "AliL3RootTypes.h" 

struct AliL3DigitData
{
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliL3DigitData AliL3DigitData;

struct AliL3DigitRowData
{
  UInt_t fNDigit;
  UInt_t fRow;
  AliL3DigitData fDigitData[0];
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
