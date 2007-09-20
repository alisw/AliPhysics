// @(#) $Id$
// Original: AliHLTDigitData.h,v 1.5 2004/05/12 11:51:27 loizides 

#ifndef _ALIHLTTPCDIGITDATA_H_
#define _ALIHLTTPCDIGITDATA_H_

#include "AliHLTTPCRootTypes.h" 

/**
 * @struct AliHLTTPCDigitData
 * Raw data structure for one row of TPC data. A digit corresponds to one
 * measured signal in 2 coordinates: pad and time <br>
 * The hight of the signal is given by the charge. <br>
 * This structure is only used while runnig HLT analysis in the offline
 * framework. Thats why it comes with a MC track id by default.
 *
 * The exact meaning of the 3 track ID fields is currently not known to me.
 * (Matthias 18.09.2007) 
 * @ingroup alihlt-tpc-datastructs
 */
struct AliHLTTPCDigitData
{
  UShort_t fCharge;
  UChar_t  fPad;
  UShort_t fTime;
  Int_t    fTrackID[3];
};
typedef struct AliHLTTPCDigitData AliHLTTPCDigitData;

/**
 * @struct AliHLTTPCDigitRowData
 * A container for TPC raw data organized in rows.
 * This is the 3rd coordinate which is missing in @ref AliHLTTPCDigitData.
 * @ingroup alihlt-tpc-datastructs
 */
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

/**
 * @struct AliHLTTPCRandomDigitData
 * Raw data structure for TPC data. A digit corresponds to one
 * measured signal in r coordinates: row, pad and time <br>
 * The hight of the signal is given by the charge. <br>
 *
 * The structure is not used for data exchange between components,
 * it's here for legacy reasons.
 * @ingroup alihlt-tpc-datastructs
 */
struct AliHLTTPCRandomDigitData{
  UChar_t fRow;
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
};
typedef struct AliHLTTPCRandomDigitData AliHLTTPCRandomDigitData;

/**
 * @struct AliHLTTPCPackedRawData
 * Container structure for TPC data.
 * It contains an array of TPC data objects, organized by pad rows.
 * @ingroup alihlt-tpc-datastructs
 */
struct AliHLTTPCUnpackedRawData
{
#ifndef __SUNPRO_CC
  AliHLTTPCDigitRowData fDigits[];
#else
  AliHLTTPCDigitRowData fDigits[1];
#endif
};

#endif /* _ALIHLTTPCDIGITDATA_H_ */
