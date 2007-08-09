#ifndef _ALIHLTTPCRAWDATAFORMAT_HPP_
#define _ALIHLTTPCRAWDATAFORMAT_HPP_

#include "Rtypes.h"
#include "AliHLTTPCDigitData.h"

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCRawData
 */

struct AliHLTTPCPackedRawData
    {
#ifndef __SUNPRO_CC
	UInt_t fPackedDigits[];
#else
	UInt_t fPackedDigits[1];
#endif
    };

struct AliHLTTPCUnpackedRawData
    {
#ifndef __SUNPRO_CC
	AliHLTTPCDigitRowData fDigits[];
#else
	AliHLTTPCDigitRowData fDigits[1];
#endif
    };

#endif // _ALIHLTTPCRAWDATAFORMAT_HPP_
