#ifndef _ALIHLTTPCRAWDATAFORMAT_HPP_
#define _ALIHLTTPCRAWDATAFORMAT_HPP_

#include "Rtypes.h"
#include "AliL3DigitData.h"

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCRawData
 */

struct AliHLTTPCPackedRawData
    {
	UInt_t fPackedDigits[];
    };

struct AliHLTTPCUnpackedRawData
    {
	AliL3DigitRowData fDigits[];
    };

#endif // _ALIHLTTPCRAWDATAFORMAT_HPP_
