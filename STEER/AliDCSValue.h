#ifndef ALI_DCS_VALUE_H
#define ALI_DCS_VALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class represents the main value structure
// which forms so called 'historical data' in any SCADA system.
//

#include "AliSimpleValue.h"

class AliDCSValue: public TObject {
public:

	AliDCSValue();
	AliDCSValue(const AliSimpleValue& value, UInt_t timeStamp);

	AliSimpleValue& GetSimpleValue() {return fValue;};
	const AliSimpleValue& GetSimpleValue() const {return fValue;};	
	void SetSimpleValue(const AliSimpleValue& value) {fValue = value;};

	UInt_t GetTimeStamp() const {return fTimeStamp;};
	void SetTimeStamp(UInt_t timeStamp) {fTimeStamp = timeStamp;};

	Int_t GetSize() const {return fValue.GetSize() + sizeof(UInt_t);};

	TString ToString() const;

private:

        AliSimpleValue fValue;

        UInt_t fTimeStamp;


	ClassDef(AliDCSValue, 1);
};

#endif
