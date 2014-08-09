#ifndef ALICUTVALUERANGE_H
#define ALICUTVALUERANGE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TObject.h>

namespace EMCalTriggerPtAnalysis{

template<typename t>
class AliCutValueRange : public TObject {
public:
	AliCutValueRange();
	AliCutValueRange(t min, t max);
	AliCutValueRange(t limit, bool isUpper);
	~AliCutValueRange() {}

	void SetLimits(t min, t max){
		fLimits[0] = min;
		fLimits[1] = max;
		fHasLimit[0] = fHasLimit[1] = true;
	}
	void UnsetLimits(){ fHasLimit[0] = fHasLimit[1] = false; }
	void SetLimit(t value, bool isUpper){
		int bin = isUpper ? 1 : 0;
		fLimits[bin] = value;
		fHasLimit[bin] = true;
	}
	void UnsetLimit(bool isUpper){
		int bin = isUpper ? 1 : 0;
		fHasLimit[bin] = false;
	}
	void Negate() { fNegate = true; }
	void SetPositive() { fNegate = false; }
	bool IsInRange(t value) const;
private:
	t       fLimits[2];
	bool    fHasLimit[2];
	bool    fNegate;

	ClassDef(AliCutValueRange, 1);     // Value range for cuts
};

}

#endif
