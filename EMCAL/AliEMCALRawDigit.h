#ifndef ALIEMCALRAWDIGIT_H
#define ALIEMCALRAWDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALRawDigit.h 17335 2007-03-10 03:40:17Z mvl $ */

// --- ROOT system ---

#include "TObject.h" 

// --- Standard library ---

// --- AliRoot header files ---
//#include "AliDigitNew.h"

class AliEMCALRawDigit : public TObject 
{
//	friend ostream& operator<<(ostream& , const AliEMCALRawDigit&);

public:
	
	AliEMCALRawDigit();
	AliEMCALRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples);
	AliEMCALRawDigit(const AliEMCALRawDigit& digit);
	
	virtual ~AliEMCALRawDigit();

	Bool_t operator==(const AliEMCALRawDigit &rValue) const;
	const AliEMCALRawDigit& operator = (const AliEMCALRawDigit&) {return *this;}

	Int_t   Compare(const TObject* obj) const;
	Bool_t  IsSortable() const {return kTRUE;}
	void    SetId(Int_t id) {fId = id;}
    Int_t   GetId() const {return fId;}
	
	Int_t   GetNSamples() const {return fNSamples;}
	Bool_t  GetTimeSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
	virtual void Print(const Option_t* opt) const;
	
private: 
 
	Int_t   fId;            //Absolute id
	Int_t   fNSamples;      //Number of time samples
	Int_t*  fSamples;	    //[fNSamples]
	
	ClassDef(AliEMCALRawDigit,1)   // Digit in EMCAL 
};
#endif

