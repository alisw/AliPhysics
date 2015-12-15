#ifndef ALIEMCALTRIGGERRAWDIGIT_H
#define ALIEMCALTRIGGERRAWDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
 Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALRawDigit.h" 
#include "AliEMCALTriggerTypes.h" 
#include "AliLog.h"

class AliEMCALTriggerRawDigit : public AliEMCALRawDigit 
{
public:
	
	AliEMCALTriggerRawDigit();
	AliEMCALTriggerRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples);
	
	virtual ~AliEMCALTriggerRawDigit();
	
	void    SetTriggerBit(const int type, const Int_t mode) {fTriggerBits = (fTriggerBits | (1 << (type + kTriggerTypeEnd * mode)));}
	
	Bool_t  SetL0Time(   Int_t i);
	
	Int_t   GetTriggerBit(const TriggerType_t type, const Int_t mode) const;

	Int_t   GetTriggerBits() const {return fTriggerBits;}
	
	Bool_t  GetL0Time(const Int_t i, Int_t& time) const;
	Bool_t  GetL0Times(Int_t times[]            ) const;
	Int_t   GetNL0Times(                        ) const {return fNL0Times;}
	
	Int_t   GetL0TimeSum(const Int_t time) const;
	
	void    SetL1TimeSum(Int_t ts) {if (fL1TimeSum >= 0) AliWarning("You're overwriting digit time sum! Please check"); fL1TimeSum = ts;}
	Int_t   GetL1TimeSum(        ) const {return fL1TimeSum;}

	void    SetL1SubRegion(Int_t sr) {if (fL1SubRegion >= 0) AliWarning("You're overwriting digit subregion! Please check"); fL1SubRegion = sr;}
        Int_t   GetL1SubRegion(        ) const {return fL1SubRegion;}
	
	virtual void Print(const Option_t* opt) const;
	
private: 
 
	AliEMCALTriggerRawDigit(const AliEMCALTriggerRawDigit &cd);            // Not implemented
	AliEMCALTriggerRawDigit &operator=(const AliEMCALTriggerRawDigit &cd); // Not implemented

	Int_t   fTriggerBits; // Trigger bits
	Int_t   fNL0Times;    // N L0 times
	Int_t   fL0Times[10]; // L0 times
	
	Int_t   fL1TimeSum;   // L1 time sum
        Int_t   fL1SubRegion; // Subregion
	
	ClassDef(AliEMCALTriggerRawDigit,2)
};
#endif

