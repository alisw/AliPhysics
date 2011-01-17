#ifndef ALIESDCALOTRIGGER_H
#define ALIESDCALOTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 


Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <TNamed.h>

class TArrayI;

class AliESDCaloTrigger : public TNamed 
{
public:
	         AliESDCaloTrigger();
	         AliESDCaloTrigger(const AliESDCaloTrigger& ctrig);
	virtual ~AliESDCaloTrigger();
	
	AliESDCaloTrigger& operator=(const AliESDCaloTrigger& ctrig);
	
	Bool_t  IsEmpty() {return (fNEntries == 0);}

	virtual void Reset() {fCurrent = -1;}

	void    Allocate(Int_t size);
	void    DeAllocate(        ); 
	
	Bool_t  Add(Int_t col, Int_t row, Float_t amp, Float_t time, Int_t trgtimes[], Int_t ntrgtimes, Int_t trgts, Int_t trgbits);
	
	void    SetL1Threshold(Int_t i, Int_t thr) {fL1Threshold[i] = thr;}
//	void    SetTriggerBits(Int_t col, Int_t row, Int_t i, Int_t j);

	void    GetPosition(     Int_t& col, Int_t& row           ) const;
	
	void    GetAmplitude(  Float_t& amp                       ) const;
	void    GetTime(       Float_t& time                      ) const;
	
	void    GetTriggerBits(  Int_t& bits                      ) const;
	void    GetNL0Times(     Int_t& ntimes                    ) const;
	void    GetL0Times(      Int_t  times[]                   ) const;
	Int_t   GetEntries(                                       ) const {return fNEntries;}

	void    GetL1TimeSum(    Int_t& timesum                   ) const;

	Int_t   GetL1Threshold(  Int_t  i                         ) const {return fL1Threshold[i];}

	virtual Bool_t Next();

	virtual void Copy(TObject& obj) const;
	
	virtual void Print(const Option_t* opt) const;
	
private:

	Int_t    fNEntries;
    Int_t    fCurrent;

	Int_t*   fColumn;         // [fNEntries]
	Int_t*   fRow;            // [fNEntries]
	Float_t* fAmplitude;      // [fNEntries]
	Float_t* fTime;           // [fNEntries]
	Int_t*   fNL0Times;       // [fNEntries]
	TArrayI* fL0Times;        //
	Int_t*   fL1TimeSum;      // [fNEntries]
	Int_t*   fTriggerBits;    // [fNEntries]
	
	Int_t    fL1Threshold[2]; // L1 thresholds from raw data
	
	ClassDef(AliESDCaloTrigger, 3)
};
#endif

