#ifndef ALIEMCALTRIGGERPATCH_H
#define ALIEMCALTRIGGERPATCH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#include "TVector2.h"

class TArrayI;

class AliEMCALTriggerPatch : public TObject {

public:
	          AliEMCALTriggerPatch();                                  // default ctor
	          AliEMCALTriggerPatch(const AliEMCALTriggerPatch& other); // copy ctor
              AliEMCALTriggerPatch(Int_t i, Int_t j, Int_t e = 0, Int_t t = 0);
	 virtual ~AliEMCALTriggerPatch();

	void      SetPosition(Int_t px, Int_t py)  {fPosition->Set(float(px), float(py));}
	void      SetPosition(const TVector2& pos) {*fPosition = pos;}
	void      SetSum(Int_t sum) {fSum = sum;}
	void      SetTime(Int_t time) {fTime = time;}
	void      SetPeak(Int_t x, Int_t y, Int_t sizeX, Int_t sizeY);

	void      Position(TVector2& pos       ) const {pos = *fPosition;}
	void      Position(Int_t& px, Int_t& py) const {px = (Int_t)fPosition->X(); py = (Int_t)fPosition->Y();}
	TVector2* Position(                    ) const {return fPosition;}
	Int_t     Sum()   const {return fSum;} // in ADC counts
	Int_t     Time()  const {return fTime;}
	Int_t     Peaks() const {return fPeaks;}
	
	void      Print(const Option_t*) const;
	
private:
	
	AliEMCALTriggerPatch& operator=(const AliEMCALTriggerPatch& other); // Not implemented
	
	TVector2*         fPosition; // Position
	Int_t             fSum;      // Amplitude
	Int_t             fTime;     // Time
	Int_t             fPeaks;    // Peaks (L0 only)
	
	ClassDef(AliEMCALTriggerPatch,1)
};

#endif
