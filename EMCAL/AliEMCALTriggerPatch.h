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
              AliEMCALTriggerPatch(Int_t i, Int_t j, Int_t e);
	 virtual ~AliEMCALTriggerPatch();

	void      Position(TVector2& pos) const {pos = *fPosition;}
	TVector2* Position(             ) const {return fPosition;}
	Int_t     Sum() const {return fSum;} // in ADC counts
	void      Print(const Option_t*) const;
	void      GetAbsCellIdsFromPatchPosition(TVector2& psize, TVector2& ssize, TArrayI& absid);
	
private:
	
	AliEMCALTriggerPatch& operator=(const AliEMCALTriggerPatch& other); // Not implemented
	
	TVector2*         fPosition;
	Int_t             fSum;
	
	ClassDef(AliEMCALTriggerPatch,1)
};

#endif
