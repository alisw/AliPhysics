/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBRunRange					   //
//  defines the run validity range of the object:		   //
//  [fFirstRun, fLastRun] 					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBRunRange.h"

#include "AliLog.h"

ClassImp(AliCDBRunRange)

//___________________________________________________________________________
AliCDBRunRange::AliCDBRunRange():
fFirstRun(-1), 
fLastRun(-1)
{
// constructor

}

//___________________________________________________________________________
AliCDBRunRange::AliCDBRunRange(Int_t firstRun, Int_t lastRun):
fFirstRun(firstRun), 
fLastRun(lastRun)
{
// constructor

}

//___________________________________________________________________________
AliCDBRunRange::~AliCDBRunRange() {
// destructor

}

//___________________________________________________________________________
Bool_t AliCDBRunRange::Overlaps(const AliCDBRunRange& other) const {
// check if this runRange overlaps other runRange
	
	if (!(IsValid() && other.IsValid())) {
		AliError("Comparing invalid run ranges!");
		return kFALSE;
	}	
	
	if (IsAnyRange() || other.IsAnyRange()) {
		AliError("Comparing unspecified ranges!");
		return kFALSE;
	}
	
	return ((fFirstRun < other.fFirstRun && other.fFirstRun <= fLastRun)
	    || (other.fFirstRun <= fFirstRun && fFirstRun <= other.fLastRun));
}

//___________________________________________________________________________
Bool_t AliCDBRunRange::Comprises(const AliCDBRunRange& other) const {
// check if this runRange contains other runRange

	if (!(IsValid() && other.IsValid())) {
		AliError("Comparing invalid run ranges!");
		return kFALSE;
	}	
	
	if (IsAnyRange()) {
		return kTRUE;
	}

	return fFirstRun <= other.fFirstRun && other.fFirstRun <= fLastRun
		&& fFirstRun <= other.fLastRun && other.fLastRun <= fLastRun;
}

//___________________________________________________________________________
Bool_t AliCDBRunRange::IsEqual(const TObject* obj) const {
// check if this runRange is equal to other runRange
	
        if (this == obj) {
                return kTRUE;
        }

        if (AliCDBRunRange::Class() != obj->IsA()) {
                return kFALSE;
        }
        AliCDBRunRange* other = (AliCDBRunRange*) obj;
	return fFirstRun == other->fFirstRun && fLastRun == other->fLastRun;
}

//___________________________________________________________________________
Bool_t AliCDBRunRange::IsValid() const {
// validity check

	if (fFirstRun < 0 && fLastRun < 0) {
		return kTRUE;
	}

	if (fFirstRun >= 0 && fLastRun >= fFirstRun) {
		return kTRUE;
	}

	return kFALSE;
}


