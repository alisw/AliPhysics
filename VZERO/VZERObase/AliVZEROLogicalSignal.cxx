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

// 
// Class AliVZEROLogicalSignal
// ---------------------------
// Describes a logical signal in the electronics. 
// Use it to generate observation windows
// which are used by AliVZEROTriggerSimulator class
// 

#include "AliLog.h"
#include "AliVZEROLogicalSignal.h"

ClassImp(AliVZEROLogicalSignal)

//_____________________________________________________________________________
AliVZEROLogicalSignal::AliVZEROLogicalSignal() : TObject(), fStart(0.), fStop(0.)
{
	// Default constructor
}
//_____________________________________________________________________________
AliVZEROLogicalSignal::AliVZEROLogicalSignal(UShort_t profilClock, UInt_t delay, Bool_t run2) : TObject(), fStart(0.), fStop(0.)
{
	// Constructor using the profilClock and delay parameters comming from the FEE
	
	Bool_t word;
	Bool_t up=kFALSE;
	Bool_t down=kFALSE;

	if (!run2) {
	  for(int i=0 ; i<5 ; i++) {
	    Int_t shift = (i<4) ? (3-i) : 4;
	    word = (profilClock >> shift) & 0x1;
	    if(word&&!up) {
	      fStart = 5. * (i + 1);
	      up = kTRUE;
	    }
	    if(!word&&up&&!down) {
	      fStop = 5. * (i + 1);
	      down = kTRUE;
	    }		
	  }
	  if(!down) fStop = 30.;
	}
	else {
	  for(int i=0 ; i<5 ; i++) {
	    Int_t shift = (i<2) ? (1-i) : (6-i);
	    word = (profilClock >> shift) & 0x1;
	    if(word&&!up) {
	      fStart = 5. * (i + 3);
	      up = kTRUE;
	    }
	    if(!word&&up&&!down) {
	      fStop = 5. * (i + 3);
	      down = kTRUE;
	    }		
	  }
	  if(!down) fStop = 40.;
	}
	
	fStart += delay*1e-2; // Add 10 ps par register unit
	fStop  += delay*1e-2; 
}
//_____________________________________________________________________________
AliVZEROLogicalSignal::AliVZEROLogicalSignal(const AliVZEROLogicalSignal &signal) : 
	TObject(), fStart(signal.fStart), 
	fStop(signal.fStop)
{
	// Copy constructor
}

//_____________________________________________________________________________
AliVZEROLogicalSignal::~AliVZEROLogicalSignal(){
	// Destructor
}

//_____________________________________________________________________________
AliVZEROLogicalSignal& AliVZEROLogicalSignal::operator = 
(const AliVZEROLogicalSignal& signal)
{
	// Operator =
        if(&signal == this) return *this;
	fStart = signal.fStart;
	fStop  = signal.fStop;
	return *this;
}

//_____________________________________________________________________________
AliVZEROLogicalSignal AliVZEROLogicalSignal::operator|(const AliVZEROLogicalSignal& signal) const 
{
	// Perform the Logical OR of two signals: C = A or B
	if((fStart>signal.fStop) || (signal.fStart>fStop))
		AliError(Form("Both signal do not superpose in time.\n  Start(A) = %f Stop(A) = %f\n   Start(B) = %f Stop(B) = %f",fStart, fStop, signal.fStart,signal.fStop));
	
	AliVZEROLogicalSignal result;
	if(fStart<signal.fStart) result.fStart = fStart;
	else result.fStart = signal.fStart;
	
	if(fStop>signal.fStop) result.fStop = fStop;
	else result.fStop = signal.fStop;
		
	return result;
}
//_____________________________________________________________________________
AliVZEROLogicalSignal AliVZEROLogicalSignal::operator&(const AliVZEROLogicalSignal& signal) const
{
	// Perform the Logical AND of two signals: C = A and B
	if((fStart>signal.fStop) || (signal.fStart>fStop))
		AliError(Form("Both signal do not superpose in time.\n  Start(A) = %f Stop(A) = %f\n   Start(B) = %f Stop(B) = %f",fStart, fStop, signal.fStart,signal.fStop));
	
	AliVZEROLogicalSignal result;
	if(fStart>signal.fStart) result.fStart = fStart;
	else result.fStart = signal.fStart;
	
	if(fStop<signal.fStop) result.fStop = fStop;
	else result.fStop = signal.fStop;
	
	return result;
}

//_____________________________________________________________________________
Bool_t AliVZEROLogicalSignal::IsInCoincidence(Float_t time) const
{
	// Check if a signal arriving at the time "time" is in coincidence with the logical signal
	Bool_t result = kFALSE;
	if((time>fStart) && (time<fStop)) result = kTRUE;
	return result;
}

