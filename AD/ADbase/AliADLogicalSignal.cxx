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
// Class AliADLogicalSignal
// ---------------------------
// Describes a logical signal in the electronics. 
// Use it to generate observation windows
// which are used by AliADTriggerSimulator class
// 
#include <iostream>
#include <bitset>

#include "AliLog.h"
#include "AliADLogicalSignal.h"

ClassImp(AliADLogicalSignal)

//_____________________________________________________________________________
AliADLogicalSignal::AliADLogicalSignal() : TObject(), fStart(0.), fStop(0.)
{
	// Default constructor
}
//_____________________________________________________________________________
AliADLogicalSignal::AliADLogicalSignal(UShort_t profilClock, UInt_t delay, UInt_t latch, UInt_t reset) : TObject(), fStart(0.), fStop(0.)
{
	/*/
	std::cout << "P " << std::bitset<5>(profilClock)<< std::endl;
	std::cout << "L " << std::bitset<5>(latch)<< std::endl;
	std::cout << "R " << std::bitset<5>(reset)<< std::endl;
	std::cout << "Delay " <<delay<< std::endl;
	/*/
	
	// Constructor using the profilClock and delay parameters comming from the FEE
	Bool_t fClock[11];
	Int_t fTimes[11];
	Bool_t risingFound = kFALSE;
	
	for(Int_t i=0; i<5; i++) {
		fClock[i+1] = (profilClock >> 4-i) & 0x1;
		fClock[i+6] = (profilClock >> 4-i) & 0x1;
		}
	fClock[0] = (profilClock >> 0) & 0x1;
	 
	if(reset>latch) for(Int_t i=0; i<11; i++) fTimes[i] = -5+5*i;
	if(reset<latch) for(Int_t i=0; i<11; i++) fTimes[i] = -30+5*i;
	
	/*/
	for(Int_t i=0; i<10; i++)std::cout<<fTimes[i]<<" ";
	std::cout<<std::endl;
	for(Int_t i=0; i<10; i++)std::cout<<fClock[i]<<" ";
	std::cout<<std::endl;
	/*/
	
	for(Int_t i=1; i<11; i++){
		if(!risingFound && !fClock[i-1] && fClock[i]){
			risingFound = kTRUE;
			fStart = fTimes[i];
			continue;
			}
		else if(risingFound && fClock[i-1] && !fClock[i]) {
			fStop = fTimes[i];
			break;
			}	
		}
	
	fStart += delay*1e-2; // Add 10 ps par register unit
	fStop  += delay*1e-2; 
}
//_____________________________________________________________________________
AliADLogicalSignal::AliADLogicalSignal(const AliADLogicalSignal &signal) : 
	TObject(), fStart(signal.fStart), 
	fStop(signal.fStop)
{
	// Copy constructor
}

//_____________________________________________________________________________
AliADLogicalSignal::~AliADLogicalSignal(){
	// Destructor
}

//_____________________________________________________________________________
AliADLogicalSignal& AliADLogicalSignal::operator = 
(const AliADLogicalSignal& signal)
{
	// Operator =
        if(&signal == this) return *this;
	fStart = signal.fStart;
	fStop  = signal.fStop;
	return *this;
}

//_____________________________________________________________________________
AliADLogicalSignal AliADLogicalSignal::operator|(const AliADLogicalSignal& signal) const 
{
	// Perform the Logical OR of two signals: C = A or B
	if((fStart>signal.fStop) || (signal.fStart>fStop))
		AliError(Form("Both signal do not superpose in time.\n  Start(A) = %f Stop(A) = %f\n   Start(B) = %f Stop(B) = %f",fStart, fStop, signal.fStart,signal.fStop));
	
	AliADLogicalSignal result;
	if(fStart<signal.fStart) result.fStart = fStart;
	else result.fStart = signal.fStart;
	
	if(fStop>signal.fStop) result.fStop = fStop;
	else result.fStop = signal.fStop;
		
	return result;
}
//_____________________________________________________________________________
AliADLogicalSignal AliADLogicalSignal::operator&(const AliADLogicalSignal& signal) const
{
	// Perform the Logical AND of two signals: C = A and B
	if((fStart>signal.fStop) || (signal.fStart>fStop))
		AliError(Form("Both signal do not superpose in time.\n  Start(A) = %f Stop(A) = %f\n   Start(B) = %f Stop(B) = %f",fStart, fStop, signal.fStart,signal.fStop));
	
	AliADLogicalSignal result;
	if(fStart>signal.fStart) result.fStart = fStart;
	else result.fStart = signal.fStart;
	
	if(fStop<signal.fStop) result.fStop = fStop;
	else result.fStop = signal.fStop;
	
	return result;
}

//_____________________________________________________________________________
Bool_t AliADLogicalSignal::IsInCoincidence(Float_t time) const
{
	// Check if a signal arriving at the time "time" is in coincidence with the logical signal
	Bool_t result = kFALSE;
	if((time>fStart) && (time<fStop)) result = kTRUE;
	return result;
}

