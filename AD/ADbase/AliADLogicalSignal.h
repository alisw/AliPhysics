#ifndef ALIADLOGICALSIGNAL_H
#define ALIADLOGICALSIGNAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */


// 
// Class AliADLogicalSignal
// ---------------------------
// Describes a logical signal in the electronics. 
// Use it to generate observation windows
// which are used by AliADTriggerSimulator class
// 


#include <TObject.h>
#include <AliLog.h>

class AliADLogicalSignal  : public TObject {
public:
	AliADLogicalSignal();
	AliADLogicalSignal(UShort_t profilClock, UInt_t delay, UInt_t latch, UInt_t reset);
	virtual ~AliADLogicalSignal();
	AliADLogicalSignal(const AliADLogicalSignal &signal);
	AliADLogicalSignal& operator= (const AliADLogicalSignal &signal);
	AliADLogicalSignal operator& (const AliADLogicalSignal &signal) const;
	AliADLogicalSignal operator| (const AliADLogicalSignal &signal) const;
	// Print method
	virtual void Print(Option_t* opt="") const { AliInfo(Form("\t%s -> Start %f Stop %f\n ",opt,fStart,fStop));}
	
	Float_t GetStartTime() const {return fStart;};
	Float_t GetStopTime() const {return fStop;};
	Float_t GetWidth() const {return (fStop - fStart);};
	
	void SetStartTime(Float_t time){fStart = time;};
	void SetStopTime(Float_t time){fStop = time;};
	
	Bool_t IsInCoincidence(Float_t time) const;
	
private:
	
	Float_t fStart; // Start Time of the signal with respect to the LHC Clock
	Float_t fStop;  // Stop Time of the signal with respect to the LHC Clock
	
	
	ClassDef( AliADLogicalSignal, 1 )  
	
};

#endif // ALIADLOGICALSIGNAL_H


