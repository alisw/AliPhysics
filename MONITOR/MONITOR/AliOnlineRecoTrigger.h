#ifndef ALIONLINERECOTRIGGER_H
#define ALIONLINERECOTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TSysEvtHandler.h>
#include <TMutex.h>
#include <TCondition.h>


class AliOnlineRecoTrigger;
class SORNotifier;

class TerminateSignalHandler: public TSignalHandler {
	
public:
	TerminateSignalHandler(): TSignalHandler((ESignals) 0,0), fTrigger(0) { }
	TerminateSignalHandler(AliOnlineRecoTrigger* trigger, ESignals signal):
	TSignalHandler(signal, kFALSE), fTrigger(trigger) {}

	virtual ~TerminateSignalHandler() { }
	virtual Bool_t Notify();

private:
	TerminateSignalHandler(const TerminateSignalHandler& other);
	TerminateSignalHandler& operator= (const TerminateSignalHandler& other);
  
	AliOnlineRecoTrigger* fTrigger;  // pointer to the current AliOnlineRecoTrigger

	ClassDef(TerminateSignalHandler, 0)
};

class AliOnlineRecoTrigger: public TObject {

public:
	AliOnlineRecoTrigger();
	virtual ~AliOnlineRecoTrigger();

	virtual Bool_t Notify();
	void Terminate();

	Int_t Run();

private:
	AliOnlineRecoTrigger(const AliOnlineRecoTrigger& other);
	AliOnlineRecoTrigger& operator= (const AliOnlineRecoTrigger& other);

	Bool_t fNotified;  		// Notified flag
	Bool_t fTerminate; 		// Terminate flag

	TMutex fMutex;  		// Mutex
	TCondition fCondition;  	// Condition 

	TerminateSignalHandler* fQuitSignalHandler; 		// Quit signal
	TerminateSignalHandler* fInterruptSignalHandler;  	// Interrupt signal

	ClassDef(AliOnlineRecoTrigger, 0)
};

#endif
