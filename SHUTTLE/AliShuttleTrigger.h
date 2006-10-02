#ifndef ALI_SHUTTLE_TRIGGER_H
#define ALI_SHUTTLE_TRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is to deal with DAQ LogBook and DAQ end run notification.
// It executes AliShuttle for retrieval of conditions data.
// For more info see AliShuttleTrigger.cxx
//

#include <TObject.h>
#include <TSysEvtHandler.h>
#include <TMutex.h>
#include <TCondition.h>


class AliShuttle;
class AliShuttleConfig;

class AliShuttleTrigger;
class DATENotifier;

class TerminateSignalHandler: public TSignalHandler {
	
public:
	TerminateSignalHandler(): TSignalHandler((ESignals) 0,0), fTrigger(0) { }
	TerminateSignalHandler(AliShuttleTrigger* trigger, ESignals signal):
		TSignalHandler(signal, kFALSE), fTrigger(trigger) {}

	virtual ~TerminateSignalHandler() { }
	virtual Bool_t Notify();

private:

	TerminateSignalHandler(const TerminateSignalHandler& other); 	
	TerminateSignalHandler& operator= (const TerminateSignalHandler& other); 	

	AliShuttleTrigger* fTrigger;  // pointer to the current AliShuttleTrigger

	ClassDef(TerminateSignalHandler, 0)
};

class AliShuttleTrigger: public TObject {

public:
	AliShuttleTrigger(const AliShuttleConfig* config, UInt_t timeout = 5000, Int_t retries = 5);
	~AliShuttleTrigger();

	AliShuttle* GetShuttle() {return fShuttle;}

	Bool_t Collect(Int_t run);
	Bool_t CollectNew();
	Bool_t CollectAll();

	virtual Bool_t Notify();
	void Terminate();

	void Run();

private:
	AliShuttleTrigger(const AliShuttleTrigger& other);
	AliShuttleTrigger& operator= (const AliShuttleTrigger& other);

	const AliShuttleConfig* fConfig;

	AliShuttle* fShuttle; 		// Pointer to the actual Shuttle instance

	Bool_t fNotified;  		// Notified flag
	Bool_t fTerminate; 		// Terminate flag

	TMutex fMutex;  		// Mutex
	TCondition fCondition;  	// Condition 

	TerminateSignalHandler fQuitSignalHandler; 		// Quit signal 
	TerminateSignalHandler fInterruptSignalHandler;  	// Interrupt signal


	ClassDef(AliShuttleTrigger, 0)
};

#endif
