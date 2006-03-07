#ifndef ALI_SHUTTLE_TRIGGER_H
#define ALI_SHUTTLE_TRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is to deal with DAQ LogBook and DAQ end run notification.
// It executes AliShuttle for retrieval of conditions data.
//

#include <TObject.h>
#include <TSysEvtHandler.h>
#include <TMutex.h>
#include <TCondition.h>


class AliCDBStorage;
class AliShuttle;
class AliShuttleConfig;

class AliShuttleTrigger;
class DATENotifier;

class TerminateSignalHandler: public TSignalHandler {
	
	AliShuttleTrigger* fTrigger;
public:
	TerminateSignalHandler(AliShuttleTrigger* trigger, ESignals signal):
		TSignalHandler(signal, kFALSE), fTrigger(trigger) {}

	virtual Bool_t Notify();

	ClassDef(TerminateSignalHandler, 0)
};

class AliShuttleTrigger: public TObject {
public:
	AliShuttleTrigger(const AliShuttleConfig* config,
			AliCDBStorage* storage, UInt_t timeout = 5000,
			Int_t retries = 5);
	~AliShuttleTrigger();

	AliShuttle* GetShuttle() {return fShuttle;}

	Bool_t Collect(Int_t run);
	Bool_t CollectNew();
	Bool_t CollectAll();
	
	virtual Bool_t Notify();
	void Terminate();

	void Run();

private:

	class DATEEntry: public TObject {
		Int_t fRun;
		UInt_t fStartTime;
		UInt_t fEndTime;
	public:
		DATEEntry(Int_t run, UInt_t startTime, UInt_t endTime):
			fRun(run), fStartTime(startTime), fEndTime(endTime) {}

		Int_t GetRun() {return fRun;}
		UInt_t GetStartTime() {return fStartTime;}
		UInt_t GetEndTime() {return fEndTime;}

		ClassDef(DATEEntry, 0)
	};

	Bool_t RetrieveDATEEntries(const char* whereClause, TObjArray& entries,
			Int_t& lastRun);
	Bool_t RetrieveConditionsData(const TObjArray& dateEntries);

	const AliShuttleConfig* fConfig;
	AliCDBStorage* fStorage;

	AliShuttle* fShuttle;

	Bool_t fNotified;
	Bool_t fTerminate;

	TMutex fMutex;
	TCondition fCondition;

	TerminateSignalHandler fQuitSignalHandler;
	TerminateSignalHandler fInterruptSignalHandler;


	ClassDef(AliShuttleTrigger, 0)
};

#endif
