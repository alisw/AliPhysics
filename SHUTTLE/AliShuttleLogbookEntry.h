#ifndef ALI_SHUTTLE_LOGBOOK_ENTRY_H
#define ALI_SHUTTLE_LOGBOOK_ENTRY_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is a container for the data queried from DAQ's logbook and logbook_shuttle tables.
//

#include <TObject.h>
#include <TString.h>
#include <TMap.h>

#include "AliShuttleInterface.h"

class AliShuttleLogbookEntry : public TObject {

public:
	enum Status {
		kUnprocessed = 0,
		kInactive,
		kFailed,  // final
		kDone // final
	};

	AliShuttleLogbookEntry();
	AliShuttleLogbookEntry(Int_t run, Status* status=0);
	~AliShuttleLogbookEntry();

	AliShuttleLogbookEntry& operator=(const AliShuttleLogbookEntry& c);
	AliShuttleLogbookEntry(const AliShuttleLogbookEntry& c);
	virtual void Copy(TObject& c) const;

	Int_t GetRun() const {return fRun;}
	UInt_t GetStartTime() const  {TString tmp(GetRunParameter("DAQ_time_start")); return tmp.Atoi();}
	UInt_t GetEndTime() const {TString tmp(GetRunParameter("DAQ_time_end")); return tmp.Atoi();}
	UInt_t GetTimeCreated() const {TString tmp(GetRunParameter("time_created")); return tmp.Atoi();}
	Bool_t GetECSSuccess() const {TString tmp(GetRunParameter("ecs_success")); return (Bool_t) tmp.Atoi();}
	Bool_t GetDATestMode() const {return fDATestMode;}

//	void SetRun(Int_t run) {fRun=run;}

	void SetRunParameter(const char* key, const char* value);
	const char* GetRunParameter(const char* key) const;

	Status GetDetectorStatus(const char* detCode) const
			{return GetDetectorStatus(AliShuttleInterface::GetDetPos(detCode));}
	Status GetDetectorStatus(Int_t detPos) const;
	Status* GetDetectorStatus() const {return (Status*) fDetectorStatus;}

	void SetDetectorStatus(const char* detCode, Status status);
	void SetDetectorStatus(Status* status);
	void SetDetectorStatus(UInt_t detPos, Status status);
	void SetDetectorStatus(const char* detCode, const char* statusName);
	void SetDetectorStatus(UInt_t detPos, const char* statusName);
	void SetDATestMode(Bool_t daTestMode) {fDATestMode = daTestMode;}

	const char* GetRunType() const { return GetRunParameter("run_type"); }

	Bool_t IsDone() const;

	static const char* GetDetectorStatusName(Status status);
        void Print(Option_t *option) const;

private:

	Int_t fRun;   			// Run number
	TMap fRunParameters;		// run parameters written in DAQ logbook
	Status fDetectorStatus[AliShuttleInterface::kNDetectors]; 	// Detector status array
	Bool_t fDATestMode;             // flag to set whether we are in the DA Test Mode for the current run

	ClassDef(AliShuttleLogbookEntry, 0)
};

#endif

