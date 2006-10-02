#ifndef ALI_SHUTTLE_LOGBOOK_ENTRY_H
#define ALI_SHUTTLE_LOGBOOK_ENTRY_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is a container for the data queried from DAQ's logbook and logbook_shuttle tables.
//

#include <TObject.h>
#include "AliShuttle.h"

class TSQLServer;

class AliShuttleLogbookEntry: public TObject {

public:
	enum Status {
		kUnprocessed = 0,
		kInactive,
		kFailed,  // final
		kDone // final
	};

	AliShuttleLogbookEntry();
	AliShuttleLogbookEntry(Int_t run, UInt_t startTime, UInt_t endTime, Status* status=0);
	~AliShuttleLogbookEntry();

	AliShuttleLogbookEntry& operator=(const AliShuttleLogbookEntry& c);
	AliShuttleLogbookEntry(const AliShuttleLogbookEntry& c);
	virtual void Copy(TObject& c) const;

	Int_t GetRun() const {return fRun;}
	UInt_t GetStartTime() const  {return fStartTime;}
	UInt_t GetEndTime() const {return fEndTime;}

	void SetRun(Int_t run) {fRun=run;}
	void SetStartTime(UInt_t startTime) {fStartTime=startTime;}
	void SetEndTime(UInt_t endTime) {fEndTime=endTime;}

	Status GetDetectorStatus(const char* detCode) const;
	Status GetDetectorStatus(Int_t detPos) const;
	Status* GetDetectorStatus() const {return (Status*) fDetectorStatus;}

	void SetDetectorStatus(const char* detCode, Status status);
	void SetDetectorStatus(Status* status);
	void SetDetectorStatus(UInt_t detPos, Status status);

	Bool_t IsDone() const;

	static const char* GetDetectorStatusName(Status status);
        void Print(Option_t *option) const;

	// TODO Test only, remove later!
	Bool_t Connect();
	Bool_t QueryShuttleLogbook(Int_t runNumber=-1);
	Bool_t UpdateShuttleLogbook();
	Bool_t UpdateShuttleLogbook(const char* detCode, Status status);
	Bool_t InsertNewRun(Int_t runNumber=-1);


private:

	Int_t fRun;   			// Run number
	UInt_t fStartTime; 		// Run start time
	UInt_t fEndTime; 		// Run end time
	Status fDetectorStatus[AliShuttle::kNDetectors]; 	// Detector status array

	// TODO Test only, remove later!
	TSQLServer* fServer;	  	// pointer to the MySQLServer which handles the DAQ logbook

	ClassDef(AliShuttleLogbookEntry, 0)
};

#endif

