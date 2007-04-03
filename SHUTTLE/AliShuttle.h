#ifndef ALI_SHUTTLE_H
#define ALI_SHUTTLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//
// This class is the main manager for 
// AliShuttle. It organizes the data retrieval
// from DCS and call the interface methods of
// AliPreprocessor.
//

#include <TMap.h>
#include <TString.h>
#include <TList.h>

#include "AliShuttleInterface.h"
#include "AliShuttleStatus.h"

class TObject;
class AliShuttleConfig;
class AliShuttleLogbookEntry;
class AliPreprocessor;
class TSQLServer;
class TMutex;
class TMonaLisaWriter;

class AliShuttle: public AliShuttleInterface {
public:
	enum DCSType {kAlias=0, kDP};
	enum returnCodes {kDCSErrCode=20051975, kStorErrCode=20051976};

	AliShuttle(const AliShuttleConfig* config, UInt_t timeout = 5000, Int_t retries = 5);
	virtual ~AliShuttle();

	virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);

	Bool_t Collect(Int_t run = -1);

	Bool_t Process(AliShuttleLogbookEntry* entry);

	// monitoring functions
	ULong_t GetTimeOfLastAction() const;
	const TString GetLastAction() const;

	Int_t GetCurrentRun() const;
	UInt_t GetCurrentStartTime() const;
	UInt_t GetCurrentEndTime() const;

	virtual Bool_t Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
			Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
	virtual Bool_t StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData);
	virtual const char* GetFile(Int_t system, const char* detector,
		const char* id, const char* source);
	virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
	virtual const char* GetRunParameter(const char* lbEntry);
	virtual AliCDBEntry* GetFromOCDB(const AliCDBPath& path);
	virtual const char* GetRunType(const char* detCode);
	virtual void Log(const char* detector, const char* message);

	//TODO Test only, remove later !
	void SetProcessDCS(Bool_t process) {fgkProcessDCS = process;}
	void SetLogbookEntry(AliShuttleLogbookEntry* entry) {fLogbookEntry=entry;}

	Bool_t Connect(Int_t system);

	static void SetMainCDB (TString mainCDB) {fgkMainCDB = mainCDB;}
	static void SetLocalCDB (TString localCDB) {fgkLocalCDB = localCDB;}

	static void SetMainRefStorage (TString mainRefStorage) {fgkMainRefStorage = mainRefStorage;}
	static void SetLocalRefStorage (TString localRefStorage) {fgkLocalRefStorage = localRefStorage;}

	static void SetShuttleTempDir (const char* tmpDir);
	static void SetShuttleLogDir (const char* logDir);

private:
	AliShuttle(const AliShuttle& other);
	AliShuttle& operator= (const AliShuttle& other);

	UInt_t ProcessCurrentDetector();

	AliShuttleLogbookEntry* QueryRunParameters(Int_t run);
	Bool_t QueryShuttleLogbook(const char* whereClause, TObjArray& entries);
	Bool_t RetrieveConditionsData(const TObjArray& shuttleLogbookEntries);

	Bool_t GetValueSet(const char* host, Int_t port, const char* alias,
				TObjArray* result, DCSType type);

	Bool_t RetrieveFile(UInt_t system, const char* daqFileName, const char* localFileName);

	Bool_t UpdateTable();
	Bool_t UpdateTableFailCase();

	Bool_t StoreLocally(const TString& localUri, const AliCDBPath& path, TObject* object,
				AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);

	Bool_t StoreOCDB();
	Bool_t StoreOCDB(const TString& uri);
	void CleanLocalStorage(const TString& uri);
	void RemoveFile(const char* filename);

	AliShuttleStatus* ReadShuttleStatus();
	Bool_t WriteShuttleStatus(AliShuttleStatus* status);
	Bool_t ContinueProcessing();
	void UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount = kFALSE);
	Bool_t UpdateShuttleLogbook(const char* detector, const char* status=0);
	Bool_t SendMail();

	void SetLastAction(const char* action);
	
	void SendMLInfo();

	void SetRunType();

	const AliShuttleConfig* fConfig; 	// pointer to configuration object

	UInt_t fTimeout; 	// DCS server connection timeout parameter
	Int_t fRetries; 	// Number of DCS server connection retries

	TMap fPreprocessorMap; 	// list of detector Preprocessors ("DET", "Preprocessor")

	AliShuttleLogbookEntry* fLogbookEntry;   //! current Shuttle logbook entry
	TString fCurrentDetector; // current detector

	TSQLServer *fServer[4]; 	// pointer to the three FXS + Run & Shuttle logbook servers
	Bool_t fFXSCalled[3];		// FXS call status
	TList  fFXSlist[3];		// List of files retrieved from each FXS

	AliCDBEntry* fStatusEntry; // last CDB entry containing a AliShuttleStatus retrieved

	TMutex* fMonitoringMutex;   // mutex to lock the monitoring class members
	UInt_t fLastActionTime;    // time of last action for monitoring
	TString fLastAction;       // string description for last action

	Bool_t fFirstUnprocessed[AliShuttleInterface::kNDetectors];       // array of flags for first unprocessed dets

	TMonaLisaWriter* fMonaLisa;  // ML instance that sends the processing information

	//TODO Test only, remove later !
	static Bool_t fgkProcessDCS; // flag to enable DCS archive data processing

	ClassDef(AliShuttle, 0);
};

#endif
