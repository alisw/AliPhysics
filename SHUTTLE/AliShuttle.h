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
	enum TestMode { kNone = 0, kSkipDCS = 1, kErrorDCS = 2, kErrorFXSSources = 4, kErrorFXSFiles = 8, kErrorOCDB = 16, kErrorStorage = 32, kErrorGrid = 64 };
	enum EMailTarget { kDCSEMail = 0, kFXSEMail, kPPEMail };

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
	UInt_t GetCurrentYear() const;
	
	const char* GetLHCPeriod() const;

	virtual Bool_t Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
			Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
	virtual Bool_t StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData);
	virtual Bool_t StoreReferenceFile(const char* detector, const char* localFile, const char* gridFileName);
	virtual Bool_t StoreRunMetadataFile(const char* localFile, const char* gridFileName);
	virtual const char* GetFile(Int_t system, const char* detector,
		const char* id, const char* source);
	virtual TList* GetFileSources(Int_t system, const char* detector, const char* id = 0);
	virtual TList* GetFileIDs(Int_t system, const char* detector, const char* source);
	virtual const char* GetRunParameter(const char* lbEntry);
	virtual UInt_t GetStartTimeDCSQuery();
	virtual UInt_t GetEndTimeDCSQuery();
	virtual AliCDBEntry* GetFromOCDB(const char* detector, const AliCDBPath& path);
	virtual const char* GetRunType();
    	virtual Bool_t GetHLTStatus();
	virtual const char* GetTriggerConfiguration(); 
	virtual const char* GetCTPTimeParams(); 
	virtual const char* GetTriggerDetectorMask(); 
	virtual void Log(const char* detector, const char* message, UInt_t level=3);

	void SetLogbookEntry(AliShuttleLogbookEntry* entry) {fLogbookEntry=entry;}
	
	void SetTestMode(TestMode testMode) { fTestMode = testMode; }
	void SetReadTestModeFromLog(Bool_t flag) { fReadTestMode = flag; }

	Bool_t Connect(Int_t system);

	static void SetMainCDB (TString mainCDB) {fgkMainCDB = mainCDB;}
	static void SetLocalCDB (TString localCDB) {fgkLocalCDB = localCDB;}

	static void SetMainRefStorage (TString mainRefStorage) {fgkMainRefStorage = mainRefStorage;}
	static void SetLocalRefStorage (TString localRefStorage) {fgkLocalRefStorage = localRefStorage;}

	static void SetShuttleTempDir (const char* tmpDir);
	static void SetShuttleLogDir (const char* logDir);

	virtual void SendMLFromDet(const char* value);
	virtual TString* GetLTUConfig(const char* det);

private:
	AliShuttle(const AliShuttle& other);
	AliShuttle& operator= (const AliShuttle& other);

	Int_t ProcessCurrentDetector();

	AliShuttleLogbookEntry* QueryRunParameters(Int_t run);
	Bool_t QueryShuttleLogbook(const char* whereClause, TObjArray& entries);
	void CountOpenRuns();
	Bool_t RetrieveConditionsData(const TObjArray& shuttleLogbookEntries);

	TMap* GetValueSet(const char* host, Int_t port, const TSeqCollection* entries,
			      DCSType type, Int_t valueSet);

	Bool_t RetrieveFile(UInt_t system, const char* daqFileName, const char* localFileName);

	Bool_t UpdateTable();
	Bool_t UpdateTableSkippedCase(const char* detector="ALL");
	Bool_t UpdateTableFailCase();

	Bool_t StoreLocally(const TString& localUri, const AliCDBPath& path, TObject* object,
				AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);

	Bool_t StoreOCDB();
	Int_t StoreOCDB(const TString& uri);
	Bool_t CopyFileLocally(const char* localFile, const TString& target);
	Bool_t CopyFilesToGrid(const char* type);
	void CleanLocalStorage(const TString& uri);
	Bool_t CleanReferenceStorage(const char* detector);
	void RemoveFile(const char* filename);
	const char* GetRefFilePrefix(const char* base, const char* detector);

	AliShuttleStatus* ReadShuttleStatus();
	Bool_t WriteShuttleStatus(AliShuttleStatus* status);
	Bool_t ContinueProcessing();
	void UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount = kFALSE);
	Bool_t UpdateShuttleLogbook(const char* detector, const char* status=0);
	Bool_t SendMail(EMailTarget target, Int_t system = -1);
	Int_t GetMem(Int_t pid);
	
	TString GetLogFileName(const char* detector) const;

	void SetLastAction(const char* action);
	
	void SendAlive();
	void SendMLDetInfo();
	void SendMLRunInfo(const char* status);
	virtual Bool_t TouchFile();

	const AliShuttleConfig* fConfig; 	// pointer to configuration object

	UInt_t fTimeout; 	// DCS server connection timeout parameter
	Int_t fRetries; 	// Number of DCS server connection retries

	TMap fPreprocessorMap; 	// list of detector Preprocessors ("DET", "Preprocessor")

	AliShuttleLogbookEntry* fLogbookEntry;   //! current Shuttle logbook entry
	TString fCurrentDetector; // current detector
	Bool_t fFirstProcessing;  // processing this detector the first time in this run

	TSQLServer *fServer[5]; 	// pointer to the four FXS + Run & Shuttle logbook servers
	Bool_t fFXSCalled[4];		// FXS call status
	TList  fFXSlist[4];		// List of files retrieved from each FXS
	Int_t  fFXSError;		// Variable to keep track of any FXS errors; contains -1 for no error, kDAQ, kDCS, kHLT otherwise

	AliCDBEntry* fStatusEntry; // last CDB entry containing a AliShuttleStatus retrieved

	TMutex* fMonitoringMutex;   // mutex to lock the monitoring class members
	UInt_t fLastActionTime;    // time of last action for monitoring
	TString fLastAction;       // string description for last action

	Bool_t fFirstUnprocessed[AliShuttleInterface::kNDetectors];       // array of flags for first unprocessed dets

	TMonaLisaWriter* fMonaLisa;  // ML instance that sends the processing information

	TestMode fTestMode;          // sets test mode flags, that e.g. simulate a dcs error etc.
	Bool_t fReadTestMode;        // Reads the test mode from the log entry of the given run (only for test)
	
	Bool_t fOutputRedirected;    // is the output redirected to a file

	ClassDef(AliShuttle, 0);
};

#endif
