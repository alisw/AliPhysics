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
class AliCDBMetaData;
class TSQLServer;
class AliCDBEntry;
class AliCDBPath;

class AliShuttle: public AliShuttleInterface {
public:
	enum { kNDetectors=17 }; // number of subdetectors in ALICE

	AliShuttle(const AliShuttleConfig* config, UInt_t timeout = 5000, Int_t retries = 5);
	virtual ~AliShuttle();

	virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);

	Bool_t Collect(Int_t run);
	Bool_t CollectNew();
	Bool_t CollectAll();

	Bool_t Process(AliShuttleLogbookEntry* entry);

	Int_t GetCurrentRun() const;
	UInt_t GetCurrentStartTime() const;
	UInt_t GetCurrentEndTime() const;

	virtual UInt_t Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
			Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
	virtual UInt_t StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData);
	virtual const char* GetFile(Int_t system, const char* detector,
		const char* id, const char* source);
	virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
	virtual void Log(const char* detector, const char* message);

	static TString GetMainCDB () {return fgkMainCDB;}
	static void SetMainCDB (TString mainCDB) {fgkMainCDB = mainCDB;}
	static TString GetLocalCDB () {return fgkLocalCDB;}
	static void SetLocalCDB (TString localCDB) {fgkLocalCDB = localCDB;}

	static TString GetMainRefStorage() {return fgkMainRefStorage;}
	static void SetMainRefStorage (TString mainRefStorage) {fgkMainRefStorage = mainRefStorage;}
	static TString GetLocalRefStorage() {return fgkLocalRefStorage;}
	static void SetLocalRefStorage (TString localRefStorage) {fgkLocalRefStorage = localRefStorage;}

	//TODO Test only, remove later !
	void SetProcessDCS(Bool_t process) {fgkProcessDCS = process;}

	static const char* GetDetCode(const char* detector);
	static const char* GetDetCode(UInt_t detPos);
	static const Int_t GetDetPos(const char* detCode);
	static const UInt_t NDetectors() {return kNDetectors;}
	static const char* GetShuttleTempDir() {return fgkShuttleTempDir;}

	Bool_t Connect(Int_t system);

private:
	AliShuttle(const AliShuttle& other);
	AliShuttle& operator= (const AliShuttle& other);

	UInt_t ProcessCurrentDetector();

	Bool_t QueryRunParameters(Int_t& run, UInt_t& startTime, UInt_t& endTime);
	Bool_t QueryShuttleLogbook(const char* whereClause, TObjArray& entries);
	Bool_t RetrieveConditionsData(const TObjArray& shuttleLogbookEntries);

	Bool_t GetValueSet(const char* host, Int_t port, const char* alias, TObjArray* result);

	const char* GetDAQFileName(const char* detector, const char* id, const char* source);
	Bool_t RetrieveDAQFile(const char* daqFileName, const char* localFileName);
	TList* GetDAQFileSources(const char* detector, const char* id);
	Bool_t UpdateDAQTable();

	const char* GetDCSFileName(const char* detector, const char* id, const char* source);
//	Bool_t RetrieveDCSFile(const char* daqFileName const char* localFileName);
	TList* GetDCSFileSources(const char* detector, const char* id);

	const char* GetHLTFileName(const char* detector, const char* id, const char* source);
//	Bool_t RetrieveHLTFile(const char* daqFileName, const char* localFileName;
	TList* GetHLTFileSources(const char* detector, const char* id);

	UInt_t WriteToCDB(const char* mainUri, const char* localUri,
				const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
				Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);

	Bool_t TryToStoreAgain();
	Bool_t TryToStoreAgain(TString& storageType);

  	AliShuttleStatus* ReadShuttleStatus();
  	Bool_t WriteShuttleStatus(AliShuttleStatus* status);
  	Bool_t ContinueProcessing();
  	void UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount = kFALSE);
  	Bool_t UpdateShuttleLogbook(const char* detector, const char* status=0);

	const AliShuttleConfig* fConfig; 	// pointer to configuration object

//	static const UInt_t fgkNDetectors = 17;		   	//! number of detectors
	static const char*  fgkDetectorName[kNDetectors]; 	//! names of detectors
	static const char*  fgkDetectorCode[kNDetectors]; 	//! codes of detectors
	static TString 	    fgkMainCDB;		// URI of the main (Grid) CDB storage
	static TString 	    fgkLocalCDB;		//! URI of the local backup CDB storage
	static TString 	    fgkMainRefStorage;	// URI of the main (Grid) REFERENCE storage
	static TString 	    fgkLocalRefStorage;	// URI of the local REFERENCE storage
	static const char*  fgkShuttleTempDir;	// base path of SHUTTLE temp folder
	static const char*  fgkShuttleLogDir;	// path of SHUTTLE log folder

	UInt_t fTimeout; 	// DCS server connection timeout parameter
	Int_t fRetries; 	// Number of DCS server connection retries

	TMap fPreprocessorMap; 	// list of detector Preprocessors ("DET", "Preprocessor")

	AliShuttleLogbookEntry* fLogbookEntry;   //! current Shuttle logbook entry
	TString fCurrentDetector; // current detector

	TSQLServer *fServer[3]; 	// pointer to the three FS logbook servers
	Bool_t fFESCalled[3];		// FES call status
	TList  fFESlist[3];		// List of files retrieved from each FES

	AliCDBEntry* fStatusEntry; // last CDB entry containing a AliShuttleStatus retrieved
	Bool_t fGridError; 	   // Grid storage error flag

	//TODO Test only, remove later !
	static Bool_t fgkProcessDCS; // flag to enable DCS archive data processing

	ClassDef(AliShuttle, 0);
};

#endif
