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
class AliPreprocessor;
class AliCDBMetaData;
class TSQLServer;
class AliCDBEntry;

class AliShuttle: public AliShuttleInterface {
public:
	AliShuttle(const AliShuttleConfig* config, UInt_t timeout = 5000, Int_t retries = 5);
	virtual ~AliShuttle();

	virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);

	Bool_t Process(Int_t run, UInt_t startTime, UInt_t endTime);
	Bool_t Process();

	Int_t GetCurrentRun() const {return fCurrentRun;};
	UInt_t GetCurrentStartTime() const {return fCurrentStartTime;};
	UInt_t GetCurrentEndTime() const {return fCurrentEndTime;};

	virtual UInt_t Store(const char* detector, TObject* object, AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
	virtual const char* GetFile(Int_t system, const char* detector,
		const char* id, const char* source);
	virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
	virtual void Log(const char* detector, const char* message);

	static TString GetLocalURI () {return fgkLocalUri;}
	static void SetLocalURI (TString localUri) {fgkLocalUri = localUri;}

	// TODO Test only, remove later!
	void SetCurrentRun(int run) {fCurrentRun=run;}

	static const char* GetDetCode(const char* detector);
	static const char* GetShuttleTempDir() {return fgkShuttleTempDir;}

	Bool_t Connect(Int_t system);


private:
	AliShuttle(const AliShuttle& other);
	AliShuttle& operator= (const AliShuttle& other);

	Bool_t GetValueSet(const char* host, Int_t port, const char* alias,
			TObjArray& result);

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

  AliShuttleStatus* ReadShuttleStatus();
  Bool_t WriteShuttleStatus(AliShuttleStatus* status);
  Bool_t ContinueProcessing();
  void UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount = kFALSE);

	const AliShuttleConfig* fConfig; 	//! pointer to configuration object

	static const Int_t fgkNDetectors = 17;		   	//! number of detectors
	static const char* fgkDetectorName[fgkNDetectors]; 	//! names of detectors
	static const char* fgkDetectorCode[fgkNDetectors]; 	//! codes of detectors
	static TString 	   fgkLocalUri;		//! URI of the local backup storage location
	static const char* fgkShuttleTempDir;	//! base path of SHUTTLE temp folder
	static const char* fgkShuttleLogDir;	//! path of SHUTTLE log folder

	UInt_t fTimeout; 	//! DCS server connection timeout parameter
	Int_t fRetries; 	//! Number of DCS server connection retries

	TMap fPreprocessorMap; 	//! list of detector Preprocessors ("DET", "Preprocessor")

	Int_t fCurrentRun;  		//! run currenty processed
	UInt_t fCurrentStartTime; 	//! Run Start time
	UInt_t fCurrentEndTime; 	//! Run end time

  TString fCurrentDetector; // current detector

	TSQLServer *fServer[3]; 	//! pointer to the three FS logbook servers

	Bool_t fFESCalled[3];		//! FES call status
	TList  fFESlist[3];		//! List of files retrieved from each FES

  AliCDBEntry* fStatusEntry; //! last CDB entry containing a AliShuttleStatus retrieved

	ClassDef(AliShuttle, 0);
};

#endif
