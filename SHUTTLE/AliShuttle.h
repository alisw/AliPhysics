#ifndef ALI_SHUTTLE_H
#define ALI_SHUTTLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//
// This class is the main manager for 
// AliShuttle. It organizes the data retrieval
// from DCS and call the interface methods of
// AliCDBPreProcessor.
//

#include <TObject.h>
#include <TMap.h>
#include <TString.h>

#include "AliShuttleInterface.h"

class AliShuttleConfig;
class AliPreprocessor;
class AliCDBMetaData;

class AliShuttle: public AliShuttleInterface {
public:
	AliShuttle(const AliShuttleConfig* config, UInt_t timeout = 5000, Int_t retries = 5);
	virtual ~AliShuttle();

	virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);

	Bool_t Process(Int_t run, UInt_t startTime, UInt_t endTime);
	Bool_t Process(Int_t run, UInt_t startTime, UInt_t endTime,
		const char* detector);

	Int_t GetCurrentRun() const {return fCurrentRun;};
	UInt_t GetCurrentStartTime() const {return fCurrentStartTime;};
	UInt_t GetCurrentEndTime() const {return fCurrentEndTime;};

	virtual UInt_t Store(const char* detector, TObject* object, AliCDBMetaData* metaData);
	virtual const char* GetFile(Int_t system, const char* detector,
		const char* id, const char* source);
	virtual TList* GetFileSources(Int_t system, const char* detector, const char* id);
	virtual void Log(const char* detector, const char* message);

	static TString GetLocalURI () {return fgkLocalUri;}
	static void SetLocalURI (TString localUri) {fgkLocalUri = localUri;}

private:

	static TString fgkLocalUri;

	void ClearLog() {fLog = "";}
	void StoreLog(Int_t run);
  	const AliShuttleConfig* fConfig;

//	AliCDBStorage* fLocalStorage;

	UInt_t fTimeout;
	Int_t fRetries;

	TMap fPreprocessorMap;

	Int_t fCurrentRun;
	UInt_t fCurrentStartTime;
	UInt_t fCurrentEndTime;

	TString fLog;

	Bool_t GetValueSet(const char* host, Int_t port, const char* alias,
			TObjArray& result);
	
	ClassDef(AliShuttle, 0);
};

#endif
