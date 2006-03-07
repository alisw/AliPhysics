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

class AliShuttleConfig;
class AliCDBStorage;
class AliCDBMetaData;
class AliCDBParam;
class AliCDBPreProcessor;

class AliShuttle: public TObject {
public:
	AliShuttle(const AliShuttleConfig* config, AliCDBStorage* cdbStorage,
		UInt_t timeout = 5000, Int_t retries = 5);
	virtual ~AliShuttle();

	void RegisterCDBPreProcessor(AliCDBPreProcessor* processor);
	
	Bool_t Process(Int_t run, UInt_t startTime, UInt_t endTime);
	Bool_t Process(Int_t run, UInt_t startTime, UInt_t endTime,
		const char* detector);

	Int_t GetCurrentRun() const {return fCurrentRun;};
	UInt_t GetCurrentStartTime() const {return fCurrentStartTime;};
	UInt_t GetCurrentEndTime() const {return fCurrentEndTime;};

	Bool_t Store(const char* detector, const char* detSpec,
			TObject* object, AliCDBMetaData* metaData);

private:
	const AliShuttleConfig* fConfig;
	AliCDBStorage* fStorage;
	UInt_t fTimeout;
	Int_t fRetries;
	
	TMap fPreProcessorMap;	

	Int_t fCurrentRun;
	UInt_t fCurrentStartTime;
	UInt_t fCurrentEndTime;

	Bool_t GetValueSet(const char* host, Int_t port, const char* alias,
			TObjArray& result);
	
	ClassDef(AliShuttle, 0);
};

#endif
