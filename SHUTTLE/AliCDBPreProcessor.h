#ifndef ALI_CDB_PRE_PROCESSOR_H
#define ALI_CDB_PRE_PROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is the CDBPreProcessor interface,
// supposed to be implemented by any detector
// interested in immediate processing of data 
// which is retrieved from DCS.
//

#include <TNamed.h>

class AliShuttle;
class AliCDBMetaData;

class AliCDBPreProcessor: public TNamed {

	friend class AliShuttle;

public:
	AliCDBPreProcessor(const char* detector);
	virtual ~AliCDBPreProcessor();

	Int_t GetRun() const;
	UInt_t GetStartTime() const;
	UInt_t GetEndTime() const;

	Bool_t Store(const char* specType, TObject* object, 
		AliCDBMetaData* metaData);

	void SetShuttle(AliShuttle* shuttle) {fShuttle = shuttle;};
	AliShuttle* GetShuttle() const {return fShuttle;};

protected:

	virtual void Initialize(Int_t /*run*/, UInt_t /*startTime*/, 
			UInt_t /*endTime*/) {};

	virtual void Finalize() {};

	virtual void Process(const char* alias, TList& valueSet, 
			Bool_t hasError) = 0;

private:

	AliShuttle* fShuttle;

	ClassDef(AliCDBPreProcessor, 0);
};

#endif
