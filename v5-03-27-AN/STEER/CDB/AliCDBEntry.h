#ifndef ALI_CDB_ENTRY_H
#define ALI_CDB_ENTRY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBEntry						   //
//  container for an object, it identity (AliCDBId)  		   //
//  and its metaData (AliCDBMetaData) 				   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBId.h"
#include "AliCDBMetaData.h"

class AliCDBEntry: public TObject {

public:
	AliCDBEntry();

	AliCDBEntry(TObject* object, const AliCDBId& id,  
			AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, const AliCDBRunRange& runRange,
			AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, const AliCDBRunRange& runRange,
			Int_t version, AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, const AliCDBRunRange& runRange,
			Int_t version, Int_t subVersion, 
			AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, Int_t firstRun, Int_t lastRun,
			AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, Int_t firstRun, Int_t lastRun,
			Int_t version, AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	AliCDBEntry(TObject* object, const AliCDBPath& path, Int_t firstRun, Int_t lastRun,
			Int_t version, Int_t subVersion, 
			AliCDBMetaData* metaData, Bool_t owner = kFALSE);

	virtual ~AliCDBEntry();


	void 		SetId(const AliCDBId& id) {fId = id;};
	AliCDBId& 	GetId() {return fId;};
	const AliCDBId& GetId() const {return fId;};
	void 		PrintId() const;
	
	void 		SetObject(TObject* object) {fObject = object;};
	TObject* 	GetObject() {return fObject;};
	const TObject* 	GetObject() const {return fObject;};	

	void 			SetMetaData(AliCDBMetaData* metaData) {fMetaData = metaData;};
	AliCDBMetaData* 	GetMetaData() {return fMetaData;};
	const AliCDBMetaData* 	GetMetaData() const {return fMetaData;};
	void 			PrintMetaData() const {fMetaData->PrintMetaData();}

	void 	SetOwner(Bool_t owner) {fIsOwner = owner;};
	Bool_t 	IsOwner() const {return fIsOwner;};
	
  	void 	SetVersion(Int_t version) {fId.SetVersion(version);}
  	void 	SetSubVersion(Int_t subVersion) {fId.SetSubVersion(subVersion);}
	
	const TString 	GetLastStorage() const {return fId.GetLastStorage();};
	void  		SetLastStorage(TString lastStorage) {fId.SetLastStorage(lastStorage);};

private:
	
	AliCDBEntry(const AliCDBEntry& other); // no copy ctor
	void operator= (const AliCDBEntry& other); // no assignment op

	TObject* fObject;		// object
	AliCDBId fId;			// entry ID
	AliCDBMetaData* fMetaData; 	// metaData
	Bool_t fIsOwner; 		// ownership flag
	
	ClassDef(AliCDBEntry, 1);
};

#endif
