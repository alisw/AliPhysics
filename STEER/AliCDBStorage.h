#ifndef ALI_CDB_STORAGE_H
#define ALI_CDB_STORAGE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBStorage						   //
//  interface to specific storage classes                          //
//  (AliCDBGrid, AliCDBLocal, AliCDBDump)			   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBId.h"
#include "AliCDBMetaData.h"

#include <TList.h>

class AliCDBEntry;

class AliCDBStorage: public TObject {

public:
	AliCDBStorage();
	virtual ~AliCDBStorage();	


	void AddSelection(const AliCDBId& selection);

	void AddSelection(const AliCDBPath& path, 
			const AliCDBRunRange& runRange,
			Int_t version,
			Int_t subVersion = -1);
	
	void AddSelection(const AliCDBPath& path,
			Int_t firstRun,
			Int_t lastRun,
			Int_t version,
			Int_t subVersion = -1);
			
	void RemoveSelection(const AliCDBId& selection);

	void RemoveSelection(const AliCDBPath& path,
			const AliCDBRunRange& runRange);
	
	void RemoveSelection(const AliCDBPath& path,
			Int_t firstRun = -1,
			Int_t lastRun = -1);

	void RemoveSelection(int position);
	void RemoveAllSelections();
		
	AliCDBEntry* Get(const AliCDBId& query);
	AliCDBEntry* Get(const AliCDBPath& path, Int_t runNumber, 
				Int_t version = -1, Int_t subVersion = -1);
	AliCDBEntry* Get(const AliCDBPath& path, const AliCDBRunRange& runRange,
				 Int_t version = -1, Int_t subVersion = -1);

	TList* GetAll(const AliCDBId& query);
	TList* GetAll(const AliCDBPath& path, Int_t runNumber, 
				Int_t version = -1, Int_t subVersion = -1);
	TList* GetAll(const AliCDBPath& path, const AliCDBRunRange& runRange,
				 Int_t version = -1, Int_t subVersion = -1); 
	
	void PrintSelectionList();


	Bool_t Put(TObject* object, AliCDBId& id,  AliCDBMetaData* metaData);
	Bool_t Put(AliCDBEntry* entry);


	virtual Bool_t IsReadOnly() = 0;
	virtual Bool_t HasSubVersion() = 0;

protected:
		
	AliCDBId GetSelection(const AliCDBId& id);
	virtual AliCDBEntry* GetEntry(const AliCDBId& query) = 0;
	virtual TList* GetEntries(const AliCDBId& query) = 0;
	virtual Bool_t PutEntry(AliCDBEntry* entry) = 0; 

private:

	TList fSelections; // list of selection criteria

	ClassDef(AliCDBStorage, 0);
};

#endif
