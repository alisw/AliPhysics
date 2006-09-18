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
class AliCDBPath;

class AliCDBStorage: public TObject {

public:
	AliCDBStorage();

	void SetURI(const TString& uri) {fURI = uri;}
	const TString& GetURI() const {return fURI;}
	const TString& GetType() const {return fType;}
	const TString& GetBaseFolder() const {return fBaseFolder;}


	void ReadSelectionFromFile(const char *fileName);
	
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
		
	void PrintSelectionList();

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
	
	Bool_t Put(TObject* object, AliCDBId& id,  AliCDBMetaData* metaData);
	Bool_t Put(AliCDBEntry* entry);


	virtual Bool_t IsReadOnly() const = 0;
	virtual Bool_t HasSubVersion() const = 0;
	virtual Bool_t Contains(const char* path) const = 0;
	virtual Bool_t IdToFilename(const AliCDBId& id, TString& filename) const = 0;

	void QueryCDB(Int_t run, const char* pathFilter="*",
			Int_t version=-1, AliCDBMetaData *mdFilter=0);
	void PrintQueryCDB();
	TList* GetQueryCDBList() {return &fValidFileIds;}

	virtual Int_t GetLatestVersion(const char* path, Int_t run)=0;
	virtual Int_t GetLatestSubVersion(const char* path, Int_t run, Int_t version=-1)=0;

protected:

	virtual ~AliCDBStorage();
	void    GetSelection(/*const*/ AliCDBId* id);
	virtual AliCDBEntry* GetEntry(const AliCDBId& query) = 0;
	virtual TList* GetEntries(const AliCDBId& query) = 0;
	virtual Bool_t PutEntry(AliCDBEntry* entry) = 0;
	virtual TList *GetIdListFromFile(const char* fileName)=0;
	virtual void   QueryValidFiles() = 0;

	TList fValidFileIds; 	// list of Id's of the files valid for a given run (cached as fRun)
	Int_t fRun;		// run number, used to manage list of valid files
	AliCDBPath fPathFilter;	// path filter, used to manage list of valid files
	Int_t fVersion;		// version, used to manage list of valid files
	AliCDBMetaData* fMetaDataFilter; // metadata, used to manage list of valid files

	TList fSelections; 	// list of selection criteria
	TString fURI;		// storage URI;
	TString fType;    //! Local, Grid: base folder name - Dump: file name
	TString fBaseFolder;    //! Local, Grid: base folder name - Dump: file name

private:
	AliCDBStorage(const AliCDBStorage & source);
	AliCDBStorage & operator=(const AliCDBStorage & source);

	ClassDef(AliCDBStorage, 0);
};

#endif
