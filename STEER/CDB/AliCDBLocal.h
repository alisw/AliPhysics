#ifndef ALI_CDB_LOCAL_H
#define ALI_CDB_LOCAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBLocal						   //
//  access class to a DataBase in a local storage                  //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBStorage.h"
#include "AliCDBManager.h"

class AliCDBLocal: public AliCDBStorage {
	friend class AliCDBLocalFactory;

public:

	virtual Bool_t IsReadOnly() const {return kFALSE;};
	virtual Bool_t HasSubVersion() const {return kTRUE;};
	virtual Bool_t Contains(const char* path) const;
	virtual Int_t  GetLatestVersion(const char* path, Int_t run);
	virtual Int_t  GetLatestSubVersion(const char* path, Int_t run, Int_t version=-1);
	virtual Bool_t IdToFilename(const AliCDBId& id, TString& filename) const;
	virtual void SetRetry(Int_t /* nretry */, Int_t /* initsec */);

protected:

	virtual AliCDBEntry*    GetEntry(const AliCDBId& queryId);
	virtual AliCDBId* 	GetEntryId(const AliCDBId& queryId);
        virtual TList* 		GetEntries(const AliCDBId& queryId);
        virtual Bool_t 		PutEntry(AliCDBEntry* entry);
	virtual TList* 		GetIdListFromFile(const char* fileName);

private:

	AliCDBLocal(const AliCDBLocal & source);
	AliCDBLocal & operator=(const AliCDBLocal & source);
	AliCDBLocal(const char* baseDir);
	virtual ~AliCDBLocal();
	
	Bool_t FilenameToId(const char* filename, AliCDBRunRange& runRange, 
			Int_t& version, Int_t& subVersion);

	Bool_t PrepareId(AliCDBId& id);
//	Bool_t GetId(const AliCDBId& query, AliCDBId& result);
	AliCDBId* GetId(const AliCDBId& query);

	virtual void QueryValidFiles();

	void GetEntriesForLevel0(const char* level0, const AliCDBId& query, TList* result);
	void GetEntriesForLevel1(const char* level0, const char* Level1,
			const AliCDBId& query, TList* result);

	TString fBaseDirectory; // path of the DB folder

	ClassDef(AliCDBLocal, 0); // access class to a DataBase in a local storage
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBLocalFactory					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBLocalFactory: public AliCDBStorageFactory {

public:

	virtual Bool_t Validate(const char* dbString);
        virtual AliCDBParam* CreateParameter(const char* dbString);

protected:
        virtual AliCDBStorage* Create(const AliCDBParam* param);

        ClassDef(AliCDBLocalFactory, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBLocalParam					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBLocalParam: public AliCDBParam {
	
public:
	AliCDBLocalParam();
	AliCDBLocalParam(const char* dbPath);
	AliCDBLocalParam(const char* dbPath, const char* uri);
	
	virtual ~AliCDBLocalParam();

	const TString& GetPath() const {return fDBPath;};

	virtual AliCDBParam* CloneParam() const;

        virtual ULong_t Hash() const;
        virtual Bool_t IsEqual(const TObject* obj) const;

private:

	TString fDBPath; // path of the DB folder

	ClassDef(AliCDBLocalParam, 0);
};

#endif
