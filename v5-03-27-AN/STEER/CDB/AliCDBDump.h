#ifndef ALI_CDB_DUMP_H
#define ALI_CDB_DUMP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBDump						   //
//  access class to a DataBase in a dump storage (single file)     //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBStorage.h"
#include "AliCDBManager.h"

class TDirectory;
class TFile;

class AliCDBDump: public AliCDBStorage {
	friend class AliCDBDumpFactory;

public:

	virtual Bool_t IsReadOnly() const {return fReadOnly;};
	virtual Bool_t HasSubVersion() const {return kFALSE;};
	virtual Bool_t Contains(const char* path) const;
	virtual Int_t GetLatestVersion(const char* path, Int_t run);
	virtual Int_t GetLatestSubVersion(const char* path, Int_t run, Int_t version=-1);
	virtual Bool_t IdToFilename(const AliCDBId& id, TString& filename) const;
	virtual void SetRetry(Int_t /* nretry */, Int_t /* initsec */);

protected:

	virtual AliCDBEntry* 	GetEntry(const AliCDBId& query);
	virtual AliCDBId* 	GetEntryId(const AliCDBId& query);
        virtual TList* 		GetEntries(const AliCDBId& query);
        virtual Bool_t 		PutEntry(AliCDBEntry* entry);
	virtual TList* 		GetIdListFromFile(const char* fileName);

private:

	AliCDBDump(const AliCDBDump & source);
	AliCDBDump & operator=(const AliCDBDump & source);
	AliCDBDump(const char* dbFile, Bool_t readOnly);
	virtual ~AliCDBDump();	

	Bool_t KeyNameToId(const char* keyname, AliCDBRunRange& runRange,
			Int_t& version, Int_t& subVersion);
	Bool_t IdToKeyName(const AliCDBRunRange& runRange, Int_t version,
		        Int_t subVersion, TString& keyname); 	

	Bool_t MkDir(const TString& dir);


	Bool_t PrepareId(AliCDBId& id);
//	Bool_t GetId(const AliCDBId& query, AliCDBId& result);
	AliCDBId* GetId(const AliCDBId& query);

	virtual void QueryValidFiles();

	void GetEntriesForLevel0(const AliCDBId& query, TList* result);
	void GetEntriesForLevel1(const AliCDBId& query, TList* result);

	TFile* fFile;		// Dump file
	Bool_t fReadOnly;	// ReadOnly flag

	ClassDef(AliCDBDump, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBDumpFactory					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBDumpFactory: public AliCDBStorageFactory {

public:

        virtual Bool_t Validate(const char* dbString);
        virtual AliCDBParam* CreateParameter(const char* dbString);

protected:
        virtual AliCDBStorage* Create(const AliCDBParam* param);

        ClassDef(AliCDBDumpFactory, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBDumpParam					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBDumpParam: public AliCDBParam {

public:
        AliCDBDumpParam();
        AliCDBDumpParam(const char* dbPath, Bool_t readOnly = kFALSE);
        
        virtual ~AliCDBDumpParam();

        const TString& GetPath() const {return fDBPath;};
	Bool_t IsReadOnly() const {return fReadOnly;};

	virtual AliCDBParam* CloneParam() const;

	virtual ULong_t Hash() const;
	virtual Bool_t IsEqual(const TObject* obj) const;
	
private:

        TString fDBPath;	// Dump file path name
	Bool_t fReadOnly;	// ReadOnly flag

	ClassDef(AliCDBDumpParam, 0);
};

#endif
