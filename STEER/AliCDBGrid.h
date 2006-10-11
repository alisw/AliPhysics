#ifndef ALICDBGRID_H
#define ALICDBGRID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBGrid						   //
//  access class to a DataBase in an AliEn storage  		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"

class AliCDBGrid: public AliCDBStorage {
	friend class AliCDBGridFactory;

public:
		  
	virtual Bool_t IsReadOnly() const {return kFALSE;}
	virtual Bool_t HasSubVersion() const {return kFALSE;}
	virtual Bool_t Contains(const char* path) const;
	virtual Int_t  GetLatestVersion(const char* path, Int_t run);
	virtual Int_t  GetLatestSubVersion(const char* path, Int_t run, Int_t version);
	virtual Bool_t IdToFilename(const AliCDBId& id, TString& filename) const;

protected:

	virtual AliCDBEntry*	GetEntry(const AliCDBId& queryId);
	virtual TList* 		GetEntries(const AliCDBId& queryId);
	virtual Bool_t          PutEntry(AliCDBEntry* entry);
	virtual TList* 		GetIdListFromFile(const char* fileName);

private:
 
	AliCDBGrid(const char *gridUrl, const char *user, const char* dbFolder, const char *se);

	virtual ~AliCDBGrid();

	AliCDBGrid(const AliCDBGrid& db);
	AliCDBGrid& operator = (const AliCDBGrid& db);

	Bool_t FilenameToId(TString& filename, AliCDBId& id);

	Bool_t PrepareId(AliCDBId& id);
	AliCDBId* GetId(const TObjArray& validFileIds, const AliCDBId& query);
	AliCDBEntry* GetEntryFromFile(TString& filename, const AliCDBId* dataId);

	Bool_t AddTag(TString& foldername, const char* tagname);
	Bool_t TagFileId(TString& filename, const AliCDBId* id);
	Bool_t TagFileMetaData(TString& filename, const AliCDBMetaData* md);

	void MakeQueryFilter(Int_t firstRun, Int_t lastRun, const AliCDBMetaData* md, TString& result) const;

	virtual void QueryValidFiles();

	TString    fGridUrl;	// Grid Url ("alien://aliendb4.cern.ch:9000")
	TString    fUser;	// User
	TString    fDBFolder;   // path of the DB folder
	TString    fSE;	  	// Storage Element

ClassDef(AliCDBGrid, 0)      // access class to a DataBase in an AliEn storage 
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBGridFactory					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBGridFactory: public AliCDBStorageFactory {

public:

	virtual Bool_t Validate(const char* gridString);
        virtual AliCDBParam* CreateParameter(const char* gridString);
	virtual ~AliCDBGridFactory(){}

protected:
        virtual AliCDBStorage* Create(const AliCDBParam* param);

        ClassDef(AliCDBGridFactory, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBGridParam					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBGridParam: public AliCDBParam {
	
public:
	AliCDBGridParam();
	AliCDBGridParam(const char* gridUrl, const char* user,
			const char* dbFolder, const char* se);
	
	virtual ~AliCDBGridParam();

	const TString& GridUrl() const {return fGridUrl;};
	const TString& GetUser() const {return fUser;};
	const TString& GetDBFolder() const {return fDBFolder;};
	const TString& GetSE() 	 const {return fSE;};

	virtual AliCDBParam* CloneParam() const;

        virtual ULong_t Hash() const;
        virtual Bool_t IsEqual(const TObject* obj) const;

private:
	TString fGridUrl;    // Grid url "Host:port"
	TString fUser;	     // User
	TString fDBFolder;   // path of the DB folder
	TString fSE;	     // Storage Element 

	ClassDef(AliCDBGridParam, 0);
};


#endif
