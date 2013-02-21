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
	virtual Bool_t IdToFilename(const AliCDBId& id, TString& filename) const;
	virtual void SetRetry(Int_t nretry, Int_t initsec);
	virtual void SetMirrorSEs(const char* mirrors) {fMirrorSEs=mirrors;}
	virtual const char* GetMirrorSEs() const {return fMirrorSEs;}


protected:

	virtual AliCDBEntry*	GetEntry(const AliCDBId& queryId);
	virtual AliCDBId*	GetEntryId(const AliCDBId& queryId);
	virtual TList* 		GetEntries(const AliCDBId& queryId);
	virtual Bool_t          PutEntry(AliCDBEntry* entry, const char* mirrors="");
	virtual TList* 		GetIdListFromFile(const char* fileName);

private:

	AliCDBGrid(const char *gridUrl, const char *user, const char* dbFolder,
	           const char *se, const char* cacheFolder, Bool_t operateDisconnected,
		   Long64_t cacheSize, Long_t cleanupInterval);

	virtual ~AliCDBGrid();

	AliCDBGrid(const AliCDBGrid& db);
	AliCDBGrid& operator = (const AliCDBGrid& db);

	Bool_t FilenameToId(TString& filename, AliCDBId& id);

	Bool_t PrepareId(AliCDBId& id);
	AliCDBId* GetId(const TObjArray& validFileIds, const AliCDBId& query);
	AliCDBEntry* GetEntryFromFile(TString& filename, AliCDBId* dataId);

	// TODO  use AliEnTag classes!
	Bool_t AddTag(TString& foldername, const char* tagname);
	Bool_t TagFileId(TString& filename, const AliCDBId* id);
	Bool_t TagFileMetaData(TString& filename, const AliCDBMetaData* md);
	Bool_t TagShortLived(TString& filename, Bool_t value);

	void MakeQueryFilter(Int_t firstRun, Int_t lastRun, const AliCDBMetaData* md, TString& result) const;

	virtual void QueryValidFiles();

	TString    fGridUrl;	 // Grid Url ("alien://aliendb4.cern.ch:9000")
	TString    fUser;	 // User
	TString    fDBFolder;    // path of the DB folder
	TString    fSE;	  	 // Storage Element
	TString    fMirrorSEs;	 // Mirror Storage Elements
	TString    fCacheFolder; // local cache folder
	Bool_t     fOperateDisconnected; // Operate disconnected flag
	Long64_t   fCacheSize;           // local cache size (in bytes)
	Long_t     fCleanupInterval;     // local cache cleanup interval

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
			const char* dbFolder, const char* se,
			const char* cacheFolder, Bool_t operateDisconnected,
			Long64_t cacheSize, Long_t cleanupInterval);
	
	virtual ~AliCDBGridParam();

	const TString& GridUrl() const {return fGridUrl;}
	const TString& GetUser() const {return fUser;}
	const TString& GetDBFolder() const {return fDBFolder;}
	const TString& GetSE() 	 const {return fSE;}
	const TString& GetCacheFolder() const {return fCacheFolder;}
	Bool_t  GetOperateDisconnected() const {return fOperateDisconnected;}
	Long64_t  GetCacheSize() const {return fCacheSize;}
	Long_t  GetCleanupInterval() const {return fCleanupInterval;}

	virtual AliCDBParam* CloneParam() const;

        virtual ULong_t Hash() const;
        virtual Bool_t IsEqual(const TObject* obj) const;

private:
	TString  fGridUrl;     // Grid url "Host:port"
	TString  fUser;	      // User
	TString  fDBFolder;    // path of the DB folder
	TString  fSE;	      // Storage Element
	TString  fCacheFolder; // Cache folder
	Bool_t   fOperateDisconnected; // Operate disconnected flag
	Long64_t fCacheSize;           // local cache size (in bytes)
	Long_t   fCleanupInterval;     // local cache cleanup interval

	ClassDef(AliCDBGridParam, 0);
};


#endif
