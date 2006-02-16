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

class AliCDBGrid: public AliCDBStorage {
	friend class AliCDBGridFactory;

public:
		  
	virtual Bool_t IsReadOnly() const {return kFALSE;};
	virtual Bool_t HasSubVersion() const {return kFALSE;};
	virtual Bool_t Contains(const char* path) const;
  
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

	Bool_t FilenameToId(const char* filename, AliCDBRunRange& runRange, Int_t& version);
	Bool_t IdToFilename(const AliCDBRunRange& runRange, Int_t version, TString& filename);

	Bool_t PrepareId(AliCDBId& id);
	Bool_t GetId(const AliCDBId& query, AliCDBId& result);


	void GetEntriesForLevel0(const char* level0, const AliCDBId& query, TList* result);
	void GetEntriesForLevel1(const char* level0, const char* level1, 
				 const AliCDBId& query, TList* result);

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
