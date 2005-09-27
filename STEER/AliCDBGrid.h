#ifndef ALICDBGRID_H
#define ALICDBGRID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBGrid						   //
//  access class to a DataBase in an AliEn storage 		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBStorage.h"
#include "AliCDBManager.h"

class AliCDBGrid: public AliCDBStorage {
	friend class AliCDBGridFactory;

public:
		  
	virtual Bool_t IsReadOnly() {return kFALSE;};
	virtual Bool_t HasSubVersion() {return kFALSE;};
  
protected:

	virtual AliCDBEntry*	GetEntry(const AliCDBId& queryId);
	virtual TList* 		GetEntries(const AliCDBId& queryId);
	virtual Bool_t          PutEntry(AliCDBEntry* entry);

private:
 
	AliCDBGrid(const char *host="aliendb4.cern.ch", 
		const Int_t port = 9000, 
		const char *user="colla", 
	        const char* dbPath = "/alice/cern.ch/user/c/colla/DB", 
		const char *SE="ALICE::CERN::Server");

	virtual ~AliCDBGrid();

	AliCDBGrid(const AliCDBGrid& db);
	AliCDBGrid& operator = (const AliCDBGrid& db);

	Bool_t FilenameToId(const char* filename, AliCDBRunRange& runRange, Int_t& version);
	Bool_t IdToFilename(const AliCDBRunRange& runRange, Int_t version, TString& filename);

	Bool_t PrepareId(AliCDBId& id);
	AliCDBId GetId(const AliCDBId& query);


	void GetEntriesForLevel0(const char* level0, const AliCDBId& query, TList* result);
	void GetEntriesForLevel1(const char* level0, const char* level1, 
				 const AliCDBId& query, TList* result);

	TString    fHost;	// Grid host
	Int_t      fPort;	// port
	TString    fUser;	// User
	TString    fDBPath;     // path of the DB folder
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
	AliCDBGridParam(const char* host, const Int_t port, const char* user, 
			const char* dbPath, const char* se);
	
	virtual ~AliCDBGridParam();

	const TString& GetHost() const {return fHost;};
	const Int_t&   GetPort() const {return fPort;};
	const TString& GetUser() const {return fUser;};
	const TString& GetPath() const {return fDBPath;};
	const TString& GetSE() 	 const {return fSE;};

	virtual AliCDBParam* CloneParam() const;

        virtual ULong_t Hash() const;
        virtual Bool_t IsEqual(const TObject* obj) const;

private:
	TString fHost;	     // Grid host
	Int_t 	fPort;	     // port
	TString fUser;	     // User
	TString fDBPath;     // path of the DB folder
	TString fSE;	     // Storage Element 

	ClassDef(AliCDBGridParam, 0);
};


#endif
