#ifndef ALI_CDB_MANAGER_H
#define ALI_CDB_MANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBManager                                            //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TList.h>
#include <TMap.h>

#include "AliCDBEntry.h"

class AliCDBStorage;
class AliCDBStorageFactory;
class AliCDBParam;

class AliCDBManager: public TObject {

 public:

	void RegisterFactory(AliCDBStorageFactory* factory);

	Bool_t HasStorage(const char* dbString);

	AliCDBParam* CreateParameter(const char* dbString);

	AliCDBStorage* GetStorage(const char* dbString);
	AliCDBStorage* GetStorage(const AliCDBParam* param);
	
	TList* GetActiveStorages();

	void SetDefaultStorage(const char* dbString);
	void SetDefaultStorage(const AliCDBParam* param);
	void SetDefaultStorage(AliCDBStorage *storage);

	Bool_t IsDefaultStorageSet() {return fDefaultStorage != 0;}
	
	AliCDBStorage* GetDefaultStorage() {return fDefaultStorage;}

	void RemoveDefaultStorage();

	void SetDrain(const char* dbString);
	void SetDrain(const AliCDBParam* param);
	void SetDrain(AliCDBStorage *storage);

	Bool_t IsDrainSet() {return fDrainStorage != 0;}

	Bool_t Drain(AliCDBEntry* entry);

	void RemoveDrain();

	void DestroyActiveStorages();
	void DestroyActiveStorage(AliCDBStorage* storage);
	
	static void Destroy();
	~AliCDBManager();

	static AliCDBManager* Instance();

 private:
		
	AliCDBManager();
	static AliCDBManager* fgInstance;
	
	AliCDBStorage* GetActiveStorage(const AliCDBParam* param);
	void PutActiveStorage(AliCDBParam* param, AliCDBStorage* storage);

	void Init();
	
	TList fFactories; 		// list of registered storage factories
	TMap fActiveStorages;		// list of active storages
	AliCDBStorage *fDefaultStorage;	// pointer to default storage
	AliCDBStorage *fDrainStorage;	// pointer to drain storage

	ClassDef(AliCDBManager, 0);
};


/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBStorageFactory                                     //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBParam;
class AliCDBStorageFactory: public TObject {
	friend class AliCDBManager;
	
public:
	virtual Bool_t Validate(const char* dbString) = 0;
	virtual AliCDBParam* CreateParameter(const char* dbString) = 0;	

protected:
	virtual AliCDBStorage* Create(const AliCDBParam* param) = 0;

	ClassDef(AliCDBStorageFactory, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBParam                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBParam: public TObject {

public:

	AliCDBParam();
	virtual ~AliCDBParam();

	const TString& GetType() const {return fType;};
	const TString& GetURI() const {return fURI;};

	virtual AliCDBParam* CloneParam() const = 0;

protected:

	void SetType(const char* type) {fType = type;};
	void SetURI(const char* uri) {fURI = uri;};

private:

	TString fType;
	TString fURI;

	ClassDef(AliCDBParam, 0);
};

#endif
