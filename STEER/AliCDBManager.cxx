/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliCDBDump.h"
#include "AliCDBLocal.h"
#include "AliCDBGrid.h"

#include <TObjString.h>
#include <TSystem.h>

ClassImp(AliCDBParam)

ClassImp(AliCDBManager)

AliCDBManager* AliCDBManager::fgInstance = 0x0;

//_____________________________________________________________________________
AliCDBManager* AliCDBManager::Instance() {
// returns AliCDBManager instance (singleton)

	if (!fgInstance) {
		fgInstance = new AliCDBManager();
		fgInstance->Init();
	}

	return fgInstance;
}

//_____________________________________________________________________________
void AliCDBManager::Init() {
// factory registering

	RegisterFactory(new AliCDBDumpFactory());
	RegisterFactory(new AliCDBLocalFactory()); 
	// AliCDBGridFactory is registered only if AliEn libraries are enabled in Root
	if(!gSystem->Exec("root-config --has-alien |grep yes 2>&1 > /dev/null")){ // returns 0 if yes
		AliInfo("AliEn classes enabled in Root. AliCDBGrid factory registered.");
		RegisterFactory(new AliCDBGridFactory());
	}
}
//_____________________________________________________________________________
void AliCDBManager::Destroy() {
// delete ALCDBManager instance and active storages

	if (fgInstance) {
		delete fgInstance;
		fgInstance = 0x0;
	}
}

//_____________________________________________________________________________
AliCDBManager::AliCDBManager():
	fDefaultStorage(NULL),
	fDrainStorage(NULL)
{
// default constuctor
	fFactories.SetOwner();
}

//_____________________________________________________________________________
AliCDBManager::~AliCDBManager() {
// destructor
	DestroyActiveStorages();
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetActiveStorage(const AliCDBParam* param) {
// get a storage object from the list of active storages 

        return (AliCDBStorage*) fActiveStorages.GetValue(param);
}

//_____________________________________________________________________________
void AliCDBManager::PutActiveStorage(AliCDBParam* param, AliCDBStorage* storage){
// put a storage object into the list of active storages

	fActiveStorages.Add(param, storage);
	AliDebug(1, Form("Active storages: %d", fActiveStorages.GetEntries()));
}

//_____________________________________________________________________________
void AliCDBManager::RegisterFactory(AliCDBStorageFactory* factory) {
// add a storage factory to the list of registerd factories
 
	if (!fFactories.Contains(factory)) {
		fFactories.Add(factory);
	}
}

//_____________________________________________________________________________
Bool_t AliCDBManager::HasStorage(const char* dbString) {
// check if dbString is a URI valid for one of the registered factories 

	TIter iter(&fFactories);

	AliCDBStorageFactory* factory;
	while ((factory = (AliCDBStorageFactory*) iter.Next())) {

		if (factory->Validate(dbString)) {
			return kTRUE;
		}	
	}

	return kFALSE;
}

//_____________________________________________________________________________
AliCDBParam* AliCDBManager::CreateParameter(const char* dbString) {
// create AliCDBParam object from URI string

	TIter iter(&fFactories);

        AliCDBStorageFactory* factory;
        while ((factory = (AliCDBStorageFactory*) iter.Next())) {

		AliCDBParam* param = factory->CreateParameter(dbString);
		if (param) {
			return param;
		}
        }

        return NULL;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetStorage(const char* dbString) {
// get storage object from URI string
	
	AliCDBParam* param = CreateParameter(dbString);
	if (!param) {
		return NULL;
	}	

	AliCDBStorage* aStorage = GetStorage(param);

	delete param;
	
	return aStorage;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetStorage(const AliCDBParam* param) {
// get storage object from AliCDBParam object

	// if the list of active storages already contains 
	// the requested storage, return it
	AliCDBStorage* aStorage = GetActiveStorage(param);
	if (aStorage) {
		return aStorage;
	}

	TIter iter(&fFactories);

        AliCDBStorageFactory* factory;
        
	// loop on the list of registered factories
	while ((factory = (AliCDBStorageFactory*) iter.Next())) {

		// each factory tries to create its storage from the parameter
		aStorage = factory->Create(param);
		if (aStorage) {
			PutActiveStorage(param->CloneParam(), aStorage);
			// if default storage is not set, set to this storage
			if(!fDefaultStorage){
				fDefaultStorage=aStorage;
				AliInfo(Form("Default storage set to: %s",(param->GetURI()).Data()));
			}
			return aStorage;
		}
        }

        return NULL;
}

//_____________________________________________________________________________
TList* AliCDBManager::GetActiveStorages() {
// return list of active storages

	TList* result = new TList();

	TIter iter(fActiveStorages.GetTable());	
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		result->Add(aPair->Value());
	}

	return result;
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(const char* dbString) {
// set drain storage from URI string

	fDrainStorage = GetStorage(dbString);	
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(const AliCDBParam* param) {
// set drain storage from AliCDBParam
	
	fDrainStorage = GetStorage(param);
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(AliCDBStorage* storage) {
// set drain storage from another active storage
	
	fDrainStorage = storage;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::Drain(AliCDBEntry *entry) {
// drain retrieved object to drain storage

	AliInfo("Draining into drain storage...");
	return fDrainStorage->Put(entry);
}

//_____________________________________________________________________________
void AliCDBManager::RemoveDrain() {
// remove drain storage

	fDrainStorage=0;
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const char* dbString) {
// sets default storage from URI string

	if(fDefaultStorage) fDefaultStorage = 0;
	fDefaultStorage = GetStorage(dbString);	
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const AliCDBParam* param) {
// set default storage from AliCDBParam object
	
	if(fDefaultStorage) fDefaultStorage = 0;
	fDrainStorage = GetStorage(param);
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(AliCDBStorage* storage) {
// set default storage from another active storage
	
	if(fDefaultStorage) fDefaultStorage = 0;
	fDefaultStorage = storage;
}

//_____________________________________________________________________________
void AliCDBManager::RemoveDefaultStorage() {
// remove default storage

	fDefaultStorage = 0;
}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorages() {
// delete list of active storages

	fActiveStorages.DeleteAll();
}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorage(AliCDBStorage* /*storage*/) {
// destroys active storage (not implemented)

}

///////////////////////////////////////////////////////////
// AliCDBManager Parameter class                         //
// interface to specific AliCDBParameter class           //
// (AliCDBGridParam, AliCDBLocalParam, AliCDBDumpParam)  //
///////////////////////////////////////////////////////////

AliCDBParam::AliCDBParam() {
// constructor

}

//_____________________________________________________________________________
AliCDBParam::~AliCDBParam() {
// destructor

}

