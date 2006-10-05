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
//-------------------------------------------------------------------------
//   Implementation of AliCDBManager and AliCDBParam classe
//   Author: Alberto Colla 
//   e-mail: Alberto.Colla@cern.ch
//-------------------------------------------------------------------------

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliCDBDump.h"
#include "AliCDBLocal.h"
#include "AliCDBGrid.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"

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
		//fgInstance->Delete();
		delete fgInstance;
		fgInstance = 0x0;
	}
}

//_____________________________________________________________________________
AliCDBManager::AliCDBManager():
  TObject(),
  fFactories(),
  fActiveStorages(),
  fSpecificStorages(),
  fDefaultStorage(NULL),
  fDrainStorage(NULL),
  fEntryCache(),
  fCache(kTRUE),
  fRun(-1)
{
// default constuctor
	fFactories.SetOwner(1);
	fActiveStorages.SetOwner(1);
	fSpecificStorages.SetOwner(1);
	fEntryCache.SetOwner(1);
}

//_____________________________________________________________________________
AliCDBManager::~AliCDBManager() {
// destructor
	ClearCache();
	DestroyActiveStorages();
	fFactories.Delete();
	fDrainStorage = 0x0;
	fDefaultStorage = 0x0;
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
Bool_t AliCDBManager::HasStorage(const char* dbString) const {
// check if dbString is a URI valid for one of the registered factories 

	TIter iter(&fFactories);

	AliCDBStorageFactory* factory=0;
	while ((factory = (AliCDBStorageFactory*) iter.Next())) {

		if (factory->Validate(dbString)) {
			return kTRUE;
		}	
	}

	return kFALSE;
}

//_____________________________________________________________________________
AliCDBParam* AliCDBManager::CreateParameter(const char* dbString) const {
// create AliCDBParam object from URI string

	TIter iter(&fFactories);

        AliCDBStorageFactory* factory=0;
        while ((factory = (AliCDBStorageFactory*) iter.Next())) {
		AliCDBParam* param = factory->CreateParameter(dbString);
		if(param) return param;
        }

        return NULL;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetStorage(const char* dbString) {
// get storage object from URI string
	
	AliCDBParam* param = CreateParameter(dbString);
	if (!param) {
		AliError(Form("Failed to activate requested storage! Check URI: %s", dbString));
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

        AliCDBStorageFactory* factory=0;

	// loop on the list of registered factories
	while ((factory = (AliCDBStorageFactory*) iter.Next())) {

		// each factory tries to create its storage from the parameter
		aStorage = factory->Create(param);
		if (aStorage) {
			PutActiveStorage(param->CloneParam(), aStorage);
			aStorage->SetURI(param->GetURI());
			if(fRun >= 0) {
				aStorage->QueryCDB(fRun);
			}
			return aStorage;
		}
        }

	AliError(Form("Failed to activate requested storage! Check URI: %s", param->GetURI().Data()));

        return NULL;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetActiveStorage(const AliCDBParam* param) {
// get a storage object from the list of active storages

        return dynamic_cast<AliCDBStorage*> (fActiveStorages.GetValue(param));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetActiveStorages() {
// return list of active storages
// user has responsibility to delete returned object

	TList* result = new TList();

	TIter iter(fActiveStorages.GetTable());
	TPair* aPair=0;
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
void AliCDBManager::SetDefaultStorage(const char* dbString) {
// sets default storage from URI string

	AliInfo(Form("Setting Default storage to: %s",dbString));
	fDefaultStorage = GetStorage(dbString);
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const AliCDBParam* param) {
// set default storage from AliCDBParam object
	
	fDrainStorage = GetStorage(param);
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(AliCDBStorage* storage) {
// set default storage from another active storage
	
	fDefaultStorage = storage;
}

//_____________________________________________________________________________
void AliCDBManager::SetSpecificStorage(const char* calibType, const char* dbString) {
// sets storage specific for detector or calibration type (works with AliCDBManager::Get(...))

	AliCDBParam *aPar = CreateParameter(dbString);
	if(!aPar) return;
	SetSpecificStorage(calibType, aPar);
	delete aPar;
}

//_____________________________________________________________________________
void AliCDBManager::SetSpecificStorage(const char* calibType, AliCDBParam* param) {
// sets storage specific for detector or calibration type (works with AliCDBManager::Get(...))
// Default storage should be defined prior to any specific storages, e.g.:
// AliCDBManager::instance()->SetDefaultStorage("alien://");
// AliCDBManager::instance()->SetSpecificStorage("TPC/*","local://DB_TPC");
// AliCDBManager::instance()->SetSpecificStorage("*/Align/*","local://DB_TPCAlign");
// calibType must be a valid CDB path! (3 level folder structure)

	if(!fDefaultStorage) {
		AliError("Please activate a default storage first!");
		return;
	}
	
	AliCDBPath aPath(calibType);
	if(!aPath.IsValid()){
		AliError(Form("Not a valid path: %s", calibType));
		return;
	}

	TObjString *objCalibType = new TObjString(aPath.GetPath());
	if(fSpecificStorages.Contains(objCalibType)){
		AliWarning(Form("Storage \"%s\" already activated! It will be replaced by the new one",
					calibType));
		fSpecificStorages.Remove(objCalibType);
	}
	GetStorage(param);
	fSpecificStorages.Add(objCalibType, param->CloneParam());
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetSpecificStorage(const char* calibType) {
// get storage specific for detector or calibration type 

	AliCDBParam *checkPar = (AliCDBParam*) fSpecificStorages.GetValue(calibType);
	if(!checkPar){
		AliError(Form("%s storage not found!",calibType));
		return NULL;
	} else {
		return GetStorage(checkPar);
	}
	
}

//_____________________________________________________________________________
AliCDBParam* AliCDBManager::SelectSpecificStorage(const TString& path) {
// select storage valid for path from the list of specific storages

	AliCDBPath aPath(path);
	if(!aPath.IsValid()){
		AliError(Form("Not a valid path: %s", path.Data()));
		return NULL;
	}

	TIter iter(&fSpecificStorages);
	TObjString *aCalibType=0;
	AliCDBParam* aPar=0;
	while((aCalibType = (TObjString*) iter.Next())){
		AliCDBPath calibTypePath(aCalibType->GetName());
		if(calibTypePath.Comprises(aPath)) {
		 	aPar = (AliCDBParam*) fSpecificStorages.GetValue(aCalibType);
			break;
		}
	}
	return aPar;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBPath& path, Int_t runNumber, 
	Int_t version, Int_t subVersion) {
// get an AliCDBEntry object from the database

	if(runNumber < 0){
		// RunNumber is not specified. Try with fRun
  		if (fRun < 0){
   	 		AliError("Run number neither specified in query nor set in AliCDBManager! Use AliCDBManager::SetRun.");
    			return NULL;
  		}
		runNumber = fRun;
	}

	return Get(AliCDBId(path, runNumber, runNumber, version, subVersion));
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBPath& path, 
	const AliCDBRunRange& runRange, Int_t version,
	Int_t subVersion) {
// get an AliCDBEntry object from the database!

	return Get(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBId& query) {	
// get an AliCDBEntry object from the database
	
	if(!fDefaultStorage) {
		AliError("No storage set!");
		return NULL;
	}

	// check if query's path and runRange are valid
	// query is invalid also if version is not specified and subversion is!
	if (!query.IsValid()) {
		AliError(Form("Invalid query: %s", query.ToString().Data()));
		return NULL;
	}
	
	// query is not specified if path contains wildcard or run range= [-1,-1]
 	if (!query.IsSpecified()) {
		AliError(Form("Unspecified query: %s", 
				query.ToString().Data()));
                return NULL;
	}

	if(fCache && query.GetFirstRun() != fRun)
		AliWarning("Run number explicitly set in query: CDB cache temporarily disabled!");

	
  	AliCDBEntry *entry=0;
	
  	// first look into map of cached objects
  	if(fCache && query.GetFirstRun() == fRun)
		entry = (AliCDBEntry*) fEntryCache.GetValue(query.GetPath());

  	if(entry) {
		AliInfo(Form("Object %s retrieved from cache !!",query.GetPath().Data()));
		return entry;
	}
  	
	// Entry is not in cache -> retrieve it from CDB and cache it!!   
	AliCDBStorage *aStorage=0;
	AliCDBParam *aPar=SelectSpecificStorage(query.GetPath());

	if(aPar) {
		aStorage=GetStorage(aPar);
		TString str = aPar->GetURI();
		AliDebug(2,Form("Looking into storage: %s",str.Data()));
		
	} else {
		aStorage=GetDefaultStorage();
		AliDebug(2,"Looking into default storage");	
	}
			
	entry = aStorage->Get(query);
  	if (!entry) return NULL;

 	if(fCache && (query.GetFirstRun() == fRun)){
		AliDebug(2,Form("Caching entry %s !",query.GetPath().Data()));
		CacheEntry(query.GetPath(), entry);
	}
  
  	return entry;
		
}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBPath& path, Int_t runNumber,
	Int_t version, Int_t subVersion) {
// get multiple AliCDBEntry objects from the database

	if(runNumber < 0){
		// RunNumber is not specified. Try with fRun
  		if (fRun < 0){
   	 		AliError("Run number neither specified in query nor set in AliCDBManager! Use AliCDBManager::SetRun.");
    			return NULL;
  		}
		runNumber = fRun;
	}

	return GetAll(AliCDBId(path, runNumber, runNumber, version, 	
			subVersion));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBPath& path, 
	const AliCDBRunRange& runRange, Int_t version, Int_t subVersion) {
// get multiple AliCDBEntry objects from the database

	return GetAll(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBId& query) {
// get multiple AliCDBEntry objects from the database
// Warning: this method works correctly only for queries of the type "Detector/*"
// 		and not for more specific queries e.g. "Detector/Calib/*" !

	if(!fDefaultStorage) {
		AliError("No storage set!");
		return NULL;
	}

	if (!query.IsValid()) {
                AliError(Form("Invalid query: %s", query.ToString().Data()));
                return NULL;
        }

	if((fSpecificStorages.GetEntries()>0) && query.GetPath().BeginsWith('*')){
                // if specific storages are active a query with "*" is ambiguous
		AliError("Query too generic in this context!");
                return NULL;		
	}

	if (query.IsAnyRange()) {
		AliError(Form("Unspecified run or runrange: %s",
				query.ToString().Data())); 	
		return NULL;
	}	
        
	TObjString objStrLev0(query.GetLevel0());
	//AliCDBParam *aPar = (AliCDBParam*) fSpecificStorages.GetValue(&objStrLev0);
	AliCDBParam *aPar=SelectSpecificStorage(query.GetPath());

	AliCDBStorage *aStorage;	
	if(aPar) {
		aStorage=GetStorage(aPar);
		TString str = aPar->GetURI();
		AliDebug(2,Form("Looking into storage: %s",str.Data()));
		
	} else {
		aStorage=GetDefaultStorage();
		AliDebug(2,"Looking into default storage");	
	}

	TList *result = aStorage->GetAll(query);

        return result;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::Put(TObject* object, AliCDBId& id,  AliCDBMetaData* metaData){
// store an AliCDBEntry object into the database

	AliCDBEntry anEntry(object, id, metaData);
	return Put(&anEntry);

}


//_____________________________________________________________________________
Bool_t AliCDBManager::Put(AliCDBEntry* entry){
// store an AliCDBEntry object into the database

	if(!fDefaultStorage) {
		AliError("No storage set!");
		return kFALSE;
	}

	if (!entry){
		AliError("No entry!");
		return kFALSE;
	}

	if (!entry->GetId().IsValid()) {
		AliError(Form("Invalid entry ID: %s", 
			entry->GetId().ToString().Data()));
		return kFALSE;
	}	

	if (!entry->GetId().IsSpecified()) {
		AliError(Form("Unspecified entry ID: %s", 
			entry->GetId().ToString().Data()));
		return kFALSE;
	}

	AliCDBId id = entry->GetId();
	AliCDBParam *aPar = SelectSpecificStorage(id.GetPath());

	AliCDBStorage *aStorage;
	
	if(aPar) {
		aStorage=GetStorage(aPar);
		TString str = aPar->GetURI();
		AliDebug(2,Form("Storing object into storage: %s",str.Data()));
		
	} else {
		aStorage=GetDefaultStorage();
		AliDebug(2,"Storing object into default storage");	
	}

	return aStorage->Put(entry);


}

//_____________________________________________________________________________
void AliCDBManager::CacheEntry(const char* path, AliCDBEntry* entry)
{
// cache AliCDBEntry. Cache is valid until run number is changed.

	AliDebug(2,Form("Filling cache with entry %s",path));
	fEntryCache.Add(new TObjString(path), entry);
	AliDebug(2,Form("Cache entries: %d",fEntryCache.GetEntries()));

}

//_____________________________________________________________________________
void AliCDBManager::Print(Option_t* /*option*/) const
{
// Print list of active storages and their URIs
	AliInfo(Form("Run number: %d\n",fRun));

	TString output;
	output = "Cache is ";
	if(!fCache) output += "NOT ";
	output += "ACTIVE\n";
	AliInfo(output.Data());

	if(fDefaultStorage) {
		AliInfo("*** Default Storage: ***");
		output = Form("%s\n",fDefaultStorage->GetURI().Data());
		AliInfo(output.Data());
	}
	if(fSpecificStorages.GetEntries()>0) {
		AliInfo("*** Specific Storages: ***");
		TIter iter(fSpecificStorages.GetTable());
		TPair *aPair=0;
		while((aPair = (TPair*) iter.Next())){
			output = Form("Key: %s - Storage: %s",
				((TObjString*) aPair->Key())->GetName(),
				((AliCDBParam*) aPair->Value())->GetURI().Data());
			AliInfo(output.Data());
		}
		printf("\n");
	}
	if(fDrainStorage) {
		AliInfo("*** Drain Storage: ***");
		output = Form("%s\n",fDrainStorage->GetURI().Data());
		AliInfo(output.Data());
	}
	AliInfo(Form("Total number of active storages: %d",fActiveStorages.GetEntries()));

}

//_____________________________________________________________________________
void AliCDBManager::SetRun(Int_t run)
{
// Sets current run number.
// When the run number changes the caching is cleared.
  
	if (fRun == run)
		return;
  
	fRun = run;
	ClearCache();
	QueryCDB();
}

//_____________________________________________________________________________
void AliCDBManager::ClearCache(){
// clear AliCDBEntry cache

	AliDebug(2,Form("Clearing cache!"));
	fEntryCache.DeleteAll();
	AliDebug(2,Form("Cache entries: %d",fEntryCache.GetEntries()));

}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorages() {
// delete list of active storages

	fActiveStorages.DeleteAll();
	fSpecificStorages.DeleteAll();
}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorage(AliCDBStorage* /*storage*/) {
// destroys active storage

/*
	TIter iter(fActiveStorages.GetTable());
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		if(storage == (AliCDBStorage*) aPair->Value())
			delete fActiveStorages.Remove(aPair->Key());
			storage->Delete(); storage=0x0;
	}
*/

}

//_____________________________________________________________________________
void AliCDBManager::QueryCDB() {
// query default and specific storages for files valid for fRun. Every storage loads the Ids into its list.

	if (fRun < 0){
		AliError("Run number not yet set! Use AliCDBManager::SetRun.");
	return;
	}
	if (!fDefaultStorage){
		AliError("Default storage is not set! Use AliCDBManager::SetDefaultStorage");
	return;
	}
	if(fDefaultStorage->GetType() == "alien"){
		fDefaultStorage->QueryCDB(fRun);
	} else {
		AliDebug(2,"Skipping query for valid files, it used only in grid...");
	}

	TIter iter(&fSpecificStorages);
	TObjString *aCalibType=0;
	AliCDBParam* aPar=0;
	while((aCalibType = dynamic_cast<TObjString*> (iter.Next()))){
		aPar = (AliCDBParam*) fSpecificStorages.GetValue(aCalibType);
		if(aPar) {
			AliDebug(2,Form("Querying specific storage %s",aCalibType->GetName()));
			AliCDBStorage *aStorage = GetStorage(aPar);
			if(aStorage->GetType() == "alien"){
				aStorage->QueryCDB(fRun,aCalibType->GetName());
			} else {
				AliDebug(2,
					"Skipping query for valid files, it is used only in grid...");
			}
		}
	}
}


///////////////////////////////////////////////////////////
// AliCDBManager Parameter class                         //
// interface to specific AliCDBParameter class           //
// (AliCDBGridParam, AliCDBLocalParam, AliCDBDumpParam)  //
///////////////////////////////////////////////////////////

AliCDBParam::AliCDBParam():
  fType(),
  fURI()
{
// constructor

}

//_____________________________________________________________________________
AliCDBParam::~AliCDBParam() {
// destructor

}

