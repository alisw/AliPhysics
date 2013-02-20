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

#include <TKey.h>
#include <TH1.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TFile.h>
#include "AliCDBStorage.h"
#include "AliCDBGrid.h"

#include "AliCDBEntry.h"
#include "AliLog.h"

ClassImp(AliCDBStorage)

//_____________________________________________________________________________
AliCDBStorage::AliCDBStorage():
fValidFileIds(),
fRun(-1),
fPathFilter(),
fVersion(-1),
fMetaDataFilter(0),
fSelections(),
fURI(),
fType(),
fBaseFolder(),
fNretry(0),
fInitRetrySeconds(0)
{
// constructor

	fValidFileIds.SetOwner(1);
	fSelections.SetOwner(1);
}

//_____________________________________________________________________________
AliCDBStorage::~AliCDBStorage() {
// destructor

	RemoveAllSelections();
	fValidFileIds.Clear();
	delete fMetaDataFilter;

}

//_____________________________________________________________________________
void AliCDBStorage::GetSelection(/*const*/ AliCDBId* id) {
// return required version and subversion from the list of selection criteria
	
	TIter iter(&fSelections);
	AliCDBId* aSelection;
        
	// loop on the list of selection criteria
	while ((aSelection = (AliCDBId*) iter.Next())) {
		// check if selection element contains id's path and run (range) 
		if (aSelection->Comprises(*id)) {
			AliDebug(2,Form("Using selection criterion: %s ", aSelection->ToString().Data()));
			// return required version and subversion
			
			id->SetVersion(aSelection->GetVersion());
			id->SetSubVersion(aSelection->GetSubVersion());
			return;  
		}
	}
	
	// no valid element is found in the list of selection criteria -> return
	AliDebug(2,"Looking for objects with most recent version");
	return;
}

//_____________________________________________________________________________
void AliCDBStorage::ReadSelectionFromFile(const char *fileName){
// read selection criteria list from file
	
	RemoveAllSelections();
	
	TList *list = GetIdListFromFile(fileName);
	if(!list) return;
	
	list->SetOwner();	
	Int_t nId = list->GetEntries();
	AliCDBId *id;
	TKey *key;
	
	for(int i=nId-1;i>=0;i--){
		key = (TKey*) list->At(i);
		id = (AliCDBId*) key->ReadObj();
		if(id) AddSelection(*id);
	}
	delete list;
	AliInfo(Form("Selection criteria list filled with %d entries",fSelections.GetEntries()));
	PrintSelectionList();
	
}

//_____________________________________________________________________________
void AliCDBStorage::AddSelection(const AliCDBId& selection) {
// add a selection criterion

	AliCDBPath path = selection.GetPath();
	if(!path.IsValid()) return;
	
	TIter iter(&fSelections);
	const AliCDBId *anId;
	while((anId = (AliCDBId*) iter.Next())){
		if(selection.Comprises(*anId)){
			AliWarning("This selection is more general than a previous one and will hide it!");
			AliWarning(Form("%s", (anId->ToString()).Data()));
			fSelections.AddBefore(anId, new AliCDBId(selection));
			return;
		}
	
	}
	fSelections.AddFirst(new AliCDBId(selection));
}

//_____________________________________________________________________________
void AliCDBStorage::AddSelection(const AliCDBPath& path,
	const AliCDBRunRange& runRange, Int_t version, Int_t subVersion){
// add a selection criterion

	AddSelection(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
void AliCDBStorage::AddSelection(const AliCDBPath& path,
	Int_t firstRun, Int_t lastRun, Int_t version, Int_t subVersion){
// add a selection criterion
	
	AddSelection(AliCDBId(path, firstRun, lastRun, version, subVersion));
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveSelection(const AliCDBId& selection) {
// remove a selection criterion

	TIter iter(&fSelections);
	AliCDBId* aSelection;

	while ((aSelection = (AliCDBId*) iter.Next())) {
		if (selection.Comprises(*aSelection)) {
			fSelections.Remove(aSelection);
		}
	}
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveSelection(const AliCDBPath& path, 
	const AliCDBRunRange& runRange){
// remove a selection criterion

	RemoveSelection(AliCDBId(path, runRange, -1, -1));
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveSelection(const AliCDBPath& path,
	Int_t firstRun, Int_t lastRun){
// remove a selection criterion

	RemoveSelection(AliCDBId(path, firstRun, lastRun, -1, -1));
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveSelection(int position){
// remove a selection criterion from its position in the list

	delete fSelections.RemoveAt(position);
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveAllSelections(){
// remove all selection criteria

	fSelections.Clear();
}

//_____________________________________________________________________________
void AliCDBStorage::PrintSelectionList(){
// prints the list of selection criteria

	TIter iter(&fSelections);
	AliCDBId* aSelection;

	// loop on the list of selection criteria
	int index=0;
	while ((aSelection = (AliCDBId*) iter.Next())) {
		AliInfo(Form("index %d -> selection: %s",index++, aSelection->ToString().Data()));
	}

}

//_____________________________________________________________________________
AliCDBEntry* AliCDBStorage::Get(const AliCDBId& query) {
// get an AliCDBEntry object from the database

	// check if query's path and runRange are valid
	// query is invalid also if version is not specified and subversion is!
	if (!query.IsValid()) {
		AliError(Form("Invalid query: %s", query.ToString().Data()));
		return NULL;
	}

	// query is not specified if path contains wildcard or runrange = [-1,-1]
	if (!query.IsSpecified()) {
		AliError(Form("Unspecified query: %s",
				query.ToString().Data()));
                return NULL;
	}

	// This is needed otherwise TH1  objects (histos, TTree's) are lost when file is closed!
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);

	AliCDBEntry* entry = GetEntry(query);

	if (oldStatus != kFALSE)
  		TH1::AddDirectory(kTRUE);

	// if drain storage is set, drain entry into drain storage
	if(entry && (AliCDBManager::Instance())->IsDrainSet())
		AliCDBManager::Instance()->Drain(entry);

	return entry;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBStorage::Get(const AliCDBPath& path, Int_t runNumber,
	Int_t version, Int_t subVersion) {
// get an AliCDBEntry object from the database

	return Get(AliCDBId(path, runNumber, runNumber, version, subVersion));
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBStorage::Get(const AliCDBPath& path,
	const AliCDBRunRange& runRange, Int_t version,
	Int_t subVersion) {
// get an AliCDBEntry object from the database

	return Get(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
TList* AliCDBStorage::GetAll(const AliCDBId& query) {
// get multiple AliCDBEntry objects from the database


	if (!query.IsValid()) {
                AliError(Form("Invalid query: %s", query.ToString().Data()));
                return NULL;
        }

	if (query.IsAnyRange()) {
		AliError(Form("Unspecified run or runrange: %s",
				query.ToString().Data())); 	
		return NULL;
	}	
        
	// This is needed otherwise TH1  objects (histos, TTree's) are lost when file is closed!
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);

	TList *result = GetEntries(query);

	if (oldStatus != kFALSE)
  		TH1::AddDirectory(kTRUE);

 	Int_t nEntries = result->GetEntries();

 	AliInfo(Form("%d objects retrieved. Request was: %s",
		 	nEntries, query.ToString().Data()));
	for(int i=0; i<nEntries;i++){
		AliCDBEntry *entry = (AliCDBEntry*) result->At(i);
		AliInfo(Form("%s",entry->GetId().ToString().Data()));
	}


	// if drain storage is set, drain entries into drain storage
	if((AliCDBManager::Instance())->IsDrainSet()){
		for(int i = 0; i<result->GetEntries(); i++){
			AliCDBEntry* entry = (AliCDBEntry*) result->At(i);
			AliCDBManager::Instance()->Drain(entry);
		}
	}
	

        return result;
}

//_____________________________________________________________________________
TList* AliCDBStorage::GetAll(const AliCDBPath& path, Int_t runNumber, 
	Int_t version, Int_t subVersion) {
// get multiple AliCDBEntry objects from the database

	return GetAll(AliCDBId(path, runNumber, runNumber, version, 	
			subVersion));
}

//_____________________________________________________________________________
TList* AliCDBStorage::GetAll(const AliCDBPath& path, 
	const AliCDBRunRange& runRange, Int_t version, Int_t subVersion) {
// get multiple AliCDBEntry objects from the database

	return GetAll(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
AliCDBId* AliCDBStorage::GetId(const AliCDBId& query) {
// get the Id of the valid object from the database (does not open the file)

	// check if query's path and runRange are valid
	// query is invalid also if version is not specified and subversion is!
	if (!query.IsValid()) {
		AliError(Form("Invalid query: %s", query.ToString().Data()));
		return NULL;
	}

	// query is not specified if path contains wildcard or runrange = [-1,-1]
	if (!query.IsSpecified()) {
		AliError(Form("Unspecified query: %s",
				query.ToString().Data()));
                return NULL;
	}

	AliCDBId* id = GetEntryId(query);

	return id;
}

//_____________________________________________________________________________
AliCDBId* AliCDBStorage::GetId(const AliCDBPath& path, Int_t runNumber,
	Int_t version, Int_t subVersion) {
// get the Id of the valid object from the database (does not open the file)

	return GetId(AliCDBId(path, runNumber, runNumber, version, subVersion));
}

//_____________________________________________________________________________
AliCDBId* AliCDBStorage::GetId(const AliCDBPath& path,
	const AliCDBRunRange& runRange, Int_t version,
	Int_t subVersion) {
// get the Id of the valid object from the database (does not open the file)

	return GetId(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
Bool_t AliCDBStorage::Put(TObject* object, AliCDBId& id, AliCDBMetaData* metaData, const char* mirrors, AliCDBManager::DataType type) {
// store an AliCDBEntry object into the database
	
	if (object==0x0) {
		AliError("Null Entry! No storage will be done!");
		return kFALSE;
	} 

	AliCDBEntry anEntry(object, id, metaData);

	return Put(&anEntry, mirrors, type);
}

//_____________________________________________________________________________
Bool_t AliCDBStorage::Put(AliCDBEntry* entry, const char* mirrors, AliCDBManager::DataType type) {
// store an AliCDBEntry object into the database

	if (!entry){
		AliError("No entry!");
		return kFALSE;
	}
	
	if (entry->GetObject()==0x0){
		AliError("No valid object in CDB entry!");
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

	AliCDBManager::DataType expectedType = GetDataType();

	if(expectedType != AliCDBManager::kPrivate && type != expectedType) {
		AliError(Form("It is forbidden to store %s data into a folder of type %s!",
			AliCDBManager::GetDataTypeName(type),
			AliCDBManager::GetDataTypeName(expectedType)));
			return 0;
	}

	
	TString strMirrors(mirrors);
	if(!strMirrors.IsNull() && !strMirrors.IsWhitespace())
	    return PutEntry(entry, mirrors);
	else
	    return PutEntry(entry);
}

//_____________________________________________________________________________
void AliCDBStorage::QueryCDB(Int_t run, const char* pathFilter,
				Int_t version, AliCDBMetaData* md){
// query CDB for files valid for given run, and fill list fValidFileIds
// Actual query is done in virtual function QueryValidFiles()

	fRun = run;

	fPathFilter = pathFilter;
	if(!fPathFilter.IsValid()) {
		AliError(Form("Filter not valid: %s", pathFilter));
		fPathFilter = "*";
		return;
	}

	fVersion = version;

	AliInfo(Form("Querying files valid for run %d and path \"%s\" into CDB storage  \"%s://%s\"",
				fRun, pathFilter, fType.Data(), fBaseFolder.Data()));

	// Clear fValidFileIds list (it cannot be filled twice!
	AliDebug(2, "Clearing list of CDB Id's previously loaded");
	fValidFileIds.Clear();

	if(fMetaDataFilter) {delete fMetaDataFilter; fMetaDataFilter=0;}
	if(md) fMetaDataFilter = dynamic_cast<AliCDBMetaData*> (md->Clone());

	QueryValidFiles();

	AliInfo(Form("%d valid files found!", fValidFileIds.GetEntriesFast()));

}

//_____________________________________________________________________________
void AliCDBStorage::PrintQueryCDB(){
// print parameters used to load list of CDB Id's (fRun, fPathFilter, fVersion)

	AliCDBId paramId(fPathFilter, fRun, fRun, fVersion);
	AliInfo(Form("**** QueryCDB Parameters **** \n\t<%s>\n",
				paramId.ToString().Data()));

	if(fMetaDataFilter) fMetaDataFilter->PrintMetaData();


	TString message = "**** Id's of valid objects found *****\n";
	TIter iter(&fValidFileIds);
	AliCDBId* anId=0;

	// loop on the list of selection criteria
	while ((anId = dynamic_cast<AliCDBId*>(iter.Next()))) {
		message += Form("\t%s\n", anId->ToString().Data());
	}
	message += Form("\n\tTotal: %d objects found\n", fValidFileIds.GetEntriesFast());
	AliInfo(Form("%s", message.Data()));
}

//_____________________________________________________________________________
AliCDBManager::DataType AliCDBStorage::GetDataType() const {
// returns the type of the data that should be stored into this storage:
// kConditions: conditions data; kReference: reference data; kPrivate: private (user-defined) data type

	if(GetType() != "alien") return AliCDBManager::kPrivate;

	TString condFolder = ((AliCDBGridParam*) AliCDBManager::Instance()->GetCondParam())->GetDBFolder();
	TString refFolder = ((AliCDBGridParam*) AliCDBManager::Instance()->GetRefParam())->GetDBFolder();

	if(GetBaseFolder().Contains(condFolder)) return AliCDBManager::kCondition;
	if(GetBaseFolder().Contains(refFolder)) return AliCDBManager::kReference;

	return AliCDBManager::kPrivate;
}

//_____________________________________________________________________________
void AliCDBStorage::SetMirrorSEs(const char* mirrors) {
// if the current storage is not of "alien" type, just issue a warning
// AliCDBGrid implements its own SetMirrorSEs method, classes for other storage types do not

        TString storageType = GetType();
	if(storageType != "alien"){
	    AliWarning(Form("The current storage is of type \"%s\". Setting of SEs to \"%s\" skipped!",storageType.Data(),mirrors));
	    return;
	}
	AliError("We should never get here!! AliCDBGrid must have masked this virtual method!");
	return;
}

//_____________________________________________________________________________
const char* AliCDBStorage::GetMirrorSEs() const {
// if the current storage is not of "alien" type, just issue a warning
// AliCDBGrid implements its own GetMirrorSEs method, classes for other storage types do not

        TString storageType = GetType();
	if(storageType != "alien"){
	    AliWarning(Form("The current storage is of type \"%s\" and cannot handle SEs. Returning empty string!",storageType.Data()));
	    return "";
	}
	AliError("We should never get here!! AliCDBGrid must have masked this virtual method!");
	return "";
}

//_____________________________________________________________________________
void AliCDBStorage::LoadTreeFromFile(AliCDBEntry *entry) const {
// Checks whether entry contains a TTree and in case loads it into memory

	TObject *obj = (TObject*) entry->GetObject();
	if (!obj) {
	  AliError("Cannot retrieve the object:");
	  entry->PrintMetaData();
	  return;
	}

	if (!strcmp(obj->ClassName(),TTree::Class_Name())) {

		AliWarning("Entry contains a TTree! Loading baskets...");

		TTree* tree = dynamic_cast<TTree*> (obj);

		if(!tree) return;

		tree->LoadBaskets();
		tree->SetDirectory(0);
	}
	else if (!strcmp(obj->ClassName(),TNtuple::Class_Name())){

		AliWarning("Entry contains a TNtuple! Loading baskets...");

		TNtuple* ntu = dynamic_cast<TNtuple*> (obj);

		if(!ntu) return;

		ntu->LoadBaskets();
		ntu->SetDirectory(0);
	}

	return;
}

// //_____________________________________________________________________________
// void AliCDBStorage::SetTreeToFile(AliCDBEntry *entry, TFile* file) const {
// // Checks whether entry contains a TTree and in case assigns it to memory
// 
// 	AliCDBMetaData *md = dynamic_cast<AliCDBMetaData*> (entry->GetMetaData());
// 	if(!md) return;
// 	TString objStr = md->GetObjectClassName();
// 	if(objStr != "TTree") return;
// 	AliWarning("Entry contains a TTree! Setting file...");
// 
// 	TTree* tree = dynamic_cast<TTree*> (entry->GetObject());
// 
// 	if(!tree) return;
// 
// //	tree->SetDirectory(file);
// 	tree->SetDirectory(0);
// 
// 	return;
// }
