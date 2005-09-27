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

#include "AliCDBEntry.h"
#include "AliLog.h"

ClassImp(AliCDBStorage)

//_____________________________________________________________________________
AliCDBStorage::AliCDBStorage() {
// constructor

	fSelections.SetOwner(1);
}

//_____________________________________________________________________________
AliCDBStorage::~AliCDBStorage() {
// destructor

}

//_____________________________________________________________________________
AliCDBId AliCDBStorage::GetSelection(const AliCDBId& id) {
// return required version and subversion from the list of selection criteria 
	
	TIter iter(&fSelections);
	AliCDBId* aSelection;
        
	// loop on the list of selection criteria
	while ((aSelection = (AliCDBId*) iter.Next())) {
		// check if selection element contains id's path and run (range) 
		if (aSelection->Comprises(id)) {
			AliInfo(Form("Using selection criterion: %s ", aSelection->ToString().Data()));
			// return required version and subversion
			return AliCDBId(id.GetAliCDBPath(), 
				id.GetAliCDBRunRange(),
				aSelection->GetVersion(), 
				aSelection->GetSubVersion());  
		}
	}
	
	// no valid element is found in the list of selection criteria -> return
	AliInfo("No matching selection criteria: highest version will be seeked!");
	return AliCDBId(id.GetAliCDBPath(), id.GetAliCDBRunRange());
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
void AliCDBStorage::RemoveSelection(const int position){
// remove a selection criterion from its position in the list

	fSelections.RemoveAt(position);
}

//_____________________________________________________________________________
void AliCDBStorage::RemoveAllSelections(){
// remove all selection criteria

	fSelections.RemoveAll();
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

	AliCDBEntry* entry = GetEntry(query);
		
  	if (entry) {
    		AliInfo(Form("Valid AliCDBEntry object found! %s", entry->GetId().ToString().Data()));
  	} else {
    		AliInfo(Form("Sorry, found no object valid for: name = <%s>, run = %d", 
		        (query.GetPath()).Data(), query.GetFirstRun()));
  	}
	
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
        
	TList *result = GetEntries(query);

 	Int_t nEntries = result->GetEntries();
 	if (nEntries) {
 		 AliInfo(Form("%d AliCDBEntry objects found!",nEntries));
		 for(int i=0; i<nEntries;i++){
		 	AliCDBEntry *entry = (AliCDBEntry*) result->At(i);
			AliInfo(Form("%s",entry->GetId().ToString().Data()));
		 
		 }
 	} else {
     		 AliInfo(Form("Sorry, found no object valid for: name = <%s>, run = %d", 
		        (query.GetPath()).Data(), query.GetFirstRun()));
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
Bool_t AliCDBStorage::Put(TObject* object, AliCDBId& id, AliCDBMetaData* metaData) {
// put an AliCDBEntry object from the database
	
	AliCDBEntry anEntry(object, id, metaData);

	return Put(&anEntry);
} 

//_____________________________________________________________________________
Bool_t AliCDBStorage::Put(AliCDBEntry* entry) {
// put an AliCDBEntry object from the database

	if (!entry->GetId().IsValid()) {
		AliWarning(Form("Invalid entry ID: %s", 
			entry->GetId().ToString().Data()));
		return kFALSE;
	}	

	if (!entry->GetId().IsSpecified()) {
		AliError(Form("Unspecified entry ID: %s", 
			entry->GetId().ToString().Data()));
		return kFALSE;
	}

	return PutEntry(entry);
}

