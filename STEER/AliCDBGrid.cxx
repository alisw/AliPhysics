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

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid										       //
// access class to a DataBase in an AliEn storage  			                       //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include <TGrid.h>
#include <TGridResult.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>

#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBGrid.h"


ClassImp(AliCDBGrid)

//_____________________________________________________________________________
AliCDBGrid::AliCDBGrid(const char *gridUrl, const char *user, const char *dbFolder, const char *se) :
AliCDBStorage(),
fGridUrl(gridUrl),
fUser(user),
fDBFolder(dbFolder),
fSE(se)
{
// constructor //

	// if the same Grid is alreay active, skip connection
	if (!gGrid || fGridUrl != gGrid->GridUrl()  
	     || (( fUser != "" ) && ( fUser != gGrid->GetUser() )) ) {
   		// connection to the Grid
		AliInfo("Connection to the Grid...");
		if(gGrid){
			AliInfo(Form("gGrid = %x; fGridUrl = %s; gGrid->GridUrl() = %s",gGrid,fGridUrl.Data(), gGrid->GridUrl()));
			AliInfo(Form("fUser = %s; gGrid->GetUser() = %s",fUser.Data(), gGrid->GetUser()));
		}
		TGrid::Connect(fGridUrl.Data(),fUser.Data());
	}

	if(!gGrid) {
		AliError("Connection failed!");
		return;
	}

	TString initDir(gGrid->Pwd(0));
	if (fDBFolder[0] != '/') {
		fDBFolder.Prepend(initDir);
	}

	// check DBFolder: trying to cd to DBFolder; if it does not exist, create it
	if(!gGrid->Cd(fDBFolder.Data(),0)){
		AliDebug(2,Form("Creating new folder <%s> ...",fDBFolder.Data()));
		if(!gGrid->Mkdir(fDBFolder.Data(),"",0)){
			AliError(Form("Cannot create folder <%s> !",fDBFolder.Data())); 
		}
	} else {
		AliDebug(2,Form("Folder <%s> found",fDBFolder.Data()));
	}

	// removes any '/' at the end of path, then append one '/'
	while(fDBFolder.EndsWith("/")) fDBFolder.Remove(fDBFolder.Last('/')); 
	fDBFolder+="/";

	fType="alien";
	fBaseFolder = fDBFolder;

	// return to the initial directory
	gGrid->Cd(initDir.Data(),0);
}

//_____________________________________________________________________________
AliCDBGrid::~AliCDBGrid()
{
// destructor
	delete gGrid; gGrid=0;

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::FilenameToId(TString& filename, AliCDBId& id) {
// build AliCDBId from full path filename (fDBFolder/path/Run#x_#y_v#z.root)

	if(filename.Contains(fDBFolder)){
		filename = filename(fDBFolder.Length(),filename.Length()-fDBFolder.Length());
	}

	TString idPath = filename(0,filename.Last('/'));
	id.SetPath(idPath);
	if(!id.IsValid()) return kFALSE;

	filename=filename(idPath.Length()+1,filename.Length()-idPath.Length());

        Ssiz_t mSize;
	// valid filename: Run#firstRun_#lastRun_v#version.root
        TRegexp keyPattern("^Run[0-9]+_[0-9]+_v[0-9]+.root$");
        keyPattern.Index(filename, &mSize);
        if (!mSize) {
		AliDebug(2,Form("Bad filename <%s>.", filename.Data()));
                return kFALSE;
        }

	filename.Resize(filename.Length() - sizeof(".root") + 1);

        TObjArray* strArray = (TObjArray*) filename.Tokenize("_");

	TString firstRunString(((TObjString*) strArray->At(0))->GetString());
	id.SetFirstRun(atoi(firstRunString.Data() + 3));
	id.SetLastRun(atoi(((TObjString*) strArray->At(1))->GetString()));

	TString verString(((TObjString*) strArray->At(2))->GetString());
	id.SetVersion(atoi(verString.Data() + 1));

        delete strArray;

        return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBGrid::IdToFilename(const AliCDBId& id, TString& filename) const {
// build file name from AliCDBId (path, run range, version) and fDBFolder

	if (!id.GetAliCDBRunRange().IsValid()) {
		AliDebug(2,Form("Invalid run range <%d, %d>.",
			id.GetFirstRun(), id.GetLastRun()));
		return kFALSE;
	}

	if (id.GetVersion() < 0) {
		AliDebug(2,Form("Invalid version <%d>.", id.GetVersion()));
                return kFALSE;
	}

	filename = Form("Run%d_%d_v%d.root",
				id.GetFirstRun(),
				id.GetLastRun(),
				id.GetVersion());

	filename.Prepend(fDBFolder + id.GetPath() + '/');

        return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBGrid::PrepareId(AliCDBId& id) {
// prepare id (version) of the object that will be stored (called by PutEntry)

	TString initDir(gGrid->Pwd(0));
	TString pathName= id.GetPath();

	TString dirName(fDBFolder);

	Bool_t dirExist=kFALSE;
 
	// go to the path; if directory does not exist, create it
	TObjArray *arrName=pathName.Tokenize("/");
	for(int i=0;i<arrName->GetEntries();i++){
		TString buffer((arrName->At(i))->GetName());
 		dirName+=buffer; dirName+="/";
		dirExist=gGrid->Cd(dirName,0);
		if (!dirExist) {
			AliDebug(2,Form("Creating new folder <%s> ...",dirName.Data()));
			if(!gGrid->Mkdir(dirName,"",0)){
				AliError(Form("Cannot create directory <%s> !",dirName.Data()));
				gGrid->Cd(initDir.Data());
			return kFALSE;
			}
		}
	}
	delete arrName;
	gGrid->Cd(initDir,0);

	TString filename;
	AliCDBId anId; // the id got from filename
	AliCDBRunRange lastRunRange(-1,-1); // highest runRange found
	Int_t lastVersion=0; // highest version found

	TGridResult *res = gGrid->Ls(dirName);

	//loop on the files in the directory, look for highest version
	for(int i=0; i < res->GetEntries(); i++){
		filename=res->GetFileNamePath(i);
		if (!FilenameToId(filename, anId)) continue;
		if (anId.GetAliCDBRunRange().Overlaps(id.GetAliCDBRunRange()) && anId.GetVersion() > lastVersion) {
			lastVersion = anId.GetVersion();
			lastRunRange = anId.GetAliCDBRunRange();
		}

	}
	delete res; 
 
	id.SetVersion(lastVersion + 1);

	TString lastStorage = id.GetLastStorage();
	if(lastStorage.Contains(TString("new"), TString::kIgnoreCase) && id.GetVersion() > 1 ){
		AliDebug(2, Form("A NEW object is being stored with version %d",
					id.GetVersion()));
		AliDebug(2, Form("and it will hide previously stored object with version %d!",
					id.GetVersion()-1));
	}

	if(!lastRunRange.IsAnyRange() && !(lastRunRange.IsEqual(&id.GetAliCDBRunRange()))) 
    		AliWarning(Form("Run range modified w.r.t. previous version (Run%d_%d_v%d)",
    		     	lastRunRange.GetFirstRun(), lastRunRange.GetLastRun(), id.GetVersion()));
    
	return kTRUE;
}

//_____________________________________________________________________________
AliCDBId* AliCDBGrid::GetId(const TList& validFileIds, const AliCDBId& query) {
// look for the Id that matches query's requests (highest or exact version)

	if(validFileIds.GetEntries() < 1) {
		return NULL;
	} else if (validFileIds.GetEntries() == 1) {
		return dynamic_cast<AliCDBId*> (validFileIds.At(0));
	}

	TIter iter(&validFileIds);

	AliCDBId *anIdPtr=0;
	AliCDBId* result=0;

	while((anIdPtr = dynamic_cast<AliCDBId*> (iter.Next()))){
		if(anIdPtr->GetPath() != query.GetPath()) continue;

		//if(!CheckVersion(query, anIdPtr, result)) return NULL;

		if (!query.HasVersion()){ // look for highest version
			if(result && result->GetVersion() > anIdPtr->GetVersion()) continue;
			if(result && result->GetVersion() == anIdPtr->GetVersion()) {
				AliDebug(2,Form("More than one object valid for run %d, version %d!",
					query.GetFirstRun(), anIdPtr->GetVersion()));
				return NULL;
			}
			result = anIdPtr;
		} else { // look for specified version
			if(query.GetVersion() != anIdPtr->GetVersion()) continue;
			if(result && result->GetVersion() == anIdPtr->GetVersion()){
				AliDebug(2,Form("More than one object valid for run %d, version %d!",
					query.GetFirstRun(), anIdPtr->GetVersion()));
				return NULL;
			}
			result = anIdPtr;
		}

	}


	return result;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntry(const AliCDBId& queryId) {
// get AliCDBEntry from the database

	AliCDBId* dataId=0;

	AliCDBId selectedId(queryId);
	if (!selectedId.HasVersion()) {
		// if version is not specified, first check the selection criteria list
		GetSelection(&selectedId);
	}

	TList validFileIds;
	validFileIds.SetOwner(1);

	// look for file matching query requests (path, runRange, version)
	if(selectedId.GetFirstRun() == fRun &&
			fPathFilter.Comprises(selectedId.GetAliCDBPath()) && fVersion < 0 && !fMetaDataFilter){
		// look into list of valid files previously loaded with AliCDBStorage::FillValidFileIds()
		AliDebug(2, Form("List of files valid for run %d and for path %s was loaded. Looking there!",
					selectedId.GetFirstRun(), selectedId.GetPath().Data()));
		dataId = GetId(fValidFileIds, selectedId);

	} else {
		// List of files valid for reqested run was not loaded. Looking directly into CDB
		AliDebug(2, Form("List of files valid for run %d and for path %s was not loaded. Looking directly into CDB!",
					selectedId.GetFirstRun(), selectedId.GetPath().Data()));

		TString filter;
		MakeQueryFilter(selectedId.GetFirstRun(), selectedId.GetLastRun(),
				selectedId.GetAliCDBPath(), selectedId.GetVersion(), 0, filter);

		TGridResult *res = gGrid->Query(fDBFolder, "Run*.root", filter, "");
		AliCDBId validFileId;
		for(int i=0; i<res->GetEntries(); i++){
			TString filename = res->GetKey(i, "lfn");
			if(filename == "") continue;
			if(FilenameToId(filename, validFileId))
					validFileIds.AddLast(validFileId.Clone());
		}
		delete res;
		dataId = GetId(validFileIds, selectedId);
	}

	if (!dataId) return NULL;

	TString filename;
	if (!IdToFilename(*dataId, filename)) {
		AliDebug(2,Form("Bad data ID encountered! Subnormal error!"));
		return NULL;
	}

	AliCDBEntry* anEntry = GetEntryFromFile(filename, dataId);

	return anEntry;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntryFromFile(TString& filename, const AliCDBId* dataId){
// Get AliCBEntry object from file "filename"

	AliDebug(2,Form("Opening file: %s",filename.Data()));

	filename.Prepend("/alien");
	TFile *file = TFile::Open(filename);
	if (!file) {
		AliDebug(2,Form("Can't open file <%s>!", filename.Data()));
		return NULL;
	}

	// get the only AliCDBEntry object from the file
	// the object in the file is an AliCDBEntry entry named "AliCDBEntry"

	AliCDBEntry* anEntry = dynamic_cast<AliCDBEntry*> (file->Get("AliCDBEntry"));

	if (!anEntry) {
		AliDebug(2,Form("Bad storage data: file does not contain an AliCDBEntry object!"));
		file->Close();
		return NULL;
	}

	// The object's Id is not reset during storage
	// If object's Id runRange or version do not match with filename,
	// it means that someone renamed file by hand. In this case a warning msg is issued.

	if(anEntry){
		AliCDBId entryId = anEntry->GetId();
		if(!((entryId.GetAliCDBRunRange()).IsEqual(&(dataId->GetAliCDBRunRange()))) ||
			entryId.GetVersion() != dataId->GetVersion()){
			AliWarning(Form("Either RunRange or gridVersion in the object's metadata"));
			AliWarning(Form("do noth match with fileName numbers:"));
			AliWarning(Form("someone renamed file by hand!"));
		}
	}

	anEntry->SetLastStorage("grid");

	// close file, return retieved entry
	file->Close(); delete file; file=0;

	return anEntry;
}

//_____________________________________________________________________________
TList* AliCDBGrid::GetEntries(const AliCDBId& queryId) {
// multiple request (AliCDBStorage::GetAll)

	TList* result = new TList();
	result->SetOwner();

	TList validFileIds;
	validFileIds.SetOwner(1);

	Bool_t alreadyLoaded = kFALSE;

	// look for file matching query requests (path, runRange)
	if(queryId.GetFirstRun() == fRun &&
			fPathFilter.Comprises(queryId.GetAliCDBPath()) && fVersion < 0 && !fMetaDataFilter){
		// look into list of valid files previously loaded with AliCDBStorage::FillValidFileIds()
		AliDebug(2,Form("List of files valid for run %d and for path %s was loaded. Looking there!",
					queryId.GetFirstRun(), queryId.GetPath().Data()));

		alreadyLoaded = kTRUE;

	} else {
		// List of files valid for reqested run was not loaded. Looking directly into CDB
		AliDebug(2,Form("List of files valid for run %d and for path %s was not loaded. Looking directly into CDB!",
					queryId.GetFirstRun(), queryId.GetPath().Data()));

		TString filter;
		MakeQueryFilter(queryId.GetFirstRun(), queryId.GetLastRun(),
				queryId.GetAliCDBPath(), queryId.GetVersion(), 0, filter);

		TGridResult *res = gGrid->Query(fDBFolder, "Run*.root", filter, "");
		AliCDBId validFileId;
		for(int i=0; i<res->GetEntries(); i++){
			TString filename = res->GetKey(i, "lfn");
			if(filename == "") continue;
			if(FilenameToId(filename, validFileId))
					validFileIds.AddLast(validFileId.Clone());
		}
		delete res;
	}

	TIter *iter=0;
	if(alreadyLoaded){
		iter = new TIter(&fValidFileIds);
	} else {
		iter = new TIter(&validFileIds);
	}

	TList selectedIds;
	selectedIds.SetOwner(1);

	// loop on list of valid Ids to select the right version to get.
	// According to query and to the selection criteria list, version can be the highest or exact
	AliCDBPath pathCopy;
	AliCDBId* anIdPtr=0;
	AliCDBId* dataId=0;
	AliCDBPath queryPath = queryId.GetAliCDBPath();
	while((anIdPtr = dynamic_cast<AliCDBId*> (iter->Next()))){
		AliCDBPath thisCDBPath = anIdPtr->GetAliCDBPath();
		if(!(queryPath.Comprises(thisCDBPath)) || pathCopy.GetPath() == thisCDBPath.GetPath()) continue;
		pathCopy = thisCDBPath;

		// check the selection criteria list for this query
		AliCDBId thisId(*anIdPtr);
		thisId.SetVersion(queryId.GetVersion());
		if(!thisId.HasVersion()) GetSelection(&thisId);

		if(alreadyLoaded){
			dataId = GetId(fValidFileIds, thisId);
		} else {
			dataId = GetId(validFileIds, thisId);
		}
		if(dataId) selectedIds.Add(dataId->Clone());
	}

	delete iter; iter=0;

	// selectedIds contains the Ids of the files matching all requests of query!
	// All the objects are now ready to be retrieved
	iter = new TIter(&selectedIds);
	while((anIdPtr = dynamic_cast<AliCDBId*> (iter->Next()))){
		TString filename;
		if (!IdToFilename(*anIdPtr, filename)) {
			AliDebug(2,Form("Bad data ID encountered! Subnormal error!"));
			continue;
		}

		AliCDBEntry* anEntry = GetEntryFromFile(filename, anIdPtr);

		if(anEntry) result->Add(anEntry);

	}
	delete iter; iter=0;

	return result;
}

//_____________________________________________________________________________
Bool_t AliCDBGrid::PutEntry(AliCDBEntry* entry) {
// put an AliCDBEntry object into the database

	AliCDBId& id = entry->GetId();

	// set version for the entry to be stored
	if (!PrepareId(id)) return kFALSE;

	// build filename from entry's id
	TString filename;
	if (!IdToFilename(id, filename)) {
		AliError("Bad ID encountered! Subnormal error!");
		return kFALSE;
	}

	TString folderToTag = Form("%s%s",
					fDBFolder.Data(),
					id.GetPath().Data());

	// add CDB and CDB_MD tag to folder
	// TODO how to check that folder has already tags?
	AddTag(folderToTag,"CDB");
	AddTag(folderToTag,"CDB_MD");

	TDirectory* saveDir = gDirectory;

	// specify SE to filename
	TString fullFilename = Form("/alien%s?se=%s", filename.Data(), fSE.Data());

	// open file
	TFile *file = TFile::Open(fullFilename,"CREATE");
	if(!file || !file->IsWritable()){
		AliError(Form("Can't open file <%s>!", filename.Data()));
		if(file && !file->IsWritable()) file->Close(); delete file; file=0;
		return kFALSE;
	}

	file->cd();

	entry->SetVersion(id.GetVersion());

	// write object (key name: "AliCDBEntry")
	Bool_t result = (entry->Write("AliCDBEntry") != 0); 
	if (!result) AliError(Form("Can't write entry to file <%s>!", filename.Data()));


	if (saveDir) saveDir->cd(); else gROOT->cd();
	file->Close(); delete file; file=0;

	if(result) {
		AliInfo(Form("CDB object stored into file %s", filename.Data()));
		AliInfo(Form("using S.E. %s", fSE.Data()));

		if(!TagFileId(filename, &id)){
			AliInfo(Form("CDB tagging failed. Deleting file %s!",filename.Data()));
			if(!gGrid->Rm(filename.Data()))
				AliError("Can't delete file!");
			return kFALSE;
		}

		TagFileMetaData(filename, entry->GetMetaData());
	}

	return result;
}
//_____________________________________________________________________________
Bool_t AliCDBGrid::AddTag(TString& folderToTag, const char* tagname){
// add "tagname" tag (CDB or CDB_MD) to folder where object will be stored

	Bool_t result = kTRUE;
	AliDebug(2, Form("adding %s tag to folder %s", tagname, folderToTag.Data()));
	TString addTag = Form("addTag %s %s", folderToTag.Data(), tagname);
	TGridResult *gridres = gGrid->Command(addTag.Data());
	const char* resCode = gridres->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliError(Form("Couldn't add %s tags to folder %s !",
						tagname, folderToTag.Data()));
		result = kFALSE;
	}
	delete gridres;
	return result;
}

//_____________________________________________________________________________
Bool_t AliCDBGrid::TagFileId(TString& filename, const AliCDBId* id){
// tag stored object in CDB table using object Id's parameters

	TString addTagValue_1 = Form("addTagValue %s CDB ", filename.Data());
	TString addTagValue_2 = Form("first_run=%d last_run=%d version=%d ",
					id->GetFirstRun(),
					id->GetLastRun(),
					id->GetVersion());
	TString addTagValue_3 = Form("path_level_0=\"%s\" path_level_1=\"%s\" path_level_2=\"%s\"",
					id->GetLevel0().Data(),
					id->GetLevel1().Data(),
					id->GetLevel2().Data());
	TString addTagValue = Form("%s%s%s",
					addTagValue_1.Data(),
					addTagValue_2.Data(),
					addTagValue_3.Data());

	Bool_t result = kFALSE;
	AliDebug(2, Form("Tagging file. Tag command: %s", addTagValue.Data()));
	TGridResult* res = gGrid->Command(addTagValue.Data());
	const char* resCode = res->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliError(Form("Couldn't add CDB tag value to file %s !",
						filename.Data()));
		result = kFALSE;
	} else {
		AliInfo("Object successfully tagged.");
		result = kTRUE;
	}
	delete res;
	return result;

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::TagFileMetaData(TString& filename, const AliCDBMetaData* md){
// tag stored object in CDB table using object Id's parameters

	TString addTagValue_1 = Form("addTagValue %s CDB_MD ", filename.Data());
	TString addTagValue_2 = Form("object_classname=\"%s\" responsible=\"%s\" beam_period=%d ",
					md->GetObjectClassName(),
					md->GetResponsible(),
					md->GetBeamPeriod());
	TString addTagValue_3 = Form("aliroot_version=\"%s\" comment=\"%s\"",
					md->GetAliRootVersion(),
					md->GetComment());
	TString addTagValue = Form("%s%s%s",
					addTagValue_1.Data(),
					addTagValue_2.Data(),
					addTagValue_3.Data());

	Bool_t result = kFALSE;
	AliDebug(2, Form("Tagging file. Tag command: %s", addTagValue.Data()));
	TGridResult* res = gGrid->Command(addTagValue.Data());
	const char* resCode = res->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliWarning(Form("Couldn't add CDB_MD tag value to file %s !",
						filename.Data()));
		result = kFALSE;
	} else {
		AliInfo("Object successfully tagged.");
		result = kTRUE;
	}
	return result;
}

//_____________________________________________________________________________
TList* AliCDBGrid::GetIdListFromFile(const char* fileName){

	TString turl(fileName);
	turl.Prepend("/alien" + fDBFolder);
	turl += "?se="; turl += fSE.Data();
	TFile *file = TFile::Open(turl);
	if (!file) {
		AliError(Form("Can't open selection file <%s>!", turl.Data()));
		return NULL;
	}

	TList *list = new TList();
	list->SetOwner();
	int i=0;
	TString keycycle;

	AliCDBId *id;
	while(1){
		i++;
		keycycle = "AliCDBId;";
		keycycle+=i;
		
		id = (AliCDBId*) file->Get(keycycle);
		if(!id) break;
		list->AddFirst(id);
	}
	file->Close(); delete file; file=0;
	
	return list;


}

//_____________________________________________________________________________
Bool_t AliCDBGrid::Contains(const char* path) const{
// check for path in storage's DBFolder

	TString initDir(gGrid->Pwd(0));
	TString dirName(fDBFolder);
	dirName += path; // dirName = fDBFolder/path
	Bool_t result=kFALSE;
	if (gGrid->Cd(dirName,0)) result=kTRUE;
	gGrid->Cd(initDir.Data(),0);
	return result;
}

//_____________________________________________________________________________
void AliCDBGrid::QueryValidFiles()
{
// Query the CDB for files valid for AliCDBStorage::fRun
// fills list fValidFileIds with AliCDBId objects created from file name

	TString filter;
	MakeQueryFilter(fRun, fRun, fPathFilter, fVersion, fMetaDataFilter, filter);

	TGridResult *res = gGrid->Query(fDBFolder, "Run*.root", filter, "");
	AliCDBId validFileId;
	for(int i=0; i<res->GetEntries(); i++){
		TString filename = res->GetKey(i, "lfn");
		if(filename == "") continue;
		AliDebug(2,Form("Found valid file: %s", filename.Data()));
		Bool_t result = FilenameToId(filename, validFileId);
		if(result) {
			fValidFileIds.AddLast(validFileId.Clone());
		}
	}
	delete res;

}

//_____________________________________________________________________________
void AliCDBGrid::MakeQueryFilter(Int_t firstRun, Int_t lastRun,
					const AliCDBPath& pathFilter, Int_t version,
					const AliCDBMetaData* md, TString& result) const
{
// create filter for file query

	result = Form("CDB:first_run<=%d and CDB:last_run>=%d", firstRun, lastRun);

	if(version >= 0) {
		result += Form(" and CDB:version=%d", version);
	}
	if(pathFilter.GetLevel0() != "*") {
		result += Form(" and CDB:path_level_0=\"%s\"", pathFilter.GetLevel0().Data());
	}
	if(pathFilter.GetLevel1() != "*") {
		result += Form(" and CDB:path_level_1=\"%s\"", pathFilter.GetLevel1().Data());
	}
	if(pathFilter.GetLevel2() != "*") {
		result += Form(" and CDB:path_level_2=\"%s\"", pathFilter.GetLevel2().Data());
	}

	if(md){
		if(md->GetObjectClassName()[0] != '\0') {
			result += Form(" and CDB_MD:object_classname=\"%s\"", md->GetObjectClassName());
		}
		if(md->GetResponsible()[0] != '\0') {
			result += Form(" and CDB_MD:responsible=\"%s\"", md->GetResponsible());
		}
		if(md->GetBeamPeriod() != 0) {
			result += Form(" and CDB_MD:beam_period=%d", md->GetBeamPeriod());
		}
		if(md->GetAliRootVersion()[0] != '\0') {
			result += Form(" and CDB_MD:aliroot_version=\"%s\"", md->GetAliRootVersion());
		}
		if(md->GetComment()[0] != '\0') {
			result += Form(" and CDB_MD:comment=\"%s\"", md->GetComment());
		}
	}
	AliDebug(2, Form("filter: %s",result.Data()));

}

//_____________________________________________________________________________
Int_t AliCDBGrid::GetLatestVersion(const char* path, Int_t run){
// get last version found in the database valid for run and path

	TList validFileIds;
	validFileIds.SetOwner(1);

	AliCDBPath aCDBPath(path);
	if(!aCDBPath.IsValid() || aCDBPath.IsWildcard()) {
		AliError(Form("Invalid path in request: %s", path));
		return -1;
	}
	AliCDBId query(path, run, run, -1, -1);
	AliCDBId* dataId = 0;

	// look for file matching query requests (path, runRange, version)
	if(run == fRun &&
			fPathFilter.Comprises(aCDBPath) && fVersion < 0){
		// look into list of valid files previously loaded with AliCDBStorage::FillValidFileIds()
		AliDebug(2, Form("List of files valid for run %d and for path %s was loaded. Looking there!",
					run, path));
		dataId = GetId(fValidFileIds, query);
		if (!dataId) return -1;
		return dataId->GetVersion();

	}
	// List of files valid for reqested run was not loaded. Looking directly into CDB
	AliDebug(2, Form("List of files valid for run %d and for path %s was not loaded. Looking directly into CDB!",
				run, path));

	TString filter;
	MakeQueryFilter(run, run, aCDBPath, -1, 0, filter);

	TGridResult *res = gGrid->Query(fDBFolder, "Run*.root", filter, "");
	AliCDBId validFileId;
	for(int i=0; i<res->GetEntries(); i++){
		TString filename = res->GetKey(i, "lfn");
		if(filename == "") continue;
		if(FilenameToId(filename, validFileId))
				validFileIds.AddLast(validFileId.Clone());
	}
	delete res;

	dataId = GetId(validFileIds, query);
	if (!dataId) return -1;

	return dataId->GetVersion();

}

//_____________________________________________________________________________
Int_t AliCDBGrid::GetLatestSubVersion(const char* /*path*/, Int_t /*run*/, Int_t /*version*/){
// get last subversion found in the database valid for run and path
	AliError("Objects in GRID storage have no sub version!");
	return -1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid factory  			                                                       //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBGridFactory)

//_____________________________________________________________________________
Bool_t AliCDBGridFactory::Validate(const char* gridString) {
// check if the string is valid Grid URI

	// pattern: alien://hostName:Port;user;dbPath;SE
	// example of a valid pattern:
	// "alien://aliendb4.cern.ch:9000;colla;DBTest;ALICE::CERN::Server"
//        TRegexp gridPattern("^alien://.+:[0-9]+;[a-zA-Z0-9_-.]+;.+;.+$");
        TRegexp gridPattern("^alien://.+$");

        return TString(gridString).Contains(gridPattern);
}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridFactory::CreateParameter(const char* gridString) {
// create AliCDBGridParam class from the URI string

	if (!Validate(gridString)) {
		return NULL;
	}
	//TString buffer(gridString + sizeof("alien://") - 1);
	TString buffer(gridString);

 	TString gridUrl 	= "alien://";
	TString user 		= "";
	TString dbFolder 	= "DBGrid";
	TString se		= "ALICE::CERN::se01";

	TObjArray *arr = buffer.Tokenize('?');
	TIter iter(arr);
	TObjString *str = 0;
	
	while((str = (TObjString*) iter.Next())){
		TString entry(str->String());
		Int_t indeq = entry.Index('=');
		if(indeq == -1) {
			if(entry.BeginsWith("alien://")) { // maybe it's a gridUrl!
				gridUrl = entry;
				continue;
			} else {
				AliError(Form("Invalid entry! %s",entry.Data()));
				continue;
			}
		}
		
		TString key = entry(0,indeq);
		TString value = entry(indeq+1,entry.Length()-indeq);

		if(key.Contains("grid",TString::kIgnoreCase)) {
			gridUrl += value;
		} 
		else if (key.Contains("user",TString::kIgnoreCase)){
			user = value;
		}
		else if (key.Contains("folder",TString::kIgnoreCase)){
			dbFolder = value;
		}
		else if (key.Contains("se",TString::kIgnoreCase)){
			se = value;
		}
		else{
			AliError(Form("Invalid entry! %s",entry.Data()));
		}
	}
	delete arr; arr=0;
		
	AliDebug(2, Form("gridUrl:	%s",gridUrl.Data()));
	AliDebug(2, Form("user:	%s",user.Data()));
	AliDebug(2, Form("dbFolder:	%s",dbFolder.Data()));
	AliDebug(2, Form("s.e.:	%s",se.Data()));

	return new AliCDBGridParam(gridUrl.Data(), user.Data(), dbFolder.Data(), se.Data());
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBGridFactory::Create(const AliCDBParam* param) {
// create AliCDBGrid storage instance from parameters
	
	AliCDBGrid *grid = 0;
	if (AliCDBGridParam::Class() == param->IsA()) {

		const AliCDBGridParam* gridParam = (const AliCDBGridParam*) param;
		grid = new AliCDBGrid(gridParam->GridUrl().Data(),
				      gridParam->GetUser().Data(),
				      gridParam->GetDBFolder().Data(),
				      gridParam->GetSE().Data());

	}

	if(!gGrid && grid) {
		delete grid; grid=0;
	}

	return grid;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid Parameter class  			                                               //                                          //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBGridParam)

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam():
 AliCDBParam(),
 fGridUrl(),
 fUser(),
 fDBFolder(),
 fSE()
 {
// default constructor

}

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam(const char* gridUrl, 
				const char* user,
			        const char* dbFolder, 
				const char* se):
 AliCDBParam(),
 fGridUrl(gridUrl),
 fUser(user),
 fDBFolder(dbFolder),
 fSE(se)
{
// constructor
	
	SetType("alien");

	TString uri = Form("%s?User=%s?DBFolder=%s?SE=%s",
			fGridUrl.Data(), fUser.Data(),
			fDBFolder.Data(), fSE.Data());

	SetURI(uri.Data());
}

//_____________________________________________________________________________
AliCDBGridParam::~AliCDBGridParam() {
// destructor

}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridParam::CloneParam() const {
// clone parameter

        return new AliCDBGridParam(fGridUrl.Data(), fUser.Data(),
					fDBFolder.Data(), fSE.Data());
}

//_____________________________________________________________________________
ULong_t AliCDBGridParam::Hash() const {
// return Hash function

        return fGridUrl.Hash()+fUser.Hash()+fDBFolder.Hash()+fSE.Hash();
}

//_____________________________________________________________________________
Bool_t AliCDBGridParam::IsEqual(const TObject* obj) const {
// check if this object is equal to AliCDBParam obj

        if (this == obj) {
                return kTRUE;
        }

        if (AliCDBGridParam::Class() != obj->IsA()) {
                return kFALSE;
        }

        AliCDBGridParam* other = (AliCDBGridParam*) obj;

        if(fGridUrl != other->fGridUrl) return kFALSE;
        if(fUser != other->fUser) return kFALSE;
        if(fDBFolder != other->fDBFolder) return kFALSE;
        if(fSE != other->fSE) return kFALSE;
	return kTRUE;
}

