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

#include <cstdlib>
#include <TGrid.h>
#include <TGridResult.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>
#include <TRegexp.h>

#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBGrid.h"
#include "AliCDBManager.h"


ClassImp(AliCDBGrid)

//_____________________________________________________________________________
AliCDBGrid::AliCDBGrid(const char *gridUrl, const char *user, const char *dbFolder,
                       const char *se, const char* cacheFolder, Bool_t operateDisconnected,
		       Long64_t cacheSize, Long_t cleanupInterval) :
AliCDBStorage(),
fGridUrl(gridUrl),
fUser(user),
fDBFolder(dbFolder),
fSE(se),
fMirrorSEs(""),
fCacheFolder(cacheFolder),
fOperateDisconnected(operateDisconnected),
fCacheSize(cacheSize),
fCleanupInterval(cleanupInterval)
{
// constructor //

	// if the same Grid is alreay active, skip connection
	if (!gGrid || fGridUrl != gGrid->GridUrl()
	     || (( fUser != "" ) && ( fUser != gGrid->GetUser() )) ) {
   		// connection to the Grid
		AliInfo("Connection to the Grid...");
		if(gGrid){
			AliInfo(Form("gGrid = %p; fGridUrl = %s; gGrid->GridUrl() = %s",gGrid,fGridUrl.Data(), gGrid->GridUrl()));
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
		TGridResult* res = gGrid->Command(Form("mkdir -p %s",fDBFolder.Data()));
		TString result = res->GetKey(0,"__result__");
		if(result == "0"){
			AliFatal(Form("Cannot create folder <%s> !",fDBFolder.Data()));
			return;
		}
	} else {
		AliDebug(2,Form("Folder <%s> found",fDBFolder.Data()));
	}

	// removes any '/' at the end of path, then append one '/'
	while(fDBFolder.EndsWith("/")) fDBFolder.Remove(fDBFolder.Last('/')); 
	fDBFolder+="/";

	fType="alien";
	fBaseFolder = fDBFolder;

	// Setting the cache

	// Check if local cache folder is already defined
	TString origCache(TFile::GetCacheFileDir());
	if(fCacheFolder.Length() > 0) {
		if(origCache.Length() == 0) {
			AliInfo(Form("Setting local cache to: %s", fCacheFolder.Data()));
		} else if(fCacheFolder != origCache) {
			AliWarning(Form("Local cache folder was already defined, changing it to: %s",
					fCacheFolder.Data()));
		}

		// default settings are: operateDisconnected=kTRUE, forceCacheread = kFALSE
		if(!TFile::SetCacheFileDir(fCacheFolder.Data(), fOperateDisconnected)) {
			AliError(Form("Could not set cache folder %s !", fCacheFolder.Data()));
			fCacheFolder = "";
		} else {
			// reset fCacheFolder because the function may have
			// slightly changed the folder name (e.g. '/' added)
			fCacheFolder = TFile::GetCacheFileDir();
		}

		// default settings are: cacheSize=1GB, cleanupInterval = 0
		if(!TFile::ShrinkCacheFileDir(fCacheSize, fCleanupInterval)) {
			AliError(Form("Could not set following values "
				"to ShrinkCacheFileDir: cacheSize = %lld, cleanupInterval = %ld !",
				fCacheSize, fCleanupInterval));
		}
	}

	// return to the initial directory
	gGrid->Cd(initDir.Data(),0);

	fNretry = 3;  // default
	fInitRetrySeconds = 5;   // default
}

//_____________________________________________________________________________
AliCDBGrid::~AliCDBGrid()
{
// destructor
	delete gGrid; gGrid=0;

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::FilenameToId(TString& filename, AliCDBId& id) {
// build AliCDBId from full path filename (fDBFolder/path/Run#x_#y_v#z_s0.root)

	if(filename.Contains(fDBFolder)){
		filename = filename(fDBFolder.Length(),filename.Length()-fDBFolder.Length());
	}

	TString idPath = filename(0,filename.Last('/'));
	id.SetPath(idPath);
	if(!id.IsValid()) return kFALSE;

	filename=filename(idPath.Length()+1,filename.Length()-idPath.Length());

        Ssiz_t mSize;
	// valid filename: Run#firstRun_#lastRun_v#version_s0.root
        TRegexp keyPattern("^Run[0-9]+_[0-9]+_v[0-9]+_s0.root$");
        keyPattern.Index(filename, &mSize);
        if (!mSize) {

		// TODO backward compatibility ... maybe remove later!
        	Ssiz_t oldmSize;
        	TRegexp oldKeyPattern("^Run[0-9]+_[0-9]+_v[0-9]+.root$");
        	oldKeyPattern.Index(filename, &oldmSize);
		if(!oldmSize) {
			AliDebug(2,Form("Bad filename <%s>.", filename.Data()));
                	return kFALSE;
		} else {
			AliDebug(2,Form("Old filename format <%s>.", filename.Data()));
			id.SetSubVersion(-11); // TODO trick to ensure backward compatibility
		}

        } else {
		id.SetSubVersion(-1); // TODO trick to ensure backward compatibility
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

	filename = Form("Run%d_%d_v%d",
				id.GetFirstRun(),
				id.GetLastRun(),
				id.GetVersion());

	if (id.GetSubVersion() != -11) filename += "_s0"; // TODO to ensure backward compatibility
	filename += ".root";

	filename.Prepend(fDBFolder + id.GetPath() + '/');

        return kTRUE;
}

//_____________________________________________________________________________
void AliCDBGrid::SetRetry(Int_t nretry, Int_t initsec) {

	// Function to set the exponential retry for putting entries in the OCDB

	AliWarning("WARNING!!! You are changing the exponential retry times and delay: this function should be used by experts!"); 
	fNretry = nretry;
	fInitRetrySeconds = initsec;
	AliDebug(2,Form("fNretry = %d, fInitRetrySeconds = %d", fNretry, fInitRetrySeconds));
} 


//_____________________________________________________________________________
Bool_t AliCDBGrid::PrepareId(AliCDBId& id) {
// prepare id (version) of the object that will be stored (called by PutEntry)

	TString initDir(gGrid->Pwd(0));

	TString dirName(fDBFolder);

	Bool_t dirExist=kFALSE;



	// go to the path; if directory does not exist, create it
	for(int i=0;i<3;i++){
	        //TString cmd("find -d ");
		//cmd += Form("%s ",dirName);
		//cmd += 
		//gGrid->Command(cmd.Data());
 		dirName+=Form("%s/",id.GetPathLevel(i).Data());
		dirExist=gGrid->Cd(dirName,0);
		if (!dirExist) {
			AliDebug(2,Form("Creating new folder <%s> ...",dirName.Data()));
			if(!gGrid->Mkdir(dirName,"",0)){
				AliError(Form("Cannot create directory <%s> !",dirName.Data()));
				gGrid->Cd(initDir.Data());
			return kFALSE;
			}

			// if folders are new add tags to them
			if(i == 1) {
				// TODO Currently disabled
				// add short lived tag!
				// AliInfo("Tagging level 1 folder with \"ShortLived\" tag");
				// if(!AddTag(dirName,"ShortLived_try")){
				//	AliError(Form("Could not tag folder %s !", dirName.Data()));
				//	if(!gGrid->Rmdir(dirName.Data())){
				//		AliError(Form("Unexpected: could not remove %s directory!", dirName.Data()));
				//	}
				//	return 0;
				//}

			} else if(i == 2) {
				AliDebug(2,"Tagging level 2 folder with \"CDB\" and \"CDB_MD\" tag");
				if(!AddTag(dirName,"CDB")){
					AliError(Form("Could not tag folder %s !", dirName.Data()));
					if(!gGrid->Rmdir(dirName.Data())){
						AliError(Form("Unexpected: could not remove %s directory!", dirName.Data()));
					}
					return 0;
				}
				if(!AddTag(dirName,"CDB_MD")){
					AliError(Form("Could not tag folder %s !", dirName.Data()));
					if(!gGrid->Rmdir(dirName.Data())){
						AliError(Form("Unexpected: could not remove %s directory!", dirName.Data()));
					}
					return 0;
				}

				// TODO Currently disabled
				// add short lived tag!
				// TString path=id.GetPath();
				// if(AliCDBManager::Instance()->IsShortLived(path.Data())) {
				//	AliInfo(Form("Tagging %s as short lived", dirName.Data()));
				//	if(!TagShortLived(dirName, kTRUE)){
				//		AliError(Form("Could not tag folder %s !", dirName.Data()));
				//		if(!gGrid->Rmdir(dirName.Data())){
				//			AliError(Form("Unexpected: could not remove %s directory!", dirName.Data()));
				//		}
				//		return 0;
				//	}
				// } else {
				//	AliInfo(Form("Tagging %s as long lived", dirName.Data()));
				//	if(!TagShortLived(dirName, kFALSE)){
				//		AliError(Form("Could not tag folder %s !", dirName.Data()));
				//		if(!gGrid->Rmdir(dirName.Data())){
				//			AliError(Form("Unexpected: could not remove %s directory!", dirName.Data()));
				//		}
				//		return 0;
				//	}
				// }
			}
		}
	}
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

	// GRP entries with explicitly set version escape default incremental versioning
	if(id.GetPath().Contains("GRP") && id.HasVersion() && lastVersion!=0)
	{
		AliDebug(5,Form("Entry %s won't be put in the destination OCDB", id.ToString().Data()));
		return kFALSE;
	}

	id.SetVersion(lastVersion + 1);
	id.SetSubVersion(0);

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
AliCDBId* AliCDBGrid::GetId(const TObjArray& validFileIds, const AliCDBId& query) {
// look for the Id that matches query's requests (highest or exact version)

	if(validFileIds.GetEntriesFast() < 1)
		return NULL;

	TIter iter(&validFileIds);

	AliCDBId *anIdPtr=0;
	AliCDBId* result=0;

	while((anIdPtr = dynamic_cast<AliCDBId*> (iter.Next()))){
		if(anIdPtr->GetPath() != query.GetPath()) continue;

		//if(!CheckVersion(query, anIdPtr, result)) return NULL;

		if (!query.HasVersion()){ // look for highest version
			if(result && result->GetVersion() > anIdPtr->GetVersion()) continue;
			if(result && result->GetVersion() == anIdPtr->GetVersion()) {
				AliError(Form("More than one object valid for run %d, version %d!",
					query.GetFirstRun(), anIdPtr->GetVersion()));
				return NULL;
			}
			result = anIdPtr;
		} else { // look for specified version
			if(query.GetVersion() != anIdPtr->GetVersion()) continue;
			if(result && result->GetVersion() == anIdPtr->GetVersion()){
				AliError(Form("More than one object valid for run %d, version %d!",
					query.GetFirstRun(), anIdPtr->GetVersion()));
				return NULL;
			}
			result = anIdPtr;
		}

	}
	
	if (!result) return NULL;

	return dynamic_cast<AliCDBId*> (result->Clone());
}

//_____________________________________________________________________________
AliCDBId* AliCDBGrid::GetEntryId(const AliCDBId& queryId) {
// get AliCDBId from the database
// User must delete returned object

	AliCDBId* dataId=0;

	AliCDBId selectedId(queryId);
	if (!selectedId.HasVersion()) {
		// if version is not specified, first check the selection criteria list
		GetSelection(&selectedId);
	}

	TObjArray validFileIds;
	validFileIds.SetOwner(1);

	// look for file matching query requests (path, runRange, version)
	if(selectedId.GetFirstRun() == fRun && fPathFilter.Comprises(selectedId.GetAliCDBPath()) &&
			fVersion == selectedId.GetVersion() && !fMetaDataFilter){
		// look into list of valid files previously loaded with AliCDBStorage::FillValidFileIds()
		AliDebug(2, Form("List of files valid for run %d was loaded. Looking there for fileids valid for path %s!",
					selectedId.GetFirstRun(), selectedId.GetPath().Data()));
		dataId = GetId(fValidFileIds, selectedId);

	} else {
		// List of files valid for reqested run was not loaded. Looking directly into CDB
		AliDebug(2, Form("List of files valid for run %d and version %d was not loaded. Looking directly into CDB for fileids valid for path %s!",
					selectedId.GetFirstRun(), selectedId.GetVersion(), selectedId.GetPath().Data()));

		TString filter;
		MakeQueryFilter(selectedId.GetFirstRun(), selectedId.GetLastRun(), 0, filter);

		TString pattern = ".root";
		TString optionQuery = "-y -m";
		if(selectedId.GetVersion() >= 0) {
			pattern.Prepend(Form("_v%d_s0",selectedId.GetVersion()));
			optionQuery = "";
		}

		TString folderCopy(Form("%s%s/Run",fDBFolder.Data(),selectedId.GetPath().Data()));

		if (optionQuery.Contains("-y")){
			AliInfo("Only latest version will be returned");
		}

		AliDebug(2,Form("** fDBFolder = %s, pattern = %s, filter = %s",folderCopy.Data(), pattern.Data(), filter.Data()));
		TGridResult *res = gGrid->Query(folderCopy, pattern, filter, optionQuery.Data());
		if (res) {
			AliCDBId validFileId;
			for(int i=0; i<res->GetEntries(); i++){
				TString filename = res->GetKey(i, "lfn");
				if(filename == "") continue;
				if(FilenameToId(filename, validFileId))
						validFileIds.AddLast(validFileId.Clone());
			}
			delete res;
		}else{
		    return 0; // this should be only in case of file catalogue glitch
		}

		dataId = GetId(validFileIds, selectedId);
	}

	return dataId;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntry(const AliCDBId& queryId) {
// get AliCDBEntry from the database

        Printf("Entering AliCDBGrid::GetEntry");

	AliCDBId* dataId = GetEntryId(queryId);

	if (!dataId){
                AliFatal(TString::Format("No valid CDB object found! request was: %s", queryId.ToString().Data()));
                return NULL;
        }

	TString filename;
	if (!IdToFilename(*dataId, filename)) {
		AliDebug(2,Form("Bad data ID encountered! Subnormal error!"));
		delete dataId;
                AliFatal(TString::Format("No valid CDB object found! request was: %s", queryId.ToString().Data()));
		return NULL;
	}

	AliCDBEntry* anEntry = GetEntryFromFile(filename, dataId);

	delete dataId;
	if(!anEntry){
                Printf("follows an alifatal");
                AliFatal(TString::Format("No valid CDB object found! request was: %s", queryId.ToString().Data()));
        }

	return anEntry;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntryFromFile(TString& filename, AliCDBId* dataId){
// Get AliCBEntry object from file "filename"

	AliDebug(2,Form("Opening file: %s",filename.Data()));

	filename.Prepend("/alien");

	// if option="CACHEREAD" TFile will use the local caching facility!
	TString option="READ";
	if(fCacheFolder != ""){

		// Check if local cache folder was changed in the meanwhile
		TString origCache(TFile::GetCacheFileDir());
		if(fCacheFolder != origCache) {
			AliWarning(Form("Local cache folder has been overwritten!! fCacheFolder = %s origCache = %s",
					fCacheFolder.Data(), origCache.Data()));
			TFile::SetCacheFileDir(fCacheFolder.Data(), fOperateDisconnected);
			TFile::ShrinkCacheFileDir(fCacheSize, fCleanupInterval);
		}

		option.Prepend("CACHE");
        }

	AliDebug(2, Form("Option: %s", option.Data()));

	TFile *file = TFile::Open(filename, option);
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
		Int_t tmpSubVersion = dataId->GetSubVersion();
		dataId->SetSubVersion(entryId.GetSubVersion()); // otherwise filename and id may mismatch
		if(!entryId.IsEqual(dataId)){
			AliWarning(Form("Mismatch between file name and object's Id!"));
			AliWarning(Form("File name: %s", dataId->ToString().Data()));
			AliWarning(Form("Object's Id: %s", entryId.ToString().Data()));
		}
		dataId->SetSubVersion(tmpSubVersion);
	}

	anEntry->SetLastStorage("grid");

	// Check whether entry contains a TTree. In case load the tree in memory!
	LoadTreeFromFile(anEntry);

	// close file, return retieved entry
	file->Close(); delete file; file=0;

	return anEntry;
}

//_____________________________________________________________________________
TList* AliCDBGrid::GetEntries(const AliCDBId& queryId) {
// multiple request (AliCDBStorage::GetAll)

	TList* result = new TList();
	result->SetOwner();

	TObjArray validFileIds;
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
		MakeQueryFilter(queryId.GetFirstRun(), queryId.GetLastRun(), 0, filter);

		TString path = queryId.GetPath();

		TString pattern = "Run*.root";
		TString optionQuery = "-y";

		TString addFolder = "";
		if (!path.Contains("*")){
		    if (!path.BeginsWith("/")) addFolder += "/";
		    addFolder += path;
		}
		else{
		    if (path.BeginsWith("/")) path.Remove(0,1);
		    if (path.EndsWith("/")) path.Remove(path.Length()-1,1);	
		    TObjArray* tokenArr = path.Tokenize("/");
		    if (tokenArr->GetEntries() != 3) {
			AliError("Not a 3 level path! Keeping old query...");
			pattern.Prepend(path+"/");
		    }
		    else{
			TString str0 = ((TObjString*)tokenArr->At(0))->String();
			TString str1 = ((TObjString*)tokenArr->At(1))->String();
			TString str2 = ((TObjString*)tokenArr->At(2))->String();
			if (str0 != "*" && str1 != "*" && str2 == "*"){
			    // e.g. "ITS/Calib/*"
			    addFolder = "/"+str0+"/"+str1;
			}
			else if (str0 != "*" && str1 == "*" && str2 == "*"){	
			    // e.g. "ITS/*/*"
			    addFolder = "/"+str0;
			}
			else if (str0 == "*" && str1 == "*" && str2 == "*"){	
			    // e.g. "*/*/*"
			    // do nothing: addFolder is already an empty string;
			}
			else{
			    // e.g. "ITS/*/RecoParam"
			    pattern.Prepend(path+"/");
			}
		    }
		    delete tokenArr; tokenArr=0;
		}

		TString folderCopy(Form("%s%s",fDBFolder.Data(),addFolder.Data()));

		AliDebug(2,Form("fDBFolder = %s, pattern = %s, filter = %s",folderCopy.Data(), pattern.Data(), filter.Data()));

		TGridResult *res = gGrid->Query(folderCopy, pattern, filter, optionQuery.Data());  

		if (!res) {
		    AliError("Grid query failed");
		    return 0;
		}

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

	TObjArray selectedIds;
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
		if(dataId) selectedIds.Add(dataId);
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
Bool_t AliCDBGrid::PutEntry(AliCDBEntry* entry, const char* mirrors) {
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

	TDirectory* saveDir = gDirectory;

	TString fullFilename = Form("/alien%s", filename.Data());
	TString seMirrors(mirrors);
	if(seMirrors.IsNull() || seMirrors.IsWhitespace()) seMirrors=GetMirrorSEs();
	// specify SE to filename
	// if a list of SEs was passed to this method or set via SetMirrorSEs, set the first as SE for opening the file.
	// The other SEs will be used in cascade in case of failure in opening the file.
	// The remaining SEs will be used to create replicas.
	TObjArray *arraySEs = seMirrors.Tokenize(',');
	Int_t nSEs = arraySEs->GetEntries();
	Int_t remainingSEs = 1;
	if(nSEs == 0){
		if (fSE != "default") fullFilename += Form("?se=%s",fSE.Data());
	}else{
		remainingSEs = nSEs;
	}

	// open file
	TFile *file=0;
	AliDebug(2, Form("fNretry = %d, fInitRetrySeconds = %d",fNretry,fInitRetrySeconds));
	TString targetSE("");

	Bool_t result = kFALSE;
	Bool_t reOpenResult = kFALSE;
	Int_t reOpenAttempts=0;
	while( !reOpenResult && reOpenAttempts<2){ //loop to check the file after closing it, to catch the unlikely but possible case when the file
		// is cleaned up by alien just before closing as a consequence of a network disconnection while writing

		while( !file && remainingSEs>0){
			if(nSEs!=0){
				TObjString *target = (TObjString*) arraySEs->At(nSEs-remainingSEs);
				targetSE=target->String();
				if ( !(targetSE.BeginsWith("ALICE::") && targetSE.CountChar(':')==4) ) {
					AliError(Form("\"%s\" is an invalid storage element identifier.",targetSE.Data()));
					continue;
				}
				if(fullFilename.Contains('?')) fullFilename.Remove(fullFilename.Last('?'));
				fullFilename += Form("?se=%s",targetSE.Data());
			}
			Int_t remainingAttempts=fNretry;
			Int_t nsleep = fInitRetrySeconds; // number of seconds between attempts. We let it increase exponentially
			AliDebug(2, Form("Uploading file into SE #%d: %s",nSEs-remainingSEs+1,targetSE.Data()));
			while(remainingAttempts > 0) {
				AliDebug(2, Form("Uploading file into OCDB at %s - Attempt #%d",targetSE.Data(),fNretry-remainingAttempts+1));
				remainingAttempts--;
				file = TFile::Open(fullFilename,"CREATE");
				if(!file || !file->IsWritable()){
					if(file) file->Close(); delete file; file=0; // file is not writable
					TString message(TString::Format("Attempt %d failed.",fNretry-remainingAttempts));
					if(remainingAttempts>0) {
						message += " Sleeping for "; message += nsleep; message += " seconds";
					}else{
						if(remainingSEs>0) message += " Trying to upload at next SE";
					}
					AliDebug(2, message.Data());
					if(remainingAttempts>0) sleep(nsleep);
				}else{
					remainingAttempts=0;
				}
				nsleep*=fInitRetrySeconds;
			}
			remainingSEs--;
		}
		if(!file){
			AliError(Form("All %d attempts have failed on all %d SEs. Returning...",fNretry,nSEs));
			return kFALSE;
		}

		file->cd();

		//SetTreeToFile(entry, file);
		entry->SetVersion(id.GetVersion());

		// write object (key name: "AliCDBEntry")
		result = (file->WriteTObject(entry, "AliCDBEntry") != 0);
		if (!result) AliError(Form("Can't write entry to file <%s>!", filename.Data()));
		file->Close();

		if(result)
		{
			AliDebug(2, Form("Reopening file %s for checking its correctness",fullFilename.Data()));
			TFile* ffile = TFile::Open(fullFilename.Data(),"READ");
			if(!ffile){
				reOpenResult = kFALSE;
				AliInfo(Form("The file %s was closed successfully but cannot be reopened. Trying now to regenerate it (regeneration attempt number %d)",
							fullFilename.Data(),++reOpenAttempts));
				delete file; file=0;
				AliDebug(2, Form("Removing file %s", filename.Data()));
				if(!gGrid->Rm(filename.Data()))
					AliError("Can't delete file!");
				remainingSEs++;
			}else{
				reOpenResult = kTRUE;
				ffile->Close();
			}
			delete ffile; ffile=0;
		}
	}

	if (saveDir) saveDir->cd(); else gROOT->cd();
	delete file; file=0;

	if(result && reOpenResult) {

		if(!TagFileId(filename, &id)){
			AliInfo(Form("CDB tagging failed. Deleting file %s!",filename.Data()));
			if(!gGrid->Rm(filename.Data()))
				AliError("Can't delete file!");
			return kFALSE;
		}

		TagFileMetaData(filename, entry->GetMetaData());
	}else{
		AliError("The file could not be opend or the object could not be written");
		if(!gGrid->Rm(filename.Data()))
			AliError("Can't delete file!");
		return kFALSE;
	}

	AliInfo(Form("CDB object stored into file %s", filename.Data()));
	if(nSEs==0)
		AliInfo(Form("Storage Element: %s", fSE.Data()));
	else
		AliInfo(Form("Storage Element: %s", targetSE.Data()));

	//In case of other SEs specified by the user, mirror the file to the remaining SEs
	for(Int_t i=0; i<nSEs; i++){
		if(i==nSEs-remainingSEs-1) continue; // skip mirroring to the SE where the file was saved
		TString mirrorCmd("mirror ");
		mirrorCmd += filename;
		mirrorCmd += " ";
		TObjString *target = (TObjString*) arraySEs->At(i);
		TString mirrorSE(target->String());
		mirrorCmd += mirrorSE;
		AliDebug(5,Form("mirror command: \"%s\"",mirrorCmd.Data()));
		AliInfo(Form("Mirroring to storage element: %s", mirrorSE.Data()));
		gGrid->Command(mirrorCmd.Data());
	}
	arraySEs->Delete(); arraySEs=0;

	return kTRUE;
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


        TString dirname(filename);
	Int_t dirNumber = gGrid->Mkdir(dirname.Remove(dirname.Last('/')),"-d");
	
	TString addTagValue1 = Form("addTagValue %s CDB ", filename.Data());
	TString addTagValue2 = Form("first_run=%d last_run=%d version=%d ",
					id->GetFirstRun(),
					id->GetLastRun(),
					id->GetVersion());
	TString addTagValue3 = Form("path_level_0=\"%s\" path_level_1=\"%s\" path_level_2=\"%s\" ",
					id->GetPathLevel(0).Data(),
					id->GetPathLevel(1).Data(),
					id->GetPathLevel(2).Data());
	//TString addTagValue4 = Form("version_path=\"%s\" dir_number=%d",Form("%d_%s",id->GetVersion(),filename.Data()),dirNumber); 
	TString addTagValue4 = Form("version_path=\"%09d%s\" dir_number=%d",id->GetVersion(),filename.Data(),dirNumber); 
	TString addTagValue = Form("%s%s%s%s",
					addTagValue1.Data(),
					addTagValue2.Data(),
					addTagValue3.Data(),
					addTagValue4.Data());

	Bool_t result = kFALSE;
	AliDebug(2, Form("Tagging file. Tag command: %s", addTagValue.Data()));
	TGridResult* res = gGrid->Command(addTagValue.Data());
	const char* resCode = res->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliError(Form("Couldn't add CDB tag value to file %s !",
						filename.Data()));
		result = kFALSE;
	} else {
		AliDebug(2, "Object successfully tagged.");
		result = kTRUE;
	}
	delete res;
	return result;

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::TagShortLived(TString& filename, Bool_t value){
// tag folder with ShortLived tag

	TString addTagValue = Form("addTagValue %s ShortLived_try value=%d", filename.Data(), value);

	Bool_t result = kFALSE;
	AliDebug(2, Form("Tagging file. Tag command: %s", addTagValue.Data()));
	TGridResult* res = gGrid->Command(addTagValue.Data());
	const char* resCode = res->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliError(Form("Couldn't add ShortLived tag value to file %s !", filename.Data()));
		result = kFALSE;
	} else {
		AliDebug(2,"Object successfully tagged.");
		result = kTRUE;
	}
	delete res;
	return result;

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::TagFileMetaData(TString& filename, const AliCDBMetaData* md){
// tag stored object in CDB table using object Id's parameters

	TString addTagValue1 = Form("addTagValue %s CDB_MD ", filename.Data());
	TString addTagValue2 = Form("object_classname=\"%s\" responsible=\"%s\" beam_period=%d ",
					md->GetObjectClassName(),
					md->GetResponsible(),
					md->GetBeamPeriod());
	TString addTagValue3 = Form("aliroot_version=\"%s\" comment=\"%s\"",
					md->GetAliRootVersion(),
					md->GetComment());
	TString addTagValue = Form("%s%s%s",
					addTagValue1.Data(),
					addTagValue2.Data(),
					addTagValue3.Data());

	Bool_t result = kFALSE;
	AliDebug(2, Form("Tagging file. Tag command: %s", addTagValue.Data()));
	TGridResult* res = gGrid->Command(addTagValue.Data());
	const char* resCode = res->GetKey(0,"__result__"); // '1' if success
	if(resCode[0] != '1') {
		AliWarning(Form("Couldn't add CDB_MD tag value to file %s !",
						filename.Data()));
		result = kFALSE;
	} else {
		AliDebug(2,"Object successfully tagged.");
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
	MakeQueryFilter(fRun, fRun, fMetaDataFilter, filter);

	TString path = fPathFilter.GetPath();

	TString pattern = "Run*";
	TString optionQuery = "-y";
	if(fVersion >= 0) {
		pattern += Form("_v%d_s0", fVersion);
		optionQuery = "";
	}
	pattern += ".root";
	AliDebug(2,Form("pattern: %s", pattern.Data()));

	TString addFolder = "";
	if (!path.Contains("*")){
		if (!path.BeginsWith("/")) addFolder += "/";
		addFolder += path;
	}
	else{
		if (path.BeginsWith("/")) path.Remove(0,1);
		if (path.EndsWith("/")) path.Remove(path.Length()-1,1);	
		TObjArray* tokenArr = path.Tokenize("/");
		if (tokenArr->GetEntries() != 3) {
			AliError("Not a 3 level path! Keeping old query...");
			pattern.Prepend(path+"/");
		}
		else{
			TString str0 = ((TObjString*)tokenArr->At(0))->String();
			TString str1 = ((TObjString*)tokenArr->At(1))->String();
			TString str2 = ((TObjString*)tokenArr->At(2))->String();
			if (str0 != "*" && str1 != "*" && str2 == "*"){
				// e.g. "ITS/Calib/*"
				addFolder = "/"+str0+"/"+str1;
			}
			else if (str0 != "*" && str1 == "*" && str2 == "*"){	
				// e.g. "ITS/*/*"
				addFolder = "/"+str0;
			}
			else if (str0 == "*" && str1 == "*" && str2 == "*"){	
				// e.g. "*/*/*"
				// do nothing: addFolder is already an empty string;
			}
			else{
				// e.g. "ITS/*/RecoParam"
				pattern.Prepend(path+"/");
			}
		}
		delete tokenArr; tokenArr=0;
	}

	TString folderCopy(Form("%s%s",fDBFolder.Data(),addFolder.Data()));

	AliDebug(2,Form("fDBFolder = %s, pattern = %s, filter = %s",folderCopy.Data(), pattern.Data(), filter.Data()));

	if (optionQuery == "-y"){
		AliInfo("Only latest version will be returned");
	} 

	TGridResult *res = gGrid->Query(folderCopy, pattern, filter, optionQuery.Data());  

	if (!res) {
		AliError("Grid query failed");
		return;
	}

	TIter next(res);
	TMap *map;
	while ((map = (TMap*)next())) {
	  TObjString *entry;
	  if ((entry = (TObjString *) ((TMap *)map)->GetValue("lfn"))) {
	    TString& filename = entry->String();
	    if(filename.IsNull()) continue;
	    AliDebug(2,Form("Found valid file: %s", filename.Data()));
	    AliCDBId *validFileId = new AliCDBId;
	    Bool_t result = FilenameToId(filename, *validFileId);
	    if(result) {
	      fValidFileIds.AddLast(validFileId);
	    }
	    else {
	      delete validFileId;
	    }
	  }
	}
	delete res;

}

//_____________________________________________________________________________
void AliCDBGrid::MakeQueryFilter(Int_t firstRun, Int_t lastRun,
					const AliCDBMetaData* md, TString& result) const
{
// create filter for file query

	result = Form("CDB:first_run<=%d and CDB:last_run>=%d", firstRun, lastRun);

//	if(version >= 0) {
//		result += Form(" and CDB:version=%d", version);
//	}
//	if(pathFilter.GetLevel0() != "*") {
//		result += Form(" and CDB:path_level_0=\"%s\"", pathFilter.GetLevel0().Data());
//	}
//	if(pathFilter.GetLevel1() != "*") {
//		result += Form(" and CDB:path_level_1=\"%s\"", pathFilter.GetLevel1().Data());
//	}
//	if(pathFilter.GetLevel2() != "*") {
//		result += Form(" and CDB:path_level_2=\"%s\"", pathFilter.GetLevel2().Data());
//	}

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


/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid factory  			                                                       //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBGridFactory)

//_____________________________________________________________________________
Bool_t AliCDBGridFactory::Validate(const char* gridString) {
// check if the string is valid Grid URI

        TRegexp gridPattern("^alien://.+$");

        return TString(gridString).Contains(gridPattern);
}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridFactory::CreateParameter(const char* gridString) {
// create AliCDBGridParam class from the URI string

	if (!Validate(gridString)) {
		return NULL;
	}

	TString buffer(gridString);

 	TString gridUrl 	= "alien://";
	TString user 		= "";
	TString dbFolder 	= "";
	TString se		= "default";
	TString cacheFolder	= "";
	Bool_t  operateDisconnected = kTRUE;
	Long64_t cacheSize          = (UInt_t) 1024*1024*1024; // 1GB
	Long_t cleanupInterval      = 0;

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
		else if (key.Contains("se",TString::kIgnoreCase)){
			se = value;
		}
		else if (key.Contains("cacheF",TString::kIgnoreCase)){
			cacheFolder = value;
			if (!cacheFolder.IsNull() && !cacheFolder.EndsWith("/"))
      				cacheFolder += "/";
		}
		else if (key.Contains("folder",TString::kIgnoreCase)){
			dbFolder = value;
		}
		else if (key.Contains("operateDisc",TString::kIgnoreCase)){
			if(value == "kTRUE") {
				operateDisconnected = kTRUE;
			} else if (value == "kFALSE") {
   				operateDisconnected = kFALSE;
			} else if (value == "0" || value == "1") {
				operateDisconnected = (Bool_t) value.Atoi();
			} else {
				AliError(Form("Invalid entry! %s",entry.Data()));
				return NULL;
			}
		}
		else if (key.Contains("cacheS",TString::kIgnoreCase)){
			if(value.IsDigit()) {
				cacheSize = value.Atoi();
			} else {
				AliError(Form("Invalid entry! %s",entry.Data()));
				return NULL;
			}
		}
		else if (key.Contains("cleanupInt",TString::kIgnoreCase)){
			if(value.IsDigit()) {
				cleanupInterval = value.Atoi();
			} else {
				AliError(Form("Invalid entry! %s",entry.Data()));
				return NULL;
			}
		}
		else{
			AliError(Form("Invalid entry! %s",entry.Data()));
			return NULL;
		}
	}
	delete arr; arr=0;

	AliDebug(2, Form("gridUrl:	%s", gridUrl.Data()));
	AliDebug(2, Form("user:	%s", user.Data()));
	AliDebug(2, Form("dbFolder:	%s", dbFolder.Data()));
	AliDebug(2, Form("s.e.:	%s", se.Data()));
	AliDebug(2, Form("local cache folder: %s", cacheFolder.Data()));
	AliDebug(2, Form("local cache operate disconnected: %d", operateDisconnected));
	AliDebug(2, Form("local cache size: %lld", cacheSize));
	AliDebug(2, Form("local cache cleanup interval: %ld", cleanupInterval));

	if(dbFolder == ""){
		AliError("Base folder must be specified!");
		return NULL;
	}

	return new AliCDBGridParam(gridUrl.Data(), user.Data(),
	                  dbFolder.Data(), se.Data(), cacheFolder.Data(),
			  operateDisconnected, cacheSize, cleanupInterval);
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
				      gridParam->GetSE().Data(),
				      gridParam->GetCacheFolder().Data(),
				      gridParam->GetOperateDisconnected(),
				      gridParam->GetCacheSize(),
				      gridParam->GetCleanupInterval());

	}

	if(!gGrid && grid) {
		delete grid; grid=0;
	}

	return grid;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid Parameter class  			                                               //                                         //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBGridParam)

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam():
 AliCDBParam(),
 fGridUrl(),
 fUser(),
 fDBFolder(),
 fSE(),
 fCacheFolder(),
 fOperateDisconnected(),
 fCacheSize(),
 fCleanupInterval()

 {
// default constructor

}

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam(const char* gridUrl, const char* user, const char* dbFolder,
				const char* se, const char* cacheFolder, Bool_t operateDisconnected,
				Long64_t cacheSize, Long_t cleanupInterval):
 AliCDBParam(),
 fGridUrl(gridUrl),
 fUser(user),
 fDBFolder(dbFolder),
 fSE(se),
 fCacheFolder(cacheFolder),
 fOperateDisconnected(operateDisconnected),
 fCacheSize(cacheSize),
 fCleanupInterval(cleanupInterval)
{
// constructor

	SetType("alien");

	TString uri = Form("%s?User=%s?DBFolder=%s?SE=%s?CacheFolder=%s"
			"?OperateDisconnected=%d?CacheSize=%lld?CleanupInterval=%ld",
			fGridUrl.Data(), fUser.Data(),
			fDBFolder.Data(), fSE.Data(), fCacheFolder.Data(),
			fOperateDisconnected, fCacheSize, fCleanupInterval);

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
					fDBFolder.Data(), fSE.Data(), fCacheFolder.Data(),
					fOperateDisconnected, fCacheSize, fCleanupInterval);
}

//_____________________________________________________________________________
ULong_t AliCDBGridParam::Hash() const {
// return Hash function

        return fGridUrl.Hash()+fUser.Hash()+fDBFolder.Hash()+fSE.Hash()+fCacheFolder.Hash();
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
        if(fCacheFolder != other->fCacheFolder) return kFALSE;
        if(fOperateDisconnected != other->fOperateDisconnected) return kFALSE;
        if(fCacheSize != other->fCacheSize) return kFALSE;
        if(fCleanupInterval != other->fCleanupInterval) return kFALSE;
	return kTRUE;
}

