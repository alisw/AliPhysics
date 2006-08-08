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
		AliInfo("Connection to the Grid!!!!");
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

	// return to the initial directory
	gGrid->Cd(initDir.Data(),0);
}

//_____________________________________________________________________________
AliCDBGrid::~AliCDBGrid()
{
// destructor

}

//_____________________________________________________________________________
Bool_t AliCDBGrid::FilenameToId(const char* filename, AliCDBRunRange& runRange, 
				Int_t& gridVersion) {
// build AliCDBId from filename numbers

        Ssiz_t mSize;
	
	// valid filename: Run#firstRun_#lastRun_v#version.root
        TRegexp keyPattern("^Run[0-9]+_[0-9]+_v[0-9]+.root$");
        keyPattern.Index(filename, &mSize);
        if (!mSize) {
		AliDebug(2,Form("Bad filename <%s>.", filename));
                return kFALSE;
        }

	TString idString(filename);
	idString.Resize(idString.Length() - sizeof(".root") + 1);

        TObjArray* strArray = (TObjArray*) idString.Tokenize("_");

	TString firstRunString(((TObjString*) strArray->At(0))->GetString());
	runRange.SetFirstRun(atoi(firstRunString.Data() + 3));
	runRange.SetLastRun(atoi(((TObjString*) strArray->At(1))->GetString()));
	
	TString verString(((TObjString*) strArray->At(2))->GetString());
	gridVersion = atoi(verString.Data() + 1);

        delete strArray;

        return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBGrid::IdToFilename(const AliCDBRunRange& runRange, Int_t gridVersion, 
				 TString& filename) {
// build file name from AliCDBId data (run range, version)

	if (!runRange.IsValid()) {
		AliDebug(2,Form("Invalid run range <%d, %d>.", 
			runRange.GetFirstRun(), runRange.GetLastRun()));
		return kFALSE;
	}

	if (gridVersion < 0) {
		AliDebug(2,Form("Invalid version <%d>.", gridVersion));
                return kFALSE;
	}
 
        filename += "Run";
	filename += runRange.GetFirstRun();
        filename += "_";
        filename += runRange.GetLastRun();
	filename += "_v";
	filename += gridVersion;
	filename += ".root";

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

	const char* filename;
	AliCDBRunRange aRunRange; // the runRange got from filename 
	AliCDBRunRange lastRunRange(-1,-1); // highest runRange found
	Int_t aVersion; // the version got from filename 
	Int_t lastVersion=0; // highest version found

	TGridResult *res = gGrid->Ls(dirName);

	//loop on the files in the directory, look for highest version
	for(int i=0; i < res->GetEntries(); i++){
		filename=res->GetFileName(i);
		if (!FilenameToId(filename, aRunRange, aVersion)) continue;
		if (aRunRange.Overlaps(id.GetAliCDBRunRange()) && aVersion > lastVersion) {
			lastVersion = aVersion;
			lastRunRange = aRunRange;
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
Bool_t AliCDBGrid::GetId(const AliCDBId& query, AliCDBId& result) {
// look for filename matching query (called by GetEntry)

	TString initDir(gGrid->Pwd(0));

	TString dirName(fDBFolder);
	dirName += query.GetPath(); // dirName = fDBFolder/idPath

	if (!gGrid->Cd(dirName,0)) {
		AliDebug(2,Form("Directory <%s> not found", (query.GetPath()).Data()));
		AliDebug(2,Form("in DB folder %s", fDBFolder.Data()));
		return kFALSE;
	}
 
	TGridResult *res = gGrid->Ls(dirName);

	const char* filename;
	AliCDBRunRange aRunRange; // the runRange got from filename
	Int_t aVersion; // the version got from filename
 
	for(int i=0; i < res->GetEntries(); i++){
		filename=res->GetFileName(i);
		if (!FilenameToId(filename, aRunRange, aVersion)) continue;
                // aRunRange and aVersion filled from filename
		
		if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue; 
		// aRunRange contains requested run!

		if (!query.HasVersion()){ // look for highest version
			if(result.GetVersion() > aVersion) continue;
			if(result.GetVersion() == aVersion) {
				AliDebug(2,Form("More than one object valid for run %d, version %d!", 
					query.GetFirstRun(), aVersion));
			return kFALSE; 
			}
			result.SetVersion(aVersion);
			result.SetFirstRun(aRunRange.GetFirstRun());
			result.SetLastRun(aRunRange.GetLastRun());

		} else { // look for specified version
			if(query.GetVersion() != aVersion) continue;
			if(result.GetVersion() == aVersion){
				AliDebug(2,Form("More than one object valid for run %d, version %d!", 
					query.GetFirstRun(), aVersion));
					return kFALSE; 
			}
			result.SetVersion(aVersion);
			result.SetFirstRun(aRunRange.GetFirstRun());
			result.SetLastRun(aRunRange.GetLastRun());
		}
	} // end loop on filenames
	delete res;
	
	gGrid->Cd(initDir.Data(),0);
 
	return kTRUE;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntry(const AliCDBId& queryId) {
// get AliCDBEntry from the database

	AliCDBId dataId(queryId.GetAliCDBPath(), -1, -1, -1, -1);
        Bool_t result;
	
	// look for a filename matching query requests (path, runRange, version, subVersion)
	if (!queryId.HasVersion()) {
		// if version is not specified, first check the selection criteria list
		AliCDBId selectedId(queryId);
		GetSelection(&selectedId);
		result = GetId(selectedId,dataId);
	} else {
		result = GetId(queryId,dataId);
	}

	if (!result || !dataId.IsSpecified()) return NULL;

	TString filename;
	if (!IdToFilename(dataId.GetAliCDBRunRange(), dataId.GetVersion(),filename)) {
		AliDebug(2,Form("Bad data ID encountered! Subnormal error!"));
		return NULL;
	}

	filename.Prepend("/alien" + fDBFolder + queryId.GetPath() + '/');
	filename += "?se="; filename += fSE.Data(); 

	AliDebug(2,Form("Opening file: %s",filename.Data()));
	TFile *file = TFile::Open(filename);
	if (!file) {
		AliDebug(2,Form("Can't open file <%s>!", filename.Data()));
		return NULL;
	}

	// get the only AliCDBEntry object from the file
	// the object in the file is an AliCDBEntry entry named "AliCDBEntry"

	TObject* anObject = file->Get("AliCDBEntry");

	if (!anObject) {
		AliDebug(2,Form("Bad storage data: NULL entry object!"));
 		return NULL;
	}

	if (AliCDBEntry::Class() != anObject->IsA()) {
		AliDebug(2,Form("Bad storage data: Invalid entry object!"));
		return NULL;
	}

	AliCDBId entryId = ((AliCDBEntry* ) anObject)->GetId();
 
	// The object's Id is not reset during storage
	// If object's Id runRange or version do not match with filename,
	// it means that someone renamed file by hand. In this case a warning msg is issued.
 
	((AliCDBEntry*) anObject)->SetLastStorage("grid");
 
	if(!((entryId.GetAliCDBRunRange()).IsEqual(&dataId.GetAliCDBRunRange())) || 
		entryId.GetVersion() != dataId.GetVersion()){
		AliWarning(Form("Either RunRange or gridVersion in the object's metadata do noth match with fileName numbers:"));
		AliWarning(Form("someone renamed file by hand!"));
	}

	// close file, return retieved entry
	file->Close(); delete file; file=0;
	return (AliCDBEntry*) anObject;
}

//_____________________________________________________________________________
void AliCDBGrid::GetEntriesForLevel0(const char* level0,
	const AliCDBId& queryId, TList* result) {
// multiple request (AliCDBStorage::GetAll)

	TString level0Dir=fDBFolder;
	level0Dir += level0;

	if (!gGrid->Cd(level0Dir,0)) {
		AliDebug(2,Form("Level0 directory <%s> not found", level0Dir.Data()));
		return;
	}

	TGridResult *res = gGrid->Ls(level0Dir);
	TString level1;
	for(int i=0; i < res->GetEntries(); i++){
		level1=res->GetFileName(i);
	if (queryId.GetAliCDBPath().Level1Comprises(level1))
		GetEntriesForLevel1(level0, level1, queryId, result);
	}
	delete res;
}

//_____________________________________________________________________________
void AliCDBGrid::GetEntriesForLevel1(const char* level0, const char* level1,
	const AliCDBId& queryId, TList* result) {
// multiple request (AliCDBStorage::GetAll)

	TString level1Dir=fDBFolder;
	level1Dir += level0;
	level1Dir += '/';
	level1Dir += level1;

	if (!gGrid->Cd(level1Dir,0)) {
		AliDebug(2,Form("Level1 directory <%s> not found", level1Dir.Data()));
		return;
	}

	TGridResult *res = gGrid->Ls(level1Dir);
	TString level2;
	for(int i=0; i < res->GetEntries(); i++){
		level2=res->GetFileName(i);
		if (queryId.GetAliCDBPath().Level2Comprises(level2)){
			AliCDBPath entryPath(level0, level1, level2);
			AliCDBId entryId(entryPath, 
				queryId.GetAliCDBRunRange(),
				queryId.GetVersion(), 
				queryId.GetSubVersion());

			AliCDBEntry* anEntry = GetEntry(entryId);
			if (anEntry) result->Add(anEntry);

		}
	}
	delete res;
}

//_____________________________________________________________________________
TList* AliCDBGrid::GetEntries(const AliCDBId& queryId) {
// multiple request (AliCDBStorage::GetAll)

	TList* result = new TList();
	result->SetOwner();

	TString initDir(gGrid->Pwd(0));

	TGridResult *res = gGrid->Ls(fDBFolder);
	TString level0;

	for(int i=0; i < res->GetEntries(); i++){
		level0=res->GetFileName(i);
		if (queryId.GetAliCDBPath().Level0Comprises(level0)) 
			GetEntriesForLevel0(level0, queryId, result);				    
	}	 
	delete res;
 
	gGrid->Cd(initDir.Data(),0);
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
	if (!IdToFilename(id.GetAliCDBRunRange(), id.GetVersion(), filename)) {
		AliError("Bad ID encountered! Subnormal error!");
		return kFALSE;
	} 

	filename.Prepend("/alien" + fDBFolder + id.GetPath() + '/');
	TString filenameCopy(filename);
	filename += "?se="; filename += fSE.Data(); 
 
	TDirectory* saveDir = gDirectory;

	// open file
	TFile *file = TFile::Open(filename,"CREATE");
	if(!file || !file->IsWritable()){
		AliError(Form("Can't open file <%s>!", filename.Data()));    
		if(file && !file->IsWritable()) file->Close(); delete file; file=0;
		return kFALSE;
	}
  
	file->cd(); 

	entry->SetVersion(id.GetVersion());

	// write object (key name: "AliCDBEntry")
	Bool_t result = (entry->Write("AliCDBEntry") != 0); 
	if (!result) AliError(Form("Can't write entry to file <%s>!",filename.Data()));


	if (saveDir) saveDir->cd(); else gROOT->cd();
	file->Close(); delete file; file=0;
	if(result) {
		AliInfo(Form("CDB object stored into file %s",filenameCopy.Data()));
		AliInfo(Form("using S.E. %s", fSE.Data()));
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
	TObjString *str;
	
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
		
	AliInfo(Form("gridUrl:	%s",gridUrl.Data()));
	AliInfo(Form("user:	%s",user.Data()));
	AliInfo(Form("dbFolder:	%s",dbFolder.Data()));
	AliInfo(Form("s.e.:	%s",se.Data()));

	return new AliCDBGridParam(gridUrl, user, dbFolder, se);       
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBGridFactory::Create(const AliCDBParam* param) {
// create AliCDBGrid storage instance from parameters
	
	if (AliCDBGridParam::Class() == param->IsA()) {
		
		const AliCDBGridParam* gridParam = (const AliCDBGridParam*) param;
		AliCDBGrid *grid = new AliCDBGrid(gridParam->GridUrl(), 
				      gridParam->GetUser(), gridParam->GetDBFolder(),
				      gridParam->GetSE()); 

		if(gGrid) return grid;
	}

	return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBGrid Parameter class  			                                               //                                          //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBGridParam)

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam() {
// default constructor

}

//_____________________________________________________________________________
AliCDBGridParam::AliCDBGridParam(const char* gridUrl, 
				const char* user, 
			        const char* dbFolder, 
				const char* se):
 fGridUrl(gridUrl),
 fUser(user),
 fDBFolder(dbFolder),
 fSE(se)
{
// constructor
	
	SetType("alien");

	TString uri;
	uri+=fGridUrl; uri+="?";
	uri+="User="; 	uri+=fUser; uri+="?"; 
	uri+="DBFolder="; uri+=fDBFolder; uri+="?";
	uri+="SE="; 	uri+=fSE;
	
	SetURI(uri);
}

//_____________________________________________________________________________
AliCDBGridParam::~AliCDBGridParam() {
// destructor

}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridParam::CloneParam() const {
// clone parameter

        return new AliCDBGridParam(fGridUrl, fUser, fDBFolder, fSE);
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

