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
#include "AliCDBGrid.h"


ClassImp(AliCDBGrid)

//_____________________________________________________________________________
AliCDBGrid::AliCDBGrid(const char *host, const Int_t port, 
                       const char *user, const char *dbPath, const char *se) :
AliCDBStorage(),
fHost(host),
fPort(port),
fUser(user),
fDBPath(dbPath),
fSE(se)
{
// constructor //

	TString grid="alien://";
	grid+=host; grid+=":"; grid+=port;

	// if the same Grid is alreay active, skip connection
	if (!gGrid || TString(host) != gGrid->GetHost() || 
	    port != gGrid->GetPort() || TString(user) != gGrid->GetUser()) {
   		// connection to the Grid
		TGrid::Connect(grid.Data(),fUser.Data());
	}

	if(!gGrid) {
		AliError("Connection failed!");
		return;
	}

	TString initDir(gGrid->Pwd(0));
	if (fDBPath[0] != '/') {
		fDBPath.Prepend(initDir);
	}

	// check DBFolder: trying to cd to DBFolder; if it does not exist, create it
	if(!gGrid->Cd(fDBPath.Data(),0)){
		AliDebug(2,Form("Creating new folder <%s> ...",fDBPath.Data()));
		if(!gGrid->Mkdir(fDBPath.Data(),"",0)){
			AliError(Form("Cannot create folder <%s> !",fDBPath.Data())); 
		}
	} else {
		AliDebug(2,Form("Folder <%s> found",fDBPath.Data()));
	}

	// removes any '/' at the end of path, then append one '/'
	while(fDBPath.EndsWith("/")) fDBPath.Remove(fDBPath.Last('/')); 
	fDBPath+="/";

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
		AliWarning(Form("Invalid run range <%d, %d>.", 
			runRange.GetFirstRun(), runRange.GetLastRun()));
		return kFALSE;
	}

	if (gridVersion < 0) {
		AliWarning(Form("Invalid version <%d>.", gridVersion));
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

	TString dirName(fDBPath);

	Bool_t dirExist=kFALSE;
 
	// go to the path; if directory does not exist, create it
	TObjArray *arrName=pathName.Tokenize("/");
	for(int i=0;i<arrName->GetEntries();i++){
		TString buffer((arrName->At(i))->GetName());
 		dirName+=buffer; dirName+="/";
		dirExist=gGrid->Cd(dirName,0);
		if (!dirExist) {
			AliInfo(Form("Creating new folder <%s> ...",dirName.Data()));
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
		AliWarning(Form("*** WARNING! a NEW object is being stored with version %d",
					id.GetVersion()));
		AliWarning(Form("and it will hide previously stored object with version %d!",
					id.GetVersion()-1));
	}

	if(!lastRunRange.IsAnyRange() && !(lastRunRange.IsEqual(&id.GetAliCDBRunRange()))) 
    		AliWarning(Form("Run range modified w.r.t. previous version (Run%d_%d_v%d)",
    		     	lastRunRange.GetFirstRun(), lastRunRange.GetLastRun(), id.GetVersion()));
    
	return kTRUE;
}

//_____________________________________________________________________________
AliCDBId AliCDBGrid::GetId(const AliCDBId& query) {
// look for filename matching query (called by GetEntry)

	TString initDir(gGrid->Pwd(0));

	AliCDBId result(query.GetAliCDBPath(), -1, -1, -1, -1);

	TString dirName(fDBPath);
	dirName += query.GetPath(); // dirName = fDBPath/idPath

	if (!gGrid->Cd(dirName,0)) {
		AliError(Form("Directory <%s> not found", (query.GetPath()).Data()));
		AliError(Form("in DB folder %s", fDBPath.Data()));
		return result;
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
				AliError(Form("More than one object valid for run %d, version %d!", 
					query.GetFirstRun(), aVersion));
				result.SetRunRange(-1,-1); result.SetVersion(-1);
			return result; 
			}
			result.SetVersion(aVersion);
			result.SetFirstRun(aRunRange.GetFirstRun());
			result.SetLastRun(aRunRange.GetLastRun());

		} else { // look for specified version
			if(query.GetVersion() != aVersion) continue;
			if(result.GetVersion() == aVersion){
				AliError(Form("More than one object valid for run %d, version %d!", 
					query.GetFirstRun(), aVersion));
					result.SetRunRange(-1,-1); result.SetVersion(-1); 
					return result; 
			}
			result.SetVersion(aVersion);
			result.SetFirstRun(aRunRange.GetFirstRun());
			result.SetLastRun(aRunRange.GetLastRun());
		}
	} // end loop on filenames
	delete res;
	
	gGrid->Cd(initDir.Data(),0);
 
	return result;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBGrid::GetEntry(const AliCDBId& queryId) {
// get AliCDBEntry from the database

	AliCDBId dataId;
	
	// look for a filename matching query requests (path, runRange, version, subVersion)
	if (!queryId.HasVersion()) {
		// if version is not specified, first check the selection criteria list
		dataId = GetId(GetSelection(queryId));
	} else {
		dataId = GetId(queryId);
	}

	if (!dataId.IsSpecified()) return NULL;

	TString filename;
	if (!IdToFilename(dataId.GetAliCDBRunRange(), dataId.GetVersion(),filename)) {
		AliError("Bad data ID encountered! Subnormal error!");
		return NULL;
	}

	filename.Prepend("/alien" + fDBPath + queryId.GetPath() + '/');
	filename += "?se="; filename += fSE.Data(); 

	AliInfo(Form("Opening file: %s",filename.Data()));
	TFile *file = TFile::Open(filename);
	if (!file) {
		AliError(Form("Can't open file <%s>!", filename.Data()));
		return NULL;
	}

	// get the only AliCDBEntry object from the file
	// the object in the file is an AliCDBEntry entry named "AliCDBEntry"

	TObject* anObject = file->Get("AliCDBEntry");

	if (!anObject) {
		AliError("Bad storage data: NULL entry object!");
 		return NULL;
	}

	if (AliCDBEntry::Class() != anObject->IsA()) {
		AliError("Bad storage data: Invalid entry object!");
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

	TString level0Dir=fDBPath;
	level0Dir += level0;

	if (!gGrid->Cd(level0Dir,0)) {
		AliError(Form("Level0 directory <%s> not found", level0Dir.Data()));
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

	TString level1Dir=fDBPath;
	level1Dir += level0;
	level1Dir += '/';
	level1Dir += level1;

	if (!gGrid->Cd(level1Dir,0)) {
		AliError(Form("Level1 directory <%s> not found", level1Dir.Data()));
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

	TGridResult *res = gGrid->Ls(fDBPath);
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

	filename.Prepend("/alien" + fDBPath + id.GetPath() + '/');
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
		AliInfo(Form("AliCDBEntry stored into file %s",filenameCopy.Data()));
		AliInfo(Form("using S.E. %s", fSE.Data()));
	}
 
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
        TRegexp gridPattern("^alien://.+:[0-9]+;[a-zA-Z0-9_-.]+;.+;.+$");

        return TString(gridString).Contains(gridPattern);
}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridFactory::CreateParameter(const char* gridString) {
// create AliCDBGridParam class from the URI string

	if (!Validate(gridString)) {
		return NULL;
	}
	TString buffer(gridString + sizeof("alien://") - 1);
 	TString host = buffer(0,buffer.First(':')); // host (ex. aliendb4.cern.ch)
	buffer = buffer(host.Sizeof(),buffer.Sizeof());
	TString strPort = buffer(0, buffer.First(';'));
	Int_t port = atoi(strPort.Data());	// port number (ex. 9000)
	buffer = buffer(strPort.Sizeof(),buffer.Sizeof());
	TString user = buffer(0,buffer.First(';')); // user (ex. colla)
	buffer = buffer(user.Sizeof(),buffer.Sizeof());
	TString dbPath = buffer(0,buffer.First(';')); // DB path (ex. /alice/cern.ch/user/c/colla/DB)
	TString se = buffer(dbPath.Sizeof(),buffer.Sizeof()); // storage element (ex. ALICE::CERN::Server)
	
	AliInfo(Form("host: %s",host.Data()));
	AliInfo(Form("port: %d",port));
	AliInfo(Form("user: %s",user.Data()));
	AliInfo(Form("dbPath: %s",dbPath.Data()));
	AliInfo(Form("s.e.: %s",se.Data()));

	return new AliCDBGridParam(host, port, user, dbPath, se);       
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBGridFactory::Create(const AliCDBParam* param) {
// create AliCDBGrid storage instance from parameters
	
	if (AliCDBGridParam::Class() == param->IsA()) {
		
		const AliCDBGridParam* gridParam = (const AliCDBGridParam*) param;
		AliCDBGrid *grid = new AliCDBGrid(gridParam->GetHost(), gridParam->GetPort(), 
				      gridParam->GetUser(), gridParam->GetPath(),
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
AliCDBGridParam::AliCDBGridParam(const char* host, 
				const Int_t port, 
				const char* user, 
			        const char* dbPath, 
				const char* se):
 fHost(host),
 fPort(port),
 fUser(user),
 fDBPath(dbPath),
 fSE(se)
{
// constructor
	
	SetType("alien");

	TString uri=("alien://");
	uri+=host; uri+=":"; uri+=port; uri+=";";
	uri+=user; uri+=";"; uri+=dbPath; uri+=";";
	uri+=se;
	
	SetURI(uri);
}

//_____________________________________________________________________________
AliCDBGridParam::~AliCDBGridParam() {
// destructor

}

//_____________________________________________________________________________
AliCDBParam* AliCDBGridParam::CloneParam() const {
// clone parameter

        return new AliCDBGridParam(fHost, fPort, fUser, fDBPath, fSE);
}

//_____________________________________________________________________________
ULong_t AliCDBGridParam::Hash() const {
// return Hash function

        return fHost.Hash()+fPort+fUser.Hash()+fDBPath.Hash()+fSE.Hash();
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

        if(fHost != other->fHost) return kFALSE;
        if(fPort != other->fPort) return kFALSE;
        if(fUser != other->fUser) return kFALSE;
        if(fDBPath != other->fDBPath) return kFALSE;
        if(fSE != other->fSE) return kFALSE;
	return kTRUE;
}

