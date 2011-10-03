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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBDump						   //
//  access class to a DataBase in a dump storage (single file)     //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <TSystem.h>
#include <TKey.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TObjString.h>
#include <TList.h>

#include "AliCDBDump.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

ClassImp(AliCDBDump)

//_____________________________________________________________________________
AliCDBDump::AliCDBDump(const char* dbFile, Bool_t readOnly):
fFile(NULL), fReadOnly(readOnly) {
// constructor

	// opening file
	fFile = TFile::Open(dbFile, fReadOnly ? "READ" : "UPDATE");	
	if (!fFile) {
		AliError(Form("Can't open file <%s>!" , dbFile));
	} else {
		AliDebug(2,Form("File <%s> opened",dbFile));
		if(fReadOnly) AliDebug(2,Form("in read-only mode"));
	}

	fType="dump";
	fBaseFolder = dbFile;
}

//_____________________________________________________________________________
AliCDBDump::~AliCDBDump() {
// destructor

	if (fFile) {
		fFile->Close();
		delete fFile;
	}
}


//_____________________________________________________________________________
Bool_t AliCDBDump::KeyNameToId(const char* keyname, AliCDBRunRange& runRange,
	Int_t& version, Int_t& subVersion) {
// build AliCDBId from keyname numbers

        Ssiz_t mSize;

	// valid keyname: Run#firstRun_#lastRun_v#version_s#subVersion.root
        TRegexp keyPattern("^Run[0-9]+_[0-9]+_v[0-9]+_s[0-9]+$");
        keyPattern.Index(keyname, &mSize);
        if (!mSize) {
                AliDebug(2,Form("Bad keyname <%s>.", keyname));
                return kFALSE;
        }

        TObjArray* strArray = (TObjArray*) TString(keyname).Tokenize("_");

	TString firstRunString(((TObjString*) strArray->At(0))->GetString());
	runRange.SetFirstRun(atoi(firstRunString.Data() + 3));
	runRange.SetLastRun(atoi(((TObjString*) strArray->At(1))->GetString()));
	
	TString verString(((TObjString*) strArray->At(2))->GetString());
	version = atoi(verString.Data() + 1);

	TString subVerString(((TObjString*) strArray->At(3))->GetString());
        subVersion = atoi(subVerString.Data() + 1);

        delete strArray;

        return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBDump::IdToKeyName(const AliCDBRunRange& runRange, Int_t version,
        Int_t subVersion, TString& keyname) {
// build key name from AliCDBId data (run range, version, subVersion)

        if (!runRange.IsValid()) {
                AliDebug(2,Form("Invalid run range <%d, %d>.",
                        runRange.GetFirstRun(), runRange.GetLastRun()));
                return kFALSE;
        }

        if (version < 0) {
                AliDebug(2,Form("Invalid version <%d>.", version));
                return kFALSE;
        }

	if (subVersion < 0) {
		AliDebug(2,Form("Invalid subversion <%d>.", subVersion));
		return kFALSE;
	}
    
        keyname += "Run";
	keyname += runRange.GetFirstRun();
        keyname += "_";
        keyname += runRange.GetLastRun();
        keyname += "_v";
        keyname += version;
	keyname += "_s";
	keyname += subVersion;

        return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBDump::MkDir(const TString& path) {
// descend into TDirectory, making TDirectories if they don't exist 
	TObjArray* strArray = (TObjArray*) path.Tokenize("/");
	
	TIter iter(strArray);
	TObjString* str;

	while ((str = (TObjString*) iter.Next())) {
		
		TString dirName(str->GetString());
		if (!dirName.Length()) {
			continue;
		}

		if (gDirectory->cd(dirName)) {
			continue;
		}

		TDirectory* aDir = gDirectory->mkdir(dirName, "");
		if (!aDir) {
			AliError(Form("Can't create directory <%s>!", 
					dirName.Data()));
			delete strArray;

			return kFALSE;
		}

		aDir->cd();
	}	
	
	delete strArray;

	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBDump::PrepareId(AliCDBId& id) {
// prepare id (version, subVersion) of the object that will be stored (called by PutEntry)

	AliCDBRunRange aRunRange; // the runRange got from filename 
 	AliCDBRunRange lastRunRange(-1,-1); // highest runRange found
	Int_t aVersion, aSubVersion; // the version subVersion got from filename
	Int_t lastVersion = 0, lastSubVersion = -1; // highest version and subVersion found

	
	TIter iter(gDirectory->GetListOfKeys());
	TKey* key;

	if (!id.HasVersion()) {	// version not specified: look for highest version & subVersion	
				
		while ((key = (TKey*) iter.Next())) { // loop on keys

			const char* keyName = key->GetName();
	
			if (!KeyNameToId(keyName, aRunRange, aVersion, 
			   aSubVersion)) {
				AliDebug(2,Form(
					"Bad keyname <%s>!I'll skip it.", keyName));
				continue;
			}
			
			if (!aRunRange.Overlaps(id.GetAliCDBRunRange())) continue;
			if(aVersion < lastVersion) continue;
			if(aVersion > lastVersion) lastSubVersion = -1;
			if(aSubVersion < lastSubVersion) continue;
			lastVersion = aVersion;
			lastSubVersion = aSubVersion;
			lastRunRange = aRunRange;
		}

		id.SetVersion(lastVersion);
		id.SetSubVersion(lastSubVersion + 1);

	} else { // version specified, look for highest subVersion only
		
		while ((key = (TKey*) iter.Next())) { // loop on the keys

			const char* keyName = key->GetName();
	
			if (!KeyNameToId(keyName, aRunRange, aVersion, 
			   aSubVersion)) {
				AliDebug(2,Form(
					"Bad keyname <%s>!I'll skip it.", keyName));
				continue;
			}

			if (aRunRange.Overlaps(id.GetAliCDBRunRange()) 
				&& aVersion == id.GetVersion()
				&& aSubVersion > lastSubVersion) {
				lastSubVersion = aSubVersion;
				lastRunRange = aRunRange;
			}
	
		}
		
		id.SetSubVersion(lastSubVersion + 1);
	}

	TString lastStorage = id.GetLastStorage();
	if(lastStorage.Contains(TString("grid"), TString::kIgnoreCase) &&
	   id.GetSubVersion() > 0 ){
		AliError(Form("Grid to Dump Storage error! local object with version v%d_s%d found:",id.GetVersion(), id.GetSubVersion()-1));
		AliError(Form("This object has been already transferred from Grid (check v%d_s0)!",id.GetVersion()));
		return kFALSE;
	}

	if(lastStorage.Contains(TString("new"), TString::kIgnoreCase) &&
	   id.GetSubVersion() > 0 ){
		AliDebug(2, Form("A NEW object is being stored with version v%d_s%d",
					id.GetVersion(),id.GetSubVersion()));
		AliDebug(2, Form("and it will hide previously stored object with v%d_s%d!",
					id.GetVersion(),id.GetSubVersion()-1));
	}

 	if(!lastRunRange.IsAnyRange() && !(lastRunRange.IsEqual(& id.GetAliCDBRunRange()))) 
    		AliWarning(Form("Run range modified w.r.t. previous version (Run%d_%d_v%d_s%d)",
    		     	lastRunRange.GetFirstRun(), lastRunRange.GetLastRun(), 
			id.GetVersion(), id.GetSubVersion()-1));

	return kTRUE;

}

// //_____________________________________________________________________________
// Bool_t AliCDBDump::GetId(const AliCDBId& query, AliCDBId& result) {
// // look for filename matching query (called by GetEntry)
// 
// 
//         AliCDBRunRange aRunRange; // the runRange got from filename
//         Int_t aVersion, aSubVersion; // the version and subVersion got from filename
// 
// 	TIter iter(gDirectory->GetListOfKeys());
// 	TKey* key;
// 
// 	if (!query.HasVersion()) { // neither version and subversion specified -> look for highest version and subVersion
// 
//                 while ((key = (TKey*) iter.Next())) { // loop on the keys
// 
// 			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
//                         // aRunRange, aVersion, aSubVersion filled from filename
// 
// 			if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue;
// 			// aRunRange contains requested run!
// 
// 			if (result.GetVersion() < aVersion) {
// 				result.SetVersion(aVersion);
// 				result.SetSubVersion(aSubVersion);
// 
// 				result.SetFirstRun(
// 					aRunRange.GetFirstRun());
// 				result.SetLastRun(
// 					aRunRange.GetLastRun());
// 
// 			} else if (result.GetVersion() == aVersion
// 				&& result.GetSubVersion()
// 					< aSubVersion) {
// 
// 				result.SetSubVersion(aSubVersion);
// 
//                         	result.SetFirstRun(
// 					aRunRange.GetFirstRun());
//                         	result.SetLastRun(
// 					aRunRange.GetLastRun());
// 			} else if (result.GetVersion() == aVersion
// 				&& result.GetSubVersion() == aSubVersion){
//               			AliDebug(2,Form("More than one object valid for run %d, version %d_%d!",
// 		       			query.GetFirstRun(), aVersion, aSubVersion));
// 	      			result.SetRunRange(-1,-1); result.SetVersion(-1); result.SetSubVersion(-1);
// 	      			return kFALSE;
// 				}
// 		}
// 
// 	} else if (!query.HasSubVersion()) { // version specified but not subversion -> look for highest subVersion
// 
// 		result.SetVersion(query.GetVersion());
// 
//                 while ((key = (TKey*) iter.Next())) { // loop on the keys
// 
// 			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
//                         // aRunRange, aVersion, aSubVersion filled from filename
// 
//                        if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue;
// 			// aRunRange contains requested run!
// 			
// 			if(query.GetVersion() != aVersion) continue;
// 			// aVersion is requested version!
// 
// 	 		if(result.GetSubVersion() == aSubVersion){
//               			AliDebug(2,Form("More than one object valid for run %d, version %d_%d!",
// 		       			query.GetFirstRun(), aVersion, aSubVersion));
// 	     			result.SetRunRange(-1,-1); result.SetVersion(-1); result.SetSubVersion(-1);
// 	     			return kFALSE;
// 	 		}
// 			if( result.GetSubVersion() < aSubVersion) {
// 
//                                 result.SetSubVersion(aSubVersion);
// 
//                                 result.SetFirstRun(
// 					aRunRange.GetFirstRun());
//                                 result.SetLastRun(
// 					aRunRange.GetLastRun());
// 			}
//                 }
// 
// 	} else { // both version and subversion specified
// 
//                 while ((key = (TKey*) iter.Next())) { // loop on the keys
// 
// 			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
//                         // aRunRange, aVersion, aSubVersion filled from filename
// 
// 			if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue;
//  			// aRunRange contains requested run!
// 
// 			if(query.GetVersion() != aVersion || query.GetSubVersion() != aSubVersion) continue;
// 			// aVersion and aSubVersion are requested version and subVersion!
// 
// 			if(result.GetVersion() == aVersion && result.GetSubVersion() == aSubVersion){
//               			AliDebug(2,Form("More than one object valid for run %d, version %d_%d!",
// 		       			query.GetFirstRun(), aVersion, aSubVersion));
// 	     			result.SetRunRange(-1,-1); result.SetVersion(-1); result.SetSubVersion(-1);
// 	     			return kFALSE; 
// 			}
// 			result.SetVersion(aVersion);
// 		        result.SetSubVersion(aSubVersion);
// 			result.SetFirstRun(aRunRange.GetFirstRun());
// 			result.SetLastRun(aRunRange.GetLastRun());
// 
// 		}
// 	}
// 
// 	return kTRUE;
// }

//_____________________________________________________________________________
AliCDBId* AliCDBDump::GetId(const AliCDBId& query) {
// look for filename matching query (called by GetEntry)


        AliCDBRunRange aRunRange; // the runRange got from filename
        Int_t aVersion, aSubVersion; // the version and subVersion got from filename

	TIter iter(gDirectory->GetListOfKeys());
	TKey* key;

	AliCDBId* result = new AliCDBId();
	result->SetPath(query.GetPath());

	if (!query.HasVersion()) { // neither version and subversion specified -> look for highest version and subVersion

                while ((key = (TKey*) iter.Next())) { // loop on the keys
			
			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
                        // aRunRange, aVersion, aSubVersion filled from filename

			if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue;
			// aRunRange contains requested run!

			if (result->GetVersion() < aVersion) {
				result->SetVersion(aVersion);
				result->SetSubVersion(aSubVersion);

				result->SetFirstRun(
					aRunRange.GetFirstRun());
				result->SetLastRun(
					aRunRange.GetLastRun());

			} else if (result->GetVersion() == aVersion
				&& result->GetSubVersion()
					< aSubVersion) {

				result->SetSubVersion(aSubVersion);

                        	result->SetFirstRun(
					aRunRange.GetFirstRun());
                        	result->SetLastRun(
					aRunRange.GetLastRun());
			} else if (result->GetVersion() == aVersion
				&& result->GetSubVersion() == aSubVersion){
              			AliError(Form("More than one object valid for run %d, version %d_%d!",
		       			query.GetFirstRun(), aVersion, aSubVersion));
				delete result;
				return NULL;
				}
		}

	} else if (!query.HasSubVersion()) { // version specified but not subversion -> look for highest subVersion

		result->SetVersion(query.GetVersion());

                while ((key = (TKey*) iter.Next())) { // loop on the keys

			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
                        // aRunRange, aVersion, aSubVersion filled from filename

                       if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue; 
			// aRunRange contains requested run!
			
			if(query.GetVersion() != aVersion) continue;
			// aVersion is requested version!

	 		if(result->GetSubVersion() == aSubVersion){
              			AliError(Form("More than one object valid for run %d, version %d_%d!",
		       			query.GetFirstRun(), aVersion, aSubVersion));
	     			delete result;
	     			return NULL;
	 		}
			if( result->GetSubVersion() < aSubVersion) {

                                result->SetSubVersion(aSubVersion);

                                result->SetFirstRun(
					aRunRange.GetFirstRun());
                                result->SetLastRun(
					aRunRange.GetLastRun());
			} 
                }

	} else { // both version and subversion specified

                while ((key = (TKey*) iter.Next())) { // loop on the keys
			
			if (!KeyNameToId(key->GetName(), aRunRange, aVersion, aSubVersion)) continue;
                        // aRunRange, aVersion, aSubVersion filled from filename

			if (!aRunRange.Comprises(query.GetAliCDBRunRange())) continue;
 			// aRunRange contains requested run!

			if(query.GetVersion() != aVersion || query.GetSubVersion() != aSubVersion) continue;
			// aVersion and aSubVersion are requested version and subVersion!

			if(result->GetVersion() == aVersion && result->GetSubVersion() == aSubVersion){
              			AliError(Form("More than one object valid for run %d, version %d_%d!",
		       			query.GetFirstRun(), aVersion, aSubVersion));
	     			delete result;
	     			return NULL;
			}
			result->SetVersion(aVersion);
		        result->SetSubVersion(aSubVersion);
			result->SetFirstRun(aRunRange.GetFirstRun());
			result->SetLastRun(aRunRange.GetLastRun());
			
		}
	}

	return result;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBDump::GetEntry(const AliCDBId& queryId) {
// get AliCDBEntry from the database

	TDirectory::TContext context(gDirectory, fFile);

	if (!(fFile && fFile->IsOpen())) {
                AliError("AliCDBDump storage is not initialized properly");
                return NULL;
        }

        if (!gDirectory->cd(queryId.GetPath())) {
                return NULL;
        }

	AliCDBId *dataId = GetEntryId(queryId);

	if (!dataId || !dataId->IsSpecified()) {
		if(dataId) delete dataId;
		return NULL;
	}

	TString keyname;
	if (!IdToKeyName(dataId->GetAliCDBRunRange(), dataId->GetVersion(),
		dataId->GetSubVersion(), keyname)) {
		AliDebug(2,Form("Bad ID encountered! Subnormal error!"));
		delete dataId;
		return NULL;
	}

	// get the only AliCDBEntry object from the file
	// the object in the file is an AliCDBEntry entry named keyname
	// keyName = Run#firstRun_#lastRun_v#version_s#subVersion

	TObject* anObject = gDirectory->Get(keyname);
	if (!anObject) {
		AliDebug(2,Form("Bad storage data: NULL entry object!"));
		delete dataId;
		return NULL;
	} 

	if (AliCDBEntry::Class() != anObject->IsA()) {
		AliDebug(2,Form("Bad storage data: Invalid entry object!"));
		delete dataId;
		return NULL;
	}

	((AliCDBEntry*) anObject)->SetLastStorage("dump");

	delete dataId;
	return (AliCDBEntry*) anObject;
}

//_____________________________________________________________________________
AliCDBId* AliCDBDump::GetEntryId(const AliCDBId& queryId) {
// get AliCDBEntry from the database

	TDirectory::TContext context(gDirectory, fFile);

	if (!(fFile && fFile->IsOpen())) {
                AliError("AliCDBDump storage is not initialized properly");
                return NULL;
        }

        if (!gDirectory->cd(queryId.GetPath())) {
                return NULL;
        }

	AliCDBId* dataId = 0;

	// look for a filename matching query requests (path, runRange, version, subVersion)
	if (!queryId.HasVersion()) {
		// if version is not specified, first check the selection criteria list
		AliCDBId selectedId(queryId);
		GetSelection(&selectedId);
		dataId = GetId(queryId);
	} else {
		dataId = GetId(queryId);
	}

	if (dataId && !dataId->IsSpecified()) {
		delete dataId;
		return NULL;
	}

	return dataId;
}

//_____________________________________________________________________________
void AliCDBDump::GetEntriesForLevel0(const AliCDBId& queryId, TList* result) {
// multiple request (AliCDBStorage::GetAll)

	TDirectory* saveDir = gDirectory;

	TIter iter(gDirectory->GetListOfKeys());
	TKey* key;

	while ((key = (TKey*) iter.Next())) {

		TString keyNameStr(key->GetName());
		if (queryId.GetAliCDBPath().Level1Comprises(keyNameStr)) {
			gDirectory->cd(keyNameStr);
			GetEntriesForLevel1(queryId, result);
			
			saveDir->cd();
		}
	}
}

//_____________________________________________________________________________
void AliCDBDump::GetEntriesForLevel1(const AliCDBId& queryId, TList* result) {
// multiple request (AliCDBStorage::GetAll)

        TIter iter(gDirectory->GetListOfKeys());
        TKey* key;

	TDirectory* level0Dir = (TDirectory*) gDirectory->GetMother();

        while ((key = (TKey*) iter.Next())) {

                TString keyNameStr(key->GetName());
                if (queryId.GetAliCDBPath().Level2Comprises(keyNameStr)) {
			
			AliCDBPath aPath(level0Dir->GetName(), 
					gDirectory->GetName(), keyNameStr);
			AliCDBId anId(aPath, queryId.GetAliCDBRunRange(), 
					queryId.GetVersion(), -1);

			AliCDBEntry* anEntry = GetEntry(anId);
			if (anEntry) {
				result->Add(anEntry);
			}
	
                }
        }

}

                        
//_____________________________________________________________________________
TList* AliCDBDump::GetEntries(const AliCDBId& queryId) {
// multiple request (AliCDBStorage::GetAll)

	TDirectory::TContext context(gDirectory, fFile);

	if (!(fFile && fFile->IsOpen())) {
                AliError("AliCDBDump storage is not initialized properly");
                return NULL;
        }

	TList* result = new TList();
	result->SetOwner();

	TIter iter(gDirectory->GetListOfKeys());
	TKey* key;

	while ((key = (TKey*) iter.Next())) {
		
		TString keyNameStr(key->GetName());
		if (queryId.GetAliCDBPath().Level0Comprises(keyNameStr)) {
			gDirectory->cd(keyNameStr);
			GetEntriesForLevel0(queryId, result);

			fFile->cd();
		}
	}

	return result;
}

//_____________________________________________________________________________
Bool_t AliCDBDump::PutEntry(AliCDBEntry* entry) {
// put an AliCDBEntry object into the database

	TDirectory::TContext context(gDirectory, fFile);

	if (!(fFile && fFile->IsOpen())) {
		AliError("AliCDBDump storage is not initialized properly");
		return kFALSE;
	}

	if (fReadOnly) {
		AliError("AliCDBDump storage is read only!");
		return kFALSE;
	}

	AliCDBId& id = entry->GetId();
	
        if (!gDirectory->cd(id.GetPath())) {
                if (!MkDir(id.GetPath())) {
                        AliError(Form("Can't open directory <%s>!", 
					id.GetPath().Data()));
                        return kFALSE;
                }
        }

	// set version and subVersion for the entry to be stored
	if (!PrepareId(id)) {
		return kFALSE;		
	}

	// build keyname from entry's id
	TString keyname;
	if (!IdToKeyName(id.GetAliCDBRunRange(), id.GetVersion(), id.GetSubVersion(), keyname)) {
		AliError("Invalid ID encountered! Subnormal error!");
		return kFALSE;
	}	

	// write object (key name: Run#firstRun_#lastRun_v#version_s#subVersion)
	Bool_t result = gDirectory->WriteTObject(entry, keyname);
	if (!result) {
		AliError(Form("Can't write entry to file: %s",
				fFile->GetName()));
	}

        if(result) {
    		AliInfo(Form("CDB object stored into file %s",fFile->GetName()));
    		AliInfo(Form("TDirectory/key name: %s/%s",id.GetPath().Data(),keyname.Data()));
        }

	return result;
}
//_____________________________________________________________________________
TList* AliCDBDump::GetIdListFromFile(const char* fileName){

	TString turl(fileName);
	if (turl[0] != '/') {
		turl.Prepend(TString(gSystem->WorkingDirectory()) + '/');
	}
	TFile *file = TFile::Open(turl);
	if (!file) {
		AliError(Form("Can't open selection file <%s>!", turl.Data()));
		return NULL;
	}
	file->cd();

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
Bool_t AliCDBDump::Contains(const char* path) const{
// check for path in storage

	TDirectory::TContext context(gDirectory, fFile);
	if (!(fFile && fFile->IsOpen())) {
                AliError("AliCDBDump storage is not initialized properly");
                return kFALSE;
        }
	
	return gDirectory->cd(path);

}

//_____________________________________________________________________________
void AliCDBDump::QueryValidFiles()
{
// Query the CDB for files valid for AliCDBStorage::fRun
// fills list fValidFileIds with AliCDBId objects created from file name

	AliError("Not yet (and maybe never) implemented");
}

//_____________________________________________________________________________
Bool_t AliCDBDump::IdToFilename(const AliCDBId& /*id*/, TString& /*filename*/) const {
// build file name from AliCDBId (path, run range, version) and fDBFolder

	AliError("Not implemented");
        return kFALSE;
}

//_____________________________________________________________________________
Int_t AliCDBDump::GetLatestVersion(const char* path, Int_t run){
// get last version found in the database valid for run and path

	AliCDBPath aCDBPath(path);
	if(!aCDBPath.IsValid() || aCDBPath.IsWildcard()) {
		AliError(Form("Invalid path in request: %s", path));
		return -1;
	}

	AliCDBId query(path, run, run, -1, -1);
	AliCDBId *dataId = GetId(query);

	if(!dataId) return -1;
	Int_t version = dataId->GetVersion();

	delete dataId;
	return version;
}

//_____________________________________________________________________________
Int_t AliCDBDump::GetLatestSubVersion(const char* path, Int_t run, Int_t version){
// get last version found in the database valid for run and path

	AliCDBPath aCDBPath(path);
	if(!aCDBPath.IsValid() || aCDBPath.IsWildcard()) {
		AliError(Form("Invalid path in request: %s", path));
		return -1;
	}

	AliCDBId query(path, run, run, version, -1);
	AliCDBId *dataId = GetId(query);

	if(!dataId) return -1;

	Int_t subVersion = dataId->GetSubVersion();

	delete dataId;
	return subVersion;

}

//_____________________________________________________________________________
void AliCDBDump::SetRetry(Int_t /* nretry */, Int_t /* initsec */) {

	// Function to set the exponential retry for putting entries in the OCDB

	AliInfo("This function sets the exponential retry for putting entries in the OCDB - to be used ONLY for AliCDBGrid --> returning without doing anything");
	return;
} 

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBDump factory  			                                                       //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBDumpFactory)

//_____________________________________________________________________________
Bool_t AliCDBDumpFactory::Validate(const char* dbString) {
// check if the string is valid dump URI

        TRegexp dbPattern("^dump://.+$");

        return TString(dbString).Contains(dbPattern);
}

//_____________________________________________________________________________
AliCDBParam* AliCDBDumpFactory::CreateParameter(const char* dbString) {
// create AliCDBDumpParam class from the URI string

        if (!Validate(dbString)) {
                return NULL;
        }

        TString pathname(dbString + sizeof("dump://") - 1);

	Bool_t readOnly;

	if (pathname.Contains(TRegexp(";ReadOnly$"))) {
		pathname.Resize(pathname.Length() - sizeof(";ReadOnly") + 1);
		readOnly = kTRUE;
	} else {
		readOnly = kFALSE;
	}

	gSystem->ExpandPathName(pathname);
	
	if (pathname[0] != '/') {
		pathname.Prepend(TString(gSystem->WorkingDirectory()) + '/');
	}

        return new AliCDBDumpParam(pathname, readOnly);
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBDumpFactory::Create(const AliCDBParam* param) {
// create AliCDBDump storage instance from parameters

        if (AliCDBDumpParam::Class() == param->IsA()) {

                const AliCDBDumpParam* dumpParam = 
			(const AliCDBDumpParam*) param;

                return new AliCDBDump(dumpParam->GetPath(),
				dumpParam->IsReadOnly());
        }

        return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// AliCDBDump parameter class  			                                               //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliCDBDumpParam)

//_____________________________________________________________________________
AliCDBDumpParam::AliCDBDumpParam():
fDBPath(), fReadOnly(kFALSE)
{
// default constructor

}

//_____________________________________________________________________________
AliCDBDumpParam::AliCDBDumpParam(const char* dbPath, Bool_t readOnly):
        fDBPath(dbPath), fReadOnly(readOnly)
{
// constructor

	TString uri;
	uri += "dump://";
	uri += dbPath;

	if (fReadOnly) {
		uri += ";ReadOnly";
	}
	
        SetURI(uri);
        SetType("dump");
}

//_____________________________________________________________________________
AliCDBDumpParam::~AliCDBDumpParam() {
// destructor

}

//_____________________________________________________________________________
AliCDBParam* AliCDBDumpParam::CloneParam() const {
// clone parameter

	return new AliCDBDumpParam(fDBPath, fReadOnly);
}

//_____________________________________________________________________________
ULong_t AliCDBDumpParam::Hash() const {
// return Hash function

	return fDBPath.Hash();
}

//_____________________________________________________________________________
Bool_t AliCDBDumpParam::IsEqual(const TObject* obj) const {
// check if this object is equal to AliCDBParam obj
	
	if (this == obj) {
		return kTRUE;
	}

	if (AliCDBDumpParam::Class() != obj->IsA()) {
		return kFALSE;
	}

	AliCDBDumpParam* other = (AliCDBDumpParam*) obj;

	return fDBPath == other->fDBPath;
}
