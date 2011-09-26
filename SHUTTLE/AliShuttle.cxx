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

/* $Id$ */

//
// This class is the main manager for AliShuttle. 
// It organizes the data retrieval from DCS and call the 
// interface methods of AliPreprocessor.
// For every detector in AliShuttleConfgi (see AliShuttleConfig),
// data for its set of aliases is retrieved. If there is registered
// AliPreprocessor for this detector then it will be used
// accroding to the schema (see AliPreprocessor).
// If there isn't registered AliPreprocessor than the retrieved
// data is stored automatically to the undelying AliCDBStorage.
// For detSpec is used the alias name.
//

#include "AliShuttle.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBRunRange.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliShuttleConfig.h"
#include "DCSClient/AliDCSClient.h"
#include "AliLog.h"
#include "AliPreprocessor.h"
#include "AliShuttleStatus.h"
#include "AliShuttleLogbookEntry.h"

#include <TSystem.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMutex.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMap.h>

#include <TMonaLisaWriter.h>

#include <fstream>

#include <sys/types.h>
#include <sys/wait.h>

#include <signal.h>

ClassImp(AliShuttle)

//______________________________________________________________________________________________
AliShuttle::AliShuttle(const AliShuttleConfig* config,
		UInt_t timeout, Int_t retries):
fConfig(config),
fTimeout(timeout), fRetries(retries),
fPreprocessorMap(),
fLogbookEntry(0),
fCurrentDetector(),
fFirstProcessing(0),
fFXSError(-1),
fStatusEntry(0),
fMonitoringMutex(0),
fLastActionTime(0),
fLastAction(),
fMonaLisa(0),
fTestMode(kNone),
fReadTestMode(kFALSE),
fOutputRedirected(kFALSE)
{
	//
	// config: AliShuttleConfig used
	// timeout: timeout used for AliDCSClient connection
	// retries: the number of retries in case of connection error.
	//

	if (!fConfig->IsValid()) AliFatal("********** !!!!! Invalid configuration !!!!! **********");
	for(int iSys=0;iSys<5;iSys++) {
		fServer[iSys]=0;
		if (iSys < 4)
			fFXSlist[iSys].SetOwner(kTRUE);
	}
	fPreprocessorMap.SetOwner(kTRUE);

	for (UInt_t iDet=0; iDet<NDetectors(); iDet++)
		fFirstUnprocessed[iDet] = kFALSE;

	fMonitoringMutex = new TMutex();
}

//______________________________________________________________________________________________
AliShuttle::~AliShuttle()
{
	//
	// destructor
	//

	fPreprocessorMap.DeleteAll();
	for(int iSys=0;iSys<5;iSys++)
		if(fServer[iSys]) {
			fServer[iSys]->Close();
			delete fServer[iSys];
		        fServer[iSys] = 0;
		}

	if (fStatusEntry){
		delete fStatusEntry;
		fStatusEntry = 0;
	}
	
	if (fMonitoringMutex) 
	{
		delete fMonitoringMutex;
		fMonitoringMutex = 0;
	}
}

//______________________________________________________________________________________________
void AliShuttle::RegisterPreprocessor(AliPreprocessor* preprocessor)
{
	//
	// Registers new AliPreprocessor.
	// It uses GetName() for indentificator of the pre processor.
	// The pre processor is registered it there isn't any other
	// with the same identificator (GetName()).
	//

	const char* detName = preprocessor->GetName();
	if(GetDetPos(detName) < 0)
		AliFatal(Form("********** !!!!! Invalid detector name: %s !!!!! **********", detName));

	if (fPreprocessorMap.GetValue(detName)) {
		AliWarning(Form("AliPreprocessor %s is already registered!", detName));
		return;
	}

	fPreprocessorMap.Add(new TObjString(detName), preprocessor);
}
//______________________________________________________________________________________________
Bool_t AliShuttle::Store(const AliCDBPath& path, TObject* object,
		AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
	// Stores a CDB object in the storage for offline reconstruction. Objects that are not needed for
	// offline reconstruction, but should be stored anyway (e.g. for debugging) should NOT be stored
	// using this function. Use StoreReferenceData instead!
	// It calls StoreLocally function which temporarily stores the data locally; when the preprocessor
	// finishes the data are transferred to the main storage (Grid).

	return StoreLocally(fgkLocalCDB, path, object, metaData, validityStart, validityInfinite);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData)
{
	// Stores a CDB object in the storage for reference data. This objects will not be available during
	// offline reconstrunction. Use this function for reference data only!
	// It calls StoreLocally function which temporarily stores the data locally; when the preprocessor
	// finishes the data are transferred to the main storage (Grid).

	return StoreLocally(fgkLocalRefStorage, path, object, metaData);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::StoreLocally(const TString& localUri,
			const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
			Int_t validityStart, Bool_t validityInfinite)
{
	// Store object temporarily in local storage. Parameters are passed by Store and StoreReferenceData functions.
	// when the preprocessor finishes the data are transferred to the main storage (Grid).
	// The parameters are:
	//   1) Uri of the backup storage (Local)
	//   2) the object's path.
	//   3) the object to be stored
	//   4) the metaData to be associated with the object
	//   5) the validity start run number w.r.t. the current run,
	//      if the data is valid only for this run leave the default 0
	//   6) specifies if the calibration data is valid for infinity (this means until updated),
	//      typical for calibration runs, the default is kFALSE
	//
	// returns 0 if fail, 1 otherwise


	if (fTestMode & kErrorStorage)
	{
		Log(fCurrentDetector, "StoreLocally - In TESTMODE - Simulating error while storing locally");
		return kFALSE;
	}
	
	const char* cdbType = (localUri == fgkLocalCDB) ? "CDB" : "Reference";

	Int_t firstRun = GetCurrentRun() - validityStart;
  	if(firstRun < 0) {
		AliWarning("First valid run happens to be less than 0! Setting it to 0.");
		firstRun=0;
  	}

	Int_t lastRun = -1;
	if(validityInfinite) {
		lastRun = AliCDBRunRange::Infinity();
	} else {
		lastRun = GetCurrentRun();
	}

	// Version is set to current run, it will be used later to transfer data to Grid
	AliCDBId id(path, firstRun, lastRun, GetCurrentRun(), -1);

	if(! dynamic_cast<TObjString*> (metaData->GetProperty("RunUsed(TObjString)"))){
		TObjString runUsed = Form("%d", GetCurrentRun());
		metaData->SetProperty("RunUsed(TObjString)", runUsed.Clone());
	}

	Bool_t result = kFALSE;

	if (!(AliCDBManager::Instance()->GetStorage(localUri))) {
		Log("SHUTTLE", Form("StoreLocally - Cannot activate local %s storage", cdbType));
	} else {
	        Int_t logLevel = AliLog::GetGlobalLogLevel();
	        AliLog::SetGlobalLogLevel(AliLog::kError);
		result = AliCDBManager::Instance()->GetStorage(localUri)
					->Put(object, id, metaData);
                AliLog::SetGlobalLogLevel((AliLog::EType_t)logLevel);
	}

	if(!result) {

		Log(fCurrentDetector, Form("StoreLocally - Can't store object <%s>!", id.ToString().Data()));
	}


	return result;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::StoreOCDB()
{
	//
	// Called when preprocessor ends successfully or when previous storage attempt failed (kStoreError status)
	// Calls underlying StoreOCDB(const char*) function twice, for OCDB and Reference storage.
	// Then calls StoreRefFilesToGrid to store reference files. 
	//
	
	UpdateShuttleStatus(AliShuttleStatus::kStoreStarted);
			
	if (fTestMode & kErrorGrid)
	{
		Log("SHUTTLE", "StoreOCDB - In TESTMODE - Simulating error while storing in the Grid");
		Log(fCurrentDetector, "StoreOCDB - In TESTMODE - Simulating error while storing in the Grid");
		return kFALSE;
	}
	
	Log("SHUTTLE","StoreOCDB - Storing OCDB data ...");
	Int_t resultCDB = StoreOCDB(fgkMainCDB);

	Log("SHUTTLE","StoreOCDB - Storing reference data ...");
	Int_t resultRef = StoreOCDB(fgkMainRefStorage);
	
	Log("SHUTTLE","StoreOCDB - Storing reference files ...");
	Bool_t resultRefFiles = CopyFilesToGrid("reference");
	
	Bool_t resultMetadata = kTRUE;
	if(fCurrentDetector == "GRP") 
	{
		Log("SHUTTLE","StoreOCDB - Storing Run Metadata file ...");
		resultMetadata = CopyFilesToGrid("metadata");
	}
	
	Int_t storeResult = 0;
	
	if (resultCDB < 0 || resultRef < 0 || resultRefFiles == kFALSE || resultMetadata == kFALSE)
		storeResult = -1;
	else if (resultCDB > 0 || resultRef > 0)
		storeResult = 1;
		
	if (storeResult < 0)
	{
		Log("SHUTTLE", 
			Form("\t\t\t****** run %d - %s: STORAGE ERROR ******",
				GetCurrentRun(), fCurrentDetector.Data()));
		UpdateShuttleStatus(AliShuttleStatus::kStoreError);
	} 
	else if (storeResult > 0)
	{
		Log("SHUTTLE", 
			Form("\t\t\t****** run %d - %s: STORAGE DELAYED ******",
				GetCurrentRun(), fCurrentDetector.Data()));
		UpdateShuttleStatus(AliShuttleStatus::kStoreDelayed);
	}
	else if (storeResult == 0) 
	{
		Log("SHUTTLE", 
			Form("\t\t\t****** run %d - %s: DONE ******",
				GetCurrentRun(), fCurrentDetector.Data()));
		UpdateShuttleStatus(AliShuttleStatus::kDone);
		UpdateShuttleLogbook(fCurrentDetector, "DONE");
	}

	return (storeResult == 0);
}

//______________________________________________________________________________________________
Int_t AliShuttle::StoreOCDB(const TString& gridURI)
{
	//
	// Called by StoreOCDB(), performs actual storage to the main OCDB and reference storages (Grid)
	//
	// Return code:
	//   -2 initialization error
	//   -1 storage error
	//   0  success
	//   1  storage delayed (e.g. previous unprocessed runs)
	//

	TObjArray* gridIds=0;

	Bool_t result = kTRUE;
	Bool_t delayed = kFALSE;

	const char* type = 0;
	TString localURI;
	if(gridURI == fgkMainCDB) {
		type = "OCDB";
		localURI = fgkLocalCDB;
	} else if(gridURI == fgkMainRefStorage) {
		type = "reference";
		localURI = fgkLocalRefStorage;
	} else {
		AliError(Form("Invalid storage URI: %s", gridURI.Data()));
		return -2;
	}

	AliCDBManager* man = AliCDBManager::Instance();

	AliCDBStorage *gridSto = man->GetStorage(gridURI);
	if(!gridSto) {
		Log("SHUTTLE",
			Form("StoreOCDB - cannot activate main %s storage", type));
		return -2;
	}

	gridIds = gridSto->GetQueryCDBList();

	// get objects previously stored in local CDB
	AliCDBStorage *localSto = man->GetStorage(localURI);
	if(!localSto) {
		Log("SHUTTLE",
			Form("StoreOCDB - cannot activate local %s storage", type));
		return -2;
	}
	AliCDBPath aPath(GetOfflineDetName(fCurrentDetector.Data()),"*","*");
	// Local objects were stored with current run as Grid version!
	TList* localEntries = localSto->GetAll(aPath.GetPath(), GetCurrentRun(), GetCurrentRun());
	localEntries->SetOwner(1);

	// loop on local stored objects
	TIter localIter(localEntries);
	AliCDBEntry *aLocEntry = 0;
	while((aLocEntry = dynamic_cast<AliCDBEntry*> (localIter.Next()))){
		aLocEntry->SetOwner(1);
		AliCDBId aLocId = aLocEntry->GetId();
		aLocEntry->SetVersion(-1);
		aLocEntry->SetSubVersion(-1);
		                                
		Log(fCurrentDetector.Data(), Form("Attempting to store %s", aLocId.ToString().Data()));

		// If local object is valid up to infinity we store it only if it is
		// the first unprocessed run!
		if (aLocId.GetLastRun() == AliCDBRunRange::Infinity() &&
			!fFirstUnprocessed[GetDetPos(fCurrentDetector)])
		{
			Log("SHUTTLE", Form("StoreOCDB - %s: object %s has validity infinite but "
						"there are previous unprocessed runs!",
						fCurrentDetector.Data(), aLocId.GetPath().Data()));
			Log(fCurrentDetector.Data(), Form("StoreOCDB - %s: object %s has validity infinite but "
						"there are previous unprocessed runs!",
						fCurrentDetector.Data(), aLocId.GetPath().Data()));
			delayed = kTRUE;
			continue;
		}

		// loop on Grid valid Id's
		Bool_t store = kTRUE;
		TIter gridIter(gridIds);
		AliCDBId* aGridId = 0;
		while ((aGridId = dynamic_cast<AliCDBId*> (gridIter.Next()))) {
			if (aGridId->GetPath() != aLocId.GetPath()) 
				continue;
			// skip all objects valid up to infinity
			if (aGridId->GetLastRun() == AliCDBRunRange::Infinity()) 
				continue;
			
			// if we get here, it means there's already some more recent object stored on Grid!
			Log(fCurrentDetector.Data(),
				Form("StoreOCDB - A more recent object already exists in %s storage: <%s>",
				type, aGridId->ToString().Data()));
			
			store = kFALSE;
			break;
		}

		Bool_t storeOk = kFALSE;
		if (store)
		{
			Log(fCurrentDetector.Data(), Form("Prechecks succeeded. Ready to store %s", aLocId.ToString().Data()));
			storeOk = gridSto->Put(aLocEntry);
			if (storeOk) {
				Log("SHUTTLE",
				Form("StoreOCDB - Object <%s> successfully put into %s storage",
					aLocId.ToString().Data(), type));
				Log(fCurrentDetector.Data(),
					Form("StoreOCDB - Object <%s> successfully put into %s storage",
					aLocId.ToString().Data(), type));
			} else	{
				Log("SHUTTLE",
					Form("StoreOCDB - Grid %s storage of object <%s> failed",
					type, aLocId.ToString().Data()));
				Log(fCurrentDetector.Data(),
					Form("StoreOCDB - Grid %s storage of object <%s> failed",
					type, aLocId.ToString().Data()));
				result = kFALSE;
			}
		}
		
		if (!store || storeOk) {
			// removing local file...
			TString filename;
			localSto->IdToFilename(aLocId, filename);
			Log("SHUTTLE", Form("StoreOCDB - Removing local file %s", filename.Data()));
			RemoveFile(filename.Data());
		}
	}
	localEntries->Clear();

	Int_t returnCode = 0;
	
	if (result == kFALSE)
		returnCode = -1;
	else if (delayed != kFALSE)
		returnCode =  1;

	Log("SHUTTLE", Form("StoreOCDB - Returning with %d (result = %d, delayed = %d)", returnCode, result, delayed));
	Log(fCurrentDetector.Data(), Form("StoreOCDB - Returning with %d (result = %d, delayed = %d)", returnCode, result, delayed));
	
	return returnCode;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::CleanReferenceStorage(const char* detector)
{
	// clears the directory used to store reference files of a given subdetector
  
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
  	TString localBaseFolder = sto->GetBaseFolder();

  	TString targetDir = GetRefFilePrefix(localBaseFolder.Data(), detector);
	
  	Log("SHUTTLE", Form("CleanReferenceStorage - Cleaning %s", targetDir.Data()));

	TString begin;
	begin.Form("%d_", GetCurrentRun());
	
	TSystemDirectory* baseDir = new TSystemDirectory("/", targetDir);
	if (!baseDir)
		return kTRUE;
		
	TList* dirList = baseDir->GetListOfFiles();
	delete baseDir;
	
	if (!dirList) return kTRUE;
			
	if (dirList->GetEntries() < 3) // to be changed to 4?
	{
		delete dirList;
		return kTRUE;
	}
				
	Int_t nDirs = 0, nDel = 0;
	TIter dirIter(dirList);
	TSystemFile* entry = 0;

	Bool_t success = kTRUE;
	
	while ((entry = dynamic_cast<TSystemFile*> (dirIter.Next())))
	{					
		if (entry->IsDirectory())
			continue;
		
		TString fileName(entry->GetName());
		if (!fileName.BeginsWith(begin))
			continue;
			
		nDirs++;
						
    		// delete file
    		Int_t result = gSystem->Unlink(fileName.Data());
		
		if (result)
		{
			Log("SHUTTLE", Form("CleanReferenceStorage - Could not delete file %s!", fileName.Data()));
			success = kFALSE;
		} else {
			nDel++;
		}
	}

	if(nDirs > 0)
		Log("SHUTTLE", Form("CleanReferenceStorage - %d (over %d) reference files in folder %s were deleted.", 
			nDel, nDirs, targetDir.Data()));

		
	delete dirList;
	return success;






  Int_t result = gSystem->GetPathInfo(targetDir, 0, (Long64_t*) 0, 0, 0);
  if (result == 0)
  {
    // delete directory
    result = gSystem->Exec(Form("rm -rf %s", targetDir.Data()));
    if (result != 0)
    {  
      Log("SHUTTLE", Form("CleanReferenceStorage - Could not clean directory %s", targetDir.Data()));
      return kFALSE;
    }
  }

  result = gSystem->mkdir(targetDir, kTRUE);
  if (result != 0)
  {
    Log("SHUTTLE", Form("CleanReferenceStorage - Error creating base directory %s", targetDir.Data()));
    return kFALSE;
  }
	
  return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::StoreReferenceFile(const char* detector, const char* localFile, const char* gridFileName)
{
	//
	// Stores reference file directly (without opening it). This function stores the file locally.
	//
	// The file is stored under the following location: 
	// <base folder of local reference storage>/<DET>/<RUN#>_<gridFileName>
	// where <gridFileName> is the second parameter given to the function
	// 
	
	if (fTestMode & kErrorStorage)
	{
		Log(fCurrentDetector, "StoreReferenceFile - In TESTMODE - Simulating error while storing locally");
		return kFALSE;
	}
	
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
	
	TString localBaseFolder = sto->GetBaseFolder();
	
	TString target = GetRefFilePrefix(localBaseFolder.Data(), detector);	
	target.Append(Form("/%d_%s", GetCurrentRun(), gridFileName));
	
	return CopyFileLocally(localFile, target);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::StoreRunMetadataFile(const char* localFile, const char* gridFileName)
{
	//
	// Stores Run metadata file to the Grid, in the run folder
	//
	// Only GRP can call this function.
	
	if (fTestMode & kErrorStorage)
	{
		Log(fCurrentDetector, "StoreRunMetaDataFile - In TESTMODE - Simulating error while storing locally");
		return kFALSE;
	}
	
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
	
	TString localBaseFolder = sto->GetBaseFolder();
	
	// Build Run level folder
	// folder = /alice/data/year/lhcPeriod/runNb/raw
	
		
	TString lhcPeriod = GetLHCPeriod();	
	if (lhcPeriod.Length() == 0) 
	{
		Log("SHUTTLE","StoreRunMetaDataFile - LHCPeriod not found in logbook!");
		return 0;
	}
	
	// TODO partitions with one detector only write data into LHCperiod_DET
	TString partition = GetRunParameter("detector");
	
	if (partition.Length() > 0 && partition != "ALICE")
	{
		lhcPeriod.Append(Form("_%s", partition.Data()));
		Log(fCurrentDetector, Form("Run data tags merged file will be written in %s", 
				lhcPeriod.Data()));
	}
		
	TString target = Form("%s/GRP/RunMetadata%s%d/%s/%09d/raw/%s", 
				localBaseFolder.Data(), fConfig->GetAlienPath(), GetCurrentYear(), 
				lhcPeriod.Data(), GetCurrentRun(), gridFileName);
					
	return CopyFileLocally(localFile, target);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::CopyFileLocally(const char* localFile, const TString& target)
{
	//
	// Stores file locally. Called by StoreReferenceFile and StoreRunMetadataFile
	// Files are temporarily stored in the local reference storage. When the preprocessor 
	// finishes, the Shuttle calls CopyFilesToGrid to transfer the files to AliEn 
	// (in reference or run level folders)
	//
	
	TString targetDir(target(0, target.Last('/')));
	
	//try to open base dir folder, if it does not exist
	void* dir = gSystem->OpenDirectory(targetDir.Data());
	if (dir == NULL) {
		if (gSystem->mkdir(targetDir.Data(), kTRUE)) {
			Log("SHUTTLE", Form("CopyFileLocally - Can't open directory <%s>", targetDir.Data()));
			return kFALSE;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}
	
	Int_t result = 0;
	
	result = gSystem->GetPathInfo(localFile, 0, (Long64_t*) 0, 0, 0);
	if (result)
	{
		Log("SHUTTLE", Form("CopyFileLocally - %s does not exist", localFile));
		return kFALSE;
	}

	result = gSystem->GetPathInfo(target, 0, (Long64_t*) 0, 0, 0);
	if (!result)
	{
		Log("SHUTTLE", Form("CopyFileLocally - target file %s already exist, removing...", target.Data()));
		if (gSystem->Unlink(target.Data()))
		{
			Log("SHUTTLE", Form("CopyFileLocally - Could not remove existing target file %s!", target.Data()));
			return kFALSE;
		}
	}	
	
	result = gSystem->CopyFile(localFile, target);

	if (result == 0)
	{
		Log("SHUTTLE", Form("CopyFileLocally - File %s stored locally to %s", localFile, target.Data()));
		return kTRUE;
	}
	else
	{
		Log("SHUTTLE", Form("CopyFileLocally - Could not store file %s to %s! Error code = %d", 
				localFile, target.Data(), result));
		return kFALSE;
	}	



}

//______________________________________________________________________________________________
Bool_t AliShuttle::CopyFilesToGrid(const char* type)
{
	//
	// Transfers local files to the Grid. Local files can be reference files 
	// or run metadata file (from GRP only).
	//
	// According to the type (ref, metadata) the files are stored under the following location: 
	// ref --> <base folder of reference storage>/<DET>/<RUN#>_<gridFileName>
	// metadata --> <run data folder>/<MetadataFileName>
	//
		
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
	if (!sto)
		return kFALSE;
	TString localBaseFolder = sto->GetBaseFolder();
	
	TString dir;
	TString alienDir;
	TString begin;
	
	if (strcmp(type, "reference") == 0) 
	{
		dir = GetRefFilePrefix(localBaseFolder.Data(), fCurrentDetector.Data());
		AliCDBStorage* gridSto = man->GetStorage(fgkMainRefStorage);
		if (!gridSto)
			return kFALSE;
		TString gridBaseFolder = gridSto->GetBaseFolder();
		alienDir = GetRefFilePrefix(gridBaseFolder.Data(), fCurrentDetector.Data());
		begin = Form("%d_", GetCurrentRun());
	} 
	else if (strcmp(type, "metadata") == 0)
	{
			
		TString lhcPeriod = GetLHCPeriod();
	
		if (lhcPeriod.Length() == 0) 
		{
			Log("SHUTTLE","CopyFilesToGrid - LHCPeriod not found in logbook!");
			return 0;
		}
		
		// TODO partitions with one detector only write data into LHCperiod_DET
		TString partition = GetRunParameter("detector");
	
		if (partition.Length() > 0 && partition != "ALICE")
		{
			lhcPeriod.Append(Form("_%s", partition.Data()));
		}
		
		dir = Form("%s/GRP/RunMetadata%s%d/%s/%09d/raw", 
				localBaseFolder.Data(), fConfig->GetAlienPath(), GetCurrentYear(), 
				lhcPeriod.Data(), GetCurrentRun());
		alienDir = dir(dir.Index(fConfig->GetAlienPath()), dir.Length());
		
		begin = "";
	}
	else 
	{
		Log("SHUTTLE", "CopyFilesToGrid - Unexpected: type label must be reference or metadata!");
		return kFALSE;
	}
		
	TSystemDirectory* baseDir = new TSystemDirectory("/", dir);
	if (!baseDir)
		return kTRUE;
		
	TList* dirList = baseDir->GetListOfFiles();
	delete baseDir;
	
	if (!dirList) return kTRUE;
		
	if (dirList->GetEntries() < 3) 
	{
		delete dirList;
		return kTRUE;
	}
			
	if (!gGrid)
	{ 
		Log("SHUTTLE", "CopyFilesToGrid - Connection to Grid failed: Cannot continue!");
		delete dirList;
		return kFALSE;
	}
	
	Int_t nDirs = 0, nTransfer = 0;
	TIter dirIter(dirList);
	TSystemFile* entry = 0;

	Bool_t success = kTRUE;
	Bool_t first = kTRUE;
	
	while ((entry = dynamic_cast<TSystemFile*> (dirIter.Next())))
	{			
		if (entry->IsDirectory())
			continue;
			
		TString fileName(entry->GetName());
		if (!fileName.BeginsWith(begin))
			continue;
			
		nDirs++;
			
		if (first)
		{
			first = kFALSE;
			// check that folder exists, otherwise create it
			TGridResult* result = gGrid->Ls(alienDir.Data(), "a");
			
			if (!result)
			{
				delete dirList;
				return kFALSE;
			}
			
			if (!result->GetFileName(1)) // TODO: It looks like element 0 is always 0!!
			{
				// TODO It does not work currently! Bug in TAliEn::Mkdir
				// TODO Manually fixed in local root v5-16-00
				if (!gGrid->Mkdir(alienDir.Data(),"-p",0))
				{
					Log("SHUTTLE", Form("CopyFilesToGrid - Cannot create directory %s",
							alienDir.Data()));
					delete dirList;
					return kFALSE;
				} else {
					Log("SHUTTLE",Form("CopyFilesToGrid - Folder %s created", alienDir.Data()));
				}
				
			} else {
					Log("SHUTTLE",Form("CopyFilesToGrid - Folder %s found", alienDir.Data()));
			}
		}
			
		TString fullLocalPath;
		fullLocalPath.Form("%s/%s", dir.Data(), fileName.Data());
		
		TString fullGridPath;
		fullGridPath.Form("alien://%s/%s", alienDir.Data(), fileName.Data());

		Bool_t result = TFile::Cp(fullLocalPath, fullGridPath);
		
		if (result)
		{
			Log("SHUTTLE", Form("CopyFilesToGrid - Copying local file %s to %s succeeded!", 
						fullLocalPath.Data(), fullGridPath.Data()));
			RemoveFile(fullLocalPath);
			nTransfer++;
		}
		else
		{
			Log("SHUTTLE", Form("CopyFilesToGrid - Copying local file %s to %s FAILED!", 
						fullLocalPath.Data(), fullGridPath.Data()));
			success = kFALSE;
		}
	}

	Log("SHUTTLE", Form("CopyFilesToGrid - %d (over %d) files in folder %s copied to Grid.", 
						nTransfer, nDirs, dir.Data()));

		
	delete dirList;
	return success;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetRefFilePrefix(const char* base, const char* detector)
{
	//
	// Get folder name of reference files 
	//

	TString offDetStr(GetOfflineDetName(detector));
	static TString dir;
	if (offDetStr == "ITS" || offDetStr == "MUON" || offDetStr == "PHOS")
	{
		dir.Form("%s/%s/%s", base, offDetStr.Data(), detector);
	} else {
		dir.Form("%s/%s", base, offDetStr.Data());
	}
	
	return dir.Data();
}

//______________________________________________________________________________________________
void AliShuttle::CleanLocalStorage(const TString& uri)
{
	//
	// Called in case the preprocessor is declared failed. Remove remaining objects from the local storages.
	//

	const char* type = 0;
	if(uri == fgkLocalCDB) {
		type = "OCDB";
	} else if(uri == fgkLocalRefStorage) {
		type = "Reference";
	} else {
		AliError(Form("Invalid storage URI: %s", uri.Data()));
		return;
	}

	AliCDBManager* man = AliCDBManager::Instance();

	// open local storage
	AliCDBStorage *localSto = man->GetStorage(uri);
	if(!localSto) {
		Log("SHUTTLE",
			Form("CleanLocalStorage - cannot activate local %s storage", type));
		return;
	}

	TString filename(Form("%s/%s/*/Run*_v%d_s*.root",
		localSto->GetBaseFolder().Data(), GetOfflineDetName(fCurrentDetector.Data()), GetCurrentRun()));

	AliDebug(2, Form("filename = %s", filename.Data()));

	Log("SHUTTLE", Form("Removing remaining local files for run %d and detector %s ...",
		GetCurrentRun(), fCurrentDetector.Data()));

	RemoveFile(filename.Data());

}

//______________________________________________________________________________________________
void AliShuttle::RemoveFile(const char* filename)
{
	//
	// removes local file
	//

	TString command(Form("rm -f %s", filename));

	Int_t result = gSystem->Exec(command.Data());
	if(result != 0)
	{
		Log("SHUTTLE", Form("RemoveFile - %s: Cannot remove file %s!",
			fCurrentDetector.Data(), filename));
	}
}

//______________________________________________________________________________________________
AliShuttleStatus* AliShuttle::ReadShuttleStatus()
{
	//
	// Reads the AliShuttleStatus from the CDB
	//

	if (fStatusEntry){
		delete fStatusEntry;
		fStatusEntry = 0;
	}

	Int_t path1 = GetCurrentRun()/10000;
	fStatusEntry = AliCDBManager::Instance()->GetStorage(GetLocalCDB())
		->Get(Form("/SHUTTLE/%s/%d", fCurrentDetector.Data(), path1), GetCurrentRun());

	if (!fStatusEntry) return 0;
	fStatusEntry->SetOwner(1);

	AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (fStatusEntry->GetObject());
	if (!status) {
		AliError("Invalid object stored to CDB!");
		return 0;
	}

	return status;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::WriteShuttleStatus(AliShuttleStatus* status)
{
	//
	// writes the status for one subdetector
	//

	if (fStatusEntry){
		delete fStatusEntry;
		fStatusEntry = 0;
	}

	Int_t run = GetCurrentRun();
	Int_t path1 = run/10000;
	TString path1_string = Form("%d",path1);

	AliCDBId id(AliCDBPath("SHUTTLE", fCurrentDetector, path1_string), run, run);

	fStatusEntry = new AliCDBEntry(status, id, new AliCDBMetaData);
	fStatusEntry->SetOwner(1);

	Int_t logLevel = AliLog::GetGlobalLogLevel();
        AliLog::SetGlobalLogLevel(AliLog::kError);

	UInt_t result = AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);

        AliLog::SetGlobalLogLevel((AliLog::EType_t)logLevel);

 	if (!result) {
		Log("SHUTTLE", Form("WriteShuttleStatus - Failed for %s, run %d",
						fCurrentDetector.Data(), run));
		return kFALSE;
	}
	
	SendMLDetInfo();

	return kTRUE;
}

//______________________________________________________________________________________________
void AliShuttle::UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount)
{
	//
	// changes the AliShuttleStatus for the given detector and run to the given status
	//

	if (!fStatusEntry){
		AliError("UNEXPECTED: fStatusEntry empty");
		return;
	}

	AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (fStatusEntry->GetObject());

	if (!status){
		Log("SHUTTLE", "UpdateShuttleStatus - UNEXPECTED: status could not be read from current CDB entry");
		return;
	}

	TString actionStr = Form("UpdateShuttleStatus - %s: Changing state from %s to %s",
				fCurrentDetector.Data(),
				status->GetStatusName(),
				status->GetStatusName(newStatus));
	Log("SHUTTLE", actionStr);
	SetLastAction(actionStr);

	status->SetStatus(newStatus);
	if (increaseCount) status->IncreaseCount();

	Int_t logLevel = AliLog::GetGlobalLogLevel();
        AliLog::SetGlobalLogLevel(AliLog::kError);

	AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);

        AliLog::SetGlobalLogLevel((AliLog::EType_t)logLevel);

	SendMLDetInfo();
}

//______________________________________________________________________________________________
void AliShuttle::SendMLDetInfo()
{
	//
	// sends ML information about the current status of the current detector being processed
	//
	
	AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (fStatusEntry->GetObject());
	
	if (!status){
		Log("SHUTTLE", "SendMLDetInfo - UNEXPECTED: status could not be read from current CDB entry");
		return;
	}
	
	TMonaLisaText  mlStatus(Form("%s_status", fCurrentDetector.Data()), status->GetStatusName());
	TMonaLisaValue mlRetryCount(Form("%s_count", fCurrentDetector.Data()), status->GetCount());

	TList mlList;
	mlList.Add(&mlStatus);
	mlList.Add(&mlRetryCount);

	TString mlID;
	mlID.Form("%d", GetCurrentRun());
	fMonaLisa->SendParameters(&mlList, mlID);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::ContinueProcessing()
{
	// this function reads the AliShuttleStatus information from CDB and
	// checks if the processing should be continued
	// if yes it returns kTRUE and updates the AliShuttleStatus with nextStatus

	if (!fConfig->HostProcessDetector(fCurrentDetector))
		return kFALSE;

	AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(fCurrentDetector));
	if (!aPreprocessor)
	{
		Log("SHUTTLE", Form("ContinueProcessing - %s: no preprocessor registered", fCurrentDetector.Data()));
		return kFALSE;
	}

	AliShuttleLogbookEntry::Status entryStatus =
		fLogbookEntry->GetDetectorStatus(fCurrentDetector);

	if (entryStatus != AliShuttleLogbookEntry::kUnprocessed) {
		Log("SHUTTLE", Form("ContinueProcessing - %s is %s",
				fCurrentDetector.Data(),
				fLogbookEntry->GetDetectorStatusName(entryStatus)));
		return kFALSE;
	}

	// if we get here, according to Shuttle logbook subdetector is in UNPROCESSED state

	// check if current run is first unprocessed run for current detector
	if (fConfig->StrictRunOrder(fCurrentDetector) &&
		!fFirstUnprocessed[GetDetPos(fCurrentDetector)])
	{
		if (fTestMode == kNone)
		{
			Log("SHUTTLE", Form("ContinueProcessing - %s requires strict run ordering"
					    " but this is not the first unprocessed run!",fCurrentDetector.Data()));
			return kFALSE;
		}
		else
		{
			Log("SHUTTLE", Form("ContinueProcessing - In TESTMODE - "
					"Although %s requires strict run ordering "
					"and this is not the first unprocessed run, "
					    "the SHUTTLE continues",fCurrentDetector.Data()));
		}
	}

	// Is the subdetector processed first time for this run?
	fFirstProcessing = kFALSE;

	AliShuttleStatus* status = ReadShuttleStatus();
	if (!status) {
		// first time
		Log("SHUTTLE", Form("ContinueProcessing - %s: Processing first time",
				fCurrentDetector.Data()));
		status = new AliShuttleStatus(AliShuttleStatus::kStarted);
		fFirstProcessing = kTRUE;
		return WriteShuttleStatus(status);
	}

	// The following case shouldn't happen if Shuttle Logbook was correctly updated.
	// If it happens it may mean Logbook updating failed... let's do it now!
	if (status->GetStatus() == AliShuttleStatus::kDone ||
	    status->GetStatus() == AliShuttleStatus::kFailed ||
	    status->GetStatus() == AliShuttleStatus::kSkipped) {
		Log("SHUTTLE", Form("ContinueProcessing - %s is already %s. Updating Shuttle Logbook",
					fCurrentDetector.Data(),
					status->GetStatusName(status->GetStatus())));
		
		if (status->GetStatus() == AliShuttleStatus::kSkipped)
		{
			UpdateShuttleLogbook(fCurrentDetector.Data(), "DONE");
		}
		else
			UpdateShuttleLogbook(fCurrentDetector.Data(), status->GetStatusName(status->GetStatus()));
			
		return kFALSE;
	}

	if (status->GetStatus() == AliShuttleStatus::kStoreStarted || status->GetStatus() == AliShuttleStatus::kStoreDelayed ||status->GetStatus() == AliShuttleStatus::kStoreError) {
		Log("SHUTTLE",
			Form("ContinueProcessing - %s: Grid storage of one or more "
				"objects failed. Trying again now",
				fCurrentDetector.Data()));
		StoreOCDB();
		return kFALSE;
	}

	// if we get here, there is a restart
	Bool_t cont = kFALSE;

	// abort conditions
	if (status->GetCount() >= fConfig->GetMaxRetries()) {
		Log("SHUTTLE", Form("ContinueProcessing - %s failed %d times in status %s - "
				"Updating Shuttle Logbook", fCurrentDetector.Data(),
				status->GetCount(), status->GetStatusName()));
		UpdateShuttleLogbook(fCurrentDetector.Data(), "FAILED");
		UpdateShuttleStatus(AliShuttleStatus::kFailed);

		// there may still be objects in local OCDB and reference storage
		// and FXS databases may be not updated: do it now!
		
		// TODO Currently disabled, we want to keep files in case of failure!
		// CleanLocalStorage(fgkLocalCDB);
		// CleanLocalStorage(fgkLocalRefStorage);
		// UpdateTableFailCase();
		
		// Send mail to detector expert!
		Log("SHUTTLE", Form("ContinueProcessing - Sending mail to %s expert...", 
				    fCurrentDetector.Data()));
		// det experts in to
		TString to="";
		TIter *iterExperts = 0;
		iterExperts = new TIter(fConfig->GetResponsibles(fCurrentDetector));
		TObjString *anExpert=0;
		while ((anExpert = (TObjString*) iterExperts->Next()))
			{
				to += Form("%s, \n", anExpert->GetName());
			}
		delete iterExperts;
		
		if (to.Length() > 0)
			to.Remove(to.Length()-3);
		AliDebug(2, Form("to: %s",to.Data()));

		if (to.IsNull()) {
			Log("SHUTTLE", Form("List of %s responsibles not set!", fCurrentDetector.Data()));
			return kFALSE;
		}

		Log(fCurrentDetector.Data(), Form("ContinueProcessing - Sending mail to %s expert(s):", 
				    fCurrentDetector.Data()));
		Log(fCurrentDetector.Data(), Form("\n%s", to.Data()));
		if (!SendMail(kPPEMail))
			Log("SHUTTLE", Form("ContinueProcessing - Could not send mail to %s expert",
					    fCurrentDetector.Data()));

	} else {
		Log("SHUTTLE", Form("ContinueProcessing - %s: restarting. "
				"Aborted before with %s. Retry number %d.", fCurrentDetector.Data(),
				status->GetStatusName(), status->GetCount()));
		Bool_t increaseCount = kTRUE;
		if (status->GetStatus() == AliShuttleStatus::kDCSError || 
		    status->GetStatus() == AliShuttleStatus::kDCSStarted ||
		    status->GetStatus() == AliShuttleStatus::kFXSError ||
		    status->GetStatus() == AliShuttleStatus::kOCDBError)
				increaseCount = kFALSE;
				
		UpdateShuttleStatus(AliShuttleStatus::kStarted, increaseCount);
		cont = kTRUE;
	}

	return cont;
}

//______________________________________________________________________________________________
void AliShuttle::SendMLRunInfo(const char* status)
{
	// 
	// Send information about this run to ML
	
	TMonaLisaText  mlStatus("SHUTTLE_status", status);
	TString runType(fLogbookEntry->GetRunType());
	if (strlen(fLogbookEntry->GetRunParameter("log")) > 0){

		runType += "(";
		runType += fLogbookEntry->GetRunParameter("log");
		runType += ")";
	}
	if (fLogbookEntry->GetDATestMode()){
		runType += " (DATest)";
	}
	TMonaLisaText  mlRunType("SHUTTLE_runtype", runType);

	TList mlList;
	mlList.Add(&mlStatus);
	mlList.Add(&mlRunType);

	TString mlID;
	mlID.Form("%d", GetCurrentRun());
	fMonaLisa->SendParameters(&mlList, mlID);	
}

//______________________________________________________________________________________________
Int_t AliShuttle::GetMem(Int_t pid)
{
	// invokes ps to get the memory consumption of the process <pid>
	// returns -1 in case of error
	
	TString checkStr;
	checkStr.Form("ps -o vsize --pid %d | tail -n 1", pid);
	FILE* pipe = gSystem->OpenPipe(checkStr, "r");
	if (!pipe)
	{
		Log("SHUTTLE", Form("Process - Error: "
			"Could not open pipe to %s", checkStr.Data()));
		return -1;
	}
		
	char buffer[100];
	if (!fgets(buffer, 100, pipe))
	{
		Log("SHUTTLE", "Process - Error: ps did not return anything");
		gSystem->ClosePipe(pipe);
		return -1;
	}
	gSystem->ClosePipe(pipe);
	
	//Log("SHUTTLE", Form("ps returned %s", buffer));
	
	Int_t mem = 0;
	if ((sscanf(buffer, "%d\n", &mem) != 1) || !mem)
	{
		Log("SHUTTLE", "Process - Error: Could not parse output of ps");
		return -1;
	}
	
	return mem;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Process(AliShuttleLogbookEntry* entry)
{
	//
	// Makes data retrieval for all detectors in the configuration.
	// entry: Shuttle logbook entry, contains run paramenters and status of detectors
	// (Unprocessed, Inactive, Failed or Done).
	// Returns kFALSE in case of error occured and kTRUE otherwise
	//

	if (!entry) return kFALSE;

	fLogbookEntry = entry;

	Log("SHUTTLE", Form("\t\t\t^*^*^*^*^*^*^*^*^*^*^*^* run %d: START ^*^*^*^*^*^*^*^*^*^*^*^*",
					GetCurrentRun()));

	CountOpenRuns();
	
	// Send the information to ML
	SendMLRunInfo("Processing");

	if (fLogbookEntry->IsDone())
	{
		Log("SHUTTLE","Process - Shuttle is already DONE. Updating logbook");
		UpdateShuttleLogbook("shuttle_done");
		fLogbookEntry = 0;
		return kTRUE;
	}

	// read test mode if flag is set
	if (fReadTestMode)
	{
		fTestMode = kNone;
		TString logEntry(entry->GetRunParameter("log"));
		//printf("log entry = %s\n", logEntry.Data());
		TString searchStr("Testmode: ");
		Int_t pos = logEntry.Index(searchStr.Data());
		//printf("%d\n", pos);
		if (pos >= 0)
		{
			TSubString subStr = logEntry(pos + searchStr.Length(), logEntry.Length());
			//printf("%s\n", subStr.String().Data());
			TString newStr(subStr.Data());
			TObjArray* token = newStr.Tokenize(' ');
			if (token)
			{
				//token->Print();
				TObjString* tmpStr = dynamic_cast<TObjString*> (token->First());
				if (tmpStr)
				{
					Int_t testMode = tmpStr->String().Atoi();
					if (testMode > 0)
					{
						Log("SHUTTLE", Form("Process - Enabling test mode %d", testMode));
						SetTestMode((TestMode) testMode);
					}
				}
				delete token;	       
			}
		}
	}
		
	fLogbookEntry->Print("all");

	// Initialization
	Bool_t hasError = kFALSE;

	// Set the CDB and Reference folders according to the year

	// build cdb paths (repeat each time, run might be a DATest run)
	if (!fLogbookEntry->GetDATestMode()){
		fgkMainCDB.Form("alien://folder=%s%d/OCDB?user=alidaq?cacheFold=/tmp/OCDBCache", 
				fConfig->GetAlienPath(), GetCurrentYear());
		
		fgkMainRefStorage.Form("alien://folder=%s%d/Reference?user=alidaq?cacheFold=/tmp/OCDBCache", 
				       fConfig->GetAlienPath(), GetCurrentYear());
	}
	else {
		fgkMainCDB.Form("alien://folder=%s%d/DATest/OCDB?user=alidaq?cacheFold=/tmp/OCDBCache",
				fConfig->GetAlienPath(), GetCurrentYear());
		
		fgkMainRefStorage.Form("alien://folder=%s%d/DATest/Reference?user=alidaq?cacheFold=/tmp/OCDBCache",
				       fConfig->GetAlienPath(), GetCurrentYear());
	}

	AliDebug(2,Form("Main CDB storage = %s",fgkMainCDB.Data()));
	AliDebug(2,Form("Main Reference storage = %s",fgkMainRefStorage.Data()));

	// Loop on detectors in the configuration
	TIter iter(fConfig->GetDetectors());
	TObjString* aDetector = 0;

	Bool_t first = kTRUE;

	while ((aDetector = (TObjString*) iter.Next()))
	{
		fCurrentDetector = aDetector->String();

		if (ContinueProcessing() == kFALSE) 
			continue;
		
		if (first)
		{
		  // only read QueryCDB when needed and only once
		  AliCDBStorage *mainCDBSto = AliCDBManager::Instance()->GetStorage(fgkMainCDB);
		  if(mainCDBSto) mainCDBSto->QueryCDB(GetCurrentRun());
		  AliCDBStorage *mainRefSto = AliCDBManager::Instance()->GetStorage(fgkMainRefStorage);
		  if(mainRefSto) mainRefSto->QueryCDB(GetCurrentRun());
		  first = kFALSE;
		}

		Log("SHUTTLE", Form("\t\t\t****** run %d - %s: START  ******",
						GetCurrentRun(), aDetector->GetName()));

		for(Int_t iSys=0;iSys<3;iSys++) fFXSCalled[iSys]=kFALSE;
		
		Int_t initialMem = GetMem(getpid());
		Log("SHUTTLE", Form("Memory consumption before forking is %d", initialMem));

		Log(fCurrentDetector.Data(), "Process - Starting processing");

		Int_t pid = fork();

		if (pid < 0)
		{
			Log("SHUTTLE", "Process - ERROR: Forking failed");
		}
		else if (pid > 0)
		{
			// parent
			Log("SHUTTLE", Form("Process - In parent process of %d - %s: Starting monitoring",
							GetCurrentRun(), aDetector->GetName()));

			Long_t begin = time(0);

			int status; // to be used with waitpid, on purpose an int (not Int_t)!
			while (waitpid(pid, &status, WNOHANG) == 0)
			{
				Long_t expiredTime = time(0) - begin;

				if (expiredTime > fConfig->GetPPTimeOut())
				{
                                        TString logMsg;
					AliShuttleStatus *currentStatus = ReadShuttleStatus();
					AliShuttleStatus::Status newStatus = AliShuttleStatus::kInvalid;
					
					if (currentStatus->GetStatus() == AliShuttleStatus::kDCSStarted)
					{
						// in case the pp goes in TimeOut while retrieving the DCS DPs
						// set status to kDCSError
						
						logMsg.Form("Process - Process of %s timed out while retrieving the DCS DataPoints. Run time: %ld seconds. Killing... and setting status to DCSError.",
								fCurrentDetector.Data(), expiredTime);
						newStatus = AliShuttleStatus::kDCSError;
					}
					else if (currentStatus->GetStatus() <= AliShuttleStatus::kPPDone)
					{
						// in case pp not yet done set status to kPPTimeOut
					
						logMsg.Form("Process - Process of %s timed out. Run time: %ld seconds. Killing...",
								fCurrentDetector.Data(), expiredTime);
						newStatus = AliShuttleStatus::kPPTimeOut;
					}
					else if (currentStatus->GetStatus() == AliShuttleStatus::kStoreStarted)
					{
						// in case the pp goes in TimeOut while storing the objects in the OCDB
						// set status to kStoreError
						
						logMsg.Form("Process - Process of %s timed out while storing the OCDB object. Run time: %ld seconds. Killing... and setting status to StoreError.",
								fCurrentDetector.Data(), expiredTime);
						newStatus = AliShuttleStatus::kStoreError;
					}
					else 
					{
						// in other cases don't change the status
						
					       	logMsg.Form("Process - Process of %s timed out in status = %s. Run time: %ld seconds. Killing... without changing the status",
								fCurrentDetector.Data(), currentStatus->GetStatusName(), expiredTime);
					}
				
					Log("SHUTTLE", logMsg);
					Log(fCurrentDetector, logMsg);

					kill(pid, 9);

					if (newStatus != AliShuttleStatus::kInvalid)
						UpdateShuttleStatus(newStatus);
					hasError = kTRUE;

					gSystem->Sleep(1000);
				}
				else
				{
					gSystem->Sleep(1000);
					
					Int_t mem = GetMem(pid);

					if (mem < 0)
						continue;
						
					mem -= initialMem;
					if (mem < 0)
						mem = 0;
					
					if (expiredTime % 60 == 0)
					{
						Log("SHUTTLE", Form("Process - %s: Checking process. "
							"Run time: %ld seconds - Memory consumption: %d KB",
							fCurrentDetector.Data(), expiredTime, mem));
						SendAlive();
					}
					
					if (mem > fConfig->GetPPMaxMem())
					{
						TString tmp;
						tmp.Form("Process - Process exceeds maximum allowed memory "
							"(%d KB > %d KB). Killing...",
							mem, fConfig->GetPPMaxMem());
						Log("SHUTTLE", tmp);
						Log(fCurrentDetector, tmp);
	
						kill(pid, 9);
	
	    					UpdateShuttleStatus(AliShuttleStatus::kPPOutOfMemory);
						hasError = kTRUE;
	
						gSystem->Sleep(1000);
					}
				}
			}

			Log("SHUTTLE", Form("Process - In parent process of %d - %s: Client has terminated.",
								GetCurrentRun(), aDetector->GetName()));

			if (WIFEXITED(status))
			{
				Int_t returnCode = WEXITSTATUS(status);

				Log("SHUTTLE", Form("Process - %s: the return code is %d", fCurrentDetector.Data(),
										returnCode));

				if (returnCode == 0) hasError = kTRUE;
			}
		}
		else if (pid == 0)
		{
			// child
			Log("SHUTTLE", Form("Process - In child process of %d - %s", GetCurrentRun(),
				aDetector->GetName()));

			Log("SHUTTLE", Form("Process - Redirecting output to %s log",fCurrentDetector.Data()));

			if ((freopen(GetLogFileName(fCurrentDetector), "a", stdout)) == 0)
			{
				Log("SHUTTLE", "Process - Could not freopen stdout");
			}
			else
			{
				fOutputRedirected = kTRUE;
				if ((dup2(fileno(stdout), fileno(stderr))) < 0)
					Log("SHUTTLE", "Process - Could not redirect stderr");
				
			}

			Log("SHUTTLE", "Executing TGrid::Connect");
			TGrid::Connect("alien://");
			
			TString wd = gSystem->WorkingDirectory();
			Int_t dir_lev1 = GetCurrentRun()/10000;
			TString tmpDir = Form("%s/%d/%d/%s_process", GetShuttleTempDir(), 
				dir_lev1, GetCurrentRun(), fCurrentDetector.Data());
			
			Int_t result = gSystem->GetPathInfo(tmpDir.Data(), 0, (Long64_t*) 0, 0, 0);
			if (!result) // temp dir already exists!
			{
				Log(fCurrentDetector.Data(), 
					Form("Process - %s dir already exists! Removing...", tmpDir.Data()));
				gSystem->Exec(Form("rm -rf %s",tmpDir.Data()));		
			} 
			
			if (gSystem->mkdir(tmpDir.Data(), 1))
			{
				Log(fCurrentDetector.Data(), "Process - could not make temp directory!!");
				gSystem->Exit(1);
			}
			
			if (!gSystem->ChangeDirectory(tmpDir.Data())) 
			{
				Log(fCurrentDetector.Data(), "Process - could not change directory!!");
				gSystem->Exit(1);			
			}
			
			Int_t success = ProcessCurrentDetector();

			gSystem->ChangeDirectory(wd.Data());
						
			if (success == 1) // Preprocessor finished successfully!
			{ 
				// remove temporary folder or DCS map
   				if (!fConfig->KeepTempFolder())
				{
					gSystem->Exec(Form("rm -rf %s",tmpDir.Data()));
				} else if (!fConfig->KeepDCSMap())
				{
					gSystem->Exec(Form("rm -f %s/DCSMap.root",tmpDir.Data()));
				}
				
				// Update time_processed field in FXS DB
				if (UpdateTable() == kFALSE)
					Log("SHUTTLE", Form("Process - %s: Could not update FXS databases!", 
							fCurrentDetector.Data()));

				// Transfer the data from local storage to main storage (Grid)
				if (StoreOCDB() == kFALSE)
					success = kFALSE;
			} 
			else if (success == 0)
			{
				Log("SHUTTLE", 
					Form("\t\t\t****** run %d - %s: ERROR ******",
						GetCurrentRun(), aDetector->GetName()));
			}

			for (UInt_t iSys=0; iSys<3; iSys++)
			{
				if (fFXSCalled[iSys]) fFXSlist[iSys].Clear();
			}

			Log("SHUTTLE", Form("Process - Client process of %d - %s is exiting now with %d.",
							GetCurrentRun(), aDetector->GetName(), success));

			// the client exits here
			gSystem->Exit(success);

			AliError("We should never get here!!!");
		}
	}

	Log("SHUTTLE", Form("\t\t\t^*^*^*^*^*^*^*^*^*^*^*^* run %d: FINISH ^*^*^*^*^*^*^*^*^*^*^*^*",
							GetCurrentRun()));

	//check if shuttle is done for this run, if so update logbook
	TObjArray checkEntryArray;
	checkEntryArray.SetOwner(1);
	TString whereClause = Form("where run=%d", GetCurrentRun());
	if (!QueryShuttleLogbook(whereClause.Data(), checkEntryArray) || 
			checkEntryArray.GetEntries() == 0) {
		Log("SHUTTLE", Form("Process - Warning: Cannot check status of run %d on Shuttle logbook!",
						GetCurrentRun()));
		return hasError == kFALSE;
	}

	AliShuttleLogbookEntry* checkEntry = dynamic_cast<AliShuttleLogbookEntry*>
						(checkEntryArray.At(0));

	if (checkEntry)
	{
		if (checkEntry->IsDone())
		{
			Log("SHUTTLE","Process - Shuttle is DONE. Updating logbook");
			UpdateShuttleLogbook("shuttle_done");
		}
		else
		{
			for (UInt_t iDet=0; iDet<NDetectors(); iDet++)
			{
				if (checkEntry->GetDetectorStatus(iDet) == AliShuttleLogbookEntry::kUnprocessed)
				{
					AliDebug(2, Form("Run %d: setting %s as \"not first time unprocessed\"",
							checkEntry->GetRun(), GetDetName(iDet)));
					fFirstUnprocessed[iDet] = kFALSE;
				}
			}
			SendMLRunInfo("Pending");
		}
	}

	fLogbookEntry = 0;

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
Int_t AliShuttle::ProcessCurrentDetector()
{
	//
        // Makes data retrieval just for a specific detector (fCurrentDetector).
	// Threre should be a configuration for this detector.

	Log("SHUTTLE", Form("ProcessCurrentDetector - Retrieving values for %s, run %d", 
						fCurrentDetector.Data(), GetCurrentRun()));

	TString wd = gSystem->WorkingDirectory();
	
	if (!CleanReferenceStorage(fCurrentDetector.Data()))
		return 0;
	
	gSystem->ChangeDirectory(wd.Data());
	
	// call preprocessor
	AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(fCurrentDetector));

	// check if the preprocessor wants to process this run type
	if (aPreprocessor->ProcessRunType() == kFALSE)
	{
		UpdateShuttleStatus(AliShuttleStatus::kSkipped);
		UpdateShuttleLogbook(fCurrentDetector, "DONE");
		if (!UpdateTableSkippedCase(fCurrentDetector.Data()))
		{
			AliError(Form("Could not update FXS tables for run %d !", GetCurrentRun()));
		}
		Log(fCurrentDetector, Form("ProcessCurrentDetector - %s preprocessor is not interested in this run type", fCurrentDetector.Data()));
	
		return 2;
	}
	
	// checking if OCDB is reachable
	AliCDBEntry* testEntry = GetFromOCDB("SHUTTLE","GRP/CTP/DummyConfig");
	if (!testEntry){
		// OCDB is not accessible, going in OCDBError for current detector
		AliError("OCDB Test entry not accessible");
		UpdateShuttleStatus(AliShuttleStatus::kOCDBError);
		return 0;
	}

	TMap* dcsMap = new TMap();
	
	aPreprocessor->Initialize(GetCurrentRun(), GetCurrentStartTime(), GetCurrentEndTime());

	Bool_t processDCS = aPreprocessor->ProcessDCS();

	if (!processDCS)
	{
		Log(fCurrentDetector, "ProcessCurrentDetector -"
			" The preprocessor requested to skip the retrieval of DCS values");
	}
	else if (fTestMode & kSkipDCS)
	{
		Log(fCurrentDetector, "ProcessCurrentDetector - In TESTMODE: Skipping DCS processing");
	} 
	else if (fTestMode & kErrorDCS)
	{
		Log(fCurrentDetector, "ProcessCurrentDetector - In TESTMODE: Simulating DCS error");
		UpdateShuttleStatus(AliShuttleStatus::kDCSStarted);
		UpdateShuttleStatus(AliShuttleStatus::kDCSError);
		delete dcsMap;
		return 0;
	} else {

		UpdateShuttleStatus(AliShuttleStatus::kDCSStarted);

		// Query DCS archive
		Int_t nServers = fConfig->GetNServers(fCurrentDetector);
		
		for (int iServ=0; iServ<nServers; iServ++)
		{
		
			TString host(fConfig->GetDCSHost(fCurrentDetector, iServ));
			Int_t port = fConfig->GetDCSPort(fCurrentDetector, iServ);
			Int_t multiSplit = fConfig->GetMultiSplit(fCurrentDetector, iServ);

			Log(fCurrentDetector, Form("ProcessCurrentDetector -"
					" Querying DCS Amanda server %s:%d (%d of %d)", 
					host.Data(), port, iServ+1, nServers));
			
			TMap* aliasMap = 0;
			TMap* dpMap = 0;

			if (fConfig->GetDCSAliases(fCurrentDetector, iServ)->GetEntries() > 0)
			{
				Log(fCurrentDetector, Form("Querying %d DCS aliases", fConfig->GetDCSAliases(fCurrentDetector, iServ)->GetEntries()));
				aliasMap = GetValueSet(host, port, 
						fConfig->GetDCSAliases(fCurrentDetector, iServ), 
						kAlias, multiSplit);
				if (!aliasMap)
				{
					Log(fCurrentDetector, 
						Form("ProcessCurrentDetector -"
							" Error retrieving DCS aliases from server %s."
							" Sending mail to DCS experts!", host.Data()));
					UpdateShuttleStatus(AliShuttleStatus::kDCSError);
					
					if (!SendMail(kDCSEMail))
						Log("SHUTTLE", Form("ProcessCurrentDetector - "
								    "Could not send mail to DCS experts!"));

					delete dcsMap;
					return 0;
				}
			}
			
			if (fConfig->GetDCSDataPoints(fCurrentDetector, iServ)->GetEntries() > 0)
			{
				Log(fCurrentDetector, Form("Querying %d DCS data points", fConfig->GetDCSDataPoints(fCurrentDetector, iServ)->GetEntries()));
				dpMap = GetValueSet(host, port, 
						fConfig->GetDCSDataPoints(fCurrentDetector, iServ), 
						kDP, multiSplit);
				if (!dpMap)
				{
					Log(fCurrentDetector, 
						Form("ProcessCurrentDetector -"
							" Error retrieving DCS data points from server %s."
							" Sending mail to DCS experts!", host.Data()));
					UpdateShuttleStatus(AliShuttleStatus::kDCSError);
					
					if (!SendMail(kDCSEMail))
						Log("SHUTTLE", Form("ProcessCurrentDetector - "
								    "Could not send mail to DCS experts!"));
					
					if (aliasMap) delete aliasMap;
					delete dcsMap;
					return 0;
				}				
			}
			
			// merge aliasMap and dpMap into dcsMap
			if(aliasMap) {
				TIter iter(aliasMap);
				TObjString* key = 0;
				while ((key = (TObjString*) iter.Next()))
					dcsMap->Add(key, aliasMap->GetValue(key->String()));
				
				aliasMap->SetOwner(kFALSE);
				delete aliasMap;
			}	
			
			if(dpMap) {
				TIter iter(dpMap);
				TObjString* key = 0;
				while ((key = (TObjString*) iter.Next()))
					dcsMap->Add(key, dpMap->GetValue(key->String()));
				
				dpMap->SetOwner(kFALSE);
				delete dpMap;
			}
		}
	}
	
	// save map into file, to help debugging in case of preprocessor error
	TFile* f = TFile::Open("DCSMap.root","recreate");
	f->cd();
	dcsMap->Write("DCSMap", TObject::kSingleKey);
	f->Close();
	delete f;
	
	// DCS Archive DB processing successful. Call Preprocessor!
	UpdateShuttleStatus(AliShuttleStatus::kPPStarted);

	fFXSError = -1; // this variable is kTRUE after ::Process if an FXS error occured
	
	UInt_t returnValue = aPreprocessor->Process(dcsMap);
	
	if (fFXSError!=-1) {
		UpdateShuttleStatus(AliShuttleStatus::kFXSError);
		SendMail(kFXSEMail, fFXSError);
		dcsMap->DeleteAll();
		delete dcsMap;
		return 0;
	}

	if (returnValue > 0) // Preprocessor error!
	{
		Log(fCurrentDetector, Form("ProcessCurrentDetector - "
				"Preprocessor failed. Process returned %d.", returnValue));
		UpdateShuttleStatus(AliShuttleStatus::kPPError);
		dcsMap->DeleteAll();
		delete dcsMap;
		return 0;
	}
	
	// preprocessor ok!
	UpdateShuttleStatus(AliShuttleStatus::kPPDone);
	Log(fCurrentDetector, Form("ProcessCurrentDetector - %s preprocessor returned success",
				fCurrentDetector.Data()));

	dcsMap->DeleteAll();
	delete dcsMap;

	return 1;
}

//______________________________________________________________________________________________
void AliShuttle::CountOpenRuns()
{
	// Query DAQ's Shuttle logbook and sends the number of open runs to ML
	
	SendAlive();
	
	// check connection, in case connect
	if (!Connect(4)) 
		return;

	TString sqlQuery;
	sqlQuery = Form("select count(*) from %s where shuttle_done=0", fConfig->GetShuttlelbTable());
	
	TSQLResult* aResult = fServer[4]->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return;
	}

	AliDebug(2,Form("Query = %s", sqlQuery.Data()));
	
	if (aResult->GetRowCount() == 0) {
		AliError(Form("No result for query %s received", sqlQuery.Data()));
		return;
	}

	if (aResult->GetFieldCount() != 1) {
		AliError(Form("Invalid field count for query %s received", sqlQuery.Data()));
		return;
	}

	TSQLRow* aRow = aResult->Next();
	if (!aRow) {
		AliError(Form("Could not receive result of query %s", sqlQuery.Data()));
		return;
	}
	
	TString result(aRow->GetField(0), aRow->GetFieldLength(0));
	Int_t count = result.Atoi();
	
	Log("SHUTTLE", Form("%d unprocessed runs", count));
	
	delete aRow;
	delete aResult;

	TMonaLisaValue mlStatus("SHUTTLE_openruns", count);

	TList mlList;
	mlList.Add(&mlStatus);

	fMonaLisa->SendParameters(&mlList, "__PROCESSINGINFO__");
}

//______________________________________________________________________________________________
Bool_t AliShuttle::QueryShuttleLogbook(const char* whereClause,
		TObjArray& entries)
{
	// Query DAQ's Shuttle logbook and fills detector status object.
	// Call QueryRunParameters to query DAQ logbook for run parameters.
	//

	entries.SetOwner(1);

	// check connection, in case connect
	if (!Connect(4)) return kFALSE;

	TString sqlQuery;
	sqlQuery = Form("select * from %s %s order by run", fConfig->GetShuttlelbTable(), whereClause);

	TSQLResult* aResult = fServer[4]->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	AliDebug(2,Form("Query = %s", sqlQuery.Data()));

	if(aResult->GetRowCount() == 0) {
		Log("SHUTTLE", "No entries in Shuttle Logbook match request");
		delete aResult;
		return kTRUE;
	}

	// TODO Check field count!
	const UInt_t nCols = 26;
	if (aResult->GetFieldCount() != (Int_t) nCols) {
		Log("SHUTTLE", "Invalid SQL result field number!");
		delete aResult;
		return kFALSE;
	}

	TSQLRow* aRow;
	while ((aRow = aResult->Next())) {
		TString runString(aRow->GetField(0), aRow->GetFieldLength(0));
		Int_t run = runString.Atoi();

		AliShuttleLogbookEntry *entry = QueryRunParameters(run);
		if (!entry)
			continue;

		// DA test mode flag
		TString daTestModeString(aRow->GetField(2), aRow->GetFieldLength(2)); // field 2 = DA test mode flag 
		Bool_t daTestMode = (Bool_t)daTestModeString.Atoi();
		entry->SetDATestMode(daTestMode);

		// loop on detectors
		for(UInt_t ii = 0; ii < nCols; ii++)
			entry->SetDetectorStatus(aResult->GetFieldName(ii), aRow->GetField(ii));

		entries.AddLast(entry);
		delete aRow;
	}

	delete aResult;
	return kTRUE;
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry* AliShuttle::QueryRunParameters(Int_t run)
{
	//
	// Retrieve run parameters written in the DAQ logbook and sets them into AliShuttleLogbookEntry object
	//

	// check connection, in case connect
	if (!Connect(4))
		return 0;

	TString sqlQuery;
	sqlQuery.Form("select * from %s where run=%d", fConfig->GetDAQlbTable(), run);

	TSQLResult* aResult = fServer[4]->Query(sqlQuery);
	if (!aResult) {
		Log("SHUTTLE", Form("Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}

	if (aResult->GetRowCount() == 0) {
		Log("SHUTTLE", Form("QueryRunParameters - No entry in DAQ Logbook for run %d. Skipping", run));
		delete aResult;
		return 0;
	}

	if (aResult->GetRowCount() > 1) {
		Log("SHUTTLE", Form("QueryRunParameters - UNEXPECTED: "
				"more than one entry in DAQ Logbook for run %d!", run));
		delete aResult;
		return 0;
	}

	TSQLRow* aRow = aResult->Next();
	if (!aRow)
	{
		Log("SHUTTLE", Form("QueryRunParameters - Could not retrieve row for run %d. Skipping", run));
		delete aResult;
		return 0;
	}

	AliShuttleLogbookEntry* entry = new AliShuttleLogbookEntry(run);

	for (Int_t ii = 0; ii < aResult->GetFieldCount(); ii++)
		entry->SetRunParameter(aResult->GetFieldName(ii), aRow->GetField(ii));

	delete aRow;
	delete aResult;
	
	UInt_t startTime = entry->GetStartTime();
	UInt_t endTime = entry->GetEndTime();
	Bool_t ecsSuccess = entry->GetECSSuccess();
	TString runType = entry->GetRunType();
	TString tmpdaqstartTime = entry->GetRunParameter("DAQ_time_start");
	TString recordingFlagString = entry->GetRunParameter("GDCmStreamRecording");
	UInt_t recordingFlag = recordingFlagString.Atoi();
	UInt_t daqstartTime = tmpdaqstartTime.Atoi();
	
	UInt_t now = time(0);
	Int_t dcsDelay = fConfig->GetDCSDelay()+fConfig->GetDCSQueryOffset();

	Bool_t skip = kFALSE;
				
	// runs are processed if
	//   a) runType is PHYSICS and ecsSuccess is set
	//   b) runType is not PHYSICS and (ecsSuccess is set or DAQ_time_start is non-0)
	// effectively this means that all runs are processed that started properly (ecsSucess behaviour is different for PHYSICS and non-PHYSICS runs (check with ECS!)
	if (startTime != 0 && endTime != 0) {
		if (endTime > startTime) { 
			if (endTime >= now - dcsDelay) {
				Log("SHUTTLE", Form("Skipping run %d for now, because DCS buffer time is not yet expired", run));
			} else {
				if ((runType == "PHYSICS" || runType == "STANDALONE") && recordingFlag == 0){
					Log("SHUTTLE", Form("QueryRunParameters - Run type for run %d is %s but the recording is OFF - Skipping!", run, runType.Data()));
					skip = kTRUE;
				} 
				else {
					if (runType == "PHYSICS") {
						if (ecsSuccess) {
							return entry;
						} else {
							Log("SHUTTLE", Form("QueryRunParameters - Run type for run %d is PHYSICS but ECS success flag not set (Reason = %s) - Skipping!", run, entry->GetRunParameter("eor_reason")));
							skip = kTRUE;
						} 
					} else {
						if (ecsSuccess || daqstartTime > 0) {
							if (ecsSuccess == kFALSE)
								Log("SHUTTLE", Form("Processing run %d although in status ECS failure (Reason: %s), since run type != PHYSICS and DAQ_time_start != 0", run, entry->GetRunParameter("eor_reason")));
							return entry;
						} else {
							Log("SHUTTLE", Form("QueryRunParameters - Run type for run %d is %s, ECS success flag was not set (Reason = %s) and DAQ_time_start was NULL - Skipping!", run, runType.Data(), entry->GetRunParameter("eor_reason")));
							skip = kTRUE;
						}
					}
				}
			}
		} else {
			Log("SHUTTLE", Form("QueryRunParameters - Invalid parameters for run %d: startTime equal to endTime: %d %d - Skipping!", run, startTime, endTime));
			skip = kTRUE;
		}
	} else {
	  Log("SHUTTLE", Form("QueryRunParameters - Invalid parameters for Run %d: "
			      "startTime = %d, endTime = %d. Skipping (Shuttle won't be marked as DONE)!",
			      run, startTime, endTime));
	}
	
	if (skip)
	{
		Log("SHUTTLE", Form("Marking SHUTTLE skipped for run %d", run));
		fLogbookEntry = entry;
		if (!UpdateShuttleLogbook("shuttle_skipped"))
		{
			AliError(Form("Could not update logbook for run %d !", run));
		}
		if (!UpdateTableSkippedCase("ALL"))
		{
			AliError(Form("Could not update FXS tables for run %d !", run));
		}
		fLogbookEntry = 0;
	}
			
	delete entry;
	return 0;
}

//______________________________________________________________________________________________
TMap* AliShuttle::GetValueSet(const char* host, Int_t port, const TSeqCollection* entries,
			      DCSType type, Int_t multiSplit)
{
	// Retrieve all "entry" data points from the DCS server
	// host, port: TSocket connection parameters
	// entries: list of name of the alias or data point
	// type: kAlias or kDP
	// returns TMap of values, 0 when failure
	
	AliDCSClient client(host, port, fTimeout, fRetries, multiSplit);

	TMap* result = 0;
	if (type == kAlias)
	{
		//result = client.GetAliasValues(entries, GetCurrentStartTime()-offset, 
		//      GetCurrentEndTime()+offset);
		result = client.GetAliasValues(entries, GetStartTimeDCSQuery(), 
			GetEndTimeDCSQuery());
	} 
	else if (type == kDP)
	{
		//result = client.GetDPValues(entries, GetCurrentStartTime()-offset, 
		//	GetCurrentEndTime()+offset);
       		result = client.GetDPValues(entries, GetStartTimeDCSQuery(), 
			GetEndTimeDCSQuery());
	}

	if (result == 0)
	{
		Log(fCurrentDetector.Data(), Form("GetValueSet - Can't get entries! Reason: %s",
			client.GetErrorString(client.GetResultErrorCode())));
		if (client.GetResultErrorCode() == AliDCSClient::fgkServerError)	
			Log(fCurrentDetector.Data(), Form("GetValueSet - Server error code: %s",
				client.GetServerError().Data()));

		return 0;
	}
		
	return result;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetFile(Int_t system, const char* detector,
		const char* id, const char* source)
{
	// Get calibration file from file exchange servers
	// First queris the FXS database for the file name, using the run, detector, id and source info
	// then calls RetrieveFile(filename) for actual copy to local disk
	// run: current run being processed (given by Logbook entry fLogbookEntry)
	// detector: the Preprocessor name
	// id: provided as a parameter by the Preprocessor
	// source: provided by the Preprocessor through GetFileSources function

	// check if test mode should simulate a FXS error
	if (fTestMode & kErrorFXSFiles)
	{
		Log(detector, Form("GetFile - In TESTMODE - Simulating error while connecting to %s FXS", GetSystemName(system)));
		return 0;
	}
	
	// check connection, in case connect
	if (!Connect(system))
	{
		Log(detector, Form("GetFile - Couldn't connect to %s FXS database", GetSystemName(system)));
		fFXSError = system;
		return 0;
	}

	// Query preparation
	TString sourceName(source);
	Int_t nFields = 3;
	TString sqlQueryStart = Form("select filePath,size,fileChecksum from %s where",
								fConfig->GetFXSdbTable(system));
	TString whereClause = Form("run=%d and detector=\"%s\" and fileId=\"%s\"",
								GetCurrentRun(), detector, id);

	if (system == kDAQ || system == kDQM)
	{
		whereClause += Form(" and DAQsource=\"%s\"", source);
	}
	else if (system == kDCS)
	{
		sourceName="none";
	}
	else if (system == kHLT)
	{
		whereClause += Form(" and DDLnumbers=\"%s\"", source);
	}

	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult = 0;
	aResult = dynamic_cast<TSQLResult*> (fServer[system]->Query(sqlQuery));
	if (!aResult) {
		Log(detector, Form("GetFile - Can't execute SQL query to %s database for: id = %s, source = %s",
				GetSystemName(system), id, sourceName.Data()));
		fFXSError = system;
		return 0;
	}

	if (aResult->GetRowCount() == 0)
	{
		Log(detector,
			Form("GetFile - No entry in %s FXS db for: id = %s, source = %s",
				GetSystemName(system), id, sourceName.Data()));
		delete aResult;
		return 0;
	}

	if (aResult->GetRowCount() > 1) {
		Log(detector,
			Form("GetFile - More than one entry in %s FXS db for: id = %s, source = %s",
				GetSystemName(system), id, sourceName.Data()));
		fFXSError = system;
		delete aResult;
		return 0;
	}

	if (aResult->GetFieldCount() != nFields) {
		Log(detector,
			Form("GetFileName - Wrong field count in %s FXS db for: id = %s, source = %s",
				GetSystemName(system), id, sourceName.Data()));
		fFXSError = system;
		delete aResult;
		return 0;
	}

	TSQLRow* aRow = dynamic_cast<TSQLRow*> (aResult->Next());

	if (!aRow){
		Log(detector, Form("GetFile - Empty set result in %s FXS db from query: id = %s, source = %s",
				GetSystemName(system), id, sourceName.Data()));
		fFXSError = system;
		delete aResult;
		return 0;
	}

	TString filePath(aRow->GetField(0), aRow->GetFieldLength(0));
	TString fileSize(aRow->GetField(1), aRow->GetFieldLength(1));
	TString fileChecksum(aRow->GetField(2), aRow->GetFieldLength(2));

	delete aResult;
	delete aRow;

	AliDebug(2, Form("filePath = %s; size = %s, fileChecksum = %s",
				filePath.Data(), fileSize.Data(), fileChecksum.Data()));

	// retrieved file is renamed to make it unique
	Int_t dir_lev1 = GetCurrentRun()/10000;
	TString localFileName = Form("%s/%d/%d/%s_process/%s_%s_%d_%s_%s.shuttle",
				     GetShuttleTempDir(), dir_lev1, GetCurrentRun(), detector,
					GetSystemName(system), detector, GetCurrentRun(), 
					id, sourceName.Data());
	Log("SHUTTLE",Form("file from FXS = %s",localFileName.Data())); 


	// file retrieval from FXS
	UInt_t nRetries = 0;
	UInt_t maxRetries = 3;
	Bool_t result = kFALSE;

	// copy!! if successful TSystem::Exec returns 0
	while (nRetries++ < maxRetries) {
		AliDebug(2, Form("Trying to copy file. Retry # %d", nRetries));
		result = RetrieveFile(system, filePath.Data(), localFileName.Data());
		if (!result)
		{
			Log(detector, Form("GetFile - Copy of file %s from %s FXS failed",
					filePath.Data(), GetSystemName(system)));
			continue;
		} 

                if (fileSize.Length()>0)
                {
                        // compare filesize of local file with the one stored in the FXS DB
			Long_t size = -1;
			Int_t sizeComp = gSystem->GetPathInfo(localFileName.Data(), 0, &size, 0, 0);

			if (sizeComp != 0 || size != fileSize.Atoi())
			{
				Log(detector, Form("GetFile - size of file %s does not match with local copy!",
							filePath.Data()));
				result = kFALSE;
				continue;
			}

		} else {
			Log(fCurrentDetector, Form("GetFile - size of file %s not set in %s database, skipping comparison",
						filePath.Data(), GetSystemName(system)));
                }

		if (fileChecksum.Length()>0)
		{
			// compare md5sum of local file with the one stored in the FXS DB
			if(fileChecksum.Contains(' ')) fileChecksum.Resize(fileChecksum.First(' '));
			Int_t md5Comp = gSystem->Exec(Form("md5sum %s |grep %s > /dev/null 2> /dev/null",
						localFileName.Data(), fileChecksum.Data()));

			if (md5Comp != 0)
			{
				Log(detector, Form("GetFile - md5sum of file %s does not match with local copy!",
							filePath.Data()));
				result = kFALSE;
				continue;
			}
		} else {
			Log(fCurrentDetector, Form("GetFile - md5sum of file %s not set in %s database, skipping comparison",
							filePath.Data(), GetSystemName(system)));
		}
		if (result) break;
	}

	if (!result) 
	{
		fFXSError = system;
		return 0;
	}

	fFXSCalled[system]=kTRUE;
	TObjString *fileParams = new TObjString(Form("%s#!?!#%s", id, sourceName.Data()));
	fFXSlist[system].Add(fileParams);

	static TString staticLocalFileName;
	staticLocalFileName.Form("%s", localFileName.Data());
	
	Log(fCurrentDetector, Form("GetFile - Retrieved file with id %s and "
			"source %s from %s to %s", id, source, 
			GetSystemName(system), localFileName.Data()));
			
	return staticLocalFileName.Data();
}

//______________________________________________________________________________________________
Bool_t AliShuttle::RetrieveFile(UInt_t system, const char* fxsFileName, const char* localFileName)
{
	//
	// Copies file from FXS to local Shuttle machine
	//

	// check temp directory: trying to cd to temp; if it does not exist, create it
	AliDebug(2, Form("Copy file %s from %s FXS into %s",
			GetSystemName(system), fxsFileName, localFileName));
			
	TString tmpDir(localFileName);
	
	tmpDir = tmpDir(0,tmpDir.Last('/'));

	Int_t noDir = gSystem->GetPathInfo(tmpDir.Data(), 0, (Long64_t*) 0, 0, 0);
	if (noDir) // temp dir does not exists!
	{
		if (gSystem->mkdir(tmpDir.Data(), 1))
		{
			Log(fCurrentDetector.Data(), "RetrieveFile - could not make temp directory!!");
			return kFALSE;
		}
	}

	TString command = Form("scp -oPort=%d -2 %s@%s:%s/%s %s",
		fConfig->GetFXSPort(system),
		fConfig->GetFXSUser(system),
		fConfig->GetFXSHost(system),
		fConfig->GetFXSBaseFolder(system),
		fxsFileName,
		localFileName);

	AliDebug(2, Form("%s",command.Data()));

	Bool_t result = (gSystem->Exec(command.Data()) == 0);

	return result;
}

//______________________________________________________________________________________________
TList* AliShuttle::GetFileSources(Int_t system, const char* detector, const char* id)
{
	//
	// Get sources producing the condition file Id from file exchange servers
	// if id is NULL all sources are returned (distinct)
	//

	if (id)
	{
		Log(detector, Form("GetFileSources - Querying %s FXS for files with id %s produced by %s", GetSystemName(system), id, detector));
	} else {
		Log(detector, Form("GetFileSources - Querying %s FXS for files produced by %s", GetSystemName(system), detector));
	}
	
	// check if test mode should simulate a FXS error
	if (fTestMode & kErrorFXSSources)
	{
		Log(detector, Form("GetFileSources - In TESTMODE - Simulating error while connecting to %s FXS", GetSystemName(system)));
		return 0;
	}

	if (system == kDCS)
	{
		Log(detector, "GetFileSources - WARNING: DCS system has only one source of data!");
		TList *list = new TList();
		list->SetOwner(1);
		list->Add(new TObjString(" "));
		return list;
	}

	// check connection, in case connect
	if (!Connect(system))
	{
		Log(detector, Form("GetFileSources - Couldn't connect to %s FXS database", GetSystemName(system)));
		fFXSError = system;
		return NULL;
	}

	TString sourceName = "";
	if (system == kDAQ || system == kDQM)
	{
		sourceName = "DAQsource";
	} else if (system == kHLT)
	{
		sourceName = "DDLnumbers";
	}

	TString sqlQueryStart = Form("select distinct %s from %s where", sourceName.Data(), fConfig->GetFXSdbTable(system));
	TString whereClause = Form("run=%d and detector=\"%s\"",
				GetCurrentRun(), detector);
	if (id)
		whereClause += Form(" and fileId=\"%s\"", id);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = fServer[system]->Query(sqlQuery);
	if (!aResult) {
		Log(detector, Form("GetFileSources - Can't execute SQL query to %s database for id: %s",
				GetSystemName(system), id));
		fFXSError = system;
		return 0;
	}

	TList *list = new TList();
	list->SetOwner(1);
	
	if (aResult->GetRowCount() == 0)
	{
		Log(detector,
			Form("GetFileSources - No entry in %s FXS table for id: %s", GetSystemName(system), id));
		delete aResult;
		return list;
	}

	Log(detector, Form("GetFileSources - Found %d sources", aResult->GetRowCount()));

	TSQLRow* aRow;
	while ((aRow = aResult->Next()))
	{

		TString source(aRow->GetField(0), aRow->GetFieldLength(0));
		AliDebug(2, Form("%s = %s", sourceName.Data(), source.Data()));
		list->Add(new TObjString(source));
		delete aRow;
	}

	delete aResult;

	return list;
}

//______________________________________________________________________________________________
TList* AliShuttle::GetFileIDs(Int_t system, const char* detector, const char* source)
{
	//
	// Get all ids of condition files produced by a given source from file exchange servers
	//
	
        Log(detector, Form("GetFileIDs - Retrieving ids with source %s with %s", source, GetSystemName(system)));

	// check if test mode should simulate a FXS error
	if (fTestMode & kErrorFXSSources)
	{
		Log(detector, Form("GetFileIDs - In TESTMODE - Simulating error while connecting to %s FXS", GetSystemName(system)));
		return 0;
	}

	// check connection, in case connect
	if (!Connect(system))
	{
		Log(detector, Form("GetFileIDs - Couldn't connect to %s FXS database", GetSystemName(system)));
		return NULL;
	}

	TString sourceName = "";
	if (system == kDAQ)
	{
		sourceName = "DAQsource";
	} else if (system == kHLT)
	{
		sourceName = "DDLnumbers";
	}

	TString sqlQueryStart = Form("select fileId from %s where", fConfig->GetFXSdbTable(system));
	TString whereClause = Form("run=%d and detector=\"%s\"",
				GetCurrentRun(), detector);
	if (sourceName.Length() > 0 && source)
		whereClause += Form(" and %s=\"%s\"", sourceName.Data(), source);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = fServer[system]->Query(sqlQuery);
	if (!aResult) {
		Log(detector, Form("GetFileIDs - Can't execute SQL query to %s database for source: %s",
				GetSystemName(system), source));
		return 0;
	}

	TList *list = new TList();
	list->SetOwner(1);
	
	if (aResult->GetRowCount() == 0)
	{
		Log(detector,
			Form("GetFileIDs - No entry in %s FXS table for source: %s", GetSystemName(system), source));
		delete aResult;
		return list;
	}

        Log(detector, Form("GetFileIDs - Found %d ids", aResult->GetRowCount()));

	TSQLRow* aRow;

	while ((aRow = aResult->Next()))
	{

		TString id(aRow->GetField(0), aRow->GetFieldLength(0));
		AliDebug(2, Form("fileId = %s", id.Data()));
		list->Add(new TObjString(id));
		delete aRow;
	}

	delete aResult;

	return list;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Connect(Int_t system)
{
	// Connect to MySQL Server of the system's FXS MySQL databases
	// DAQ Logbook, Shuttle Logbook and DAQ FXS db are on the same host
	//

	// check connection: if already connected return

	if(fServer[system] && fServer[system]->IsConnected()) {
		// ping the server	       	
	        if (fServer[system]->PingVerify()==kTRUE){ // connection is still alive
			return kTRUE;
		}		
		else{
			AliWarning(Form("Connection got lost to FXS database for %s. Closing and reconnecting.",
					AliShuttleInterface::GetSystemName(system)));
			fServer[system]->Close();
			delete fServer[system];
			fServer[system] = 0x0;
		}
	}

	TString dbHost, dbUser, dbPass, dbName;

	if (system < 4) // FXS db servers
	{
		dbHost = Form("mysql://%s:%d", fConfig->GetFXSdbHost(system), fConfig->GetFXSdbPort(system));
		dbUser = fConfig->GetFXSdbUser(system);
		dbPass = fConfig->GetFXSdbPass(system);
		dbName =   fConfig->GetFXSdbName(system);
	} else { // Run & Shuttle logbook servers
	// TODO Will the Shuttle logbook server be the same as the Run logbook server ???
		dbHost = Form("mysql://%s:%d", fConfig->GetDAQlbHost(), fConfig->GetDAQlbPort());
		dbUser = fConfig->GetDAQlbUser();
		dbPass = fConfig->GetDAQlbPass();
		dbName =   fConfig->GetDAQlbDB();
	}

	fServer[system] = TSQLServer::Connect(dbHost.Data(), dbUser.Data(), dbPass.Data());
		if (!fServer[system] || !fServer[system]->IsConnected()) {
		if(system < 4)
		{
		AliError(Form("Can't establish connection to FXS database for %s",
					AliShuttleInterface::GetSystemName(system)));
		} else {
		AliError("Can't establish connection to Run logbook.");
		}
		if(fServer[system]) delete fServer[system];
		return kFALSE;
	}

	// Get tables
	TSQLResult* aResult=0;
	switch(system){
		case kDAQ:
			aResult = fServer[kDAQ]->GetTables(dbName.Data());
			break;
		case kDCS:
			aResult = fServer[kDCS]->GetTables(dbName.Data());
			break;
		case kHLT:
			aResult = fServer[kHLT]->GetTables(dbName.Data());
			break;
		case kDQM:
			aResult = fServer[kDQM]->GetTables(dbName.Data());
			break;
		default:
			aResult = fServer[4]->GetTables(dbName.Data());
			break;
	}

	delete aResult;
	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateTable()
{
	//
	// Update FXS table filling time_processed field in all rows corresponding to current run and detector
	//

	Bool_t result = kTRUE;

	for (UInt_t system=0; system<3; system++)
	{
		if(!fFXSCalled[system]) continue;

		// check connection, in case connect
		if (!Connect(system))
		{
			Log(fCurrentDetector, Form("UpdateTable - Couldn't connect to %s FXS database", GetSystemName(system)));
			result = kFALSE;
			continue;
		}

		TTimeStamp now; // now

		// Loop on FXS list entries
		TIter iter(&fFXSlist[system]);
		TObjString *aFXSentry=0;
		while ((aFXSentry = dynamic_cast<TObjString*> (iter.Next())))
		{
			TString aFXSentrystr = aFXSentry->String();
			TObjArray *aFXSarray = aFXSentrystr.Tokenize("#!?!#");
			if (!aFXSarray || aFXSarray->GetEntries() != 2 )
			{
				Log(fCurrentDetector, Form("UpdateTable - error updating %s FXS entry. Check string: <%s>",
					GetSystemName(system), aFXSentrystr.Data()));
				if(aFXSarray) delete aFXSarray;
				result = kFALSE;
				continue;
			}
			const char* fileId = ((TObjString*) aFXSarray->At(0))->GetName();
			const char* source = ((TObjString*) aFXSarray->At(1))->GetName();

			TString whereClause;
			if (system == kDAQ || system == kDQM)
			{
				whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\";",
							GetCurrentRun(), fCurrentDetector.Data(), fileId, source);
			}
			else if (system == kDCS)
			{
				whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\";",
							GetCurrentRun(), fCurrentDetector.Data(), fileId);
			}
			else if (system == kHLT)
			{
				whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\" and DDLnumbers=\"%s\";",
							GetCurrentRun(), fCurrentDetector.Data(), fileId, source);
			}

			delete aFXSarray;

			TString sqlQuery = Form("update %s set time_processed=%ld %s", fConfig->GetFXSdbTable(system),
						(ULong_t)now.GetSec(), whereClause.Data());

			AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

			// Query execution
			TSQLResult* aResult;
			aResult = dynamic_cast<TSQLResult*> (fServer[system]->Query(sqlQuery));
			if (!aResult)
			{
				Log(fCurrentDetector, Form("UpdateTable - %s db: can't execute SQL query <%s>",
								GetSystemName(system), sqlQuery.Data()));
				result = kFALSE;
				continue;
			}
			delete aResult;
		}
	}

	return result;
}

//_______________________________________________________________________________
Bool_t AliShuttle::UpdateTableSkippedCase(const char* detector)
{
	//
	// Update FXS table filling time_processed field in all rows corresponding to current run and detector
	// if detector = "ALL" update all detectors
	//

	Bool_t result = kTRUE;

	TString detName(detector);

	for (UInt_t system=0; system<3; system++)
	{

		// check connection, in case connect
		if (!Connect(system))
		{
			Log(fCurrentDetector, Form("UpdateTableSkippedCase - Couldn't connect to %s FXS database", GetSystemName(system)));
			result = kFALSE;
			continue;
		}

		TTimeStamp now; // now

		// Loop on FXS list entries
		TIter iter(&fFXSlist[system]);
			
		TString whereClause;
		if (detName == "ALL") whereClause = Form("where run=%d and time_processed IS NULL;",GetCurrentRun());
		else whereClause = Form("where run=%d and detector=\"%s\" and time_processed IS NULL;",GetCurrentRun(), detector);

		//Log("SHUTTLE",Form(" whereClause = %s ",whereClause.Data()));

		TString sqlQuery = Form("update %s set time_processed=%ld %s", fConfig->GetFXSdbTable(system),
					(ULong_t)now.GetSec(), whereClause.Data());

		AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

		// Query execution
		TSQLResult* aResult;
		aResult = dynamic_cast<TSQLResult*> (fServer[system]->Query(sqlQuery));
		if (!aResult)
		{
			Log("SHUTTLE", Form("UpdateTableSkippedCase - %s db: can't execute SQL query <%s>",
							GetSystemName(system), sqlQuery.Data()));
			result = kFALSE;
			continue;
		}
		delete aResult;
		
	}

	return result;
}
//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateTableFailCase()
{
	// Update FXS table filling time_processed field in all rows corresponding to current run and detector
	// this is called in case the preprocessor is declared failed for the current run, because
	// the fields are updated only in case of success

	Bool_t result = kTRUE;

	for (UInt_t system=0; system<3; system++)
	{
		// check connection, in case connect
		if (!Connect(system))
		{
			Log(fCurrentDetector, Form("UpdateTableFailCase - Couldn't connect to %s FXS database",
							GetSystemName(system)));
			result = kFALSE;
			continue;
		}

		TTimeStamp now; // now

		// Loop on FXS list entries

		TString whereClause = Form("where run=%d and detector=\"%s\";",
						GetCurrentRun(), fCurrentDetector.Data());


		TString sqlQuery = Form("update %s set time_processed=%ld %s", fConfig->GetFXSdbTable(system),
					(ULong_t)now.GetSec(), whereClause.Data());

		AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

		// Query execution
		TSQLResult* aResult;
		aResult = dynamic_cast<TSQLResult*> (fServer[system]->Query(sqlQuery));
		if (!aResult)
		{
			Log(fCurrentDetector, Form("UpdateTableFailCase - %s db: can't execute SQL query <%s>",
							GetSystemName(system), sqlQuery.Data()));
			result = kFALSE;
			continue;
		}
		delete aResult;
	}

	return result;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateShuttleLogbook(const char* detector, const char* status)
{
	//
	// Update Shuttle logbook filling detector or shuttle_done column
	// ex. of usage: UpdateShuttleLogbook("PHOS", "DONE") or UpdateShuttleLogbook("shuttle_done")
	//

	// check connection, in case connect
	if(!Connect(4)){
		Log("SHUTTLE", "UpdateShuttleLogbook - Couldn't connect to DAQ Logbook.");
		return kFALSE;
	}

	TString detName(detector);
	TString setClause;
	if (detName == "shuttle_done" || detName == "shuttle_skipped")
	{
		setClause = "set shuttle_done=1";
		
		if (detName == "shuttle_done")
		{
			if (TouchFile() != kTRUE)
			{
				SendMLRunInfo("Pending");
				return kFALSE;
			}
			
			SendMLRunInfo("Done");
		}
		else 
			SendMLRunInfo("Skipped");
	} 
	else {
		TString statusStr(status);
		if(statusStr.Contains("done", TString::kIgnoreCase) ||
		   statusStr.Contains("failed", TString::kIgnoreCase)){
			setClause = Form("set %s=\"%s\"", detector, status);
		} else {
			Log("SHUTTLE",
				Form("UpdateShuttleLogbook - Invalid status <%s> for detector %s",
					status, detector));
			return kFALSE;
		}
	}

	TString whereClause = Form("where run=%d", GetCurrentRun());

	TString sqlQuery = Form("update %s %s %s",
					fConfig->GetShuttlelbTable(), setClause.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = dynamic_cast<TSQLResult*> (fServer[4]->Query(sqlQuery));
	if (!aResult) {
		Log("SHUTTLE", Form("UpdateShuttleLogbook - Can't execute query <%s>", sqlQuery.Data()));
		return kFALSE;
	}
	delete aResult;

	return kTRUE;
}

//______________________________________________________________________________________________
Int_t AliShuttle::GetCurrentRun() const
{
	//
	// Get current run from logbook entry
	//

	return fLogbookEntry ? fLogbookEntry->GetRun() : -1;
}

//______________________________________________________________________________________________
UInt_t AliShuttle::GetCurrentStartTime() const
{
	//
	// get current start time
	//

	return fLogbookEntry ? fLogbookEntry->GetStartTime() : 0;
}

//______________________________________________________________________________________________
UInt_t AliShuttle::GetCurrentEndTime() const
{
	//
	// get current end time from logbook entry
	//

	return fLogbookEntry ? fLogbookEntry->GetEndTime() : 0;
}
//______________________________________________________________________________________________
UInt_t AliShuttle::GetCurrentYear() const
{
	//
	// Get current year from logbook entry
	//

	if (!fLogbookEntry) return 0;
	
	TTimeStamp startTime(GetCurrentStartTime());
	TString year =  Form("%d",startTime.GetDate());
	year = year(0,4);
	
	return year.Atoi();
}

//______________________________________________________________________________________________
const char* AliShuttle::GetLHCPeriod() const
{
	//
	// Get current LHC period from logbook entry
	//

	if (!fLogbookEntry) return 0;
		
	return fLogbookEntry->GetRunParameter("LHCperiod");
}

//______________________________________________________________________________________________
void AliShuttle::Log(const char* detector, const char* message, UInt_t level)
{
	//
	// Fill log string with a message
	//
	
	TString logRunDir = GetShuttleLogDir();
	if (GetCurrentRun() >=0) {
		Int_t logDir_lev1 = GetCurrentRun()/10000;
		logRunDir += Form("/%d/%d", logDir_lev1, GetCurrentRun());
	}		
	void* dir = gSystem->OpenDirectory(logRunDir.Data());
	if (dir == NULL) {
		if (gSystem->mkdir(logRunDir.Data(), kTRUE)) {
			AliError(Form("Can't open directory <%s>", GetShuttleLogDir()));
			return;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	TString toLog = Form("%s UTC (%d): %s - ", TTimeStamp(time(0)).AsString("s"), getpid(), detector);
	if (GetCurrentRun() >= 0) 
		toLog += Form("run %d - ", GetCurrentRun());
	toLog += Form("%s", message);

	AliLog::Message(level, toLog, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);
	
	// if we redirect the log output already to the file, leave here
	if (fOutputRedirected && strcmp(detector, "SHUTTLE") != 0)
		return;

  	TString fileName = GetLogFileName(detector);
	
  	gSystem->ExpandPathName(fileName);

  	ofstream logFile;
  	logFile.open(fileName, ofstream::out | ofstream::app);

  	if (!logFile.is_open()) {
    		AliError(Form("Could not open file %s", fileName.Data()));
    		return;
  	}

  	logFile << toLog.Data() << "\n";

  	logFile.close();
}

//______________________________________________________________________________________________
TString AliShuttle::GetLogFileName(const char* detector) const
{
	// 
	// returns the name of the log file for a given sub detector
	//
	
	TString fileName;
	
	if (GetCurrentRun() >= 0) 
	{
		Int_t logDir_lev1 = GetCurrentRun()/10000;
		fileName.Form("%s/%d/%d/%s.log", GetShuttleLogDir(), logDir_lev1, GetCurrentRun(), 
			detector);
	} else {
		fileName.Form("%s/%s.log", GetShuttleLogDir(), detector);
	}

	return fileName;
}

//______________________________________________________________________________________________
void AliShuttle::SendAlive()
{
	// sends alive message to ML
	
	TMonaLisaText mlStatus("SHUTTLE_status", "Alive");

	TList mlList;
	mlList.Add(&mlStatus);

	fMonaLisa->SendParameters(&mlList, "__PROCESSINGINFO__");
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Collect(Int_t run)
{
	//
	// Collects conditions data for all UNPROCESSED run written to DAQ LogBook in case of run = -1 (default)
	// If a dedicated run is given this run is processed
	//
	// In operational mode, this is the Shuttle function triggered by the EOR signal.
	//

	if (run == -1)
		Log("SHUTTLE","Collect - Shuttle called. Collecting conditions data for unprocessed runs");
	else
		Log("SHUTTLE", Form("Collect - Shuttle called. Collecting conditions data for run %d", run));

	SetLastAction("Starting");

	// create ML instance
	if (!fMonaLisa)
		fMonaLisa = new TMonaLisaWriter(fConfig->GetMonitorHost(), fConfig->GetMonitorTable());
		
	CountOpenRuns();

	TString whereClause("where shuttle_done=0");
	if (run != -1)
		whereClause += Form(" and run=%d", run);

	TObjArray shuttleLogbookEntries;
	if (!QueryShuttleLogbook(whereClause, shuttleLogbookEntries))
	{
		Log("SHUTTLE", "Collect - Can't retrieve entries from Shuttle logbook");
		return kFALSE;
	}

	if (shuttleLogbookEntries.GetEntries() == 0)
	{
		if (run == -1)
			Log("SHUTTLE","Collect - Found no UNPROCESSED runs in Shuttle logbook");
		else
			Log("SHUTTLE", Form("Collect - Run %d is already DONE "
						"or it does not exist in Shuttle logbook", run));
		return kTRUE;
	}

	for (UInt_t iDet=0; iDet<NDetectors(); iDet++)
		fFirstUnprocessed[iDet] = kTRUE;

	if (run != -1)
	{
		// query Shuttle logbook for earlier runs, check if some detectors are unprocessed,
		// flag them into fFirstUnprocessed array
		TString whereClauseBis(Form("where shuttle_done=0 and run < %d", run));
		TObjArray tmpLogbookEntries;
		if (!QueryShuttleLogbook(whereClauseBis, tmpLogbookEntries))
		{
			Log("SHUTTLE", "Collect - Can't retrieve entries from Shuttle logbook");
			return kFALSE;
		}

		TIter iter(&tmpLogbookEntries);
		AliShuttleLogbookEntry* anEntry = 0;
		while ((anEntry = dynamic_cast<AliShuttleLogbookEntry*> (iter.Next())))
		{
			for (UInt_t iDet=0; iDet<NDetectors(); iDet++)
			{
				if (anEntry->GetDetectorStatus(iDet) == AliShuttleLogbookEntry::kUnprocessed)
				{
					AliDebug(2, Form("Run %d: setting %s as \"not first time unprocessed\"",
							anEntry->GetRun(), GetDetName(iDet)));
					fFirstUnprocessed[iDet] = kFALSE;
				}
			}

		}

	}

	if (!RetrieveConditionsData(shuttleLogbookEntries))
	{
		Log("SHUTTLE", "Collect - Process of at least one run failed");
		CountOpenRuns();
		return kFALSE;
	}

	Log("SHUTTLE", "Collect - Requested run(s) successfully processed");
	CountOpenRuns();
	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::RetrieveConditionsData(const TObjArray& dateEntries)
{
	//
	// Retrieve conditions data for all runs that aren't processed yet
	//

	Bool_t hasError = kFALSE;

	TIter iter(&dateEntries);
	AliShuttleLogbookEntry* anEntry;

	while ((anEntry = (AliShuttleLogbookEntry*) iter.Next())){
		if (!Process(anEntry)){
			hasError = kTRUE;
		}

		// clean SHUTTLE temp directory
		//TString filename = Form("%s/*.shuttle", GetShuttleTempDir());
		//RemoveFile(filename.Data());
	}

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
ULong_t AliShuttle::GetTimeOfLastAction() const
{
	//
	// Gets time of last action
	//

	ULong_t tmp;

	fMonitoringMutex->Lock();

	tmp = fLastActionTime;

	fMonitoringMutex->UnLock();

	return tmp;
}

//______________________________________________________________________________________________
const TString AliShuttle::GetLastAction() const
{
	//
	// returns a string description of the last action
	//

	TString tmp;

	fMonitoringMutex->Lock();
	
	tmp = fLastAction;
	
	fMonitoringMutex->UnLock();

	return tmp;
}

//______________________________________________________________________________________________
void AliShuttle::SetLastAction(const char* action)
{
	//
	// updates the monitoring variables
	//

	fMonitoringMutex->Lock();

	fLastAction = action;
	fLastActionTime = time(0);
	
	fMonitoringMutex->UnLock();
}

//______________________________________________________________________________________________
const char* AliShuttle::GetRunParameter(const char* param)
{
	//
	// returns run parameter read from DAQ logbook
	//

	if(!fLogbookEntry) {
		AliError("No logbook entry!");
		return 0;
	}

	return fLogbookEntry->GetRunParameter(param);
}

//______________________________________________________________________________________________
AliCDBEntry* AliShuttle::GetFromOCDB(const char* detector, const AliCDBPath& path)
{
	//
	// returns object from OCDB valid for current run
	//

	if (fTestMode & kErrorOCDB)
	{
		Log(detector, "GetFromOCDB - In TESTMODE - Simulating error with OCDB");
		return 0;
	}
	
	AliCDBStorage *sto = AliCDBManager::Instance()->GetStorage(fgkMainCDB);
	if (!sto)
	{
		Log(detector, "GetFromOCDB - Cannot activate main OCDB for query!");
		return 0;
	}

	return dynamic_cast<AliCDBEntry*> (sto->Get(path, GetCurrentRun()));
}

//______________________________________________________________________________________________
Bool_t AliShuttle::SendMail(EMailTarget target, Int_t system)
{
	//
	// sends a mail to the subdetector expert in case of preprocessor error
	//
	
	if (fTestMode != kNone)
		return kTRUE;
		
	if (!fConfig->SendMail()) 
		return kTRUE;

	if (target == kDCSEMail || target == kFXSEMail) {
		if (!fFirstProcessing)
			return kTRUE;
	}

	Int_t runMode = (Int_t)fConfig->GetRunMode();
	TString tmpStr;
	if (runMode == 0) tmpStr = " Nightly Test:";
	else tmpStr = " Data Taking:"; 
	void* dir = gSystem->OpenDirectory(GetShuttleLogDir());
	if (dir == NULL)
	{
		if (gSystem->mkdir(GetShuttleLogDir(), kTRUE))
		{
			Log("SHUTTLE", Form("SendMail - Can't open directory <%s>", GetShuttleLogDir()));
			return kFALSE;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	// det experts in to
	TString to="";
	TIter *iterExperts = 0;
	if (target == kDCSEMail) {
		iterExperts = new TIter(fConfig->GetAdmins(AliShuttleConfig::kAmanda));
	}
	else if (target == kFXSEMail) {
		iterExperts = new TIter(fConfig->GetAdmins(system));
	}
	if (iterExperts) {
		TObjString *anExpert=0;
		while ((anExpert = (TObjString*) iterExperts->Next()))
		{
			to += Form("%s,", anExpert->GetName());
		}
		delete iterExperts;
	}

	// add subdetector experts	
	iterExperts = new TIter(fConfig->GetResponsibles(fCurrentDetector));
	TObjString *anExpert=0;
	while ((anExpert = (TObjString*) iterExperts->Next()))
	{
		to += Form("%s,", anExpert->GetName());
	}
	delete iterExperts;
	
	if (to.Length() > 0)
	  to.Remove(to.Length()-1);
	AliDebug(2, Form("to: %s",to.Data()));

	if (to.IsNull()) {
		Log("SHUTTLE", Form("List of %d responsibles not set!", (Int_t) target));
		return kFALSE;
	}

	// SHUTTLE responsibles in cc
	TString cc="";
	TIter iterAdmins(fConfig->GetAdmins(AliShuttleConfig::kGlobal));
	TObjString *anAdmin=0;
	while ((anAdmin = (TObjString*) iterAdmins.Next()))
	{
		cc += Form("%s,", anAdmin->GetName());
	}
	if (cc.Length() > 0)
	  cc.Remove(cc.Length()-1);
	AliDebug(2, Form("cc: %s",to.Data()));

	// mail body 
  	TString bodyFileName;
  	bodyFileName.Form("%s/mail.body", GetShuttleLogDir());
  	gSystem->ExpandPathName(bodyFileName);

  	ofstream mailBody;
  	mailBody.open(bodyFileName, ofstream::out);

  	if (!mailBody.is_open())
	{
    		Log("SHUTTLE", Form("Could not open mail body file %s", bodyFileName.Data()));
    		return kFALSE;
  	}


	TString subject;
	TString body;

	if (target == kDCSEMail){
		subject = Form("%s CRITICAL Retrieval of data points for %s FAILED in run %d !",
				tmpStr.Data(), fCurrentDetector.Data(), GetCurrentRun());
		AliDebug(2, Form("subject: %s", subject.Data()));
		
		body = Form("Dear DCS experts, \n\n");
		body += Form("SHUTTLE couldn\'t retrieve the data points for detector %s "
			     "in run %d!!\n\n", fCurrentDetector.Data(), GetCurrentRun());
	}
	else if (target == kFXSEMail){
		subject = Form("%s CRITICAL FXS communication for %s FAILED in run %d !",
				tmpStr.Data(), fCurrentDetector.Data(), GetCurrentRun());
		AliDebug(2, Form("subject: %s", subject.Data()));
		TString sys;
		if (system == kDAQ) sys="DAQ";
		else if (system == kDCS) sys="DCS";
		else if (system == kHLT) sys="HLT";
		else return kFALSE;
		body = Form("Dear  %s FXS experts, \n\n",sys.Data());
		body += Form("SHUTTLE couldn\'t retrieve data from the FXS for detector %s "
			     "in run %d!!\n\n", fCurrentDetector.Data(), GetCurrentRun());
		body += Form("The contacted server was:\nDB: %s\nFXS:%s\n\n", fConfig->GetFXSdbHost(system), fConfig->GetFXSHost(system));
	}
	else {
		subject = Form("%s %s Shuttle preprocessor FAILED in run %d (run type = %s)!",
				       tmpStr.Data(), fCurrentDetector.Data(), GetCurrentRun(), GetRunType());
		AliDebug(2, Form("subject: %s", subject.Data()));
	
		body = Form("Dear %s expert(s), \n\n", fCurrentDetector.Data());
		body += Form("SHUTTLE just detected that your preprocessor "
			     "failed processing run %d (run type = %s)!!\n\n", 
			     GetCurrentRun(), GetRunType());
	}

	body += Form("Please check %s status on the SHUTTLE monitoring page: \n\n", 
				fCurrentDetector.Data());
	if (fConfig->GetRunMode() == AliShuttleConfig::kTest)
	{
		body += Form("\thttp://pcalimonitor.cern.ch/shuttle.jsp?time=24 \n\n");
	} else {
		body += Form("\thttp://pcalimonitor.cern.ch/shuttle.jsp?instance=PROD&time=24 \n\n");
	}
	
	
	TString logFolder = "logs";
	if (fConfig->GetRunMode() == AliShuttleConfig::kProd) 
		logFolder += "_PROD";
	
	
	body += Form("Find the %s log for the current run on \n\n"
		"\thttp://pcalishuttle02.cern.ch/%s/%d/%d/%s.log \n\n", 
		     fCurrentDetector.Data(), logFolder.Data(), GetCurrentRun()/10000,  
				GetCurrentRun(), fCurrentDetector.Data());
	body += Form("The last 15 lines of %s log file are following:\n\n", fCurrentDetector.Data());

	AliDebug(2, Form("Body begin: %s", body.Data()));

	mailBody << body.Data();
  	mailBody.close();
  	mailBody.open(bodyFileName, ofstream::out | ofstream::app);

	TString logFileName = Form("%s/%d/%d/%s.log", GetShuttleLogDir(), 
		GetCurrentRun()/10000, GetCurrentRun(), fCurrentDetector.Data());
	TString tailCommand = Form("tail -n 15 %s >> %s", logFileName.Data(), bodyFileName.Data());
	if (gSystem->Exec(tailCommand.Data()))
	{
		mailBody << Form("%s log file not found ...\n\n", fCurrentDetector.Data());
	}

	TString endBody = Form("------------------------------------------------------\n\n");
	endBody += Form("In case of problems please contact the SHUTTLE core team.\n\n");
	endBody += "Please do not answer this message directly, it is automatically generated.\n\n";
	endBody += "Greetings,\n\n \t\t\tthe SHUTTLE\n";

	AliDebug(2, Form("Body end: %s", endBody.Data()));

	mailBody << endBody.Data();

  	mailBody.close();

	// send mail!
	TString mailCommand = Form("mail -s \"%s\" -c %s %s < %s",
						subject.Data(),
						cc.Data(),
						to.Data(),
						bodyFileName.Data());
	AliDebug(2, Form("mail command: %s", mailCommand.Data()));

	Bool_t result = gSystem->Exec(mailCommand.Data());

	return result == 0;
}
//______________________________________________________________________________________________
const char* AliShuttle::GetRunType()
{
	//
	// returns run type read from "run type" logbook
	//

	if(!fLogbookEntry) {
		AliError("No logbook entry!");
		return 0;
	}

	return fLogbookEntry->GetRunType();
}

//______________________________________________________________________________________________
Bool_t AliShuttle::GetHLTStatus()
{
  	// Return HLT status (ON=1 OFF=0)
  	// Converts the HLT status from the mode string read in the run logbook (not just a bool)

	if(!fLogbookEntry) {
		AliError("No logbook entry!");
		return 0;
	}

	// TODO implement when HLTMode is inserted in run logbook
	TString hltMode = fLogbookEntry->GetRunParameter("HLTmode");
	TSubString firstChar = hltMode(0,1);
	AliDebug(2,Form("First char = %s ",firstChar.Data())); 
	if (firstChar == "A") {
		return kFALSE;
	}
	else if ((firstChar == "B") || (firstChar == "C") || (firstChar == "D") || (firstChar == "E")) {
		return kTRUE;
	}
	else {
		Log("SHUTTLE","Unexpected HLT mode! Returning 0....");
		return kFALSE;
	}
}

//______________________________________________________________________________________________
const char* AliShuttle::GetTriggerConfiguration()
{
	// Receives the trigger configuration from the DAQ logbook for the current run
	
	// check connection, if needed reconnect
	if (!Connect(4)) 
		return 0;

	TString sqlQuery;
	sqlQuery.Form("SELECT configFile FROM logbook_trigger_config WHERE run = %d", GetCurrentRun());
	TSQLResult* result = fServer[4]->Query(sqlQuery);
	if (!result)
	{
		Log("SHUTTLE", Form("ERROR: Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}
	
	if (result->GetRowCount() == 0)
	{
		Log("SHUTTLE", "WARNING: Trigger configuration not found in logbook_trigger_config");
		delete result;
		return 0;
	}
	
	TSQLRow* row = result->Next();
	if (!row)
	{
		Log("SHUTTLE", "ERROR: Could not receive logbook_trigger_config data");
		delete result;
		return 0;
	}

	// static, so that pointer remains valid when it is returned to the calling class	
	static TString triggerConfig(row->GetField(0));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Log("SHUTTLE", Form("Found trigger configuration: %s", triggerConfig.Data()));
	
	return triggerConfig;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetCTPTimeParams()
{
	// Receives the CTP time parameters from the DAQ logbook for the current run
	
	// check connection, if needed reconnect
	if (!Connect(4)) 
		return 0;

	TString sqlQuery;
	sqlQuery.Form("SELECT alignmentFile FROM logbook_trigger_config WHERE run = %d", GetCurrentRun());
	TSQLResult* result = fServer[4]->Query(sqlQuery);
	if (!result)
	{
		Log("SHUTTLE", Form("ERROR: Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}
	
	if (result->GetRowCount() == 0)
	{
		Log("SHUTTLE", "WARNING: CTP time params not found in logbook_trigger_config");
		delete result;
		return 0;
	}
	
	TSQLRow* row = result->Next();
	if (!row)
	{
		Log("SHUTTLE", "ERROR: Could not receive logbook_trigger_config data");
		delete result;
		return 0;
	}

	// static, so that pointer remains valid when it is returned to the calling class	
	static TString triggerTimeParams(row->GetField(0));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Log("SHUTTLE", Form("Found trigger time parameters: %s", triggerTimeParams.Data()));
	
	return triggerTimeParams;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetTriggerDetectorMask()
{
	// Receives the trigger detector mask from DAQ logbook
	
	// check connection, if needed reconnect
	if (!Connect(4)) 
		return 0;

	TString sqlQuery;
	sqlQuery.Form("SELECT BIN(BIT_OR(inputDetectorMask)) from logbook_trigger_clusters WHERE run = %d;", GetCurrentRun());
	TSQLResult* result = fServer[4]->Query(sqlQuery);
	if (!result)
	{
		Log("SHUTTLE", Form("ERROR: Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}
	
	if (result->GetRowCount() == 0)
	{
		Log("SHUTTLE", "ERROR: Trigger Detector Mask not found in logbook_trigger_clusters");
		delete result;
		return 0;
	}
	
	TSQLRow* row = result->Next();
	if (!row)
	{
		Log("SHUTTLE", "ERROR: Could not receive logbook_trigger_clusters data");
		delete result;
		return 0;
	}

	// static, so that pointer remains valid when it is returned to the calling class	
	static TString triggerDetectorMask(row->GetField(0));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Log("SHUTTLE", Form("Found Trigger Detector Mask: %s", triggerDetectorMask.Data()));
	
	return triggerDetectorMask;
}

//______________________________________________________________________________________________
void AliShuttle::SetShuttleTempDir(const char* tmpDir)
{
	//
	// sets Shuttle temp directory
	//

	fgkShuttleTempDir = gSystem->ExpandPathName(tmpDir);
}

//______________________________________________________________________________________________
void AliShuttle::SetShuttleLogDir(const char* logDir)
{
	//
	// sets Shuttle log directory
	//

	fgkShuttleLogDir = gSystem->ExpandPathName(logDir);
}
//______________________________________________________________________________________________
Bool_t AliShuttle::TouchFile()
{
	//
	// touching a file on the grid if run has been DONE
	//
	
	if (!gGrid)
	{
		Log("SHUTTLE",Form("No TGrid connection estabilished!"));
		Log("SHUTTLE",Form("Could not touch file for run %i",GetCurrentRun()));
		return kFALSE;
	}

	TString dir;
	dir.Form("%s%d/SHUTTLE_DONE", fConfig->GetAlienPath(), GetCurrentYear());
	// checking whether directory for touch command exists
	TString commandLs;
	commandLs.Form("ls %s",dir.Data());
	TGridResult *resultLs = dynamic_cast<TGridResult*>(gGrid->Command(commandLs));
	if (!resultLs){
		Log("SHUTTLE",Form("No result for %s command, returning without touching",commandLs.Data()));
		return kFALSE;
	}
	TMap *mapLs = dynamic_cast<TMap*>(resultLs->At(0));
	if (!mapLs){
		Log("SHUTTLE",Form("No map for %s command, returning without touching",commandLs.Data()));
		delete resultLs;
		resultLs = 0x0;
		return kFALSE;
	}
	TObjString *valueLsPath = dynamic_cast<TObjString*>(mapLs->GetValue("path"));
	if (!valueLsPath || (valueLsPath->GetString()).CompareTo(dir)!=1){ 
		Log("SHUTTLE",Form("No directory %s found, creating it",dir.Data()));

		// creating the directory

		Bool_t boolMkdir = gGrid->Mkdir(dir.Data());
		if (!boolMkdir) {
			Log("SHUTTLE",Form("Impossible to create dir %s in alien catalogue for run %i!",dir.Data(),GetCurrentRun()));
			delete resultLs;
			resultLs = 0x0;
			return kFALSE;
		}
		Log("SHUTTLE",Form("Directory %s successfully created in alien catalogue for run %i",dir.Data(),GetCurrentRun()));
      	}
	else {
		Log("SHUTTLE",Form("Directory %s correctly found for run %i",dir.Data(),GetCurrentRun()));
	}

	delete resultLs;
	resultLs = 0x0;

	TString command;
	command.Form("touch %s/%i", dir.Data(), GetCurrentRun());
	Log("SHUTTLE", Form("Creating entry in file catalog: %s", command.Data()));
	TGridResult *resultTouch = dynamic_cast<TGridResult*>(gGrid->Command(command));
	if (!resultTouch){
		Log("SHUTTLE",Form("No result for touching command, returning without touching for run %i",GetCurrentRun()));
		return kFALSE;
	}
	TMap *mapTouch = dynamic_cast<TMap*>(resultTouch->At(0));
	if (!mapTouch){
		Log("SHUTTLE",Form("No map for touching command, returning without touching for run %i",GetCurrentRun()));
		delete resultTouch;
		resultTouch = 0x0; 
		return kFALSE;
	}
	TObjString *valueTouch = dynamic_cast<TObjString*>(mapTouch->GetValue("__result__"));
	if (!valueTouch){
		Log("SHUTTLE",Form("No value for \"__result__\" key set in the map for touching command, returning without touching for run %i",GetCurrentRun()));
		delete resultTouch;
		resultTouch = 0x0; 
		return kFALSE;
	}
	if (valueTouch->GetString()!="1"){
		Log("SHUTTLE",Form("Failing the touching command, returning without touching for run %i",GetCurrentRun()));
		delete resultTouch;
		resultTouch = 0x0; 
		return kFALSE;
	}
	delete resultTouch;
	resultTouch = 0x0; 
	Log("SHUTTLE", "Sucessfully touched the file");
	return kTRUE;
}
//______________________________________________________________________________________________
UInt_t AliShuttle::GetStartTimeDCSQuery()
{
	// Return Start Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return GetCurrentStartTime()-fConfig->GetDCSQueryOffset();
}
//______________________________________________________________________________________________
UInt_t AliShuttle::GetEndTimeDCSQuery()
{
	// Return End Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return GetCurrentEndTime()+fConfig->GetDCSQueryOffset();
}
//______________________________________________________________________________________________
void AliShuttle::SendMLFromDet(const char* value)
{
	// 
	// Sending an information coming from the current detector to ML
	//
	
	TMonaLisaText  mlText(Form("%s_RunCondition", fCurrentDetector.Data()), value);

	TList mlList;
	mlList.Add(&mlText);

	TString mlID;
	mlID.Form("%d", GetCurrentRun());
	fMonaLisa->SendParameters(&mlList, mlID);

	return;
}
//______________________________________________________________________________________________
TString* AliShuttle::GetLTUConfig(const char* det)
{
	// 
	// Getting ltuFineDelay1, ltuFineDelay2, ltuBCDelay for detector det from logbook_detectors table in logbook
	//
	
	if (!Connect(4)) 
		return 0;

	TString sqlQuery;
	sqlQuery.Form("select LTUFineDelay1, LTUFineDelay2, LTUBCDelayAdd from logbook_detectors WHERE run_number = %d and detector = \"%s\";", GetCurrentRun(),det);

	TSQLResult* result = fServer[4]->Query(sqlQuery);
	if (!result){
		Log("SHUTTLE","ERROR: No result found for the LTU configuration query");
		return 0x0;
	}
	if (result->GetRowCount() == 0){
		Log("SHUTTLE",Form("ERROR: LTU configuration not found in logbook_detectors for detector %s, returning null pointer",det));
		delete result;
		return 0x0;
	}
	if (result->GetFieldCount() != 3){
		Log("SHUTTLE",Form("ERROR: not all the required fields are there for the LTU configuration for detector %s (only %d found), returning a null pointer",det, result->GetFieldCount()));
		delete result;
		return 0x0;
	}
	TSQLRow* row = result->Next();
	if (!row){
		Printf("ERROR: Could not receive logbook_detectors data, returning null pointer");
		delete result;
		return 0x0;
	}
	TString* ltuConfigString = new TString[3];

	ltuConfigString[0] = row->GetField(0);
	ltuConfigString[1] = row->GetField(1);
	ltuConfigString[2] = row->GetField(2);

	return ltuConfigString;

}
