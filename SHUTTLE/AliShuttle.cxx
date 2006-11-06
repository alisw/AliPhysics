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

/*
$Log$
Revision 1.18  2006/10/20 15:22:59  jgrosseo
o) Adding time out to the execution of the preprocessors: The Shuttle forks and the parent process monitors the child
o) Merging Collect, CollectAll, CollectNew function
o) Removing implementation of empty copy constructors (declaration still there!)

Revision 1.17  2006/10/05 16:20:55  jgrosseo
adapting to new CDB classes

Revision 1.16  2006/10/05 15:46:26  jgrosseo
applying to the new interface

Revision 1.15  2006/10/02 16:38:39  jgrosseo
update (alberto):
fixed memory leaks
storing of objects that failed to be stored to the grid before
interfacing of shuttle status table in daq system

Revision 1.14  2006/08/29 09:16:05  jgrosseo
small update

Revision 1.13  2006/08/15 10:50:00  jgrosseo
effc++ corrections (alberto)

Revision 1.12  2006/08/08 14:19:29  jgrosseo
Update to shuttle classes (Alberto)

- Possibility to set the full object's path in the Preprocessor's and
Shuttle's  Store functions
- Possibility to extend the object's run validity in the same classes
("startValidity" and "validityInfinite" parameters)
- Implementation of the StoreReferenceData function to store reference
data in a dedicated CDB storage.

Revision 1.11  2006/07/21 07:37:20  jgrosseo
last run is stored after each run

Revision 1.10  2006/07/20 09:54:40  jgrosseo
introducing status management: The processing per subdetector is divided into several steps,
after each step the status is stored on disk. If the system crashes in any of the steps the Shuttle
can keep track of the number of failures and skips further processing after a certain threshold is
exceeded. These thresholds can be configured in LDAP.

Revision 1.9  2006/07/19 10:09:55  jgrosseo
new configuration, accesst to DAQ FES (Alberto)

Revision 1.8  2006/07/11 12:44:36  jgrosseo
adding parameters for extended validity range of data produced by preprocessor

Revision 1.7  2006/07/10 14:37:09  jgrosseo
small fix + todo comment

Revision 1.6  2006/07/10 13:01:41  jgrosseo
enhanced storing of last sucessfully processed run (alberto)

Revision 1.5  2006/07/04 14:59:57  jgrosseo
revision of AliDCSValue: Removed wrapper classes, reduced storage size per value by factor 2

Revision 1.4  2006/06/12 09:11:16  jgrosseo
coding conventions (Alberto)

Revision 1.3  2006/06/06 14:26:40  jgrosseo
o) removed files that were moved to STEER
o) shuttle updated to follow the new interface (Alberto)

Revision 1.2  2006/03/07 07:52:34  hristov
New version (B.Yordanov)

Revision 1.6  2005/11/19 17:19:14  byordano
RetrieveDATEEntries and RetrieveConditionsData added

Revision 1.5  2005/11/19 11:09:27  byordano
AliShuttle declaration added

Revision 1.4  2005/11/17 17:47:34  byordano
TList changed to TObjArray

Revision 1.3  2005/11/17 14:43:23  byordano
import to local CVS

Revision 1.1.1.1  2005/10/28 07:33:58  hristov
Initial import as subdirectory in AliRoot

Revision 1.2  2005/09/13 08:41:15  byordano
default startTime endTime added

Revision 1.4  2005/08/30 09:13:02  byordano
some docs added

Revision 1.3  2005/08/29 21:15:47  byordano
some docs added

*/

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

#include <fstream>

#include <sys/types.h>
#include <sys/wait.h>

ClassImp(AliShuttle)

TString AliShuttle::fgkMainCDB("alien://folder=ShuttleCDB");
TString AliShuttle::fgkLocalCDB("local://LocalShuttleCDB");
TString AliShuttle::fgkMainRefStorage("alien://folder=ShuttleReference");
TString AliShuttle::fgkLocalRefStorage("local://LocalReferenceStorage");

Bool_t AliShuttle::fgkProcessDCS(kTRUE); 

const char* AliShuttle::fgkShuttleTempDir = gSystem->ExpandPathName("$ALICE_ROOT/SHUTTLE/temp");
const char* AliShuttle::fgkShuttleLogDir = gSystem->ExpandPathName("$ALICE_ROOT/SHUTTLE/log");

//______________________________________________________________________________________________
AliShuttle::AliShuttle(const AliShuttleConfig* config,
		UInt_t timeout, Int_t retries):
fConfig(config),
fTimeout(timeout), fRetries(retries),
fPreprocessorMap(),
fLogbookEntry(0),
fCurrentDetector(),
fStatusEntry(0),
fGridError(kFALSE),
fMonitoringMutex(0),
fLastActionTime(0),
fLastAction()
{
	//
	// config: AliShuttleConfig used
	// timeout: timeout used for AliDCSClient connection
	// retries: the number of retries in case of connection error.
	//

	if (!fConfig->IsValid()) AliFatal("********** !!!!! Invalid configuration !!!!! **********");
	for(int iSys=0;iSys<3;iSys++) {
		fServer[iSys]=0;
		fFESlist[iSys].SetOwner(kTRUE);
	}
	fPreprocessorMap.SetOwner(kTRUE);
	
	fMonitoringMutex = new TMutex();
}

//______________________________________________________________________________________________
AliShuttle::~AliShuttle()
{
// destructor

	fPreprocessorMap.DeleteAll();
	for(int iSys=0;iSys<3;iSys++)
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
UInt_t AliShuttle::Store(const AliCDBPath& path, TObject* object,
		AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // Stores a CDB object in the storage for offline reconstruction. Objects that are not needed for
  // offline reconstruction, but should be stored anyway (e.g. for debugging) should NOT be stored
  // using this function. Use StoreReferenceData instead!
  // It calls WriteToCDB function which perform actual storage

	return WriteToCDB(fgkMainCDB, fgkLocalCDB, path, object,
				metaData, validityStart, validityInfinite);

}

//______________________________________________________________________________________________
UInt_t AliShuttle::StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData)
{
  // Stores a CDB object in the storage for reference data. This objects will not be available during
  // offline reconstrunction. Use this function for reference data only!
  // It calls WriteToCDB function which perform actual storage

	return WriteToCDB(fgkMainRefStorage, fgkLocalRefStorage, path, object, metaData);

}

//______________________________________________________________________________________________
UInt_t AliShuttle::WriteToCDB(const char* mainUri, const char* localUri,
			const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
			Int_t validityStart, Bool_t validityInfinite)
{
  // write object into the CDB. Parameters are passed by Store and StoreReferenceData functions.
  // The parameters are:
  //   1) Uri of the main storage (Grid)
  //   2) Uri of the backup storage (Local)
  //   3) the object's path.
  //   4) the object to be stored
  //   5) the metaData to be associated with the object
  //   6) the validity start run number w.r.t. the current run,
  //      if the data is valid only for this run leave the default 0
  //   7) specifies if the calibration data is valid for infinity (this means until updated),
  //      typical for calibration runs, the default is kFALSE
  //
  // returns 0 if fail
  // 	     1 if stored in main (Grid) storage
  // 	     2 if stored in backup (Local) storage

	const char* cdbType = (mainUri == fgkMainCDB) ? "CDB" : "Reference";

	Int_t firstRun = GetCurrentRun() - validityStart;
  	if(firstRun < 0) {
		AliError("First valid run happens to be less than 0! Setting it to 0.");
		firstRun=0;
  	}

	Int_t lastRun = -1;
	if(validityInfinite) {
		lastRun = AliCDBRunRange::Infinity();
	} else {
		lastRun = GetCurrentRun();
	}

	AliCDBId id(path, firstRun, lastRun, -1, -1);

	if(! dynamic_cast<TObjString*> (metaData->GetProperty("RunUsed(TObjString)"))){
		TObjString runUsed = Form("%d", GetCurrentRun());
		metaData->SetProperty("RunUsed(TObjString)",&runUsed);
	}

	UInt_t result = 0;

	if (!(AliCDBManager::Instance()->GetStorage(mainUri))) {
		AliError(Form("WriteToCDB - Cannot activate main %s storage", cdbType));
	} else {
		result = (UInt_t) AliCDBManager::Instance()->GetStorage(mainUri)
					->Put(object, id, metaData);
	}

	if(!result) {

		Log(fCurrentDetector,
			Form("WriteToCDB - Problem with main %s storage. Putting <%s> into backup storage",
				cdbType, path.GetPath().Data()));

		// Set Grid version to current run number, to ease retrieval later
		id.SetVersion(GetCurrentRun());

		result = AliCDBManager::Instance()->GetStorage(localUri)
					->Put(object, id, metaData);

		if(result) {
			result = 2;
      			fGridError = kTRUE;
		}else{
			Log(fCurrentDetector, "WriteToCDB - Can't store data!");
		}
	}

	return result;

}

//______________________________________________________________________________________________
AliShuttleStatus* AliShuttle::ReadShuttleStatus()
{
// Reads the AliShuttleStatus from the CDB

	if (fStatusEntry){
		delete fStatusEntry;
		fStatusEntry = 0;
	}

	fStatusEntry = AliCDBManager::Instance()->GetStorage(AliShuttle::GetLocalCDB())
		->Get(Form("/SHUTTLE/STATUS/%s", fCurrentDetector.Data()), GetCurrentRun());

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
// writes the status for one subdetector

	if (fStatusEntry){
		delete fStatusEntry;
		fStatusEntry = 0;
	}

	Int_t run = GetCurrentRun();

	AliCDBId id(AliCDBPath("SHUTTLE", "STATUS", fCurrentDetector), run, run);

	fStatusEntry = new AliCDBEntry(status, id, new AliCDBMetaData);
	fStatusEntry->SetOwner(1);

	UInt_t result = AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);

	if (!result) {
		AliError(Form("WriteShuttleStatus for %s, run %d failed", fCurrentDetector.Data(), run));
		return kFALSE;
	}

	return kTRUE;
}

//______________________________________________________________________________________________
void AliShuttle::UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount)
{
  // changes the AliShuttleStatus for the given detector and run to the given status

	if (!fStatusEntry){
		AliError("UNEXPECTED: fStatusEntry empty");
		return;
	}

	AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (fStatusEntry->GetObject());

	if (!status){
		AliError("UNEXPECTED: status could not be read from current CDB entry");
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

	AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);
}
//______________________________________________________________________________________________
Bool_t AliShuttle::ContinueProcessing()
{
// this function reads the AliShuttleStatus information from CDB and
// checks if the processing should be continued
// if yes it returns kTRUE and updates the AliShuttleStatus with nextStatus

	AliShuttleLogbookEntry::Status entryStatus =
		fLogbookEntry->GetDetectorStatus(fCurrentDetector);

	if(entryStatus != AliShuttleLogbookEntry::kUnprocessed) {
		Log("SHUTTLE", Form("ContinueProcessing - %s is %s",
				fCurrentDetector.Data(),
				fLogbookEntry->GetDetectorStatusName(entryStatus)));
		return kFALSE;
	}

	// if we get here, according to Shuttle logbook subdetector is in UNPROCESSED state
	AliShuttleStatus* status = ReadShuttleStatus();
	if (!status) {
		// first time
		Log("SHUTTLE", Form("ContinueProcessing - %s: Processing first time",
				fCurrentDetector.Data()));
		status = new AliShuttleStatus(AliShuttleStatus::kStarted);
		return WriteShuttleStatus(status);
	}

	// The following two cases shouldn't happen if Shuttle Logbook was correctly updated.
	// If it happens it may mean Logbook updating failed... let's do it now!
	if (status->GetStatus() == AliShuttleStatus::kDone ||
	    status->GetStatus() == AliShuttleStatus::kFailed){
		Log("SHUTTLE", Form("ContinueProcessing - %s is already %s. Updating Shuttle Logbook",
					fCurrentDetector.Data(),
					status->GetStatusName(status->GetStatus())));
		UpdateShuttleLogbook(fCurrentDetector.Data(),
					status->GetStatusName(status->GetStatus()));
		return kFALSE;
	}

	if (status->GetStatus() == AliShuttleStatus::kStoreFailed) {
		Log("SHUTTLE",
			Form("ContinueProcessing - %s: Grid storage of one or more objects failed. Trying again now",
				fCurrentDetector.Data()));
		if(TryToStoreAgain()){
			Log(fCurrentDetector.Data(), "ContinueProcessing - All objects successfully stored into OCDB");
			UpdateShuttleStatus(AliShuttleStatus::kDone);
			UpdateShuttleLogbook(fCurrentDetector.Data(), "DONE");
		} else {
			Log("SHUTTLE",
				Form("ContinueProcessing - %s: Grid storage failed again",
					fCurrentDetector.Data()));
		}
		return kFALSE;
	}

	// if we get here, there is a restart

	// abort conditions
	if (status->GetCount() >= fConfig->GetMaxRetries()) {
		Log("SHUTTLE",
			Form("ContinueProcessing - %s failed %d times in status %s - Updating Shuttle Logbook",
				fCurrentDetector.Data(),
				status->GetCount(), status->GetStatusName()));
		UpdateShuttleLogbook(fCurrentDetector.Data(), "FAILED");
		return kFALSE;
	}

	Log("SHUTTLE", Form("ContinueProcessing - %s: restarting. Aborted before with %s. Retry number %d.",
  			fCurrentDetector.Data(),
			status->GetStatusName(), status->GetCount()));

	UpdateShuttleStatus(AliShuttleStatus::kStarted, kTRUE);

	return kTRUE;
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

	if(!entry) return kFALSE;

	fLogbookEntry = entry;

	if(fLogbookEntry->IsDone()){
		Log("SHUTTLE","Process - Shuttle is already DONE. Updating logbook");
		UpdateShuttleLogbook("shuttle_done");
		fLogbookEntry = 0;
		return kTRUE;
	}


	AliInfo(Form("\n\n \t\t\t^*^*^*^*^*^*^*^*^*^*^*^* run %d: START ^*^*^*^*^*^*^*^*^*^*^*^* \n",
					GetCurrentRun()));

	fLogbookEntry->Print("all");

	// Initialization
	Bool_t hasError = kFALSE;
	for(Int_t iSys=0;iSys<3;iSys++) fFESCalled[iSys]=kFALSE;

	AliCDBStorage *mainCDBSto = AliCDBManager::Instance()->GetStorage(fgkMainCDB);
	if(mainCDBSto) mainCDBSto->QueryCDB(GetCurrentRun());
	AliCDBStorage *mainRefSto = AliCDBManager::Instance()->GetStorage(fgkMainRefStorage);
	if(mainRefSto) mainRefSto->QueryCDB(GetCurrentRun());

	// Loop on detectors in the configuration
	TIter iter(fConfig->GetDetectors());
	TObjString* aDetector = 0;

	while ((aDetector = (TObjString*) iter.Next())) {
		fCurrentDetector = aDetector->String();

		if (!fConfig->HostProcessDetector(fCurrentDetector)) continue;

		AliPreprocessor* aPreprocessor =
			dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(fCurrentDetector));
		if(!aPreprocessor){
			Log("SHUTTLE",Form("Process - %s: no preprocessor registered. Skipping", 
							fCurrentDetector.Data()));
			continue;
		}

		if (ContinueProcessing() == kFALSE) continue;

		AliInfo(Form("\n\n \t\t\t****** run %d - %s: START  ******",
						GetCurrentRun(), aDetector->GetName()));


    Int_t pid = fork();

    if (pid < 0)
    {
      Log("SHUTTLE", "ERROR: Forking failed");
    }
    else if (pid > 0)
    {
      // parent
      AliInfo(Form("In parent process of %d - %s: Starting monitoring", GetCurrentRun(), aDetector->GetName()));

      Long_t begin = time(0);

      int status; // to be used with waitpid, on purpose an int (not Int_t)!
      while (waitpid(pid, &status, WNOHANG) == 0)
      {
        Long_t expiredTime = time(0) - begin;

        if (expiredTime > fConfig->GetPPTimeOut())
        {
          Log("SHUTTLE", Form("Process time out. Run time: %d seconds. Killing...", expiredTime));

          kill(pid, 9);

          hasError = kTRUE;

          gSystem->Sleep(1000);
        }
        else
        {
          if (expiredTime % 60 == 0)
            Log("SHUTTLE", Form("Checked process. Run time: %d seconds.", expiredTime));

          gSystem->Sleep(1000);
        }
      }

      AliInfo(Form("In parent process of %d - %s: Client has terminated.", GetCurrentRun(), aDetector->GetName()));

      if (WIFEXITED(status))
      {
        Int_t returnCode = WEXITSTATUS(status);

        Log("SHUTTLE", Form("The return code is %d", returnCode));

        if (returnCode != 0)
          hasError = kTRUE;
      }
    }
    else if (pid == 0)
    {
      // client
      AliInfo(Form("In client process of %d - %s", GetCurrentRun(), aDetector->GetName()));

      UInt_t result = ProcessCurrentDetector();

      Int_t returnCode = 0; // will be set to 1 in case of an error

      if (!result) {
        returnCode = 1;
        AliInfo(Form("\n \t\t\t****** run %d - %s: PREPROCESSOR ERROR ****** \n\n",
                GetCurrentRun(), aDetector->GetName()));
      }
      else if(result == 2) {
        AliInfo(Form("\n \t\t\t****** run %d - %s: STORAGE ERROR ****** \n\n",
                GetCurrentRun(), aDetector->GetName()));
      } else {
        AliInfo(Form("\n \t\t\t****** run %d - %s: DONE ****** \n\n",
                GetCurrentRun(), aDetector->GetName()));
      }

      if (result > 0)
      {
        // Process successful: Update time_processed field in FES logbooks!
        if(fFESCalled[kDAQ]) {
          if (UpdateDAQTable() == kFALSE)
            returnCode = 1;
          fFESlist[kDAQ].Clear();
        }
        //if(fFESCalled[kDCS]) {
        //  if (UpdateDCSTable(aDetector->GetName()) == kFALSE)
        //    returnCode = 1;
        //  fFESlist[kDCS].Clear();
        //}
        //if(fFESCalled[kHLT]) {
        //  if (UpdateHLTTable(aDetector->GetName()) == kFALSE)
        //    returnCode = 1;
        //	fFESlist[kHLT].Clear();
        //}
      }

      AliInfo(Form("Client process of %d - %s is exiting now with %d.", GetCurrentRun(), aDetector->GetName(), returnCode));

      // the client exits here
      gSystem->Exit(returnCode);

      AliError("We should never get here!!!");
    }
	}

	AliInfo(Form("\n\n \t\t\t^*^*^*^*^*^*^*^*^*^*^*^* run %d: FINISH ^*^*^*^*^*^*^*^*^*^*^*^* \n",
							GetCurrentRun()));

	//check if shuttle is done for this run, if so update logbook
	TObjArray checkEntryArray;
	checkEntryArray.SetOwner(1);
	TString whereClause = Form("where run=%d",GetCurrentRun());
	if(QueryShuttleLogbook(whereClause.Data(), checkEntryArray)) {

		AliShuttleLogbookEntry* checkEntry = dynamic_cast<AliShuttleLogbookEntry*>
							(checkEntryArray.At(0));

		if(checkEntry && checkEntry->IsDone()){
			Log("SHUTTLE","Process - Shuttle is DONE. Updating logbook");
			UpdateShuttleLogbook("shuttle_done");
		}
	}

	fLogbookEntry = 0;

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
UInt_t AliShuttle::ProcessCurrentDetector()
{
	//
        // Makes data retrieval just for a specific detector (fCurrentDetector).
	// Threre should be a configuration for this detector.

	AliInfo(Form("Retrieving values for %s, run %d", fCurrentDetector.Data(), GetCurrentRun()));

	UpdateShuttleStatus(AliShuttleStatus::kDCSStarted);

	TString host(fConfig->GetDCSHost(fCurrentDetector));
	Int_t port = fConfig->GetDCSPort(fCurrentDetector);

	TIter iter(fConfig->GetDCSAliases(fCurrentDetector));
	TObjString* anAlias;
	TMap aliasMap;
	aliasMap.SetOwner(1);

	Bool_t aDCSError = kFALSE;
	fGridError = kFALSE;

	while ((anAlias = (TObjString*) iter.Next())) {
		TObjArray *valueSet = new TObjArray();
		valueSet->SetOwner(1);
		// TODO Test only... I've added a flag that allows to
		// exclude DCS archive DB query
		if(fgkProcessDCS){
			AliInfo("Querying DCS archive DB data...");
			aDCSError = (GetValueSet(host, port, anAlias->String(), valueSet) == 0);
		} else {
			AliInfo(Form("Skipping DCS processing. Port = %d",port));
			aDCSError = kFALSE;
		}
		if(!aDCSError) {
			aliasMap.Add(anAlias->Clone(), valueSet);
		}else{
			Log(fCurrentDetector, Form("ProcessCurrentDetector - Error while retrieving alias %s",
					anAlias->GetName()));
			UpdateShuttleStatus(AliShuttleStatus::kDCSError, kTRUE);
			aliasMap.DeleteAll();
			return 0;
		}
	}

	// DCS Archive DB processing successful. Call Preprocessor!
	UpdateShuttleStatus(AliShuttleStatus::kPPStarted);

	AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(fCurrentDetector));

	aPreprocessor->Initialize(GetCurrentRun(), GetCurrentStartTime(), GetCurrentEndTime());
	UInt_t aPPResult = aPreprocessor->Process(&aliasMap);

	UInt_t returnValue = 0;
	if (aPPResult == 0) { // Preprocessor error
		UpdateShuttleStatus(AliShuttleStatus::kPPError);
		returnValue = 0;
	} else if (fGridError == kFALSE) { // process and Grid storage ok!
    		UpdateShuttleStatus(AliShuttleStatus::kDone);
		UpdateShuttleLogbook(fCurrentDetector, "DONE");
		Log(fCurrentDetector.Data(),
			"ProcessCurrentDetector - Preprocessor and Grid storage ended successfully");
		returnValue = 1;
        } else { // Grid storage error (process ok, but object put in local storage)
     		UpdateShuttleStatus(AliShuttleStatus::kStoreFailed);
		returnValue = 2;
	}

	aliasMap.DeleteAll();

	return returnValue;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::QueryShuttleLogbook(const char* whereClause,
		TObjArray& entries)
{
// Query DAQ's Shuttle logbook and fills detector status object.
// Call QueryRunParameters to query DAQ logbook for run parameters.

	// check connection, in case connect
	if(!Connect(kDAQ)) return kFALSE;

	TString sqlQuery;
	sqlQuery = Form("select * from logbook_shuttle %s order by run", whereClause);

	TSQLResult* aResult = fServer[kDAQ]->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	if(aResult->GetRowCount() == 0) {
		if(sqlQuery.Contains("where shuttle_done=0")){
			Log("SHUTTLE", "QueryShuttleLogbook - All runs in Shuttle Logbook are already DONE");
			delete aResult;
			return kTRUE;
		} else {
			AliError("No entries in Shuttle Logbook match request");
			delete aResult;
			return kFALSE;
		}
	}

	// TODO Check field count!
	const UInt_t nCols = 24;
	if (aResult->GetFieldCount() != (Int_t) nCols) {
		AliError("Invalid SQL result field number!");
		delete aResult;
		return kFALSE;
	}

	entries.SetOwner(1);

	TSQLRow* aRow;
	while ((aRow = aResult->Next())) {
		TString runString(aRow->GetField(0), aRow->GetFieldLength(0));
		Int_t run = runString.Atoi();

		AliShuttleLogbookEntry *entry = QueryRunParameters(run);
		if (!entry)
			continue;

		// loop on detectors
		for(UInt_t ii = 0; ii < nCols; ii++)
			entry->SetDetectorStatus(aResult->GetFieldName(ii), aRow->GetField(ii));

		entries.AddLast(entry);
		delete aRow;
	}

	if(sqlQuery.Contains("where shuttle_done=0"))
		Log("SHUTTLE", Form("QueryShuttleLogbook - Found %d unprocessed runs in Shuttle Logbook",
							entries.GetEntriesFast()));
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
	if (!Connect(kDAQ))
		return 0;

	TString sqlQuery;
	sqlQuery.Form("select * from logbook where run=%d", run);

	TSQLResult* aResult = fServer[kDAQ]->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}

	if (aResult->GetRowCount() == 0) {
		Log("SHUTTLE", Form("QueryRunParameters - No entry in DAQ Logbook for run %d. Skipping", run));
		delete aResult;
		return 0;
	}

	if (aResult->GetRowCount() > 1) {
		AliError(Form("More than one entry in DAQ Logbook for run %d. Skipping", run));
		delete aResult;
		return 0;
	}

	TSQLRow* aRow = aResult->Next();
	if (!aRow)
	{
		AliError(Form("Could not retrieve row for run %d. Skipping", run));
		delete aResult;
		return 0;
	}

	AliShuttleLogbookEntry* entry = new AliShuttleLogbookEntry(run);

	for (Int_t ii = 0; ii < aResult->GetFieldCount(); ii++)
		entry->SetRunParameter(aResult->GetFieldName(ii), aRow->GetField(ii));

	UInt_t startTime = entry->GetStartTime();
	UInt_t endTime = entry->GetEndTime();

	if (!startTime || !endTime || startTime > endTime) {
		Log("SHUTTLE",
			Form("QueryRunParameters - Invalid parameters for Run %d: startTime = %d, endTime = %d",
				run, startTime, endTime));
		delete entry;
		delete aRow;
		delete aResult;
		return 0;
	}

	delete aRow;
	delete aResult;

	return entry;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::TryToStoreAgain()
{
  // Called in case the detector failed to store the object in Grid OCDB
  // It tries to store the object again, if it does not find more recent and overlapping objects
  // Calls underlying TryToStoreAgain(const char*) function twice, for OCDB and Reference storage.

	AliInfo("Trying to store OCDB data again...");
	Bool_t resultCDB = TryToStoreAgain(fgkMainCDB);

	AliInfo("Trying to store reference data again...");
	Bool_t resultRef = TryToStoreAgain(fgkMainRefStorage);

	return resultCDB && resultRef;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::TryToStoreAgain(TString& gridURI)
{
  // Called by TryToStoreAgain(), performs actual storage retry

	TObjArray* gridIds=0;

	Bool_t result = kTRUE;

	const char* type = 0;
	TString backupURI;
	if(gridURI == fgkMainCDB) {
		type = "OCDB";
		backupURI = fgkLocalCDB;
	} else if(gridURI == fgkMainRefStorage) {
		type = "reference";
		backupURI = fgkLocalRefStorage;
	} else {
		AliError(Form("Invalid storage URI: %s", gridURI.Data()));
		return kFALSE;
	}

	AliCDBManager* man = AliCDBManager::Instance();

	AliCDBStorage *gridSto = man->GetStorage(gridURI);
	if(!gridSto) {
		Log(fCurrentDetector.Data(),
			Form("TryToStoreAgain - cannot activate main %s storage", type));
		return kFALSE;
	}

	gridIds = gridSto->GetQueryCDBList();

	// get objects previously stored in local CDB
	AliCDBStorage *backupSto = man->GetStorage(backupURI);
	AliCDBPath aPath(GetOfflineDetName(fCurrentDetector.Data()),"*","*");
	// Local objects were stored with current run as Grid version!
	TList* localEntries = backupSto->GetAll(aPath.GetPath(), GetCurrentRun(), GetCurrentRun());
	localEntries->SetOwner(1);

	// loop on local stored objects
	TIter localIter(localEntries);
	AliCDBEntry *aLocEntry = 0;
	while((aLocEntry = dynamic_cast<AliCDBEntry*> (localIter.Next()))){
		aLocEntry->SetOwner(1);
		AliCDBId aLocId = aLocEntry->GetId();
		aLocEntry->SetVersion(-1);
		aLocEntry->SetSubVersion(-1);

		// loop on Grid valid Id's
		Bool_t store = kTRUE;
		TIter gridIter(gridIds);
		AliCDBId* aGridId = 0;
		while((aGridId = dynamic_cast<AliCDBId*> (gridIter.Next()))){
			// If local object is valid up to infinity we store it anyway
			// TODO This does not work! It may hide more recent objects...
			if(aLocId.GetLastRun() == AliCDBRunRange::Infinity()) {
				// TODO Check that it won't hide more recent files! how????
				break;
			}
			if(aGridId->GetPath() != aLocId.GetPath()) continue;
			// skip all objects valid up to infinity
			if(aGridId->GetLastRun() == AliCDBRunRange::Infinity()) continue;
			// if we get here, it means there's already some more recent object stored on Grid!
			store = kFALSE;
			break;
		}

		if(!store){
			Log(fCurrentDetector.Data(),
				Form("TryToStoreAgain - A more recent object already exists in %s storage: <%s>",
					type, aGridId->ToString().Data()));
			// removing local filename...
			// TODO maybe it's better not to remove it, it was not copied to the Grid!
			TString filename;
			backupSto->IdToFilename(aLocId, filename);
			AliInfo(Form("Removing local file %s", filename.Data()));
			gSystem->Exec(Form("rm %s",filename.Data()));
			continue;
		}

		// If we get here, the file can be stored!
		Bool_t storeOk = gridSto->Put(aLocEntry);
		if(storeOk){
			Log(fCurrentDetector.Data(),
				Form("TryToStoreAgain - Object <%s> successfully put into %s storage",
					aLocId.ToString().Data(), type));

			// removing local filename...
			TString filename;
			backupSto->IdToFilename(aLocId, filename);
			AliInfo(Form("Removing local file %s", filename.Data()));
			gSystem->Exec(Form("rm %s", filename.Data()));
			continue;
		} else	{
			Log(fCurrentDetector.Data(),
				Form("TryToStoreAgain - Grid %s storage of object <%s> failed again",
					type, aLocId.ToString().Data()));
			result = kFALSE;
		}
	}
	localEntries->Clear();

	return result;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::GetValueSet(const char* host, Int_t port, const char* alias,
				TObjArray* valueSet)
{
// Retrieve all "alias" data points from the DCS server
// host, port: TSocket connection parameters
// alias: name of the alias
// valueSet: array of retrieved AliDCSValue's

	AliDCSClient client(host, port, fTimeout, fRetries);
	if (!client.IsConnected()) {
		return kFALSE;
	}

	Int_t result = client.GetAliasValues(alias,
		GetCurrentStartTime(), GetCurrentEndTime(), valueSet);

	if (result < 0) {
		Log(fCurrentDetector.Data(), Form("GetValueSet - Can't get '%s'! Reason: %s",
			alias, AliDCSClient::GetErrorString(result)));

		if (result == AliDCSClient::fgkServerError) {
			Log(fCurrentDetector.Data(), Form("GetValueSet - Server error: %s",
				client.GetServerError().Data()));
		}

		return kFALSE;
	}

	return kTRUE;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetFile(Int_t system, const char* detector,
		const char* id, const char* source)
{
// Get calibration file from file exchange servers
// calls specific getter according to system index (kDAQ, kDCS, kHLT)

	switch(system){
		case kDAQ:
			return GetDAQFileName(detector, id, source);
			break;
		case kDCS:
			return GetDCSFileName(detector, id, source);
			break;
		case kHLT:
			return GetHLTFileName(detector, id, source);
			break;
		default:
			AliError(Form("No valid system index: %d",system));
	}

	return 0;
}

//______________________________________________________________________________________________
TList* AliShuttle::GetFileSources(Int_t system, const char* detector, const char* id)
{
// Get sources producing the condition file Id from file exchange servers
// calls specific getter according to system index (kDAQ, kDCS, kHLT)

	switch(system){
		case kDAQ:
			return GetDAQFileSources(detector, id);
			break;
		case kDCS:
			return GetDCSFileSources(detector, id);
			break;
		case kHLT:
			return GetHLTFileSources(detector, id);
			break;
		default:
			AliError(Form("No valid system index: %d",system));
	}

	return NULL;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Connect(Int_t system)
{
// Connect to MySQL Server of the system's FES logbook
// DAQ Logbook, Shuttle Logbook and DAQ FES Logbook are on the same host

	// check connection: if already connected return
	if(fServer[system] && fServer[system]->IsConnected()) return kTRUE;

	TString aFESlbHost= Form("mysql://%s", fConfig->GetFESlbHost(system));

	fServer[system] = TSQLServer::Connect(aFESlbHost,
			fConfig->GetFESlbUser(system),
			fConfig->GetFESlbPass(system));
	if (!fServer[system] || !fServer[system]->IsConnected()) {
		AliError(Form("Can't establish connection to FES logbook for %s",
					AliShuttleInterface::GetSystemName(system)));
		if(fServer[system]) delete fServer[system];
		return kFALSE;
	}

	// Get tables
	// TODO in the configuration should the table name be there too?
	TSQLResult* aResult=0;
	switch(system){
		case kDAQ:
			aResult = fServer[kDAQ]->GetTables("REFSYSLOG");
			break;
		case kDCS:
			//aResult = fServer[kDCS]->GetTables("REFSYSLOG");
			break;
		case kHLT:
			//aResult = fServer[kHLT]->GetTables("REFSYSLOG");
			break;
		default:
			break;
	}

	delete aResult;
	return kTRUE;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetDAQFileName(const char* detector, const char* id, const char* source)
{
// Retrieves a file from the DAQ FES.
// First queris the DAQ logbook_fs for the DAQ file name, using the run, detector, id and source info
// then calls RetrieveDAQFile(DAQfilename) for actual copy to local disk
// run: current run being processed (given by Logbook entry fLogbookEntry)
// detector: the Preprocessor name
// id: provided as a parameter by the Preprocessor
// source: provided by the Preprocessor through GetFileSources function

	// check connection, in case connect
	if (!Connect(kDAQ))
	{
		Log(detector, "GetDAQFileName - Couldn't connect to DAQ Logbook");
		return 0;
	}

	// Query preparation
	TString sqlQueryStart = "select filePath from logbook_fs where";
	TString whereClause = Form("run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\"",
				GetCurrentRun(), detector, id, source);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult = 0;
	aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
	if (!aResult) {
		Log(detector, Form("GetDAQFileName - Can't execute SQL query for: id = %s, source = %s",
				id, source));
		return 0;
	}

	if(aResult->GetRowCount() == 0)
	{
		Log(detector,
			Form("GetDAQFileName - No entry in FES table for: id = %s, source = %s",
				id, source));
		delete aResult;
		return 0;
	}

	if (aResult->GetRowCount() > 1) {
		Log(detector,
			Form("GetDAQFileName - More than one entry in FES table for: id = %s, source = %s",
				id, source));
		delete aResult;
		return 0;
	}

	TSQLRow* aRow = dynamic_cast<TSQLRow*> (aResult->Next());

	if (!aRow){
		Log(detector, Form("GetDAQFileName - Empty set result from query: id = %s, source = %s",
				id, source));
		delete aResult;
		return 0;
	}

	TString filePath(aRow->GetField(0), aRow->GetFieldLength(0));

	delete aResult;
	delete aRow;

	AliDebug(2, Form("filePath = %s",filePath.Data()));

	// retrieved file is renamed to make it unique
	TString localFileName = Form("%s_%d_%s_%s.shuttle",
					detector, GetCurrentRun(), id, source);

	// file retrieval from DAQ FES
	Bool_t result = RetrieveDAQFile(filePath.Data(), localFileName.Data());
	if(!result) {
		Log(detector, Form("GetDAQFileName - Copy of file %s from DAQ FES failed", filePath.Data()));
		return 0;
	} else {
		AliInfo(Form("File %s copied from DAQ FES into %s/%s",
			filePath.Data(), fgkShuttleTempDir, localFileName.Data()));
	}


	fFESCalled[kDAQ]=kTRUE;
	TObjString *fileParams = new TObjString(Form("%s_!?!_%s", id, source));
	fFESlist[kDAQ].Add(fileParams);

	return localFileName.Data();

}

//______________________________________________________________________________________________
Bool_t AliShuttle::RetrieveDAQFile(const char* daqFileName, const char* localFileName)
{

	// check temp directory: trying to cd to temp; if it does not exist, create it
	AliDebug(2, Form("Copy file %s from DAQ FES into folder %s and rename it as %s",
			daqFileName,fgkShuttleTempDir, localFileName));

	void* dir = gSystem->OpenDirectory(fgkShuttleTempDir);
	if (dir == NULL) {
		if (gSystem->mkdir(fgkShuttleTempDir, kTRUE)) {
			AliError(Form("Can't open directory <%s>", fgkShuttleTempDir));
			return kFALSE;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	TString baseDAQFESFolder = "DAQ";
	TString command = Form("scp %s@%s:%s/%s %s/%s",
		fConfig->GetFESUser(kDAQ),
		fConfig->GetFESHost(kDAQ),
		baseDAQFESFolder.Data(),
		daqFileName,
		fgkShuttleTempDir,
		localFileName);

	AliDebug(2, Form("%s",command.Data()));

	UInt_t nRetries = 0;
	UInt_t maxRetries = 3;

	// copy!! if successful TSystem::Exec returns 0
	while(nRetries++ < maxRetries) {
		AliDebug(2, Form("Trying to copy file. Retry # %d", nRetries));
		if(gSystem->Exec(command.Data()) == 0) return kTRUE;
	}

	return kFALSE;

}

//______________________________________________________________________________________________
TList* AliShuttle::GetDAQFileSources(const char* detector, const char* id)
{
// Retrieves a file from the DCS FES.

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(detector, "GetDAQFileSources - Couldn't connect to DAQ Logbook");
		return 0;
	}

	// Query preparation
	TString sqlQueryStart = "select DAQsource from logbook_fs where";
	TString whereClause = Form("run=%d and detector=\"%s\" and fileId=\"%s\"",
				GetCurrentRun(), detector, id);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = fServer[kDAQ]->Query(sqlQuery);
	if (!aResult) {
		Log(detector, Form("GetDAQFileSources - Can't execute SQL query for id: %s", id));
		return 0;
	}

	if (aResult->GetRowCount() == 0) {
		Log(detector,
			Form("GetDAQFileSources - No entry in FES table for id: %s", id));
		delete aResult;
		return 0;
	}

	TSQLRow* aRow;
	TList *list = new TList();
	list->SetOwner(1);

	while((aRow = aResult->Next())){

		TString daqSource(aRow->GetField(0), aRow->GetFieldLength(0));
		AliDebug(2, Form("daqSource = %s", daqSource.Data()));
		list->Add(new TObjString(daqSource));
		delete aRow;
	}
	delete aResult;

	return list;

}

//______________________________________________________________________________________________
const char* AliShuttle::GetDCSFileName(const char* /*detector*/, const char* /*id*/, const char* /*source*/){
// Retrieves a file from the DCS FES.

return "You're in DCS";

}

//______________________________________________________________________________________________
TList* AliShuttle::GetDCSFileSources(const char* /*detector*/, const char* /*id*/){
// Retrieves a file from the DCS FES.

return NULL;

}

//______________________________________________________________________________________________
const char* AliShuttle::GetHLTFileName(const char* /*detector*/, const char* /*id*/, const char* /*source*/){
// Retrieves a file from the HLT FES.

return "You're in HLT";

}

//______________________________________________________________________________________________
TList* AliShuttle::GetHLTFileSources(const char* /*detector*/, const char* /*id*/){
// Retrieves a file from the HLT FES.

return NULL;

}

//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateDAQTable()
{
// Update DAQ table filling time_processed field in all rows corresponding to current run and detector

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(fCurrentDetector, "UpdateDAQTable - Couldn't connect to DAQ Logbook");
		return kFALSE;
	}

	TTimeStamp now; // now

	// Loop on FES list entries
	TIter iter(&fFESlist[kDAQ]);
	TObjString *aFESentry=0;
	while((aFESentry = dynamic_cast<TObjString*> (iter.Next()))){
		TString aFESentrystr = aFESentry->String();
		TObjArray *aFESarray = aFESentrystr.Tokenize("_!?!_");
		if(!aFESarray || aFESarray->GetEntries() != 2 ) {
			Log(fCurrentDetector, Form("UpdateDAQTable - error updating FES entry. Check string: <%s>",
				aFESentrystr.Data()));
			if(aFESarray) delete aFESarray;
			return kFALSE;
		}
		const char* fileId = ((TObjString*) aFESarray->At(0))->GetName();
		const char* daqSource = ((TObjString*) aFESarray->At(1))->GetName();
		TString whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\";",
			GetCurrentRun(), fCurrentDetector.Data(), fileId, daqSource);

		delete aFESarray;

		TString sqlQuery = Form("update logbook_fs set time_processed=%d %s", now.GetSec(), whereClause.Data());

		AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

		// Query execution
		TSQLResult* aResult;
		aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
		if (!aResult) {
			Log(fCurrentDetector, Form("UpdateDAQTable - Can't execute SQL query <%s>", sqlQuery.Data()));
			return kFALSE;
		}
		delete aResult;
	}

	return kTRUE;
}


//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateShuttleLogbook(const char* detector, const char* status)
{
// Update Shuttle logbook filling detector or shuttle_done column
// ex. of usage: UpdateShuttleLogbook("PHOS", "DONE") or UpdateShuttleLogbook("shuttle_done")

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log("SHUTTLE", "UpdateShuttleLogbook - Couldn't connect to DAQ Logbook.");
		return kFALSE;
	}

	TString detName(detector);
	TString setClause;
	if(detName == "shuttle_done") {
		setClause = "set shuttle_done=1";
	} else {
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

	TString sqlQuery = Form("update logbook_shuttle %s %s",
					setClause.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
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
// Get current run from logbook entry

	return fLogbookEntry ? fLogbookEntry->GetRun() : -1;
}

//______________________________________________________________________________________________
UInt_t AliShuttle::GetCurrentStartTime() const
{
// get current start time

	return fLogbookEntry ? fLogbookEntry->GetStartTime() : 0;
}

//______________________________________________________________________________________________
UInt_t AliShuttle::GetCurrentEndTime() const
{
// get current end time from logbook entry

	return fLogbookEntry ? fLogbookEntry->GetEndTime() : 0;
}

//______________________________________________________________________________________________
void AliShuttle::Log(const char* detector, const char* message)
{
// Fill log string with a message

	void* dir = gSystem->OpenDirectory(fgkShuttleLogDir);
	if (dir == NULL) {
		if (gSystem->mkdir(fgkShuttleLogDir, kTRUE)) {
			AliError(Form("Can't open directory <%s>", fgkShuttleTempDir));
			return;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	TString toLog = Form("%s (%d): %s - ", TTimeStamp(time(0)).AsString("s"), getpid(), detector);
	if(GetCurrentRun()>=0 ) toLog += Form("run %d - ", GetCurrentRun());
	toLog += Form("%s", message);

  	AliInfo(toLog.Data());

  	TString fileName;
  	fileName.Form("%s/%s.log", fgkShuttleLogDir, detector);
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

	TString whereClause("where shuttle_done=0");
	if (run != -1)
		whereClause += Form(" and run=%d", run);

	TObjArray shuttleLogbookEntries;
	if (!QueryShuttleLogbook(whereClause, shuttleLogbookEntries)) {
		Log("SHUTTLE", "Collect - Can't retrieve entries from Shuttle logbook");
		return kFALSE;
	}

	if (!RetrieveConditionsData(shuttleLogbookEntries)) {
		Log("SHUTTLE", "Collect - Process of at least one run failed");
		return kFALSE;
	}

	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::RetrieveConditionsData(const TObjArray& dateEntries)
{
// Retrieve conditions data for all runs that aren't processed yet

	Bool_t hasError = kFALSE;

	TIter iter(&dateEntries);
	AliShuttleLogbookEntry* anEntry;

	while ((anEntry = (AliShuttleLogbookEntry*) iter.Next())){
		if (!Process(anEntry)){
			hasError = kTRUE;
		}
	}

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
ULong_t AliShuttle::GetTimeOfLastAction() const
{
	ULong_t tmp;
	
	fMonitoringMutex->Lock();
	
	tmp = fLastActionTime;
	
	fMonitoringMutex->UnLock();
	
	return tmp;
}

//______________________________________________________________________________________________
const TString AliShuttle::GetLastAction() const
{
	// returns a string description of the last action

	TString tmp;
	
	fMonitoringMutex->Lock();
	
	tmp = fLastAction;
	
	fMonitoringMutex->UnLock();

	return tmp;	
}

//______________________________________________________________________________________________
void AliShuttle::SetLastAction(const char* action)
{
	// updates the monitoring variables
	
	fMonitoringMutex->Lock();
	
	fLastAction = action;
	fLastActionTime = time(0);
	
	fMonitoringMutex->UnLock();
}

//______________________________________________________________________________________________
const char* AliShuttle::GetRunParameter(const char* param)
{
// returns run parameter read from DAQ logbook

	if(!fLogbookEntry) {
		AliError("No logbook entry!");
		return 0;
	}

	return fLogbookEntry->GetRunParameter(param);
}
