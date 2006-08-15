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
#include "AliDCSClient.h"
#include "AliLog.h"
#include "AliPreprocessor.h"
#include "AliShuttleStatus.h"

#include <TSystem.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <fstream>

ClassImp(AliShuttle)

TString AliShuttle::fgkMainCDB("alien://DBFolder=ShuttleCDB");
TString AliShuttle::fgkLocalCDB("local://LocalShuttleCDB");
TString AliShuttle::fgkMainRefStorage("alien://DBFolder=ShuttleReference");
TString AliShuttle::fgkLocalRefStorage("local://LocalReferenceStorage");

Bool_t AliShuttle::fgkProcessDCS(kTRUE); 


const char* AliShuttle::fgkShuttleTempDir = gSystem->ExpandPathName("$ALICE_ROOT/SHUTTLE/temp");
const char* AliShuttle::fgkShuttleLogDir = gSystem->ExpandPathName("$ALICE_ROOT/SHUTTLE/log");

const char* AliShuttle::fgkDetectorName[AliShuttle::fgkNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
	"PHOS", "CPV", "RICH", "EMCAL", "MUON_TRK", "MUON_TRG", "FMD", "ZDC", "PMD", "START", "VZERO"};

const char* AliShuttle::fgkDetectorCode[AliShuttle::fgkNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
	"PHS", "CPV", "HMP", "EMC", "MCH", "MTR", "FMD", "ZDC", "PMD", "T00", "V00"};

//______________________________________________________________________________________________
AliShuttle::AliShuttle(const AliShuttleConfig* config,
		UInt_t timeout, Int_t retries):
fConfig(config),
fTimeout(timeout), fRetries(retries),
fPreprocessorMap(),
fCurrentRun(-1),
fCurrentStartTime(0), fCurrentEndTime(0),
fCurrentDetector(""),
fStatusEntry(0)
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
}

//______________________________________________________________________
AliShuttle::AliShuttle(const AliShuttle& /*other*/):
AliShuttleInterface(),
fConfig(0),
fTimeout(0), fRetries(0),
fPreprocessorMap(),
fCurrentRun(-1),
fCurrentStartTime(0), fCurrentEndTime(0),
fCurrentDetector(""),
fStatusEntry(0)

{
// copy constructor (not implemented)

}

//______________________________________________________________________
AliShuttle &AliShuttle::operator=(const AliShuttle& /*other*/)
{
// assignment operator (not implemented)

return *this;
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

	if (fPreprocessorMap.GetValue(preprocessor->GetName())) {
		AliWarning(Form("AliPreprocessor %s is already registered!",
			preprocessor->GetName()));
		return;
	}

	fPreprocessorMap.Add(new TObjString(preprocessor->GetName()), preprocessor);
}

//______________________________________________________________________________________________
UInt_t AliShuttle::Store(const AliCDBPath& path, TObject* object,
		AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // Stores a CDB object in the storage for offline reconstruction. Objects that are not needed for
  // offline reconstruction, but should be stored anyway (e.g. for debugging) should NOT be stored
  // using this function. Use StoreReferenceData instead!
  //

  // The parameters are
  //   1) the object's path.
  //   2) the object to be stored
  //   3) the metaData to be associated with the object
  //   4) the validity start run number w.r.t. the current run,
  //      if the data is valid only for this run leave the default 0
  //   5) specifies if the calibration data is valid for infinity (this means until updated),
  //      typical for calibration runs, the default is kFALSE
  //
  //
  // returns 0 if fail
  // 	     1 if stored in main (Grid) CDB storage
  // 	     2 if stored in backup (Local) CDB storage


	Int_t firstRun = GetCurrentRun() - validityStart;
  	if(firstRun < 0) {
		AliError("First valid run happens to be less than 0! Setting it to 0...");
		firstRun=0;
  	}

	Int_t lastRun = -1;
	if(validityInfinite) {
		lastRun = AliCDBRunRange::Infinity();
	} else {
		lastRun = GetCurrentRun();
	}

	AliCDBId id(path, firstRun, lastRun);

	UInt_t result = 0;

	if (!(AliCDBManager::Instance()->GetStorage(fgkMainCDB))) {
		Log(fCurrentDetector, "Cannot activate main CDB storage!");
	} else {
		result = (UInt_t) AliCDBManager::Instance()->GetStorage(fgkMainCDB)
					->Put(object, id, metaData);
	}

	if(!result) {

		Log(fCurrentDetector,
			"Error while storing object in main storage: it will go to local storage!");

		result = AliCDBManager::Instance()->GetStorage(fgkLocalCDB)
					->Put(object, id, metaData);

		if(result) {
			result = 2;
		}else{
			Log(fCurrentDetector, "Can't store data!");
		}
	}
	return result;

}

//______________________________________________________________________________________________
UInt_t AliShuttle::StoreReferenceData(const AliCDBPath& path, TObject* object,
		AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // Stores a CDB object in the storage for reference data. This objects will not be available during
  // offline reconstrunction. Use this function for reference data only!
  //

  // The parameters are
  //   1) the object's path.
  //   2) the object to be stored
  //   3) the metaData to be associated with the object
  //   4) the validity start run number w.r.t. the current run,
  //      if the data is valid only for this run leave the default 0
  //   5) specifies if the calibration data is valid for infinity (this means until updated),
  //      typical for calibration runs, the default is kFALSE
  //
  //
  // returns 0 if fail
  // 	     1 if stored in main (Grid) reference storage
  // 	     2 if stored in backup (Local) reference storage

 	Int_t firstRun = GetCurrentRun() - validityStart;
  	if(firstRun < 0) {
		AliError("First valid run happens to be less than 0! Setting it to 0...");
		firstRun=0;
  	}

	Int_t lastRun = -1;
	if(validityInfinite) {
		lastRun = AliCDBRunRange::Infinity();
	} else {
		lastRun = GetCurrentRun();
	}

 	AliCDBId id(path, firstRun, lastRun);

	UInt_t result = 0;

	if (!(AliCDBManager::Instance()->GetStorage(fgkMainRefStorage))) {
		Log(fCurrentDetector, "Cannot activate main reference storage!");
	} else {
		result = (UInt_t) AliCDBManager::Instance()->GetStorage(fgkMainRefStorage)
					->Put(object, id, metaData);
	}

	if(!result) {

		Log(fCurrentDetector,
			"Error while storing object in main reference storage: it will go to local ref storage!");

		result = AliCDBManager::Instance()->GetStorage(fgkLocalRefStorage)
					->Put(object, id, metaData);

		if(result) {
			result = 2;
		}else{
			Log(fCurrentDetector, "Can't store reference data!");
		}
	}
	return result;

}

//______________________________________________________________________________________________
AliShuttleStatus* AliShuttle::ReadShuttleStatus()
{
  // Reads the AliShuttleStatus from the CDB

  if (fStatusEntry)
  {
    delete fStatusEntry;
    fStatusEntry = 0;
  }

  fStatusEntry = AliCDBManager::Instance()->GetStorage(AliShuttle::GetLocalCDB())
      ->Get(Form("/SHUTTLE/STATUS/%s", fCurrentDetector.Data()), fCurrentRun);

  if (!fStatusEntry)
    return 0;

  TObject* anObject = fStatusEntry->GetObject();
  if (anObject == NULL || anObject->IsA() != AliShuttleStatus::Class())
  {
    AliError("Invalid object stored to CDB!");
    return 0;
  }

  AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (anObject);
  return status;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::WriteShuttleStatus(AliShuttleStatus* status)
{
  // writes the status for one subdetector

  if (fStatusEntry)
  {
    delete fStatusEntry;
    fStatusEntry = 0;
  }

  AliCDBId id(AliCDBPath("SHUTTLE", "STATUS", fCurrentDetector), fCurrentRun, fCurrentRun);

  fStatusEntry = new AliCDBEntry(status, id, new AliCDBMetaData);

  UInt_t result = AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);

  if (!result)
  {
    AliError(Form("WriteShuttleStatus for %s, run %d failed", fCurrentDetector.Data(), fCurrentRun));
    return kFALSE;
  }

  return kTRUE;
}

//______________________________________________________________________________________________
void AliShuttle::UpdateShuttleStatus(AliShuttleStatus::Status newStatus, Bool_t increaseCount)
{
  // changes the AliShuttleStatus for the given detector and run to the given status

  if (!fStatusEntry)
  {
    AliError("UNEXPECTED: fStatusEntry empty");
    return;
  }

  TObject* anObject = fStatusEntry->GetObject();
  AliShuttleStatus* status = dynamic_cast<AliShuttleStatus*> (anObject);

  if (!status)
  {
    AliError("UNEXPECTED: status could not be read from current CDB entry");
    return;
  }

  Log("SHUTTLE", Form("%s: Changing state from %s to %s", fCurrentDetector.Data(),
  				status->GetStatusName(), status->GetStatusName(newStatus)));

  status->SetStatus(newStatus);
  if (increaseCount)
    status->IncreaseCount();

  AliCDBManager::Instance()->GetStorage(fgkLocalCDB)->Put(fStatusEntry);
}

//______________________________________________________________________________________________
Bool_t AliShuttle::ContinueProcessing()
{
  // this function reads the AliShuttleStatus information from CDB and
  // checks if the processing should be continued
  // if yes it returns kTRUE and updates the AliShuttleStatus with nextStatus

  AliShuttleStatus* status = ReadShuttleStatus();
  if (!status)
  {
    // first time

    Log("SHUTTLE", Form("%s: Processing first time.", fCurrentDetector.Data()));
    status = new AliShuttleStatus(AliShuttleStatus::kStarted);
    return WriteShuttleStatus(status);
  }

  if (status->GetStatus() == AliShuttleStatus::kDone)
  {
    Log("SHUTTLE", Form("%s already done for run %d", fCurrentDetector.Data(), fCurrentRun));
    return kFALSE;
  }

  if (status->GetStatus() == AliShuttleStatus::kFailed)
  {
    Log("SHUTTLE", Form("%s already in failed state for run %d", fCurrentDetector.Data(), fCurrentRun));
    return kFALSE;
  }

  // if we get here, there is a restart

  // abort conditions
  if (status->GetStatus() == AliShuttleStatus::kPPStarted && status->GetCount() >= fConfig->GetMaxPPRetries() ||
      status->GetCount() >= fConfig->GetMaxRetries())
  {
    Log("SHUTTLE", Form("%s, run %d failed too often, %d times, status %s. Skipping processing.",
    		fCurrentDetector.Data(), fCurrentRun, status->GetCount(), status->GetStatusName()));

    return kFALSE;
  }

  Log("SHUTTLE", Form("Restart of %s, run %d. Got stuck before in %s, count %d",
  		fCurrentDetector.Data(), fCurrentRun, status->GetStatusName(), status->GetCount()));

  UpdateShuttleStatus(AliShuttleStatus::kStarted, kTRUE);

  return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Process(Int_t run, UInt_t startTime, UInt_t endTime)
{
	//
	// Makes data retrieval for all detectors in the configuration.
	// run: is the run number used
	// startTime: is the run start time
	// endTime: is the run end time
	// Returns kFALSE in case of error occured and kTRUE otherwise
	//

	AliInfo(Form("\n\n ^*^*^*^*^*^* Processing run %d ^*^*^*^*^*^*", run));

	// Initialization
	Bool_t hasError = kFALSE;
	for(Int_t iSys=0;iSys<3;iSys++) fFESCalled[iSys]=kFALSE;

	fCurrentRun = run;
	fCurrentStartTime = startTime;
	fCurrentEndTime = endTime;

	// Loop on detectors in the configuration
	TIter iter(fConfig->GetDetectors());
	TObjString* aDetector;

	while ((aDetector = (TObjString*) iter.Next())) {
		fCurrentDetector = aDetector->String();

		Bool_t detectorError=kFALSE;
		if (!fConfig->HostProcessDetector(fCurrentDetector)) continue;

		if (ContinueProcessing() == kFALSE) continue;

		if(!Process()) {
			hasError = kTRUE;
			detectorError=kTRUE;
			continue;
		}
		AliInfo(Form("Process ended successfully for detector %s!",aDetector->GetName()));

		// Process successful: Update time_processed field in FES logbooks!
		if(fFESCalled[kDAQ]) {
			hasError = (UpdateDAQTable() == kFALSE);
			fFESlist[kDAQ].Clear();
		}
		//if(fFESCalled[kDCS]) {
		//	hasError = UpdateDCSTable(aDetector->GetName());
		//	fFESlist[kDCS].Clear();
		//}
		//if(fFESCalled[kHLT]) {
		//	hasError = UpdateHLTTable(aDetector->GetName());
		//	fFESlist[kHLT].Clear();
		//}

		UpdateShuttleStatus(AliShuttleStatus::kDone);
	}

	fCurrentRun = -1;
	fCurrentStartTime = 0;
	fCurrentEndTime = 0;

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Process()
{
	//
        // Makes data retrieval just for one specific detector.
	// Threre should be a configuration for this detector.
        // run: is the run number used
        // startTime: is the run start time
        // endTime: is the run end time
        // detector: detector for which the retrieval will be made
	// Returns kFALSE in case of error occured and kTRUE otherwise
	//

	AliInfo(Form("Retrieving values for %s, run %d", fCurrentDetector.Data(), fCurrentRun));

	if (!fConfig->HasDetector(fCurrentDetector)) {
		Log(fCurrentDetector, "There isn't any configuration for %s !");
		UpdateShuttleStatus(AliShuttleStatus::kFailed);
		return kFALSE;
	}

	UpdateShuttleStatus(AliShuttleStatus::kDCSStarted);

	TString host(fConfig->GetDCSHost(fCurrentDetector));
	Int_t port = fConfig->GetDCSPort(fCurrentDetector);

	TIter iter(fConfig->GetDCSAliases(fCurrentDetector));
	TObjString* anAlias;
	TMap aliasMap;

	Bool_t hasError = kFALSE;
	Bool_t result=kFALSE;

	while ((anAlias = (TObjString*) iter.Next())) {
		TObjArray valueSet;
		// TODO Test only... I've added a flag that allows to
		// exclude DCS archive DB query
		if(fgkProcessDCS){
			AliInfo("Querying DCS archive DB data...");
			result = GetValueSet(host, port, anAlias->String(), valueSet);
		} else {
			AliInfo(Form("Skipping DCS processing. Port = %d",port));
			result = kTRUE;
		}
		if(result) {
			aliasMap.Add(anAlias->Clone(), valueSet.Clone());
		}else{
			TString message = Form("Error while retrieving alias %s !",
					anAlias->GetName());
			Log(fCurrentDetector, message.Data());
			hasError = kTRUE;
			break;
		}
	}

	if (hasError)
	{
		UpdateShuttleStatus(AliShuttleStatus::kDCSError);
		return kFALSE;
	}

  UpdateShuttleStatus(AliShuttleStatus::kPPStarted);

  AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(fCurrentDetector));
	if(aPreprocessor)
	{
		aPreprocessor->Initialize(fCurrentRun, fCurrentStartTime, fCurrentEndTime);

		// TODO Think about what to do in case of "Grid storage error"
		// (-> object put in local backup storage, return 2)
		hasError = (aPreprocessor->Process(&aliasMap) == 0);
	}else{
    // TODO default behaviour?
		AliInfo(Form("No Preprocessor for %s: storing TMap of DP arrays into CDB!", fCurrentDetector.Data()));
		AliCDBMetaData metaData;
		AliDCSValue dcsValue(fCurrentStartTime, fCurrentEndTime);
		metaData.SetResponsible(Form("Duck, Donald"));
  		metaData.SetProperty("StartEndTime", &dcsValue);
  		metaData.SetComment("Automatically stored by Shuttle!");
		AliCDBPath path(fCurrentDetector,"DCS","Data");
		hasError = (Store(path, &aliasMap, &metaData) == 0);
	}

  if (hasError)
    UpdateShuttleStatus(AliShuttleStatus::kPPError);
  else
    UpdateShuttleStatus(AliShuttleStatus::kPPDone);

  	aliasMap.Delete();

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::GetValueSet(const char* host, Int_t port, const char* alias,
				TObjArray& valueSet)
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
		AliError(Form("Can't get '%s'! Reason: %s",
			alias, AliDCSClient::GetErrorString(result)));

		if (result == AliDCSClient::fgkServerError) {
			AliError(Form("Server error: %s",
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
Bool_t AliShuttle::Connect(Int_t system){
// Connect to MySQL Server of the system's FES logbook

	// check connection: if already connected return
	if(fServer[system] && fServer[system]->IsConnected()) return kTRUE;

	TString aFESlbHost= Form("mysql://%s", fConfig->GetFESlbHost(system));

	fServer[system] = TSQLServer::Connect(aFESlbHost,
			fConfig->GetFESlbUser(system),
			fConfig->GetFESlbPass(system));
	if (!fServer[system] || !fServer[system]->IsConnected()) {
		AliError(Form("Can't establish connection to FES logbook for %s !",fkSystemNames[system]));
		return kFALSE;
	}

	// Get tables
	// TODO in the configuration should the table name be there too?
	switch(system){
		case kDAQ:
			fServer[kDAQ]->GetTables("REFSYSLOG");
			break;
		case kDCS:
			//fServer[kDCS]->GetTables("REFSYSLOG");
			break;
		case kHLT:
			//fServer[kHLT]->GetTables("REFSYSLOG");
			break;
		default:
			break;
	}

	return kTRUE;
}

//______________________________________________________________________________________________
const char* AliShuttle::GetDAQFileName(const char* detector, const char* id, const char* source){
// Retrieves a file from the DAQ FES.
// First queris the DAQ logbook_fs for the DAQ file name, using the run, detector, id and source info
// then calls RetrieveDAQFile(DAQfilename) for actual copy to local disk
// run: current run being processed (fCurrentRun)
// detector: comes from the Preprocessor name (must be converted into detector code with GetDetCode)
// id: provided as a parameter by the Preprocessor
// source: provided by the Preprocessor through GetFileSources function

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(detector, "GetDAQFileName: Couldn't connect to DAQ Logbook !");
		return 0;
	}

	// Query preparation
	TString sqlQueryStart = "select filePath from logbook_fs where";
	TString whereClause = Form("run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\"",
				fCurrentRun, GetDetCode(detector), id, source);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = fServer[kDAQ]->Query(sqlQuery);
	if (!aResult) {
		Log(detector, Form("Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}

	if (aResult->GetRowCount() == 0) {
		Log(detector,
			Form("GetDAQFileName: No result from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		return 0;
	}

	if (aResult->GetRowCount() >1) {
		Log(detector,
			Form("GetDAQFileName: More than one row resulting from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		return 0;
	}

	TSQLRow* aRow = aResult->Next();

	if(!aRow){
		Log(detector, Form("GetDAQFileName: Empty set result from query <%s>!", sqlQuery.Data()));
		delete aResult;
		return 0;
	}

	TString filePath(aRow->GetField(0), aRow->GetFieldLength(0));

	delete aResult;

	AliDebug(2, Form("filePath = %s",filePath.Data()));

	// retrieved file is renamed to make it unique
	TString localFileName = Form("%s_%d_%s_%s.shuttle",
					detector, fCurrentRun, id, source);

	// file retrieval from DAQ FES
	Bool_t result = RetrieveDAQFile(filePath.Data(), localFileName.Data());
	if(!result) {
		Log(detector, Form("copying file %s from DAQ FES failed!", filePath.Data()));
		return 0;
	} else {
		AliInfo(Form("File %s copied from DAQ FES into %s/%s !",
			filePath.Data(), fgkShuttleTempDir, localFileName.Data()));
	}


	fFESCalled[kDAQ]=kTRUE;
	TObjString *fileParams = new TObjString(Form("%s_!?!_%s", id, source));
	fFESlist[kDAQ].Add(fileParams);

	return localFileName.Data();

}

//______________________________________________________________________________________________
Bool_t AliShuttle::RetrieveDAQFile(const char* daqFileName, const char* localFileName){

	// check temp directory: trying to cd to temp; if it does not exist, create it
	AliDebug(2, Form("Copy file %s from DAQ FES into folder %s and rename it as %s",
			daqFileName,fgkShuttleTempDir, localFileName));

	void* dir = gSystem->OpenDirectory(fgkShuttleTempDir);
	if (dir == NULL) {
		if (gSystem->mkdir(fgkShuttleTempDir, kTRUE)) {
			AliError(Form("Can't open directory <%s>!", fgkShuttleTempDir));
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
TList* AliShuttle::GetDAQFileSources(const char* detector, const char* id){
// Retrieves a file from the DCS FES.

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(detector, "GetDAQFileName: Couldn't connect to DAQ Logbook !");
		return 0;
	}

	// Query preparation
	TString sqlQueryStart = "select DAQsource from logbook_fs where";
	TString whereClause = Form("run=%d and detector=\"%s\" and fileId=\"%s\"",
				fCurrentRun, GetDetCode(detector), id);
	TString sqlQuery = Form("%s %s", sqlQueryStart.Data(), whereClause.Data());

	AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

	// Query execution
	TSQLResult* aResult;
	aResult = fServer[kDAQ]->Query(sqlQuery);
	if (!aResult) {
		Log(detector, Form("GetDAQFileSources: Can't execute query <%s>!", sqlQuery.Data()));
		return 0;
	}

	if (aResult->GetRowCount() == 0) {
		Log(detector,
			Form("GetDAQFileSources: No result from SQL query <%s>!", sqlQuery.Data()));
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
	}
	delete aResult;

	return list;

}

//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateDAQTable(){
// Update DAQ table filling time_processed field in all rows corresponding to current run and detector

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(fCurrentDetector, "UpdateDAQTable: Couldn't connect to DAQ Logbook !");
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
			Log(fCurrentDetector,Form("UpdateDAQTable: error updating FES entry! string = %s",
				aFESentrystr.Data()));
			if(aFESarray) delete aFESarray;
			return kFALSE;
		}
		const char* fileId = ((TObjString*) aFESarray->At(0))->GetName();
		const char* daqSource = ((TObjString*) aFESarray->At(1))->GetName();
		TString whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\";",
			fCurrentRun,GetDetCode(fCurrentDetector), fileId, daqSource);

		delete aFESarray;

		TString sqlQuery = Form("update logbook_fs set time_processed=%d %s", now.GetSec(), whereClause.Data());

		AliDebug(2, Form("SQL query: \n%s",sqlQuery.Data()));

		// Query execution
		TSQLResult* aResult;
		aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
		if (!aResult) {
			Log(fCurrentDetector, Form("UpdateDAQTable: Can't execute query <%s>!", sqlQuery.Data()));
			return kFALSE;
		}
		delete aResult;

		// check result - TODO Is it necessary?
		sqlQuery = Form("select time_processed from logbook_fs %s", whereClause.Data());
		AliDebug(2, Form(" CHECK - SQL query: \n%s",sqlQuery.Data()));

		aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
		if (!aResult) {
			AliWarning("Can't check result!");
			continue;
		}

	if (aResult->GetRowCount() == 0) {
		Log(fCurrentDetector,
			Form("GetDAQFileName: No result from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		//return 0;
	}

	if (aResult->GetRowCount() >1) {
		Log(fCurrentDetector,
			Form("GetDAQFileName: More than one row resulting from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		//return 0;
	}

		TSQLRow *row = dynamic_cast<TSQLRow*> (aResult->Next());
		TString processedTimeString(row->GetField(0), row->GetFieldLength(0));
		Int_t processedTime = processedTimeString.Atoi();
		if(processedTime != now.GetSec()){
			Log(fCurrentDetector, Form("UpdateDAQTable: Update table error: processed_time=%d, now=%d !",
				processedTime, now.GetSec()));
			delete aResult;
			return kFALSE;
		}

		delete aResult;

	}

	return kTRUE;
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
const char* AliShuttle::GetDetCode(const char* detector){
// Return detector code

	for(int iDet=0; iDet < fgkNDetectors; iDet++){
		if(!strcmp(fgkDetectorName[iDet], detector)) return fgkDetectorCode[iDet];
	}

	return 0;
}

//______________________________________________________________________________________________
void AliShuttle::Log(const char* detector, const char* message)
{
// Fill log string with a message

	void* dir = gSystem->OpenDirectory(fgkShuttleLogDir);
	if (dir == NULL) {
		if (gSystem->mkdir(fgkShuttleLogDir, kTRUE)) {
			AliError(Form("Can't open directory <%s>!", fgkShuttleTempDir));
			return;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

  	TString toLog = Form("%s: %s, run %d - %s", TTimeStamp(time(0)).AsString(),
	detector, GetCurrentRun(), message);
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
