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
#include "AliShuttleConfig.h"
#include "AliDCSClient.h"
#include "AliLog.h"
#include "AliPreprocessor.h"
#include "AliDefaultPreprocessor.h"

#include <TSystem.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

ClassImp(AliShuttle)

TString AliShuttle::fgkLocalUri("local://$ALICE_ROOT/SHUTTLE/ShuttleCDB");
const char* AliShuttle::fgkShuttleTempDir = "$ALICE_ROOT/SHUTTLE/temp";

const char* AliShuttle::fgkDetectorName[AliShuttle::fgkNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
	"PHOS", "CPV", "RICH", "EMCAL", "MUON_TRK", "MUON_TRG", "FMD", "ZDC", "PMD", "START", "VZERO"};

const char* AliShuttle::fgkDetectorCode[AliShuttle::fgkNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
	"PHS", "CPV", "HMP", "EMC", "MCH", "MTR", "FMD", "ZDC", "PMD", "T00", "V00"};

//______________________________________________________________________________________________
AliShuttle::AliShuttle(const AliShuttleConfig* config,
		UInt_t timeout, Int_t retries):
	fConfig(config),
	fTimeout(timeout),
	fRetries(retries), fCurrentRun(-1), fCurrentStartTime(0),
	fCurrentEndTime(0),
	fLog("")
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
AliShuttleInterface()
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
UInt_t AliShuttle::Store(const char* detector,
		TObject* object, AliCDBMetaData* metaData, Int_t /*validityStart*/, Bool_t /*validityInfinite*/)
{
	// store data into CDB
  //
  // validityStart is the start validity of the data, if not 0 GetCurrentRun() - validityStart is taken
  // validityInfinite defines if the data is valid until new data arrives (e.g. for calibration runs)
  //
	// returns 0 if fail
	// 	   1 if stored in main (Grid) storage
	// 	   2 if stored in backup (Local) storage

  // TODO implement use of two parameters

  // TODO shouldn't the path be given by the preprocessor???
	AliCDBId id(AliCDBPath(detector, "DCS", "Data"),
		GetCurrentRun(), GetCurrentRun());

	UInt_t result = 0;
	if (!(AliCDBManager::Instance()->IsDefaultStorageSet())) {
		Log(detector, "No CDB storage set!");
	} else {
		result = (UInt_t) AliCDBManager::Instance()->Put(object, id, metaData);
	}
	if(!result) {

		Log(detector, "Error while storing object in main storage!");
		AliError("local storage will be used!");

		AliCDBStorage *origStorage = AliCDBManager::Instance()->GetDefaultStorage();

		result = AliCDBManager::Instance()->GetStorage(fgkLocalUri)
					->Put(object, id, metaData);

		AliCDBManager::Instance()->SetDefaultStorage(origStorage);

		if(result) {
			result = 2;
		}else{
			Log(detector, "Can't store data!");
		}
	}
	return result;

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
	ClearLog();
	Bool_t hasError = kFALSE;
	for(Int_t iSys=0;iSys<3;iSys++) fFESCalled[iSys]=kFALSE;
	fCurrentRun = run;
	fCurrentStartTime = startTime;
	fCurrentEndTime = endTime;

	// Loop on detectors in the configuration
	TIter iter(fConfig->GetDetectors());
	TObjString* aDetector;

	while ((aDetector = (TObjString*) iter.Next())) {
		Bool_t detectorError=kFALSE;
		if(!fConfig->HostProcessDetector(aDetector->GetName())) continue;
		if(!Process(run, startTime, endTime, aDetector->String())) {
			hasError = kTRUE;
			detectorError=kTRUE;
			continue;
		}
		AliInfo(Form("Process ended successfully for detector %s!",aDetector->GetName()));

		// Process successful: Update time_processed field in FES logbooks!
		if(fFESCalled[kDAQ]) {
			hasError = (UpdateDAQTable(aDetector->GetName()) == kFALSE);
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
	}

	if(fLog != "") StoreLog(run);
	fCurrentRun = -1;
	fCurrentStartTime = 0;
	fCurrentEndTime = 0;

	return hasError == kFALSE;
}

//______________________________________________________________________________________________
Bool_t AliShuttle::Process(Int_t run, UInt_t startTime, UInt_t endTime,
		const char* detector)
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

	AliInfo(Form("Retrieving values for %s, run %d", detector, run));

	if (!fConfig->HasDetector(detector)) {
		Log(detector, "There isn't any configuration for %s !");
		return kFALSE;
	}

	TString host(fConfig->GetDCSHost(detector));
	Int_t port = fConfig->GetDCSPort(detector);

	TIter iter(fConfig->GetDCSAliases(detector));
	TObjString* anAlias;
	TMap aliasMap;

	Bool_t hasError = kFALSE;
	Bool_t result=kFALSE;

	while ((anAlias = (TObjString*) iter.Next())) {
		TObjArray valueSet;
		result = GetValueSet(host, port, anAlias->String(), valueSet);
		//AliInfo(Form("Port = %d",port));
		//result = kTRUE;
		if(result) {
			aliasMap.Add(anAlias->Clone(), valueSet.Clone());
		}else{
			TString message = Form("Error while retrieving alias %s !",
					anAlias->GetName());
			Log(detector, message.Data());
			hasError = kTRUE;
		}
	}

	// even if hasError is TRUE the Shuttle should keep on processing the detector (calib files!)

	if(hasError) return kFALSE;
	// TODO if(hasError) mark DCS error

	AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(detector));
	if(aPreprocessor)
	{
		aPreprocessor->Initialize(run, startTime, endTime);
		hasError = (aPreprocessor->Process(&aliasMap) == 0);
	}else{
    // TODO default behaviour?
		AliInfo(Form("No Preprocessor for %s: storing TMap of DP arrays into CDB!",detector));
		AliCDBMetaData metaData;
		AliDCSValue dcsValue(startTime, endTime);
		metaData.SetResponsible(Form("Duck, Donald"));
  		metaData.SetProperty("StartEndTime", &dcsValue);
  		metaData.SetComment("Automatically stored by Shuttle!");
		hasError = (Store(detector, &aliasMap, &metaData) == 0);
	}


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

	AliInfo(Form("SQL query: \n%s",sqlQuery.Data()));

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

	AliInfo(Form("filePath = %s",filePath.Data()));

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
	AliInfo(Form("Copy file %s from DAQ FES into folder %s and rename it as %s",
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

	AliInfo(Form("%s",command.Data()));

	UInt_t nRetries = 0;
	UInt_t maxRetries = 3;

	// copy!! if successful TSystem::Exec returns 0
	while(nRetries++ < maxRetries) {
		AliInfo(Form("Trying to copy file. Retry # %d", nRetries));
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

	AliInfo(Form("SQL query: \n%s",sqlQuery.Data()));

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
		AliInfo(Form("daqSource = %s", daqSource.Data()));
		list->Add(new TObjString(daqSource));
	}
	delete aResult;

	return list;

}

//______________________________________________________________________________________________
Bool_t AliShuttle::UpdateDAQTable(const char* detector){
// Update DAQ table filling time_processed field in all rows corresponding to current run and detector

	// check connection, in case connect
	if(!Connect(kDAQ)){
		Log(detector, "UpdateDAQTable: Couldn't connect to DAQ Logbook !");
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
			Log(detector,Form("UpdateDAQTable: error updating FES entry! string = %s",
				aFESentrystr.Data()));
			if(aFESarray) delete aFESarray;
			return kFALSE;
		}
		const char* fileId = ((TObjString*) aFESarray->At(0))->GetName();
		const char* daqSource = ((TObjString*) aFESarray->At(1))->GetName();
		TString whereClause = Form("where run=%d and detector=\"%s\" and fileId=\"%s\" and DAQsource=\"%s\";",
			fCurrentRun,GetDetCode(detector), fileId, daqSource);

		delete aFESarray;

		TString sqlQuery = Form("update logbook_fs set time_processed=%d %s", now.GetSec(), whereClause.Data());

		AliInfo(Form("SQL query: \n%s",sqlQuery.Data()));

		// Query execution
		TSQLResult* aResult;
		aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
		if (!aResult) {
			Log(detector, Form("UpdateDAQTable: Can't execute query <%s>!", sqlQuery.Data()));
			return kFALSE;
		}
		delete aResult;

		// check result - TODO Is it necessary?
		sqlQuery = Form("select time_processed from logbook_fs %s", whereClause.Data());
		AliInfo(Form(" CHECK - SQL query: \n%s",sqlQuery.Data()));

		aResult = dynamic_cast<TSQLResult*> (fServer[kDAQ]->Query(sqlQuery));
		if (!aResult) {
			AliWarning("Can't check result!");
			continue;
		}

	if (aResult->GetRowCount() == 0) {
		Log(detector,
			Form("GetDAQFileName: No result from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		//return 0;
	}

	if (aResult->GetRowCount() >1) {
		Log(detector,
			Form("GetDAQFileName: More than one row resulting from SQL query <%s>!", sqlQuery.Data()));
		delete aResult;
		//return 0;
	}

		TSQLRow *row = dynamic_cast<TSQLRow*> (aResult->Next());
		TString processedTimeString(row->GetField(0), row->GetFieldLength(0));
		Int_t processedTime = processedTimeString.Atoi();
		if(processedTime != now.GetSec()){
			Log(detector, Form("UpdateDAQTable: Update table error: processed_time=%d, now=%d !",
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

	TString toLog = Form("%s - %s", detector, message);
	AliError(toLog.Data());

	fLog += toLog;
	fLog += "\n";

}

//______________________________________________________________________________________________
void AliShuttle::StoreLog(Int_t run)
{
// store error log string to SHUTTLE/SYSTEM/ERROR (on local storage)

	AliInfo("Printing fLog...");
	AliInfo(fLog.Data());
	// Storing log string for runs with errors in "SHUTTLE/SYSTEM/ERRORLOGS"
	TObjString *logString = new TObjString(fLog);
	AliCDBId badRunId("SHUTTLE/SYSTEM/ERRORLOGS",run,run);
	AliCDBMetaData metaData;
	AliCDBManager::Instance()->GetStorage(fgkLocalUri)
					->Put(logString, badRunId,&metaData);
	delete logString;


}

