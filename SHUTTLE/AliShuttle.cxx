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

#include <TObject.h>
#include <TString.h>
#include <TObjString.h>

ClassImp(AliShuttle)

TString AliShuttle::fgkLocalUri("local://ShuttleCDB");

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
	// mainStorage: underlying AliCDBStorage
	// localStorage (local) CDB storage to be used if mainStorage is unavailable
	// timeout: timeout used for AliDCSClient connection
	// retries: the number of retries in case of connection error.
	//

	RegisterPreprocessor(new AliDefaultPreprocessor("DEFAULT", 0));

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
		TObject* object, AliCDBMetaData* metaData)
{
	// store data into CDB
	// returns 0 if fail
	// 	   1 if stored in main (Grid) storage
	// 	   2 if stored in backup (Local) storage

	if (!(AliCDBManager::Instance()->IsDefaultStorageSet())) {
		AliError("No CDB storage set!");
		return 0;
	}

	AliCDBId id(AliCDBPath(detector, "DCS", "Data"),
		GetCurrentRun(), GetCurrentRun());

	UInt_t result = (UInt_t) AliCDBManager::Instance()->Put(object, id, metaData);
	if(!result) {

		Log(detector, "Error while storing object in main storage!");
		AliError("local storage will be used!");

//		result = fLocalStorage->Put(object, id, metaData);
		result = AliCDBManager::Instance()->GetStorage(fgkLocalUri)
					->Put(object, id, metaData);

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

	Bool_t hasError = kFALSE;

	TIter iter(fConfig->GetDetectors());
	TObjString* aDetector;

	ClearLog();

	while ((aDetector = (TObjString*) iter.Next())) {
		if(!fConfig->HostProcessDetector(aDetector->GetName())) continue;
		if(!Process(run, startTime, endTime, aDetector->String())) {
			hasError = kTRUE;
		}
	}

	if(fLog != "") StoreLog(run);

	return !hasError;
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

	fCurrentRun = run;
	fCurrentStartTime = startTime;
	fCurrentEndTime = endTime;

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
		if(result) {
			aliasMap.Add(anAlias->Clone(), valueSet.Clone());
		}else{
			TString message = Form("Error while retrieving alias %s !", 
					anAlias->GetName());
			Log(detector, message.Data());
			hasError = kTRUE;
		}
	}

	AliPreprocessor* aPreprocessor =
		dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue(detector));
	if(!aPreprocessor){
		AliInfo(Form("No Preprocessor for %s: Using default Preprocessor!",detector));
		aPreprocessor = dynamic_cast<AliPreprocessor*> (fPreprocessorMap.GetValue("DEFAULT"));
	}

	aPreprocessor->Initialize(run, startTime, endTime);
	hasError = (Bool_t) !(aPreprocessor->Process(&aliasMap));

  aliasMap.Delete();

	fCurrentRun = -1;
	fCurrentStartTime = 0;
	fCurrentEndTime = 0;

	return !hasError;
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
const char* AliShuttle::GetFile(Int_t /*system*/, const char* /*detector*/,
		const char* /*id*/, const char* /*source*/)
{
// Get calibration file from DAQ transient file system

	AliInfo("You are in AliShuttle::GetFile!");
	return 0;
}


//______________________________________________________________________________________________
TList* AliShuttle::GetFileSources(Int_t /*system*/, const char* /*detector*/, const char* /*id*/)
{
// Get list of sources that provided the files to be retrieved from DAQ

	AliInfo("You are in AliShuttle::GetFileSources!");
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
