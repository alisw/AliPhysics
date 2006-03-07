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
// interface methods of AliCDBPreProcessor.
// For every detector in AliShuttleConfgi (see AliShuttleConfig),
// data for its set of aliases is retrieved. If there is registered
// AliCDBPreProcessor for this detector than it will be used
// accroding to the schema (see AliCDBPreProcessor).
// If there isn't registered AliCDBPreProcessor than the retrieved
// data is stored automatically to the undelying AliCDBStorage.
// For detSpec is used the alias name.
//

#include "AliShuttle.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPreProcessor.h"
#include "AliShuttleConfig.h"
#include "AliDCSClient.h"
#include "AliLog.h"

#include <TObjString.h>

ClassImp(AliShuttle)

AliShuttle::AliShuttle(const AliShuttleConfig* config, 
		AliCDBStorage* cdbStorage, UInt_t timeout, Int_t retries):
	fConfig(config), fStorage(cdbStorage), fTimeout(timeout), 
	fRetries(retries), fCurrentRun(-1), fCurrentStartTime(0), 
	fCurrentEndTime(0)
{
	//
	// config: AliShuttleConfig used
	// cdbStorage: underlying AliCDBStorage
	// timeout: timeout used for AliDCSClient connection
	// retries: the number of retries in case of connection error.
	//

}

AliShuttle::~AliShuttle() {
	fPreProcessorMap.DeleteAll();
}

void AliShuttle::RegisterCDBPreProcessor(AliCDBPreProcessor* processor) {
	//
	// Registers new AliCDBPreProcessor.
	// It uses GetName() for indentificator of the pre processor.
	// The pre processor is registered it there isn't any other
	// with the same identificator (GetName()).
	//

	if (fPreProcessorMap.GetValue(processor->GetName())) {
		AliWarning(Form("AliCDBPreProcessor %s is already registered!",
			processor->GetName()));
		return;
	}

	fPreProcessorMap.Add(new TObjString(processor->GetName()), processor);
	processor->SetShuttle(this);
}

Bool_t AliShuttle::Store(const char* detector, const char* specType,
		TObject* object, AliCDBMetaData* metaData)
{
	if (!fStorage) {
		AliError("Invalid storage object!");
		return kFALSE;
	}

	AliCDBId id(AliCDBPath(detector, "DCS", specType), 
		GetCurrentRun(), GetCurrentRun());
	return fStorage->Put(object, id, metaData);
}

Bool_t AliShuttle::Process(Int_t run, UInt_t startTime, UInt_t endTime) {
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
	while ((aDetector = (TObjString*) iter.Next())) {
		if(!Process(run, startTime, endTime, aDetector->String())) {
			hasError = kTRUE;
		}
	}

	return !hasError;
}

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
		AliError(Form("There isn't any configuration for %s",
				detector));
		return kFALSE;
	}

	fCurrentRun = run;
	fCurrentStartTime = startTime;
	fCurrentEndTime = endTime;

	TString host(fConfig->GetHost(detector));
	Int_t port = fConfig->GetPort(detector);

	AliCDBPreProcessor* aPreProcessor = 
		(AliCDBPreProcessor*) fPreProcessorMap.GetValue(detector);

	TIter iter(fConfig->GetAliases(detector));
	TObjString* anAlias;

	Bool_t hasError = kFALSE;

	if (aPreProcessor) {
		aPreProcessor->Initialize(run, startTime, endTime);

		TObjArray valueSet;
		while ((anAlias = (TObjString*) iter.Next())) {
			Bool_t result = GetValueSet(host, port, 
					anAlias->String(), valueSet);
			
			aPreProcessor->Process(anAlias->String(), valueSet, 
					!result);

                        valueSet.Delete();
                }

		aPreProcessor->Finalize();	

	} else {
		AliCDBMetaData metaData;
		metaData.SetProperty("StartTime", 
				new AliSimpleValue(startTime));
		metaData.SetProperty("EndTime",
				new AliSimpleValue(endTime));
		metaData.SetComment("Automatically stored by AliShuttle!");

		TObjArray valueSet;
		while ((anAlias = (TObjString*) iter.Next())) {
			if (GetValueSet(host, port, anAlias->String(), 
				valueSet)) {
				if (!Store(detector, anAlias->String(),
					&valueSet, &metaData)) {
					AliError(Form("Can't store %s for %s!",
						anAlias->String().Data(),
						detector));
					hasError = kTRUE;
				}		
			}

			valueSet.Delete();
		}
	}
	
	fCurrentRun = -1;
	fCurrentStartTime = 0;
	fCurrentEndTime = 0;

	return !hasError;
}

Bool_t AliShuttle::GetValueSet(const char* host, Int_t port, const char* alias,
				TObjArray& valueSet)
{
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
