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
Revision 1.3  2005/11/17 17:47:34  byordano
TList changed to TObjArray

Revision 1.2  2005/11/17 14:43:22  byordano
import to local CVS

Revision 1.1.1.1  2005/10/28 07:33:58  hristov
Initial import as subdirectory in AliRoot

Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.2  2005/08/29 21:15:47  byordano
some docs added

*/

// Description:
// This class is the CDBPreProcessor interface,
// supposed to be implemented by any detector
// interested in immediate processing of data 
// which is retrieved from DCS.
// For every particular run set of aliases and
// their corespoding value sets are returned.
// Usage schema:
//	1) virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime)	
//	This method is called at the begining of data retrieval.
//	run: run number
//	startTime: when the run started
//	endTime: when the run finished	
//
//	2) virtual void Process(const char* alias, TObjArray& valueSet,
//                       Bool_t hasError)	
//	
//	This method is called for every particular alias which the detector
//	is registered in the configuration (see AliShuttleConfig).
//	alias: alias name
//	valueSet: collection of AliDCSValue
//	hasError: flag indicating if some error has happened during
//		the retrieaval of the data for this alias.
//
//	3) virtual void Finalize()
//	This method is called after the last alias has been processed.
//


#include "AliCDBPreProcessor.h"

#include "AliShuttle.h"
#include "AliLog.h"

ClassImp(AliCDBPreProcessor)

AliCDBPreProcessor::AliCDBPreProcessor(const char* detector):
	TNamed(detector, "")
{
	SetTitle(Form("AliCDBPreProcessor for %s subdetector.", detector));
}

AliCDBPreProcessor::~AliCDBPreProcessor() {

}

Int_t AliCDBPreProcessor::GetRun() const {
	//
	// Returns current run number.
	//

	if (!fShuttle) {
		AliError(Form("Not registered AliCDBPreProcessor: %s",
				GetName()));
		return -1;
	}

	return fShuttle->GetCurrentRun();
}

UInt_t AliCDBPreProcessor::GetStartTime() const {
	//
	// Returns currernt run start time.
	//
	if (!fShuttle) {
                AliError(Form("Not registered AliCDBPreProcessor: %s",
                                GetName()));
                return 0;
        }

        return fShuttle->GetCurrentStartTime();
}

UInt_t AliCDBPreProcessor::GetEndTime() const {
	//
	// Returns current run end time.
	//
        if (!fShuttle) {
                AliError(Form("Not registered AliCDBPreProcessor: %s",
                                GetName()));
                return 0;
        }

        return fShuttle->GetCurrentEndTime();
}

Bool_t AliCDBPreProcessor::Store(const char* specType, TObject* object,
		AliCDBMetaData* metaData)
{
	//
	// Store object and metaData into the default AliCDBStorage.
	// Storage path: <detector>/DCS/<specType>
	//

	if (!fShuttle) {
                AliError(Form("Not registered AliCDBPreProcessor: %s",
                                GetName()));
                return kFALSE;
        }

	return fShuttle->Store(GetName(), specType, object, metaData);
}

