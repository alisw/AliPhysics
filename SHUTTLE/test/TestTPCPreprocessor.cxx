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

//
// Example of a Shuttle Preprocessor
//

#include "TestTPCPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TH2F.h>

ClassImp(TestTPCPreprocessor)

//________________________________________________________________________________________
TestTPCPreprocessor::TestTPCPreprocessor():
	AliPreprocessor("TPC",0)
{
// default constructor - Don't use this!

 	fData = 0;
}

//________________________________________________________________________________________
TestTPCPreprocessor::TestTPCPreprocessor(const char* detector, AliShuttleInterface* shuttle):
	AliPreprocessor(detector,shuttle)
{
// constructor - shuttle must be instantiated!

	fData = 0;

}

//________________________________________________________________________________________
TestTPCPreprocessor::~TestTPCPreprocessor()
{
// destructor

	delete fData;
	fData = 0;

}

//________________________________________________________________________________________
void TestTPCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
// Initialize preprocessor

	fRun=run;
	fStartTime = startTime;
	fEndTime = endTime;
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliTPCDataDCS(fRun, fStartTime, fEndTime);
}

//________________________________________________________________________________________
UInt_t TestTPCPreprocessor::Process(TMap* aliasMap)
{
// Process data

	fData->ProcessData(*aliasMap);
	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Alberto Colla");
	metaData.SetComment("This preprocessor fills an AliTPCDataDCS object.");

	return Store(fData, &metaData);
	delete fData;
	fData = 0;
}

