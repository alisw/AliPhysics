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
// Prototype of ZDC Preprocessor
//

#include "TestZDCPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliShuttleInterface.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSystem.h>

ClassImp(TestZDCPreprocessor)

//________________________________________________________________________________________
TestZDCPreprocessor::TestZDCPreprocessor():
	AliPreprocessor("ZDC",0)
{
// default constructor - Don't use this!

}

//________________________________________________________________________________________
TestZDCPreprocessor::TestZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle):
	AliPreprocessor(detector,shuttle)
{
// constructor - shuttle must be instantiated!

}

//________________________________________________________________________________________
void TestZDCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
// Initialize preprocessor

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;
	AliInfo("This preprocessor is to test the GetRunParameter function.");
}

//________________________________________________________________________________________
UInt_t TestZDCPreprocessor::Process(TMap* /*valueMap*/)
{
// process data retrieved by the Shuttle

	Int_t result=0;

	const char* dataRate = GetRunParameter("averageDataRate");
	if (dataRate) {
		Log(Form("Average data rate for run %d: %s",fRun, dataRate));
	} else {
		Log(Form("Average data rate not put in logbook!"));
	}

	return dataRate !=0;
}

