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

#include "TestITSPreprocessorSPD.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>

ClassImp(TestITSPreprocessorSPD)

//________________________________________________________________________________________
TestITSPreprocessorSPD::TestITSPreprocessorSPD():
	AliPreprocessor("SPD",0)
{
// default constructor - Don't use this!

}

//________________________________________________________________________________________
TestITSPreprocessorSPD::TestITSPreprocessorSPD(AliShuttleInterface* shuttle):
	AliPreprocessor("SPD", shuttle)
{
// constructor - shuttle must be instantiated!

}

//________________________________________________________________________________________
void TestITSPreprocessorSPD::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
// Initialize preprocessor

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, 
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));
}

//________________________________________________________________________________________
UInt_t TestITSPreprocessorSPD::Process(TMap* valueMap)
{
// process data retrieved by the Shuttle

	AliInfo(Form("You're in AliITSPreprocessor::Process!"));

	TIter iter(valueMap);
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		aPair->Print();
	}
	AliCDBMetaData metaData;
	metaData.SetComment("This is a test!");

	return Store("Calib", "ITSDataSPD", valueMap, &metaData, 0, kTRUE);
}

