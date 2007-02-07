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
// Prototype of HMPID Preprocessor
//

#include "TestHMPIDPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliShuttleInterface.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSystem.h>

ClassImp(TestHMPIDPreprocessor)

//________________________________________________________________________________________
TestHMPIDPreprocessor::TestHMPIDPreprocessor():
	AliPreprocessor("HMP",0)
{
// default constructor - Don't use this!

}

//________________________________________________________________________________________
TestHMPIDPreprocessor::TestHMPIDPreprocessor(AliShuttleInterface* shuttle):
	AliPreprocessor("HMP", shuttle)
{
// constructor - shuttle must be instantiated!

}

//________________________________________________________________________________________
void TestHMPIDPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime) 
{
// Initialize preprocessor

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;
}

//________________________________________________________________________________________
UInt_t TestHMPIDPreprocessor::Process(TMap* /*valueMap*/)
{
// process data retrieved by the Shuttle

	//TIter iter(valueMap);
	//TPair* aPair;
	//while ((aPair = (TPair*) iter.Next())) {
		//aPair->Print();
	//}
	//AliCDBMetaData metaData;
	//metaData.SetComment("This is a test!");

	// return Store(valueMap, &metaData);

	// DAQ
	TList* filesources = GetFileSources(AliShuttleInterface::kDAQ, "DAQFile");

	if(!filesources) {
		AliError(Form("No sources found for thresholds.txt for run %d !", fRun));
		return 0;
	}

	AliInfo("Here's the list of sources for thresholds.txt");
	filesources->Print();

	TIter iter(filesources);
	TObjString* source;
	int i=0;
	UInt_t result = 0;
	while((source=dynamic_cast<TObjString*> (iter.Next()))){
		printf("\n\n Getting file #%d\n",++i);
		//if(i==1) continue;
		TString filename = GetFile(AliShuttleInterface::kDAQ, "DAQFile", source->GetName());
		if(!filename.Length()) {
			AliError(Form("Error: retrieval of file from source %s failed!", source->GetName()));
			delete filesources;
			return 0;
		}
		TString command = Form("more %s",filename.Data());
		gSystem->Exec(command.Data());

		// STORAGE! The First file name will be stored into CDB, the second into reference storage
		TObjString filenameObj(filename);
		AliCDBMetaData metaData;
		if(i==1) result = Store("Calib", "DAQData", &filenameObj, &metaData);
		if(i==2) result = StoreReferenceData("Calib", "RefData", &filenameObj, &metaData);

	}
	delete filesources;
	
	// DCS
	TString filename = GetFile(AliShuttleInterface::kDCS, "DCSFile", 0);
	if(!filename.Length()) {
		AliError(Form("Error: retrieval of file from DCS failed!"));
		return 0;
	}
	TString command = Form("more %s", filename.Data());
	gSystem->Exec(command.Data());

	// STORAGE! The First file name will be stored into CDB, the second into reference storage
	TObjString filenameObj(filename);
	AliCDBMetaData metaData;
	result = Store("Calib", "DCSData", &filenameObj, &metaData);

	// HLT
	filesources = GetFileSources(AliShuttleInterface::kHLT, "HLTFile");

	if(!filesources) {
		Log(Form("No sources found for HLTFile for run %d !", fRun));
		return 0;
	}

	AliInfo("Here's the list of sources for HLTFile");
	filesources->Print();

	TIter iter2(filesources);
	i = 0;
	while((source=dynamic_cast<TObjString*> (iter2.Next()))){
		printf("\n\n Getting file #%d\n",++i);
		//if(i==1) continue;
		TString filename = GetFile(AliShuttleInterface::kHLT, "HLTFile", source->GetName());
		if(!filename.Length()) {
			AliError(Form("Error: retrieval of file from source %s failed!", source->GetName()));
			delete filesources;
			return 0;
		}
		TString command = Form("more %s",filename.Data());
		gSystem->Exec(command.Data());

		// STORAGE! The First file name will be stored into CDB, the second into reference storage
		TObjString filenameObj(filename);
		AliCDBMetaData metaData;
		if(i==1) result = Store("Calib", "HLTData", &filenameObj, &metaData);
		if(i==2) result = StoreReferenceData("Calib", "RefHLTData", &filenameObj, &metaData);

	}
	delete filesources;
	

	return result;
}

