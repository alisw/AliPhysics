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

/* $Id$*/

// Trigger preprocessor class. 
// According to the TriggerDetectorMask read from the logbook_trigger_clusters
// DAQ table, the triggering detectors are identified, and the 
// corresponding procedure is called.
// Data are stored in the OCDB, in /TRIGGER/<DET>/<level3>, where
// <DET> correspond to the triggering detector,
// and <level3> is defined in the detector procedure

#include "AliTRIPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliITSTriggerConditions.h"
          
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TList.h>
#include <TROOT.h>
#include <TSystem.h>

ClassImp(AliTRIPreprocessor)

// names of detectors/systems in the DETECTORS_MAP in /date/db/detCodes.h
const char* AliTRIPreprocessor::fgkDetectorsMapName[AliTRIPreprocessor::kNDetectorsMap] = {"SPD"/*0*/, "SDD"/*1*/, "SSD"/*2*/, "TPC"/*3*/, "TRD"/*4*/, 
									       "TOF"/*5*/, "HMP"/*6*/, "PHS"/*7*/, "CPV"/*8*/, "PMD"/*9*/, 
									       "MCH"/*10*/,"MTR"/*11*/,"FMD"/*12*/,"T00"/*13*/,"V00"/*14*/, 
									       "ZDC"/*15*/,"ACO"/*16*/,"TRI"/*17*/,"EMC"/*18*/,"TST"/*19*/, 
									       ""/*20*/,   ""/*21*/,   ""/*22*/,   ""/*23*/,   ""/*24*/,   
									       ""/*25*/,   ""/*26*/,   ""/*27*/,   ""/*28*/,   "GRP"/*29*/, 
									       "HLT"/*30*/};

//______________________________________________________________________________________________
AliTRIPreprocessor::AliTRIPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TRI", shuttle),
  fShuttle(shuttle)

{
	//
	// constructor
	//
	
	AddRunType("PHYSICS");
}

//______________________________________________________________________________________________
AliTRIPreprocessor::~AliTRIPreprocessor()
{
	//
	// destructor
	//
}

//______________________________________________________________________________________________
void AliTRIPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{

	//
	// Initialize preprocessor
	//

	AliPreprocessor::Initialize(run, startTime, endTime);

	Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

}

//______________________________________________________________________________________________
Bool_t AliTRIPreprocessor::ProcessDCS()
{
	//
	// DCS data are never needed
	//
	
	return kFALSE;
}

//______________________________________________________________________________________________
UInt_t AliTRIPreprocessor::Process(TMap* /*dcsAliasMap*/)
{

	// Procees function:
	// After reading the TriggerDetectorMask, the
	// corresponding triggering detector procedures to 
	// process the trigger data are called.

	typedef Short_t (AliTRIPreprocessor::*AliProcessTriggerData)();
	const AliProcessTriggerData processTriggerDataArray[AliTRIPreprocessor::kNDetectorsMap]= { 
		&AliTRIPreprocessor::ProcessSPDTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessTOFTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData,
		&AliTRIPreprocessor::ProcessEmptyTriggerData}; 


	// getting the list of triggering detectors from DAQ logbook

	TString triggerDetectorMask = (TString)GetTriggerDetectorMask();
	Int_t result=0;
	if (!triggerDetectorMask.IsNull()){
		Int_t length = triggerDetectorMask.Length();
		Log(Form("mask = %s", triggerDetectorMask.Data()));
		  for (Int_t i = 0; i<length; i++){
			AliDebug(2,Form("%d-th bit = %c in index %d",i,triggerDetectorMask[length-1-i],length-1-i));
			if (triggerDetectorMask[length-1-i] == '1'){
				Log("****************************************");
				Log(Form("Processing Trigger data for %s",fgkDetectorsMapName[i]));
				Log("****************************************");
			       
				result+=(this->*processTriggerDataArray[i])();
			}
		}
	}

	// result should be 0 to end successfully

	return result;

}
//______________________________________________________________________________________________
Short_t AliTRIPreprocessor::ProcessSPDTriggerData() 
{
	//
	// Processing SPD Trigger Data
	//
	
	Log("************** Processing SPD Trigger data... **************");

	// Read new conditions from dcs fxs
	AliITSTriggerConditions* newCond = new AliITSTriggerConditions();
	TString fxsID = "PITConditions";
	TList* list = GetFileSources(kDCS, fxsID.Data());
	if (!list) {
		AliError("FXS file not found.");
		return 1;
	}
	UInt_t nFiles = 0;
	while (list->At(nFiles)!=NULL) {
		TObjString* fileNameEntry = (TObjString*) list->At(nFiles);
		TString fileName = GetFile(kDCS, fxsID.Data(), fileNameEntry->GetString().Data());
		if (fileName.IsNull()) {
			Log(Form("GetFile failed to retrieve file %s.",fileNameEntry->GetString().Data()));
			return 1;
		}
		if (nFiles==0) newCond->ReadFromTextFile(fileName.Data());
		nFiles++;
	}
	if (nFiles!=1) {
		AliWarning(Form("Found %d files with id %s (expected exactly 1).",nFiles,fxsID.Data()));
	}
	
	// Read old conditions from ocdb
	AliITSTriggerConditions* oldCond = NULL;
	AliCDBEntry* pitCond = GetFromOCDB("SPD", "PITConditions");
	if (pitCond) {
		oldCond = (AliITSTriggerConditions*) pitCond->GetObject();
		if (!oldCond) {
			AliError("AliCDBEntry::GetObject() returned NULL.");
			return 1;
		}
	}
	else {
		Log("Old conditions not found in database.");
	}
	
	// Do we need to update db?
	Bool_t doUpdate = kTRUE;
	if (oldCond) {
		// compare to see if there were any changes...
		if (newCond->IsEqualTo(oldCond)) {
			Log("Old conditions equal to new conditions. Do nothing.");
			doUpdate = kFALSE;
		}
	}
	
	if (doUpdate) {
		// store new conditions in ocdb
		AliCDBMetaData metaData;
		metaData.SetResponsible("Henrik Tydesjo");
		metaData.SetComment("Created by Trigger PreProcessor");
		if (!Store("SPD", "PITConditions", newCond, &metaData, 0, kTRUE)) {
			Log("Failed to store conditions data.");
			return 1;
		}
		Log("Database updated.");
	}
	
	delete newCond;
	
	Log("************************* ...done.*************************");

	return 0; // 0 means success
	
}
//______________________________________________________________________________________________
Short_t AliTRIPreprocessor::ProcessTOFTriggerData() 
{
	//
	// Processing TOF Trigger Data
	//

	Log("************** Processing TOF Trigger data... **************");
	Log("************** Fake function **************");
	Log("************************* ...done.*************************");
	return 0;
}
//______________________________________________________________________________________________
Short_t AliTRIPreprocessor::ProcessEmptyTriggerData() 
{
	//
	// Processing TOF Trigger Data
	//

	Log("************** Trigger data Processing not yet implemented **************");
	return 0;
}


