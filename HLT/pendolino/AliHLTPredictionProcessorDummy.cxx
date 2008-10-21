// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sebastian Bablok <Sebastian.Bablok@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTPredictionProcessorDummy.cxx
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorDummy.h"

#include <AliCDBMetaData.h>
#include <AliCDBEntry.h>

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <AliDCSValue.h>

ClassImp(AliHLTPredictionProcessorDummy)

AliHLTPredictionProcessorDummy::AliHLTPredictionProcessorDummy(
			const char* detector, AliHLTPendolino* pendolino) :
				AliHLTPredictionProcessorInterface(detector, pendolino),
				 fPredict(false), fRun(0), fStartTime(0), fEndTime(0) {
	// C-tor of AliHLTPredictionProcessorDummy
}


AliHLTPredictionProcessorDummy::~AliHLTPredictionProcessorDummy() {
	// D-tor of AliHLTPredictionProcessorDummy

}


UInt_t AliHLTPredictionProcessorDummy::makePrediction(Bool_t doPrediction) {
	// switches prediction making on or off
	Log("Prediction switched on");
	fPredict = doPrediction;
	return 0;
}


void AliHLTPredictionProcessorDummy::Initialize(Int_t run, UInt_t startTime, 
			UInt_t endTime) {
	// initializes AliHLTPredictionProcessorDummy
	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;

	TString msg("Initialized Dummy PredictProc. Run: ");
	msg += fRun;
	msg += ", start time: ";
	msg += fStartTime;
	msg += ", end time: ";
	msg += fEndTime;
	msg += ".";	
	Log(msg.Data());

	if (fPredict) {
		Log("Dummy PredictProc has prediction switched ON.");
	} else {
		Log("Prediction is switched OFF.");
	}
}


UInt_t AliHLTPredictionProcessorDummy::Process(TMap* dcsAliasMap) {
	// processes the handed in DCS value map
	UInt_t retVal = 0;
	Int_t start = 0;
	TString path2("Dummy");
	TString path3("DCSValley");
		
	//test GetFromOCDB() calls
	AliCDBEntry* entry = GetFromOCDB("Calib", "LocalVdrift");
	if (entry != 0) {
		TString msg("AliCDBEntry -> last storage: ");
	   	msg += entry->GetLastStorage();
		Log(msg.Data());
		entry->PrintMetaData();
	} else {
		Log("Error. Cannot retrieve HCDB entry.");
	}

	entry = GetFromOCDB("Dummy", "DCSValley");
	if (entry != 0) {
		TString msg("AliCDBEntry -> last storage: ");
		msg += entry->GetLastStorage();
		Log(msg.Data());
		entry->PrintMetaData();
	} else {
		Log("Error. Cannot retrieve HCDB entry.");
	}

	// test call for non existing object
    entry = GetFromOCDB("Nix", "Nada");
    if (entry != 0) {
     	Log("Error: call should failed, but it does NOT!");
    } else {
        Log("Tested call for missing CDB entry -> succeeded!");
		if (!includeAliCDBEntryInList("TPC/Nix/Nada")) {
            Log("Adding of AliCDBEntry request for 'TPC/Nix/Nada' failed!");
        } else {
            Log("Successfully added an AliCDBEntry request (TPC/Nix/Nada).");
        }

    }
	
	
	//TODO test GetRunParameter() calls
		
	if (GetHLTStatus()) {
		Log("HLT status is: ON");
	} else {
		Log("ERROR. HLT is OFF.");
		retVal = 9;
	}

	TString runType("RunType is set to: ");
	runType += GetRunType();
	Log(runType.Data());

	//transform dcsAliasMap to ROOT object 
	TString comment("Beam it up");
	if (fPredict) {
		comment += " with PREDICTION!";
	}
	if (fRun > 3) {
		start = 3;
	}
	AliCDBMetaData meta(this->GetName(), 666, "unknownAliRoot", 
			comment.Data());

	if (Store(path2.Data(), path3.Data(), (TObject*) dcsAliasMap, &meta, start, 
			kTRUE)) {
		Log(" +++ Successfully stored object ;-)");
	} else {
		Log(" *** Storing of OBJECT failed!!");
		retVal = 7;
	}
	
	return retVal;
}

TMap* AliHLTPredictionProcessorDummy::produceTestData(TString /*aliasName*/) {
	// generates test dummy data for AliHLTPredictionProcessorDummy
	TMap* resultMap = 0;

	// here has to come real dummy data :-)
	resultMap = new TMap();
	TTimeStamp tt;
	Float_t fval = 33.3;
	TObjString* name = new TObjString("DummyData");
	AliDCSValue* val = new AliDCSValue(fval, tt.GetTime());
	TObjArray* arr = new TObjArray();
	arr->Add(val);
	resultMap->Add(name, arr);
	
	return resultMap;
}

