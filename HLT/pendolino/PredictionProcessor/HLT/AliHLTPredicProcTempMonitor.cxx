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

/** @file   AliHLTPredicProcTempMonitor.cxx
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPredicProcTempMonitor.h"

#include <AliCDBMetaData.h>
#include <AliCDBEntry.h>
#include <AliDCSValue.h>

#include <TMap.h>
#include <TObjString.h>
#include <TObjArray.h>

#include <TTimeStamp.h>


ClassImp(AliHLTPredicProcTempMonitor)


const TString AliHLTPredicProcTempMonitor::kPath2("CalibMonitor");

const TString AliHLTPredicProcTempMonitor::kPath3("DCSTempMon");

const TString AliHLTPredicProcTempMonitor::kCreator("S. Bablok (HLT)");

const TString AliHLTPredicProcTempMonitor::kComment("Calib Object for monitoring DCS temperature values in HLT.");

const TString AliHLTPredicProcTempMonitor::kAliRootVersion("");

const UInt_t AliHLTPredicProcTempMonitor::kUnableToStoreObject = 7;

//const TString AliHLTPredicProcTempMonitor::kAmandaTempSensor = "TPC_PT_%d_TEMPERATURE";

		
AliHLTPredicProcTempMonitor::AliHLTPredicProcTempMonitor(
			const char* detector, AliHLTPendolino* pendolino) :
				AliHLTPredictionProcessorInterface(detector, pendolino),
				fPredict(true), fRun(0), fStartTime(0), fEndTime(0), fBField("") {
	// C-tor for AliHLTPredicProcTempMonitor
//	fPredict = true;
//	fRun = 0;
//	fStartTime = 0;
//	fEndTime = 0;
}


AliHLTPredicProcTempMonitor::~AliHLTPredicProcTempMonitor() {
	// D-tor for AliHLTPredicProcTempMonitor

}


UInt_t AliHLTPredicProcTempMonitor::makePrediction(Bool_t doPrediction) {
	// switch for prediction making in AliHLTPredicProcTempMonitor
	Log("AliHLTPredicProcTempMonitor + B-Field extractor: prediction switched on");
	fPredict = doPrediction;
	return 0;
}


void AliHLTPredicProcTempMonitor::Initialize(Int_t run, UInt_t startTime, 
			UInt_t endTime) {
	// initializes AliHLTPredicProcTempMonitor
	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;

	TString msg("Initialized HLT PredictionProcessor; Run: ");
	msg += fRun;
	msg += ", start time: ";
	msg += fStartTime;
	msg += ", end time: ";
	msg += fEndTime;
	msg += ".";	
	Log(msg.Data());
}


UInt_t AliHLTPredicProcTempMonitor::Process(TMap* dcsAliasMap) {
	// processes the DCS value map in AliHLTPredicProcTempMonitor
	UInt_t beamPeriod = 0;
	
	UInt_t retVal = 0;
	UInt_t tempRet = 0;
	Int_t start = 0;
	TMap* tempMap = 0;
	Bool_t infiniteValid = kFALSE;
	
		
	//test GetFromOCDB() calls
	AliCDBEntry* entry = GetFromOCDB(kPath2.Data(), kPath3.Data());
	if (entry == 0) {
		// No object in HCDB -> discarding old values
		TString msg("No '" + kPath2 + "/" + kPath3 + 
				"' in HCDB, most likely first round, filling object now...");
		Log(msg.Data());
		tempMap = dcsAliasMap;
		
	} else {
		Log("Old data is already stored in HCDB, but now discarding that...");
		tempMap = dcsAliasMap; 
	
// If old data shall be included use lines below and uncomment lines above
/*
		Log("Adding new DCS value to old Temperature map...");
		// Adding new DCS values to old Temperature map 
		tempMap = (TMap*) entry->GetObject();

		TMapIter iter(dcsAliasMap);        // iterator for values in TMap
		TObject* aKey;
			   
		while ((aKey = iter.Next())) {
			tempMap->Add(aKey, dcsAliasMap->GetValue(aKey));
		}
*/
	}

	AliCDBMetaData meta(kCreator.Data(), beamPeriod, kAliRootVersion.Data(), 
			kComment.Data());

	if (Store(kPath2.Data(), kPath3.Data(), (TObject*) tempMap, &meta, start, 
			infiniteValid)) {
		TString msg(" +++ Successfully stored object '" + kPath2 + "/" + kPath3 + 
				"' in HCDB.");
		Log(msg.Data());
	} else {
		TString msg(" *** Storing of object '" + kPath2 + "/" + kPath3 + 
				"' in HCDB failed.");
		Log(msg.Data());
		retVal = kUnableToStoreObject;
	}

	// extract B-Field
	tempRet = ExtractBField(dcsAliasMap);
	retVal = tempRet & retVal; // combine retvals
	
	return retVal;
}

UInt_t AliHLTPredicProcTempMonitor::ExtractBField(TMap* dcsAliasMap) {
	// extracts the B-field value from DCS value map
  
	TString stringId = "L3Current"; // "dcs_magnet:Magnet/ALICESolenoid.Current";
  
	Float_t BField = 0; 
	Bool_t bRet = GetSensorValue(dcsAliasMap,stringId.Data(),&BField);

	if (bRet) {
		BField = BField / 60000; // If we get field, take this away and change SensorValue
		TString dummy("-solenoidBZ ");
		dummy += BField;
		TObjString dummy2(dummy.Data());
		fBField = dummy2; 
		Log(Form("BField set to %s", fBField.String().Data())); 
	} else {
		return 1;
	}

	TString path2("Config");
	TString path3("BField");
	Int_t start = 0;

	TString comment("BField");
	AliCDBMetaData meta(this->GetName(), 0, "unknownAliRoot", comment.Data());
  
	if (Store(path2.Data(), path3.Data(), (TObject*) (&fBField), &meta, start, 
			kTRUE)) {
		Log(" +++ Successfully stored object ;-)");
	} else {
		Log(" *** Storing of OBJECT failed!!");
		return 7;
	}

	return 0;
}

Bool_t AliHLTPredicProcTempMonitor::GetSensorValue(TMap* dcsAliasMap,
			const char* stringId, Float_t *value) {
	// retreives the sensor value
  // return last value read from sensor specified by stringId
  
	TObjArray* valueSet;
	TPair* pair = (TPair*) (dcsAliasMap->FindObject(stringId));
	if (pair) {
		valueSet = (TObjArray*) (pair->Value());
		if (valueSet) {
			Int_t nentriesDCS = (valueSet->GetEntriesFast()) - 1;
			if (nentriesDCS >= 0) {
				AliDCSValue* val = (AliDCSValue*) (valueSet->At(nentriesDCS));
				if (val) {
					*value = val->GetFloat();
					return kTRUE;
				}
			}
		}
	}
	return kFALSE;
}

TMap* AliHLTPredicProcTempMonitor::produceTestData(TString /*aliasName*/) {
	// produces test data for AliHLTPredicProcTempMonitor
    TMap* resultMap = 0;

    // here has to come real dummy data :-)
    resultMap = new TMap();
    TTimeStamp tt;
	Float_t fval = 33.3;
    TObjString* name = new TObjString("L3Current");
    AliDCSValue* val = new AliDCSValue(fval, tt.GetTime());
    TObjArray* arr = new TObjArray();
    arr->Add(val);
    resultMap->Add(name, arr);

    return resultMap;
}


