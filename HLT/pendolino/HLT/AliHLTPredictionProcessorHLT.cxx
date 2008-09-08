// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk                                        *
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

/** @file   AliHLTPredictionProcessorHLT.cxx
    @author Gaute Ovrebekk
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorHLT.h"

#include <AliCDBMetaData.h>
#include <AliCDBEntry.h>

// new
#include <TObjArray.h>
#include <AliDCSValue.h>

#include <TTimeStamp.h>
#include <TObjString.h>


ClassImp(AliHLTPredictionProcessorHLT)

AliHLTPredictionProcessorHLT::AliHLTPredictionProcessorHLT(
			const char* detector, AliHLTPendolino* pendolino) :
				AliHLTPredictionProcessorInterface(detector, pendolino),
				fPredict(true), fRun(0), fStartTime(0), fEndTime(0), fBField("") {
	// C-tor for AliHLTPredictionProcessorHLT
//	fPredict = false;
//	fRun = 0;
//	fStartTime = 0;
//	fEndTime = 0;
//	fBField = 0;
}


AliHLTPredictionProcessorHLT::~AliHLTPredictionProcessorHLT() {
	// D-tor for AliHLTPredictionProcessorHLT
}


UInt_t AliHLTPredictionProcessorHLT::makePrediction(Bool_t doPrediction) {
	// switch for prediction making
	Log("Prediction switched on");
	fPredict = doPrediction;
	return 0;
}


void AliHLTPredictionProcessorHLT::Initialize(Int_t run, UInt_t startTime, 
			UInt_t endTime) {
	// initializes AliHLTPredictionProcessorHLT
	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;

	TString msg("Initialized HLT PredictProc. Run: ");
	msg += fRun;
	msg += ", start time: ";
	msg += fStartTime;
	msg += ", end time: ";
	msg += fEndTime;
	msg += ".";	
	Log(msg.Data());

	if (fPredict) {
		Log("HLT PredictProc has prediction switched ON.");
	} else {
		Log("Prediction is switched OFF.");
	}
}


UInt_t AliHLTPredictionProcessorHLT::Process(TMap* dcsAliasMap) {
	// processes the DCS value map
  
  if (!dcsAliasMap) return 9;
  if (dcsAliasMap->GetEntries() == 0 ) return 9;

  UInt_t retVal = 0;
  Int_t start = 0;
  TString path2("ConfigHLT"); // "Config" 
  TString path3("SolenoidBz"); // "BField"
  
  UInt_t BFieldResult = ExtractBField(dcsAliasMap);

  if (BFieldResult != 0) {
	Log(" *** Extraction of BField failed - no entry for HCDB!!");
	return 8;
  }


  //transform dcsAliasMap to ROOT object 
  TString comment("BField");
  AliCDBMetaData meta(this->GetName(), 0, "unknownAliRoot",comment.Data());
  
  if (Store(path2.Data(),path3.Data(),(TObject*) &fBField,&meta,start,kTRUE)) {
    Log(" +++ Successfully stored object ;-)");
  } else {
    Log(" *** Storing of OBJECT failed!!");
    retVal = 7;
  }
  
  return retVal;
}

UInt_t AliHLTPredictionProcessorHLT::ExtractBField(TMap* dcsAliasMap){
  // extracts the b-field from DCS value map
//  TString stringId = "dcs_magnet:Magnet/ALICESolenoid.Current"; // old name
	TString stringId = "L3Current";
  
  Float_t BField = 0; 
  Bool_t bRet = GetSensorValue(dcsAliasMap,stringId.Data(),&BField);

  if(bRet){
	// new	 
    BField = BField/60000; // If we get field, take this away and change SensorValue
    TString dummy("-solenoidBz ");
    dummy += BField;
    TObjString dummy2(dummy.Data());
    fBField = dummy2;
	
    Log(Form("BField set to %s",fBField.String().Data()));
    return 0; 
  }
  
  return 1;

}

Bool_t AliHLTPredictionProcessorHLT::GetSensorValue(TMap* dcsAliasMap,
							  const char* stringId, Float_t *value)
{
	// extracts the sensor value
  // return last value read from sensor specified by stringId
  
  TObjArray* valueSet;
  TPair* pair = (TPair*)dcsAliasMap->FindObject(stringId);
  if (pair) {
    valueSet = (TObjArray*)pair->Value();
    Int_t nentriesDCS = valueSet->GetEntriesFast() - 1;
    if(nentriesDCS>=0){
	AliDCSValue *val = (AliDCSValue *)valueSet->At(nentriesDCS);
	//new
	*value=val->GetFloat();
	return kTRUE;
    }
  }
  return kFALSE;
}

TMap* AliHLTPredictionProcessorHLT::produceTestData(TString aliasName) {
	// produces test data for AliHLTPredictionProcessorHLT
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


