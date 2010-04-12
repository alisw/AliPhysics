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


ClassImp(AliHLTPredictionProcessorHLT)

AliHLTPredictionProcessorHLT::AliHLTPredictionProcessorHLT(
			const char* detector, AliHLTPendolino* pendolino) :
				AliHLTPredictionProcessorInterface(detector, pendolino),
				fPredict(true), fRun(0), fStartTime(0), fEndTime(0) {
	// C-tor for AliHLTPredictionProcessorHLT
//	fPredict = false;
//	fRun = 0;
//	fStartTime = 0;
//	fEndTime = 0;
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
  // there is currently no object to create
  
  return retVal;
}

TMap* AliHLTPredictionProcessorHLT::produceTestData(TString /*aliasName*/) {
	// produces test data for AliHLTPredictionProcessorHLT
    TMap* resultMap = 0;

    // here has to come real dummy data :-)
    resultMap = new TMap();

    return resultMap;
}


