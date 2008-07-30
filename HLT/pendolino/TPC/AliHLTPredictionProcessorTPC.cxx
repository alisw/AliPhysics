// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Haavard Helstrup                                      *
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

/** @file   AliHLTPredictionProcessorTPC.cxx
    @author Haavard Helstrup
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorTPC.h"

#include "AliHLTDCSArray.h"
#include "TROOT.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TString.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliDCSValue.h"

#include <TTimeStamp.h>
#include <TObjString.h>


ClassImp(AliHLTPredictionProcessorTPC)

const TString kAmandaString = "TPC_PT_%d_TEMPERATURE";

//______________________________________________________________________________________________

AliHLTPredictionProcessorTPC::AliHLTPredictionProcessorTPC(
           const char* detector, AliHLTPendolino* pendolino) :
               AliHLTPredictionProcessorInterface(detector, pendolino),
               fConfigOK(kTRUE),
	       fConfTreeTemp(0),
               fTemp(0),
	       fPredict(kFALSE),
	       fRun(0),
	       fStartTime(0),
	       fEndTime(0)
{
 // constructor
}

//______________________________________________________________________________________________

AliHLTPredictionProcessorTPC::~AliHLTPredictionProcessorTPC()
{

  // destructor
	if (fConfTreeTemp) {
		delete fConfTreeTemp;
	}
	if (fTemp) {
		fTemp->Delete();
		delete fTemp;
	}

}

//______________________________________________________________________________________________

UInt_t AliHLTPredictionProcessorTPC::makePrediction(Bool_t doPrediction) {
//
// Signal whether prediction making is required
//
   fPredict = doPrediction; // own member indicating that prediction making is required
   return 0;
}

//______________________________________________________________________________________________

void AliHLTPredictionProcessorTPC::Initialize(Int_t run, UInt_t startTime,
           UInt_t endTime) {
//
//  Initialise TPC prediction processor. Read config entry from HCDB
//
   fRun = run;               // storing locally the run number
   fStartTime = startTime;   // storing locally the start time
   fEndTime = endTime;       // storing locally the end time


   AliCDBEntry* entry = GetFromOCDB("Config", "HLTTemperature");
   if (entry) fConfTreeTemp = (TTree*) entry->GetObject();
   if ( fConfTreeTemp==0 ) {
     Log("AliHLTPredictionProcessorTPC: Temperature Config HCDB entry missing.\n");
     fConfigOK = kFALSE;
     return;
   }
}

//______________________________________________________________________________________________

UInt_t AliHLTPredictionProcessorTPC::Process(TMap* dcsAliasMap) {
//
// Event-by-event processing of the PredictionProcessorTPC
//

// test if configuration available and DCS map valid

  if (!dcsAliasMap) return 91;
  if (dcsAliasMap->GetEntries() == 0 ) return 92;
  if (!fConfigOK) return 93;

//
// Extract values from DCS maps. Store updated entry to HCDB
//
   UInt_t retVal = 0;
   Int_t start = 0;
   TString path2("Calib");        // for the storage path in HCDB
   TString path3("Temperature");    // for the storage path in HCDB

   //get data from the HCDB
   AliCDBEntry* entry = GetFromOCDB(path2.Data(), path3.Data());
   if (entry != 0) {
       entry->PrintMetaData();  // use the AliCDBEntry
       fTemp = (TObjArray*)entry->GetObject();
   } else {
       Log("Cannot retrieve HCDB entry. New TObjArray generated");
       fTemp = new TObjArray();
   }

   // process temperatures
   
   UInt_t tempResult = ExtractTemperature(dcsAliasMap); 

   //process dcsAliasMap to ROOT object
   // and create AliCDBEntry 

   if (tempResult == 0) {
		// create meta data entry for HCDB
     	TString comment("HLT temperatures");
     	AliCDBMetaData meta(this->GetName(), 0, "unknownAliRoot", comment.Data());

     	// store AliCDBEntry in HCDB
     	if (Store(path2.Data(), path3.Data(), fTemp, &meta, start, kTRUE)) {
         	Log(" +++ Successfully stored object ;-)");
     	} else {
         	Log(" *** Storing of OBJECT failed!!");
         	retVal = 7;
     	}
    } else {
      	retVal = 94;
    }
   
   return retVal;
}

//______________________________________________________________________________________________

UInt_t AliHLTPredictionProcessorTPC::ExtractTemperature(TMap* dcsAliasMap)
{
// Extract temperature values from DCS maps, according to infoamtion from configuration tree

  const Int_t error = 9; 
  Int_t nentries = fConfTreeTemp->GetEntries();
  if (nentries<1) return error;
    
  TString stringId;
  Int_t echa=0;
  fConfTreeTemp->SetBranchAddress("ECha",&echa);
  fConfTreeTemp->GetEntry(0);
  stringId = Form(kAmandaString.Data(),echa);
  UInt_t time = GetCurrentTime(dcsAliasMap,stringId.Data());   
  AliHLTDCSArray* temperatures = new AliHLTDCSArray(nentries);
  temperatures->SetTime(time);
  
  for (Int_t isensor=0; isensor<nentries; isensor++ ){
     fConfTreeTemp->GetEntry(isensor);
     stringId = Form(kAmandaString.Data(),echa);
     Float_t temp = GetSensorValue(dcsAliasMap,stringId.Data());
     temperatures->SetValue(isensor,temp);
  }
  
  fTemp->Add(temperatures);
  return 0;  
}

  
//______________________________________________________________________________________________

Float_t AliHLTPredictionProcessorTPC::GetSensorValue(TMap* dcsAliasMap,
                                     const char* stringId)
{
  // return last value read from sensor specified by stringId
  Float_t sensorValue=0;
  TObjArray* valueSet;
  TPair* pair = (TPair*)dcsAliasMap->FindObject(stringId);
  if (pair) {
     valueSet = (TObjArray*)pair->Value();
	 if (valueSet) {
     	Int_t nentriesDCS = valueSet->GetEntriesFast();
     	AliDCSValue *val = (AliDCSValue *)valueSet->At(nentriesDCS - 1);
     	if (val) {
			sensorValue=val->GetFloat();
		}
	 }
  }
  return sensorValue;
}
    	
//______________________________________________________________________________________________

UInt_t AliHLTPredictionProcessorTPC::GetCurrentTime(TMap* dcsAliasMap,
                                     const char* stringId)
{
  // return last time value read from sensor specified by stringId
  
  UInt_t time=0;
  TObjArray* valueSet;
  TPair* pair = (TPair*) (dcsAliasMap->FindObject(stringId));
  if (pair) {
     valueSet = (TObjArray*) (pair->Value());
	 if (valueSet) {
     	Int_t nentriesDCS = valueSet->GetEntriesFast();
     	AliDCSValue *val = (AliDCSValue *) (valueSet->At(nentriesDCS - 1));
	 	if (val) {
     		time = val->GetTimeStamp();
		}
	 }
  }
  return time;
}

TMap* AliHLTPredictionProcessorTPC::produceTestData(TString aliasName) {
    TMap* resultMap = 0;

    // here has to come real dummy data :-)
    resultMap = new TMap();
    TTimeStamp tt;
    TObjString* name = new TObjString("DummyData");
	Float_t fval = 33.3;
    AliDCSValue* val = new AliDCSValue(fval, tt.GetTime());
    TObjArray* arr = new TObjArray();
    arr->Add(val);
    resultMap->Add(name, arr);

    return resultMap;
}

			     

