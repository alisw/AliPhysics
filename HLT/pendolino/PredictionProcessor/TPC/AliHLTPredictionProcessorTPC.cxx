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

#include "TROOT.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TString.h"
#include "TMap.h"
#include "TRandom2.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliTPCSensorTempArray.h"

ClassImp(AliHLTPredictionProcessorTPC)

const TString kAmandaTemp = "TPC_PT_%d_TEMPERATURE";
const Int_t kValCutTemp = 100;               // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	             // discard temperature differences > 5 degrees

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
  if (fConfTreeTemp) delete fConfTreeTemp;
  if (fTemp) delete fTemp; 
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

  if (!dcsAliasMap) return 9;
  if (dcsAliasMap->GetEntries() == 0 ) return 9;
  if (!fConfigOK) return 9;

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
       fTemp = (AliTPCSensorTempArray*)entry->GetObject();
   } else {
       Log("Cannot retrieve HCDB entry. New AliTPCSensorTempArray generated");
       UInt_t startTimeLocal = fStartTime-3600;
       UInt_t endTimeLocal = fEndTime+1800;
       fTemp = new AliTPCSensorTempArray(startTimeLocal, endTimeLocal, fConfTreeTemp, kAmandaTemp);
       fTemp->SetValCut(kValCutTemp);
       fTemp->SetDiffCut(kDiffCutTemp);     
       fTemp->ExtractDCS(dcsAliasMap);
   }

   // process temperatures
   
   UInt_t tempResult = ExtractTemperature(dcsAliasMap); 

   //process dcsAliasMap to ROOT object
   // and create AliCDBEntry 

   if (tempResult==0) {
     // create meta data entry for HCDB
     TString comment("HLT temperatures");
     AliCDBMetaData meta(this->GetName(), 666, "unknownAliRoot", comment.Data());

     // store AliCDBEntry in HCDB
     if (Store(path2.Data(), path3.Data(), fTemp, &meta, start, kTRUE)) {
         Log(" +++ Successfully stored object ;-)");
     } else {
         Log(" *** Storing of OBJECT failed!!");
         retVal = 7;
     }
    } else {
      retVal = 9;
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

  TMap *map = fTemp->ExtractDCS(dcsAliasMap,kTRUE);
  if (map) {
    fTemp->MakeSplineFitAddPoints(map);
  } else {
    Log("No temperature map extracted.\n");
    return error;
  }
  return 0;
}

//______________________________________________________________________________________________
  
TMap* AliHLTPredictionProcessorTPC::produceTestData(TString aliasName) 
{
   // produce dummy values for TPC temperature sensors

   const Float_t defTemp = 22.0;
   const Float_t tempVar = 2.0;
   TMap* resultMap = 0;
   TRandom2 random;
   Float_t temp=0;
   Float_t tempDelta=0;
   TObjArray* arr = new TObjArray();

   resultMap = new TMap();
   TTimeStamp tt;

   // loop through all sensors (extracted from OCDB config entry) if no input parameter is given
   
   if (aliasName.Length() == 0 ) {
    if (!fConfTreeTemp) return 0;
    Int_t nentries = fConfTreeTemp->GetEntries();
    if (nentries<1) return 0;

    TString stringId;
    Int_t echa=0;
    fConfTreeTemp->SetBranchAddress("ECha",&echa);

    for (Int_t isensor=0; isensor<nentries; isensor++ ){
      fConfTreeTemp->GetEntry(isensor);
      stringId = Form(kAmandaString.Data(),echa);
      tempDelta = (random.Rndm()-0.5)*tempVar;
      temp = defTemp + tempDelta;
      AliDCSValue* val = new AliDCSValue(temp, tt.GetTime());
      TObjString* name = new TObjString(stringId);   
      arr->Add(val);
      resultMap->Add(name, arr);
    }
   } else {

      // simulate value for given sensor if specified as parameter
 
      tempDelta = (random.Rndm()-0.5)*tempVar;
      temp = defTemp + tempDelta;
      AliDCSValue* val = new AliDCSValue(temp, tt.GetTime());
      TObjString* name = new TObjString(aliasName.Data());   
      arr->Add(val);
      resultMap->Add(name, arr);
   }
   return resultMap;
}
