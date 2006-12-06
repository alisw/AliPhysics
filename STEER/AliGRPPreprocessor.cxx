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

//-------------------------------------------------------------------------
//                          Class AliGRPPreprocessor
//                  Global Run Parameters (GRP) preprocessor
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "AliGRPPreprocessor.h"
#include "AliGRPDCS.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TMap.h>
#include <TObjString.h>

class AliDCSValue;
class AliShuttleInterface;

#include <TH1.h>

ClassImp(AliGRPPreprocessor)

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor():
  AliPreprocessor("GRP",0) {
  // default constructor - Don't use this!
  
}

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
  AliPreprocessor("GRP",shuttle) {
  // constructor - shuttle must be instantiated!
  
}

//_______________________________________________________________
AliGRPPreprocessor::~AliGRPPreprocessor() {
  //destructor
}

//_______________________________________________________________
void AliGRPPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime) {
  // Initialize preprocessor
  
  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, TTimeStamp(startTime).AsString(), TTimeStamp(endTime).AsString()));
  
  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
  AliInfo("This preprocessor is to test the GetRunParameter function.");
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::Process(TMap* valueMap) {
  // process data retrieved by the Shuttle
  const char* timeStart = GetRunParameter("time_start");
  const char* timeEnd = GetRunParameter("time_end");
  TObjArray *alias1 = (TObjArray *)valueMap->GetValue("SFTTemp1.FloatValue");
  AliGRPDCS *dcs = new AliGRPDCS(alias1);
  TH1F *h1 = new TH1F("alias1","",100,15,25);
  TString sAlias1Mean = dcs->ProcessDCS(h1);  
  
  Int_t result=0;
  
  if (sAlias1Mean) {
    Log(Form("<alias1> for run %d: %s",fRun, sAlias1Mean.Data()));
  } else {
    Log(Form("DCSAlias1 not put in TMap!"));
  }
  if (timeStart) {
    Log(Form("Start time for run %d: %s",fRun, timeStart));
  } else {
    Log(Form("Start time not put in logbook!"));
  }
  if (timeEnd) {
    Log(Form("End time for run %d: %s",fRun, timeEnd));
  } else {
    Log(Form("End time not put in logbook!"));
  }
  
  TList *values = new TList();
  values->SetOwner(1);
  
  TMap *map1 = new TMap();
  map1->Add(new TObjString("histoDCS1"),h1);
  values->Add(map1);

  TMap *map2 = new TMap();
  map2->Add(new TObjString("DCS1"),new TObjString(sAlias1Mean));
  values->Add(map2);

  TMap *map3 = new TMap();
  map3->Add(new TObjString("fAliceStartTime"),new TObjString(timeStart));
  values->Add(map3);

  TMap *map4 = new TMap();
  map4->Add(new TObjString("fAliceStopTime"),new TObjString(timeEnd));
  values->Add(map4);

  AliCDBMetaData md;
  md.SetResponsible("Panos");
  
  result = Store("GRP", "Values", values, &md);
  
  delete values;
  
  return result;
}

