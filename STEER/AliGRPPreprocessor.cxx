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

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliGRPPreprocessor
//                  Global Run Parameters (GRP) preprocessor
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliGRPPreprocessor.h"
#include "AliGRPDCS.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

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
  UInt_t iStartTime = atoi(timeStart);
  const char* timeEnd = GetRunParameter("time_end");
  UInt_t iStopTime = atoi(timeEnd);
  const char* beamEnergy = GetRunParameter("beamEnergy");
  const char* beamType = GetRunParameter("beamType");
  const char* numberOfDetectors = GetRunParameter("numberOfDetectors");
  const char* detectorMask = GetRunParameter("detectorMask");
  const char* lhcPeriod = GetRunParameter("LHCperiod");

  TList* list = GetFileSources(kDAQ);
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next())) {
    TObjString* objStr = dynamic_cast<TObjString*> (obj);
    if (objStr) {
      Log(Form("Found source %s", objStr->String().Data()));
      TList* list2 = GetFileIDs(kDAQ, objStr->String());
      list2->Print();
      delete list2;
    }
  }
  delete iter;
  delete list;  
  
  //===========//
  //DCS data points
  //===========//
  TObjArray *aliasLHCState = (TObjArray *)valueMap->GetValue("LHCState");
  if(!aliasLHCState) {
    Log(Form("LHCState not found!!!"));
    return 1;
  }
  AliInfo(Form("==========LHCState==========="));
  AliGRPDCS *dcs1 = new AliGRPDCS(aliasLHCState,iStartTime,iStopTime);
  TString sLHCState = dcs1->ProcessDCS(3);  
  if (sLHCState) {
    Log(Form("<LHCState> for run %d: %s",fRun, sLHCState.Data()));
  } else {
    Log(Form("LHCState not put in TMap!"));
  }

  TObjArray *aliasLHCPeriod = (TObjArray *)valueMap->GetValue("LHCPeriod");
  if(!aliasLHCPeriod) {
    Log(Form("LHCPeriod not found!!!"));
    return 1;
  }
  AliInfo(Form("==========LHCPeriod==========="));
  AliGRPDCS *dcs2 = new AliGRPDCS(aliasLHCPeriod,iStartTime,iStopTime);
  TString sLHCPeriod = dcs2->ProcessDCS(3);  
  if (sLHCPeriod) {
    Log(Form("<LHCPeriod> for run %d: %s",fRun, sLHCPeriod.Data()));
  } else {
    Log(Form("LHCPeriod not put in TMap!"));
  }

  TObjArray *aliasLHCLuminosity = (TObjArray *)valueMap->GetValue("LHCLuminosity");
  if(!aliasLHCLuminosity) {
    Log(Form("LHCLuminosity not found!!!"));
    return 1;
  }
  AliInfo(Form("==========LHCLuminosity==========="));
  AliGRPDCS *dcs3 = new AliGRPDCS(aliasLHCLuminosity,iStartTime,iStopTime);
  TString sMeanLHCLuminosity = dcs3->ProcessDCS(2);  
  if (sMeanLHCLuminosity) {
    Log(Form("<LHCLuminosity> for run %d: %s",fRun, sMeanLHCLuminosity.Data()));
  } else {
    Log(Form("LHCLuminosity not put in TMap!"));
  }

  TObjArray *aliasBeamIntensity = (TObjArray *)valueMap->GetValue("BeamIntensity");
  if(!aliasBeamIntensity) {
    Log(Form("BeamIntensity not found!!!"));
    return 1;
  }
  AliInfo(Form("==========BeamIntensity==========="));
  AliGRPDCS *dcs4 = new AliGRPDCS(aliasBeamIntensity,iStartTime,iStopTime);
  TString sMeanBeamIntensity = dcs4->ProcessDCS(2);  
  if (sMeanBeamIntensity) {
    Log(Form("<BeamIntensity> for run %d: %s",fRun, sMeanBeamIntensity.Data()));
  } else {
    Log(Form("BeamIntensity not put in TMap!"));
  }

  TObjArray *aliasL3Current = (TObjArray *)valueMap->GetValue("L3Current");
  if(!aliasL3Current) {
    Log(Form("L3Current not found!!!"));
    return 1;
  }
  AliInfo(Form("==========L3Current==========="));
  AliGRPDCS *dcs5 = new AliGRPDCS(aliasL3Current,iStartTime,iStopTime);
  TString sMeanL3Current = dcs5->ProcessDCS(2);  
  if (sMeanL3Current) {
    Log(Form("<L3Current> for run %d: %s",fRun, sMeanL3Current.Data()));
  } else {
    Log(Form("L3Current not put in TMap!"));
  }

  TObjArray *aliasL3Polarity = (TObjArray *)valueMap->GetValue("L3Polarity");
  if(!aliasL3Polarity) {
    Log(Form("L3Polarity not found!!!"));
    return 1;
  }
  AliInfo(Form("==========L3Polarity==========="));
  AliGRPDCS *dcs6 = new AliGRPDCS(aliasL3Polarity,iStartTime,iStopTime);
  TString sL3Polarity = dcs6->ProcessDCS(4);  
  if (sL3Polarity) {
    Log(Form("<L3Polarity> for run %d: %s",fRun, sL3Polarity.Data()));
  } else {
    Log(Form("L3Polarity not put in TMap!"));
  }

  TObjArray *aliasDipoleCurrent = (TObjArray *)valueMap->GetValue("DipoleCurrent");
  if(!aliasDipoleCurrent) {
    Log(Form("DipoleCurrent not found!!!"));
    return 1;
  }
  AliInfo(Form("==========DipoleCurrent==========="));
  AliGRPDCS *dcs7 = new AliGRPDCS(aliasDipoleCurrent,iStartTime,iStopTime);
  TString sMeanDipoleCurrent = dcs7->ProcessDCS(2);  
  if (sMeanDipoleCurrent) {
    Log(Form("<DipoleCurrent> for run %d: %s",fRun, sMeanDipoleCurrent.Data()));
  } else {
    Log(Form("DipoleCurrent not put in TMap!"));
  }

  TObjArray *aliasDipolePolarity = (TObjArray *)valueMap->GetValue("DipolePolarity");
  if(!aliasDipolePolarity) {
    Log(Form("DipolePolarity not found!!!"));
    return 1;
  }
  AliInfo(Form("==========DipolePolarity==========="));
  AliGRPDCS *dcs8 = new AliGRPDCS(aliasDipolePolarity,iStartTime,iStopTime);
  TString sDipolePolarity = dcs8->ProcessDCS(4);  
  if (sDipolePolarity) {
    Log(Form("<DipolePolarity> for run %d: %s",fRun, sDipolePolarity.Data()));
  } else {
    Log(Form("DipolePolarity not put in TMap!"));
  }

  TObjArray *aliasCavernTemperature = (TObjArray *)valueMap->GetValue("CavernTemperature");
  if(!aliasCavernTemperature) {
    Log(Form("CavernTemperature not found!!!"));
    return 1;
  }
  AliInfo(Form("==========CavernTemperature==========="));
  AliGRPDCS *dcs9 = new AliGRPDCS(aliasCavernTemperature,iStartTime,iStopTime);
  TString sMeanCavernTemperature = dcs9->ProcessDCS(2);  
  if (sMeanCavernTemperature) {
    Log(Form("<CavernTemperature> for run %d: %s",fRun, sMeanCavernTemperature.Data()));
  } else {
    Log(Form("CavernTemperature not put in TMap!"));
  }

  TObjArray *aliasCavernPressure = (TObjArray *)valueMap->GetValue("CavernPressure");
  if(!aliasCavernPressure) {
    Log(Form("CavernPressure not found!!!"));
    return 1;
  }
  AliInfo(Form("==========CavernPressure==========="));
  AliGRPDCS *dcs10 = new AliGRPDCS(aliasCavernPressure,iStartTime,iStopTime);
  TString sMeanCavernPressure = dcs10->ProcessDCS(2);  
  if (sMeanCavernPressure) {
    Log(Form("<CavernPressure> for run %d: %s",fRun, sMeanCavernPressure.Data()));
  } else {
    Log(Form("CavernPressure not put in TMap!"));
  }

  //===========//
  //DAQ logbook
  //===========//
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
  if (beamEnergy) {
    Log(Form("Beam energy for run %d: %s",fRun, beamEnergy));
  } else {
    Log(Form("Beam energy not put in logbook!"));
  }
  if (beamType) {
    Log(Form("Beam type for run %d: %s",fRun, beamType));
  } else {
    Log(Form("Beam type not put in logbook!"));
  }
  if (numberOfDetectors) {
    Log(Form("Number of active detectors for run %d: %s",fRun, numberOfDetectors));
  } else {
    Log(Form("Number of active detectors not put in logbook!"));
  }
  if (detectorMask) {
    Log(Form("Detector mask for run %d: %s",fRun, detectorMask));
  } else {
    Log(Form("Detector mask not put in logbook!"));
  }
  if (lhcPeriod) {
    Log(Form("LHC period (DAQ) for run %d: %s",fRun, lhcPeriod));
  } else {
    Log(Form("LHCperiod not put in logbook!"));
  }

  TList *values = new TList();
  values->SetOwner(1);
  
  //DAQ logbook
  TMap *mapDAQ1 = new TMap();
  mapDAQ1->Add(new TObjString("fAliceStartTime"),new TObjString(timeStart));
  values->Add(mapDAQ1);

  TMap *mapDAQ2 = new TMap();
  mapDAQ2->Add(new TObjString("fAliceStopTime"),new TObjString(timeEnd));
  values->Add(mapDAQ2);

  TMap *mapDAQ3 = new TMap();
  mapDAQ3->Add(new TObjString("fAliceBeamEnergy"),new TObjString(beamEnergy));
  values->Add(mapDAQ3);

  TMap *mapDAQ4 = new TMap();
  mapDAQ4->Add(new TObjString("fAliceBeamType"),new TObjString(beamType));
  values->Add(mapDAQ4);

  TMap *mapDAQ5 = new TMap();
  mapDAQ5->Add(new TObjString("fNumberOfDetectors"),new TObjString(numberOfDetectors));
  values->Add(mapDAQ5);

  TMap *mapDAQ6 = new TMap();
  mapDAQ6->Add(new TObjString("fDetectorMask"),new TObjString(detectorMask));
  values->Add(mapDAQ6);

  TMap *mapDAQ7 = new TMap();
  mapDAQ7->Add(new TObjString("fLHCPeriod"),new TObjString(lhcPeriod));
  values->Add(mapDAQ7);

  //DCS dp
  TMap *mapDCS1 = new TMap();
  mapDCS1->Add(new TObjString("fLHCState"),new TObjString(sLHCState));
  values->Add(mapDCS1);

  TMap *mapDCS2 = new TMap();
  mapDCS2->Add(new TObjString("fLHCCondition"),new TObjString(sLHCPeriod));
  values->Add(mapDCS2);

  TMap *mapDCS3 = new TMap();
  mapDCS3->Add(new TObjString("fLHCLuminosity"),new TObjString(sMeanLHCLuminosity));
  values->Add(mapDCS3);

  TMap *mapDCS4 = new TMap();
  mapDCS4->Add(new TObjString("fBeamIntensity"),new TObjString(sMeanBeamIntensity));
  values->Add(mapDCS4);

  TMap *mapDCS5 = new TMap();
  mapDCS5->Add(new TObjString("fL3Current"),new TObjString(sMeanL3Current));
  values->Add(mapDCS5);

  TMap *mapDCS6 = new TMap();
  mapDCS6->Add(new TObjString("fL3Polarity"),new TObjString(sL3Polarity));
  values->Add(mapDCS6);

  TMap *mapDCS7 = new TMap();
  mapDCS7->Add(new TObjString("fDipoleCurrent"),new TObjString(sMeanDipoleCurrent));
  values->Add(mapDCS7);

  TMap *mapDCS8 = new TMap();
  mapDCS8->Add(new TObjString("fDipolePolarity"),new TObjString(sDipolePolarity));
  values->Add(mapDCS8);

  TMap *mapDCS9 = new TMap();
  mapDCS9->Add(new TObjString("fCavernTemperature"),new TObjString(sMeanCavernTemperature));
  values->Add(mapDCS9);

  TMap *mapDCS10 = new TMap();
  mapDCS10->Add(new TObjString("fCavernPressure"),new TObjString(sMeanCavernPressure));
  values->Add(mapDCS10);

  AliCDBMetaData md;
  md.SetResponsible("Panos");
  
  Bool_t result = Store("GRP", "Values", values, &md);
  
  delete values;
  
  if (result)
    return 0;
  else
    return 1;
}

