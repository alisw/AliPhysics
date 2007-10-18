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

#include <TChain.h>
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliGRPPreprocessor.h"
#include "AliGRPDCS.h"
#include "AliDCSSensorArray.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

class AliDCSValue;
class AliShuttleInterface;

#include <TH1.h>

const Double_t kFitFraction = 0.7;                 // Fraction of DCS sensor fits required

ClassImp(AliGRPPreprocessor)

//_______________________________________________________________
  const char* AliGRPPreprocessor::fgkDCSDataPoints[12] = {"LHCState","LHCPeriod","LHCLuminosity","BeamIntensity","L3Current","L3Polarity","DipoleCurrent","DipolePolarity","CavernTemperature","CavernAtmosPressure","gva_cr5AtmosphericPressure","gva_meyrinAtmosphericPressure"};

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor():
  AliPreprocessor("GRP",0), fPressure(0) {
  // default constructor - Don't use this!
  
}

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
  AliPreprocessor("GRP",shuttle), fPressure(0) {
  // constructor - shuttle must be instantiated!
  
}

//_______________________________________________________________
AliGRPPreprocessor::~AliGRPPreprocessor() {
  //destructor
  delete fPressure;
}

//_______________________________________________________________
void AliGRPPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime) {
  // Initialize preprocessor
  
  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, TTimeStamp(startTime).AsString(), TTimeStamp(endTime).AsString()));
  
  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
  AliInfo("Initialization of the GRP preprocessor.");

  TClonesArray * array = new TClonesArray("AliDCSSensor",2);
  for(Int_t j = 0; j < 2; j++) {
    AliDCSSensor * sens = new ((*array)[j])AliDCSSensor;
    sens->SetStringID(fgkDCSDataPoints[j+10]);
  }
  AliInfo(Form("Pressure Entries: %d",array->GetEntries()));

  fPressure = new AliDCSSensorArray(fStartTime, fEndTime, array);
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::Process(TMap* valueMap) {
  // process data retrieved by the Shuttle
  
  //=================//
  // DAQ logbook     //
  //=================//
  TList *daqlblist = ProcessDaqLB();
  if(!daqlblist) {
    Log(Form("Problem with the DAQ logbook parameters!!!"));
    return 0;
  }
  TMap *m1 = (TMap *)daqlblist->At(0);
  TObjString *s1 = (TObjString *)m1->GetValue("fAliceStartTime");
  UInt_t iStartTime = atoi(s1->String().Data());
  TMap *m2 = (TMap *)daqlblist->At(1);
  TObjString *s2 = (TObjString *)m2->GetValue("fAliceStopTime");
  UInt_t iStopTime = atoi(s2->String().Data());
  TMap *m3 = (TMap *)daqlblist->At(6);
  TObjString *s3 = (TObjString *)m3->GetValue("fLHCPeriod");
  TString productionYear = "";

  //=================//
  // DAQ FXS         //
  //=================//
  UInt_t iDaqFxs = ProcessDaqFxs(s3->String(),productionYear);
  if(iDaqFxs == 0) Log(Form("Raw data merged tags copied succesfully in AliEn!!!"));
 
  //=================//
  // DCS data points //
  //=================//
  TList *dcsdplist = ProcessDcsDPs(valueMap, iStartTime, iStopTime);
  if(!dcsdplist) {
    Log(Form("Problem with the DCS data points!!!"));
    return 0;
  }    
  if(dcsdplist->GetEntries() != 10) {
    Log(Form("Problem with the DCS data points!!!"));
    return 0;
  }
  //NEEDS TO BE REVISED - BREAKS!!!
//   AliDCSSensorArray *dcsSensorArray = GetPressureMap(valueMap,fPressure);
//   if(!dcsSensorArray) {
//     Log(Form("Problem with the pressure sensor values!!!"));
//     return 0;
//   }

  TList * list = new TList();
  list = GetGlobalList(daqlblist,dcsdplist);
  list->SetOwner(1);
  AliInfo(Form("Final list entries: %d",list->GetEntries()));
  
  AliCDBMetaData md;
  md.SetResponsible("Panos Christakoglou");
  md.SetComment("Output parameters from the GRP preprocessor.");
  
  Bool_t result = Store("GRP", "Data", list, &md);
  
  delete list;
  
  if (result)
    return 0;
  else
    return 1;
}

//_______________________________________________________________
TList *AliGRPPreprocessor::GetGlobalList(TList *l1, TList *l2) {
  //Getting the global output TList
  TList *list = new TList();
  TMap *map = new TMap();
  for(Int_t i = 0; i < l1->GetEntries(); i++) 
    list->AddLast(map = (TMap *)l1->At(i));
  for(Int_t i = 0; i < l2->GetEntries(); i++) 
    list->AddLast(map = (TMap *)l2->At(i));

  return list;
}

//_______________________________________________________________
TList *AliGRPPreprocessor::ProcessDaqLB() {
  //Getting the DAQ lb informnation
  const char* timeStart = GetRunParameter("time_start");
  const char* timeEnd = GetRunParameter("time_end");
  const char* beamEnergy = GetRunParameter("beamEnergy");
  const char* beamType = GetRunParameter("beamType");
  const char* numberOfDetectors = GetRunParameter("numberOfDetectors");
  const char* detectorMask = GetRunParameter("detectorMask");
  const char* lhcPeriod = GetRunParameter("LHCperiod");

  if (timeStart) {
    Log(Form("Start time for run %d: %s",fRun, timeStart));
  } else {
    Log(Form("Start time not put in logbook!"));
  }
  TMap *mapDAQ1 = new TMap();
  mapDAQ1->Add(new TObjString("fAliceStartTime"),new TObjString(timeStart));

  if (timeEnd) {
    Log(Form("End time for run %d: %s",fRun, timeEnd));
  } else {
    Log(Form("End time not put in logbook!"));
  }
  TMap *mapDAQ2 = new TMap();
  mapDAQ2->Add(new TObjString("fAliceStopTime"),new TObjString(timeEnd));

  if (beamEnergy) {
    Log(Form("Beam energy for run %d: %s",fRun, beamEnergy));
  } else {
    Log(Form("Beam energy not put in logbook!"));
  }
  TMap *mapDAQ3 = new TMap();
  mapDAQ3->Add(new TObjString("fAliceBeamEnergy"),new TObjString(beamEnergy));

  if (beamType) {
    Log(Form("Beam type for run %d: %s",fRun, beamType));
  } else {
    Log(Form("Beam type not put in logbook!"));
  }
  TMap *mapDAQ4 = new TMap();
  mapDAQ4->Add(new TObjString("fAliceBeamType"),new TObjString(beamType));

  if (numberOfDetectors) {
    Log(Form("Number of active detectors for run %d: %s",fRun, numberOfDetectors));
  } else {
    Log(Form("Number of active detectors not put in logbook!"));
  }
  TMap *mapDAQ5 = new TMap();
  mapDAQ5->Add(new TObjString("fNumberOfDetectors"),new TObjString(numberOfDetectors));

  if (detectorMask) {
    Log(Form("Detector mask for run %d: %s",fRun, detectorMask));
  } else {
    Log(Form("Detector mask not put in logbook!"));
  }
  TMap *mapDAQ6 = new TMap();
  mapDAQ6->Add(new TObjString("fDetectorMask"),new TObjString(detectorMask));

  if (lhcPeriod) {
    Log(Form("LHC period (DAQ) for run %d: %s",fRun, lhcPeriod));
  } else {
    Log(Form("LHCperiod not put in logbook!"));
  }
  TMap *mapDAQ7 = new TMap();
  mapDAQ7->Add(new TObjString("fLHCPeriod"),new TObjString(lhcPeriod));

  TList *list = new TList();
  list->Add(mapDAQ1); list->Add(mapDAQ2);
  list->Add(mapDAQ3); list->Add(mapDAQ4);
  list->Add(mapDAQ5); list->Add(mapDAQ6);
  list->Add(mapDAQ7);

  return list;
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDaqFxs(TString lhcperiod, TString productionYear) {
  //======DAQ FXS======//
  TChain *fRawTagChain = new TChain("T");
  TString fRawDataFileName;
  TList* list = GetFileSources(kDAQ);  
  if (!list) {
    Log("No raw data tag list found!!!");
    return 1;
  }
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next())) {
    TObjString* objStr = dynamic_cast<TObjString*> (obj);
    if (objStr) {
      Log(Form("Found source %s", objStr->String().Data()));
      TList* list2 = GetFileIDs(kDAQ, objStr->String());
      if (!list2) {
	Log("No list with ids from DAQ was found!!!");
	return 2;
      }
      Log(Form("Number of ids: %d",list2->GetEntries()));
      for(Int_t i = 0; i < list2->GetEntries(); i++) {
	TObjString *idStr = (TObjString *)list2->At(i);
	//Log(Form("Filename1: %s",idStr->String().Data()));
	TString fileName = GetFile(kDAQ,idStr->String().Data(),objStr->String().Data());      
	Log(Form("Adding file in the chain: %s",fileName.Data()));
	fRawTagChain->Add(fileName.Data());
	fRawDataFileName = fileName(0,fileName.First("_"));
      }
      delete list2;
    }
  }
  delete iter;
  delete list;
  fRawDataFileName += "_GRP_Merged.tag.root";
  Log(Form("Merging raw data tags into file: %s",fRawDataFileName.Data()));

  TString outputfile = "alien:///alice/data/"; 
  outputfile += productionYear.Data(); outputfile += "/";
  outputfile += lhcperiod.Data(); outputfile += "/";
  outputfile += fRun; outputfile += "/raw/"; 
  //StoreTagFiles(fRawDataFileName.Data(),outputfile.Data());

  return 0;
}

//_______________________________________________________________
TList *AliGRPPreprocessor::ProcessDcsDPs(TMap* valueMap, UInt_t iStartTime, UInt_t iStopTime) {
  //Getting the DCS dps
  //===========//
  
  TList *list = new TList();

  //DCS data points
  //===========//
  AliInfo(Form("==========LHCState==========="));
  TObjArray *aliasLHCState = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[0]);
  if(!aliasLHCState) {
    Log(Form("LHCState not found!!!"));
    return list;
  }
  AliGRPDCS *dcs1 = new AliGRPDCS(aliasLHCState,iStartTime,iStopTime);
  TString sLHCState = dcs1->ProcessDCS(3);  
  if (sLHCState) {
    Log(Form("<LHCState> for run %d: %s",fRun, sLHCState.Data()));
  } else {
    Log(Form("LHCState not put in TMap!"));
  }
  TMap *mapDCS1 = new TMap();
  mapDCS1->Add(new TObjString("fLHCState"),new TObjString(sLHCState));
  list->Add(mapDCS1);

  AliInfo(Form("==========LHCPeriod==========="));
  TObjArray *aliasLHCPeriod = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[1]);
  if(!aliasLHCPeriod) {
    Log(Form("LHCPeriod not found!!!"));
    return list;
  }
  AliGRPDCS *dcs2 = new AliGRPDCS(aliasLHCPeriod,iStartTime,iStopTime);
  TString sLHCPeriod = dcs2->ProcessDCS(3);  
  if (sLHCPeriod) {
    Log(Form("<LHCPeriod> for run %d: %s",fRun, sLHCPeriod.Data()));
  } else {
    Log(Form("LHCPeriod not put in TMap!"));
  }
  TMap *mapDCS2 = new TMap();
  mapDCS2->Add(new TObjString("fLHCCondition"),new TObjString(sLHCPeriod));
  list->Add(mapDCS2);

  AliInfo(Form("==========LHCLuminosity==========="));
  TObjArray *aliasLHCLuminosity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[2]);
  if(!aliasLHCLuminosity) {
    Log(Form("LHCLuminosity not found!!!"));
    return list;
  }
  AliGRPDCS *dcs3 = new AliGRPDCS(aliasLHCLuminosity,iStartTime,iStopTime);
  TString sMeanLHCLuminosity = dcs3->ProcessDCS(2);  
  if (sMeanLHCLuminosity) {
    Log(Form("<LHCLuminosity> for run %d: %s",fRun, sMeanLHCLuminosity.Data()));
  } else {
    Log(Form("LHCLuminosity not put in TMap!"));
  }
  TMap *mapDCS3 = new TMap();
  mapDCS3->Add(new TObjString("fLHCLuminosity"),new TObjString(sMeanLHCLuminosity));
  list->Add(mapDCS3);

  AliInfo(Form("==========BeamIntensity==========="));
  TObjArray *aliasBeamIntensity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[3]);
  if(!aliasBeamIntensity) {
    Log(Form("BeamIntensity not found!!!"));
    return list;
  }
  AliGRPDCS *dcs4 = new AliGRPDCS(aliasBeamIntensity,iStartTime,iStopTime);
  TString sMeanBeamIntensity = dcs4->ProcessDCS(2);  
  if (sMeanBeamIntensity) {
    Log(Form("<BeamIntensity> for run %d: %s",fRun, sMeanBeamIntensity.Data()));
  } else {
    Log(Form("BeamIntensity not put in TMap!"));
  }
  TMap *mapDCS4 = new TMap();
  mapDCS4->Add(new TObjString("fBeamIntensity"),new TObjString(sMeanBeamIntensity));
  list->Add(mapDCS4);

  AliInfo(Form("==========L3Current==========="));
  TObjArray *aliasL3Current = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[4]);
  if(!aliasL3Current) {
    Log(Form("L3Current not found!!!"));
    return list;
  }
  AliGRPDCS *dcs5 = new AliGRPDCS(aliasL3Current,iStartTime,iStopTime);
  TString sMeanL3Current = dcs5->ProcessDCS(2);  
  if (sMeanL3Current) {
    Log(Form("<L3Current> for run %d: %s",fRun, sMeanL3Current.Data()));
  } else {
    Log(Form("L3Current not put in TMap!"));
  }
  TMap *mapDCS5 = new TMap();
  mapDCS5->Add(new TObjString("fL3Current"),new TObjString(sMeanL3Current));
  list->Add(mapDCS5);

  AliInfo(Form("==========L3Polarity==========="));
  TObjArray *aliasL3Polarity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[5]);
  if(!aliasL3Polarity) {
    Log(Form("L3Polarity not found!!!"));
    return list;
  }
  AliGRPDCS *dcs6 = new AliGRPDCS(aliasL3Polarity,iStartTime,iStopTime);
  TString sL3Polarity = dcs6->ProcessDCS(4);  
  if (sL3Polarity) {
    Log(Form("<L3Polarity> for run %d: %s",fRun, sL3Polarity.Data()));
  } else {
    Log(Form("L3Polarity not put in TMap!"));
  }
  TMap *mapDCS6 = new TMap();
  mapDCS6->Add(new TObjString("fL3Polarity"),new TObjString(sL3Polarity));
  list->Add(mapDCS6);

  AliInfo(Form("==========DipoleCurrent==========="));
  TObjArray *aliasDipoleCurrent = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[6]);
  if(!aliasDipoleCurrent) {
    Log(Form("DipoleCurrent not found!!!"));
    return list;
  }
  AliGRPDCS *dcs7 = new AliGRPDCS(aliasDipoleCurrent,iStartTime,iStopTime);
  TString sMeanDipoleCurrent = dcs7->ProcessDCS(2);  
  if (sMeanDipoleCurrent) {
    Log(Form("<DipoleCurrent> for run %d: %s",fRun, sMeanDipoleCurrent.Data()));
  } else {
    Log(Form("DipoleCurrent not put in TMap!"));
  }
  TMap *mapDCS7 = new TMap();
  mapDCS7->Add(new TObjString("fDipoleCurrent"),new TObjString(sMeanDipoleCurrent));
  list->Add(mapDCS7);

  AliInfo(Form("==========DipolePolarity==========="));
  TObjArray *aliasDipolePolarity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[7]);
  if(!aliasDipolePolarity) {
    Log(Form("DipolePolarity not found!!!"));
    return list;
  }
  AliGRPDCS *dcs8 = new AliGRPDCS(aliasDipolePolarity,iStartTime,iStopTime);
  TString sDipolePolarity = dcs8->ProcessDCS(4);  
  if (sDipolePolarity) {
    Log(Form("<DipolePolarity> for run %d: %s",fRun, sDipolePolarity.Data()));
  } else {
    Log(Form("DipolePolarity not put in TMap!"));
  }
  TMap *mapDCS8 = new TMap();
  mapDCS8->Add(new TObjString("fDipolePolarity"),new TObjString(sDipolePolarity));
  list->Add(mapDCS8);

  AliInfo(Form("==========CavernTemperature==========="));
  TObjArray *aliasCavernTemperature = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[8]);
  if(!aliasCavernTemperature) {
    Log(Form("CavernTemperature not found!!!"));
    return list;
  }
  AliGRPDCS *dcs9 = new AliGRPDCS(aliasCavernTemperature,iStartTime,iStopTime);
  TString sMeanCavernTemperature = dcs9->ProcessDCS(2);  
  if (sMeanCavernTemperature) {
    Log(Form("<CavernTemperature> for run %d: %s",fRun, sMeanCavernTemperature.Data()));
  } else {
    Log(Form("CavernTemperature not put in TMap!"));
  }
  TMap *mapDCS9 = new TMap();
  mapDCS9->Add(new TObjString("fCavernTemperature"),new TObjString(sMeanCavernTemperature));
  list->Add(mapDCS9);

  AliInfo(Form("==========CavernPressure==========="));
  TObjArray *aliasCavernPressure = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[9]);
  if(!aliasCavernPressure) {
    Log(Form("CavernPressure not found!!!"));
    return list;
  }
  AliGRPDCS *dcs10 = new AliGRPDCS(aliasCavernPressure,iStartTime,iStopTime);
  TString sMeanCavernPressure = dcs10->ProcessDCS(2);  
  if (sMeanCavernPressure) {
    Log(Form("<CavernPressure> for run %d: %s",fRun, sMeanCavernPressure.Data()));
  } else {
    Log(Form("CavernPressure not put in TMap!"));
  }
  TMap *mapDCS10 = new TMap();
  mapDCS10->Add(new TObjString("fCavernPressure"),new TObjString(sMeanCavernPressure));
  list->Add(mapDCS10);

  return list;
}

//_______________________________________________________________
AliDCSSensorArray *AliGRPPreprocessor::GetPressureMap(TMap* dcsAliasMap, AliDCSSensorArray *fPressure) {
  // extract DCS pressure maps. Perform fits to save space
  
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    Double_t fitFraction = fPressure->NumFits()/fPressure->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Pressure values extracted, fits performed.\n"));
    } else { 
      AliInfo("Too few pressure maps fitted. \n");
    }
  } else {
    AliInfo("AliGRPDCS: no atmospheric pressure map extracted. \n");
  }
  delete map;
 
  return fPressure;
}

//_______________________________________________________________
/*UInt_t AliGRPPreprocessor::MapPressure(TMap* dcsAliasMap) {
  // extract DCS pressure maps. Perform fits to save space
  
  UInt_t result=0;
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    Double_t fitFraction = fPressure->NumFits()/fPressure->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Pressure values extracted, fits performed.\n"));
    } else { 
      Log ("Too few pressure maps fitted. \n");
      result = 9;
    }
  } else {
    Log("AliTPCPreprocsessor: no atmospheric pressure map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file
  
  if ( result == 0 ) {
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Panos Christakoglou");
    metaData.SetComment("Preprocessor AliGRP data base pressure entries.");
    
    Bool_t storeOK = Store("Calib", "Pressure", fPressure, &metaData, 0, 0);
    if ( !storeOK ) result=1; 
  }
  
  return result; 
  }*/
