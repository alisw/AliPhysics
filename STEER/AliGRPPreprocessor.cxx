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
#include <TSystem.h>
#include <TFile.h>

#include "AliGRPPreprocessor.h"
#include "AliGRPDCS.h"
#include "AliDCSSensorArray.h"

#include "AliTriggerConfiguration.h"
#include "AliTriggerRunScalers.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

class AliDCSValue;
class AliShuttleInterface;

#include <TH1.h>

// needed for ReceivePromptRecoParameters
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliCDBId.h>
#include <AliTriggerConfiguration.h>

const Double_t kFitFraction = 0.7;                 // Fraction of DCS sensor fits required

ClassImp(AliGRPPreprocessor)

//_______________________________________________________________
  const char* AliGRPPreprocessor::fgkDCSDataPoints[12] = {"LHCState","LHCPeriod","LHCLuminosity","BeamIntensity","L3Current","L3Polarity","DipoleCurrent","DipolePolarity","CavernTemperature","CavernAtmosPressure","gva_cr5AtmosphericPressure","gva_meyrinAtmosphericPressure"};

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
  AliPreprocessor("GRP",shuttle), fPressure(0) {
  // constructor - shuttle must be instantiated!

  AddRunType("PHYSICS");
}

//_______________________________________________________________
AliGRPPreprocessor::~AliGRPPreprocessor() {
  //destructor
  delete fPressure;
}

//_______________________________________________________________
void AliGRPPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime) {
  // Initialize preprocessor

  AliPreprocessor::Initialize(run, startTime, endTime);
    
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
    return 1;
  }

  //=================//
  // DAQ FXS         //
  //=================//
  UInt_t iDaqFxs = ProcessDaqFxs();
  if(iDaqFxs == 0) {
  	Log(Form("ProcessDaqFxs successful!"));
  } else {
  	Log(Form("Could not store run raw tag file!"));
	return 1;
  }
  
  //=================//
  // DCS FXS         //
  //=================//
  UInt_t iDcsFxs = ProcessDcsFxs();
  if(iDcsFxs == 0) {
  	Log(Form("ProcessDcsFxs successful!"));
  } else {
  	Log(Form("Could not store CTP run configuration and scalers!"));
	return 1;
  }
  
  //=================//
  // DCS data points //
  //=================//
  TList *dcsdplist = ProcessDcsDPs(valueMap);
  if(!dcsdplist) {
    Log(Form("Problem with the DCS data points!!!"));
    return 1; 
  }    
  if(dcsdplist->GetEntries() != 10) {
    Log(Form("Problem with the DCS data points!!!"));
    // return 1; // TODO:COMMENTED FOR TESTING PURPOSES!
  }
  //NEEDS TO BE REVISED - BREAKS!!!
//   AliDCSSensorArray *dcsSensorArray = GetPressureMap(valueMap,fPressure);
//   if(!dcsSensorArray) {
//     Log(Form("Problem with the pressure sensor values!!!"));
//     return 0;
//   }

  daqlblist->AddAll(dcsdplist);
  daqlblist->SetOwner(1);
  AliInfo(Form("Final list entries: %d",daqlblist->GetEntries()));
  
  AliCDBMetaData md;
  md.SetResponsible("Panos Christakoglou");
  md.SetComment("Output parameters from the GRP preprocessor.");
  
  Bool_t result = Store("GRP", "Data", daqlblist, &md);
  
  delete daqlblist;
  
  if (result)
    return 0;
  else
    return 1;
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
  
  TMap* mapDAQ8 = new TMap;
  mapDAQ8->Add(new TObjString("fRunType"), new TObjString(GetRunType()));
  list->Add(mapDAQ8);

  return list;
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDaqFxs() {
  //======DAQ FXS======//

  TList* list = GetFileSources(kDAQ);  
  if (!list) {
    Log("No raw data tag list: connection problems with DAQ FXS logbook!");
    return 1;
  }
  
  if (list->GetEntries() == 0)
  {
  	Log("no raw data tags in this run: nothing to merge!");
	delete  list; list=0;
	return 0;
  }

  TChain *fRawTagChain = new TChain("T");
  Int_t nFiles=0;
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next())) {
    TObjString* objStr = dynamic_cast<TObjString*> (obj);
    if (objStr) {
      Log(Form("Found source %s", objStr->String().Data()));
      TList* list2 = GetFileIDs(kDAQ, objStr->String());
      if (!list2) {
	Log("No list with ids from DAQ was found: connection problems with DAQ FXS logbook!");
	delete fRawTagChain; fRawTagChain=0;
	return 1;
      }
      Log(Form("Number of ids: %d",list2->GetEntries()));
      for(Int_t i = 0; i < list2->GetEntries(); i++) {
	TObjString *idStr = (TObjString *)list2->At(i);
	TString fileName = GetFile(kDAQ,idStr->String().Data(),objStr->String().Data());
	if (fileName.Length() > 0) 
	{      
		Log(Form("Adding file in the chain: %s",fileName.Data()));
		fRawTagChain->Add(fileName.Data());
		nFiles++;
	} else {
		Log(Form("Could not retrieve file with id %s from source %s: "
			"connection problems with DAQ FXS!",
				idStr->String().Data(),objStr->String().Data()));
		delete list; list=0;
		delete list2; list2=0;
		delete fRawTagChain; fRawTagChain=0;
		return 2;
	}
      }
      delete list2;
    }
  }
  
  TString fRawDataFileName = "GRP_Merged.tag.root";
  Log(Form("Merging %d raw data tags into file: %s", nFiles, fRawDataFileName.Data()));
  fRawTagChain->Merge(fRawDataFileName);
  
  TString outputfile = Form("Run%d.Merged.RAW.tag.root", fRun);
  Bool_t result = StoreRunMetadataFile(fRawDataFileName.Data(),outputfile.Data());
  
  if (!result)
  {
  	Log("Problem storing raw data tags in local file!!");
  } else {
  	Log("Raw data tags merged successfully!!");  
  }
  
  delete iter;
  delete list;
  delete fRawTagChain; fRawTagChain=0;
  
  if (result == kFALSE)
  {
  	return 3;
  }
  
  return 0;
  
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDcsFxs() {
  //======DCS FXS======//
  // Get the CTP run configuration
  // and scalers from DCS FXS

  {
    // Get the CTP run configuration
    TList* list = GetFileSources(kDCS,"CTP_runconfig");  
    if (!list) {
      Log("No CTP runconfig file: connection problems with DCS FXS logbook!");
      return 1;
    }
  
    if (list->GetEntries() == 0) {
      Log("No CTP runconfig file to be processed!");
    }
    else {
      TIter iter(list);
      TObjString *source;
      while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
	TString runcfgfile = GetFile(kDCS, "CTP_runconfig", source->GetName());
	if (runcfgfile.IsNull()) {
	  Log("No CTP runconfig files has been found: empty source!");
	}
	else {
	  Log(Form("File with Id CTP_runconfig found in source %s! Copied to %s",source->GetName(),runcfgfile.Data()));
	  AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfiguration(runcfgfile);
	  if (!runcfg) {
	    Log("Bad CTP run configuration file! The corresponding CDB entry will not be filled!");
	  }
	  else {
	    AliCDBMetaData metaData;
	    metaData.SetBeamPeriod(0);
	    metaData.SetResponsible("Roman Lietava");
	    metaData.SetComment("CTP run configuration");
	    if (!Store("CTP","Config", runcfg, &metaData, 0, 0)) {
	      Log("Unable to store the CTP run configuration object to OCDB!");
	    }
	  }
	}
      }
    }
    delete list;
  }

  {
    // Get the CTP counters information
    TList* list = GetFileSources(kDCS,"CTP_xcounters");  
    if (!list) {
      Log("No CTP counters file: connection problems with DAQ FXS logbook!");
      return 1;
    }
  
    if (list->GetEntries() == 0) {
      Log("No CTP counters file to be processed!");
    }
    else {
      TIter iter(list);
      TObjString *source;
      while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
	TString countersfile = GetFile(kDCS, "CTP_xcounters", source->GetName());
	if (countersfile.IsNull()) {
	  Log("No CTP counters files has been found: empty source!");
	}
	else {
	  Log(Form("File with Id CTP_xcounters found in source %s! Copied to %s",source->GetName(),countersfile.Data()));
	  AliTriggerRunScalers *scalers = AliTriggerRunScalers::ReadScalers(countersfile);
	  if (!scalers) {
	    Log("Bad CTP counters file! The corresponding CDB entry will not be filled!");
	  }
	  else {
	    AliCDBMetaData metaData;
	    metaData.SetBeamPeriod(0);
	    metaData.SetResponsible("Roman Lietava");
	    metaData.SetComment("CTP scalers");
	    if (!Store("CTP","Scalers", scalers, &metaData, 0, 0)) {
	      Log("Unable to store the CTP scalers object to OCDB!");
	    }
	  }
	}
      }
    }
    delete list;
  }

  return 0;
}

//_______________________________________________________________
TList *AliGRPPreprocessor::ProcessDcsDPs(TMap* valueMap) {
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
  AliGRPDCS *dcs1 = new AliGRPDCS(aliasLHCState,fStartTime,fEndTime);
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
  AliGRPDCS *dcs2 = new AliGRPDCS(aliasLHCPeriod,fStartTime,fEndTime);
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
  AliGRPDCS *dcs3 = new AliGRPDCS(aliasLHCLuminosity,fStartTime,fEndTime);
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
  AliGRPDCS *dcs4 = new AliGRPDCS(aliasBeamIntensity,fStartTime,fEndTime);
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
  AliGRPDCS *dcs5 = new AliGRPDCS(aliasL3Current,fStartTime,fEndTime);
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
  AliGRPDCS *dcs6 = new AliGRPDCS(aliasL3Polarity,fStartTime,fEndTime);
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
  AliGRPDCS *dcs7 = new AliGRPDCS(aliasDipoleCurrent,fStartTime,fEndTime);
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
  AliGRPDCS *dcs8 = new AliGRPDCS(aliasDipolePolarity,fStartTime,fEndTime);
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
  AliGRPDCS *dcs9 = new AliGRPDCS(aliasCavernTemperature,fStartTime,fEndTime);
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
  AliGRPDCS *dcs10 = new AliGRPDCS(aliasCavernPressure,fStartTime,fEndTime);
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

  
//_______________________________________________________________
Int_t AliGRPPreprocessor::ReceivePromptRecoParameters(UInt_t run, const char* dbHost, Int_t dbPort, const char* dbName, const char* user, const char* password, const char *cdbRoot)
{
	//
	// Retrieves logbook and trigger information from the online logbook 
	// This information is needed for prompt reconstruction
	//
	// Parameters are:
	// Run number
	// DAQ params: dbHost, dbPort, dbName, user, password, logbookTable, triggerTable
	// cdbRoot
	//
	// returns:
	//         positive on success: the return code is the run number of last run processed of the same run type already processed by the SHUTTLE
	//         0 on success and no run was found
	//         negative on error
	//
	// This function is NOT called during the preprocessor run in the Shuttle!
	//
	
	// defaults
	if (dbPort == 0)
		dbPort = 3306;
		
	// CDB connection
	AliCDBManager* cdb = AliCDBManager::Instance();
	cdb->SetDefaultStorage(cdbRoot);
	
	// SQL connection
	TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d/%s", dbHost, dbPort, dbName), user, password);
	
	if (!server)
	{
		Printf("ERROR: Could not connect to DAQ LB");
		return -1;
	}
	
	// main logbook
	TString sqlQuery;
	sqlQuery.Form("SELECT time_start, run_type, detectorMask FROM logbook WHERE run = %d", run);
	TSQLResult* result = server->Query(sqlQuery);
	if (!result) 
	{
		Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
		return -2;
	}

	if (result->GetRowCount() == 0) 
	{
		Printf("ERROR: Run %d not found", run);
		delete result;
		return -3;
	}

	TSQLRow* row = result->Next();
	if (!row)
	{
		Printf("ERROR: Could not receive data from run %d", run);
		delete result;
		return -4;
	}
	
	TString runType(row->GetField(1));
	
	TMap grpData;
	grpData.Add(new TObjString("time_start"), new TObjString(row->GetField(0)));
	grpData.Add(new TObjString("run_type"), new TObjString(runType));
	grpData.Add(new TObjString("detectorMask"), new TObjString(row->GetField(2)));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Printf("Storing GRP/GRP/Data object with the following content");
	grpData.Print();

	AliCDBMetaData metadata;
	metadata.SetResponsible("Jan Fiete Grosse-Oetringhaus");
	metadata.SetComment("GRP Output parameters received during online running");
	
	AliCDBId id("GRP/GRP/Data", run, run);
	Bool_t success = cdb->Put(&grpData, id, &metadata);
	
	grpData.DeleteAll();
	
	if (!success)
	{
		Printf("ERROR: Could not store GRP/GRP/Data into OCDB");
		return -5;
	}
	
	sqlQuery.Form("SELECT configFile FROM logbook_trigger_config WHERE run = %d", run);
	result = server->Query(sqlQuery);
	if (!result) 
	{
		Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
		return -11;
	}

	if (result->GetRowCount() == 0) 
	{
		Printf("ERROR: Run %d not found in logbook_trigger_config", run);
		delete result;
		return -12;
	}

	row = result->Next();
	if (!row)
	{
		Printf("ERROR: Could not receive logbook_trigger_config data from run %d", run);
		delete result;
		return -13;
	}
	
	TString triggerConfig(row->GetField(0));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Printf("Found trigger configuration: %s", triggerConfig.Data());
	
	// add a function that takes the configuration from a string...
	AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfigurationFromString(triggerConfig);
	if (!runcfg) 
	{
		Printf("ERROR: Could not create CTP configuration object");
		return -14;
	}
	
	metadata.SetComment("CTP run configuration received during online running");
	
	AliCDBId id2("GRP/CTP/Config", run, run);
	success = cdb->Put(runcfg, id2, &metadata);
	
	delete runcfg;
	runcfg = 0;
	
	if (!success)
	{
		Printf("ERROR: Could not store GRP/CTP/Config into OCDB");
		return -15;
	}
	
	// get last run with same run type that was already processed by the SHUTTLE
	
	sqlQuery.Form("SELECT max(logbook.run) FROM logbook LEFT JOIN logbook_shuttle ON logbook_shuttle.run = logbook.run WHERE run_type = '%s' AND shuttle_done = 1", runType.Data());
	result = server->Query(sqlQuery);
	if (!result) 
	{
		Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
		return -21;
	}

	if (result->GetRowCount() == 0) 
	{
		Printf("ERROR: No result with query <%s>", sqlQuery.Data());
		delete result;
		return -22;
	}

	row = result->Next();
	if (!row)
	{
		Printf("ERROR: Could not receive data for query <%s>", sqlQuery.Data());
		delete result;
		return -23;
	}
	
	TString lastRunStr(row->GetField(0));
	Int_t lastRun = lastRunStr.Atoi();
	
	Printf("Last run with same run type %s is %d", runType.Data(), lastRun);
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
		
	server->Close();
	delete server;
	server = 0;
	
	return lastRun;
}
