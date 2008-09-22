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
//    Modified: Ernesto.Lopez.Torres@cern.ch  CEADEN-CERN
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
  const Int_t AliGRPPreprocessor::fgknDAQLbPar = 8; // num parameters in the logbook
  const Int_t AliGRPPreprocessor::fgknDCSDP = 10;   // number of dcs dps
  const char* AliGRPPreprocessor::fgkDCSDataPoints[AliGRPPreprocessor::fgknDCSDP] = {
                   "LHCState",              // missing in DCS
                   "L3Polarity",
                   "DipolePolarity",
                   "LHCLuminosity",         // missing in DCS
                   "BeamIntensity",         // missing in DCS
                   "L3Current",
                   "DipoleCurrent",
                   "CavernTemperature",
                   "CavernAtmosPressure",
                   "SurfaceAtmosPressure"
                 };
                 
  const Short_t kSensors = 9; // start index position of sensor in DCS DPs
  const Short_t kNumSensors = 1; // Number of sensors in DCS DPs

  const char* AliGRPPreprocessor::fgkLHCState[20] = {
                   "P", "PREPARE",
                   "J", "PREINJECTION",
                   "I", "INJECTION",
                   "F", "FILLING",
                   "A", "ADJUST",
                   "U", "UNSTABLE BEAMS",
                   "S", "STABLE BEAMS",
                   "D", "BEAM DUMP",
                   "R", "RECOVER",
                   "C", "PRECYCLE"
                 };

  const char* kppError[] = {
                   "",
                   "(DAQ logbook ERROR)",
                   "(DAQ FXS ERROR)",
                   "(DCS FXS ERROR)",
                   "(DCS data points ERROR)",
                   "(Trigger Configuration ERROR)"
  };

//_______________________________________________________________
AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
	AliPreprocessor("GRP",shuttle),  fPressure(0)
{
  // constructor - shuttle must be instantiated!

  AddRunType("PHYSICS");
}

//_______________________________________________________________
AliGRPPreprocessor::~AliGRPPreprocessor()
{
  //destructor

  delete fPressure;
}

//_______________________________________________________________
void AliGRPPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  // Initialize preprocessor

  AliPreprocessor::Initialize(run, startTime, endTime);

  AliInfo("Initialization of the GRP preprocessor.");
  TClonesArray * array = new TClonesArray("AliDCSSensor",kNumSensors); 
  for(Int_t j = 0; j < kNumSensors; j++) {
    AliDCSSensor * sens = new ((*array)[j])AliDCSSensor;
    sens->SetStringID(fgkDCSDataPoints[j+kSensors]);
  }
  AliInfo(Form("Pressure Entries: %d",array->GetEntries()));

  //  fPressure = new AliDCSSensorArray(fStartTime, fEndTime, array);
  fPressure = new AliDCSSensorArray(GetStartTimeDCSQuery(), GetEndTimeDCSQuery(), array);
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::Process(TMap* valueMap)
{
  // process data retrieved by the Shuttle
  
  //=================//
  // DAQ logbook     //
  //=================//
  UInt_t error = 0;
  
  TMap *grpmap = ProcessDaqLB();
  if( grpmap->GetEntries() == fgknDAQLbPar ) {
    Log(Form("DAQ logbook, successful!"));
  } else {
    Log(Form("DAQ logbook, missing parameters!!!"));
    error |= 1;
  }
  //=================//
  // DAQ FXS         //
  //=================//
  UInt_t iDaqFxs = ProcessDaqFxs();
  if( iDaqFxs == 0 ) {
    Log(Form("DAQ FXS, successful!"));
  } else {
    Log(Form("DAQ FXS, could not store run raw tag file!!!"));
    error |= 2;
  }
  
  //=================//
  // DCS FXS         //
  //=================//
  UInt_t iDcsFxs = ProcessDcsFxs();
  if( iDcsFxs == 0 ) {
     Log(Form("DCS FXS, successful!"));
  } else {
     Log(Form("DCS FXS, Could not store CTP run configuration and scalers!!!"));
    error |= 4;
  }
  
  //=================//
  // DCS data points //
  //=================//
  Int_t entries = ProcessDcsDPs( valueMap, grpmap );
  if( entries < fgknDCSDP-3 ) { // FIXME (!= ) LHState and pressure map are not working yet...
    Log(Form("Problem with the DCS data points!!!"));
    error |= 8;
  } else  Log(Form("DCS data points, successful!"));

  //=======================//
  // Trigger Configuration //
  //=======================//
  // either from DAQ logbook.....
  const char * triggerConf = GetTriggerConfiguration();
  if (triggerConf!= NULL) {
    Log("Found trigger configuration in DAQ logbook");
    AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfigurationFromString(triggerConf);	  
    if (!runcfg) {
      Log("Bad CTP run configuration file from DAQ logbook! The corresponding CDB entry will not be filled!");
      error |= 16;
    }
    else {
      TString titleCTPcfg = Form("CTP cfg for run %i from DAQ",fRun);
      runcfg->SetTitle(titleCTPcfg);
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible("Roman Lietava");
      metaData.SetComment("CTP run configuration from DAQ logbook");
      if (!Store("CTP","Config", runcfg, &metaData, 0, 0)) {
        Log("Unable to store the CTP run configuration object to OCDB!");
	error |= 16;
      }
    }
  }
  // ...or from DCS FXS
  else{
     Log("No trigger configuration found in the DAQ logbook!! Trying reading from DCS FXS...");
     TString runcfgfile = GetFile(kDCS, "CTP_runconfig", "");
     if (runcfgfile.IsNull()) {
       Log("No CTP runconfig files has been found in DCS FXS!");
       error |= 16;
     }
     else {
       Log(Form("File with Id CTP_runconfig found! Copied to %s",runcfgfile.Data()));
       AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfiguration(runcfgfile);
       if (!runcfg) {
         Log("Bad CTP run configuration file from DCS FXS! The corresponding CDB entry will not be filled!");
         error |= 16;;
       }
       else {
	 TString titleCTPcfg = Form("CTP cfg for run %i from DCS",fRun);
         runcfg->SetTitle(titleCTPcfg);
         AliCDBMetaData metaData;
         metaData.SetBeamPeriod(0);
         metaData.SetResponsible("Roman Lietava");
         metaData.SetComment("CTP run configuration from DCS FXS");
         if (!Store("CTP","Config", runcfg, &metaData, 0, 0)) {
           Log("Unable to store the CTP run configuration object to OCDB!");
           error |= 16;
         }
       }
     }
  }

  grpmap->SetOwner(1);
  AliInfo(Form("Final list entries: %d",grpmap->GetEntries()));
  
  AliCDBMetaData md;
  md.SetResponsible("Ernesto Lopez Torres");
  md.SetComment("Output parameters from the GRP preprocessor.");
  
  Bool_t result = Store("GRP", "Data", grpmap, &md);
  
  delete grpmap;
  
  if (result && !error ) {
    Log("GRP Preprocessor Success");
    return 0;
  } else {
    Log( Form("GRP Preprocessor FAILS!!! %s%s%s%s%s",
                                 kppError[(error&1)?1:0],
                                 kppError[(error&2)?2:0],
                                 kppError[(error&4)?3:0],
	                         kppError[(error&8)?4:0],
                                 kppError[(error&16)?5:0]
                                  ));
    return error;
  }
}

//_______________________________________________________________
TMap *AliGRPPreprocessor::ProcessDaqLB()
{
  //Getting the DAQ lb information
  
  const char* timeStart         = GetRunParameter("DAQ_time_start");
  const char* timeEnd           = GetRunParameter("DAQ_time_end");
  const char* beamEnergy        = GetRunParameter("beamEnergy");
  const char* beamType          = GetRunParameter("beamType");
  const char* numberOfDetectors = GetRunParameter("numberOfDetectors");
  const char* detectorMask      = GetRunParameter("detectorMask");
  const char* lhcPeriod         = GetRunParameter("LHCperiod");
  
  TMap *mapDAQ = new TMap();
  
  if (timeStart) {
    Log(Form("Start time for run %d: %s",fRun, timeStart));
  } else {
    Log(Form("Start time not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fAliceStartTime"), new TObjString(timeStart));

  if (timeEnd) {
    Log(Form("End time for run %d: %s",fRun, timeEnd));
  } else {
    Log(Form("End time not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fAliceStopTime"), new TObjString(timeEnd));

  if (beamEnergy) {
    Log(Form("Beam energy for run %d: %s",fRun, beamEnergy));
  } else {
    Log(Form("Beam energy not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fAliceBeamEnergy"), new TObjString(beamEnergy));

  if (beamType) {
    Log(Form("Beam type for run %d: %s",fRun, beamType));
  } else {
    Log(Form("Beam type not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fAliceBeamType"), new TObjString(beamType));

  if (numberOfDetectors) {
    Log(Form("Number of active detectors for run %d: %s",fRun, numberOfDetectors));
  } else {
    Log(Form("Number of active detectors not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fNumberOfDetectors"), new TObjString(numberOfDetectors));

  if (detectorMask) {
    Log(Form("Detector mask for run %d: %s",fRun, detectorMask));
  } else {
    Log(Form("Detector mask not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fDetectorMask"), new TObjString(detectorMask));

  if (lhcPeriod) {
    Log(Form("LHC period (DAQ) for run %d: %s",fRun, lhcPeriod));
  } else {
    Log(Form("LHCperiod not put in logbook!"));
  }
  mapDAQ->Add(new TObjString("fLHCPeriod"), new TObjString(lhcPeriod));
  
  mapDAQ->Add(new TObjString("fRunType"), new TObjString(GetRunType()));
  Log( Form("Retrived %d parameters from logbook", mapDAQ->GetEntries() ) );

  return mapDAQ;
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDaqFxs()
{
  //======DAQ FXS======//

  TList* list = GetFileSources(kDAQ);  
  if (!list) {
    Log("No raw data tag list: connection problems with DAQ FXS logbook!");
    return 1;
  }

  if (list->GetEntries() == 0) {
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
        if (fileName.Length() > 0) {
           Log(Form("Adding file in the chain: %s",fileName.Data()));
           fRawTagChain->Add(fileName.Data());
           nFiles++;
        } else {
           Log(Form("Could not retrieve file with id %s from source %s: "
                    "connection problems with DAQ FXS!",
                     idStr->String().Data(), objStr->String().Data()));
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
  if( fRawTagChain->Merge(fRawDataFileName) < 1 ) {
    Log("Error merging raw data files!!!");
    return 3;
  }
  
  TString outputfile = Form("Run%d.Merged.RAW.tag.root", fRun);
  Bool_t result = StoreRunMetadataFile(fRawDataFileName.Data(),outputfile.Data());
  
  if (!result) {
    Log("Problem storing raw data tags in local file!!!");
  } else {
    Log("Raw data tags merged successfully!!");
  }
  
  delete iter;
  delete list;
  delete fRawTagChain; fRawTagChain=0;

  if (result == kFALSE) {
    return 4;
  }

  return 0;

}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDcsFxs()
{
  //======DCS FXS======//
  // Get the CTP run configuration
  // and scalers from DCS FXS

/*
  {

    // Get the CTP run configuration
    TList* list = GetFileSources(kDCS,"CTP_runconfig");  
    if (!list) {
      Log("No CTP runconfig file: connection problems with DCS FXS logbook!");
      return 1;
    }
  
    if (list->GetEntries() == 0) {
      Log("No CTP runconfig file to be processed!");
      return 1;
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
            return 1;
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
*/
  {
    // Get the CTP counters information
    TList* list = GetFileSources(kDCS,"CTP_xcounters");  
    if (!list) {
      Log("No CTP counters file: connection problems with DAQ FXS logbook!");
      return 1;
    }
  
    if (list->GetEntries() == 0) {
      Log("No CTP counters file to be processed!");
      return 1;
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
            return 1;
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
Int_t AliGRPPreprocessor::ProcessDcsDPs(TMap* valueMap, TMap* mapDCS)
{
  //Getting the DCS dps
  //===========//

  //DCS data points
  //===========//
  
  Int_t entries = 0;

  AliInfo(Form("==========LHCState==========="));
  TObjArray *aliasLHCState = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[0]);
  if(!aliasLHCState) {
    Log(Form("LHCState not found!!!"));
  } else {
    AliGRPDCS *dcs1 = new AliGRPDCS(aliasLHCState,fStartTime,fEndTime);
    TString sLHCState = dcs1->ProcessDCS(2);
    if (sLHCState) {
      for( Int_t i=0; i<20; i+=2 ) {
         if( sLHCState.CompareTo(fgkLHCState[i]) == 0 ) {
            sLHCState = fgkLHCState[i+1];
            break;
         }
      }
      Log(Form("<LHCState> for run %d: %s",fRun, sLHCState.Data()));
    } else {
      Log("LHCState not put in TMap!");
    }
    mapDCS->Add(new TObjString("fLHCState"),new TObjString(sLHCState));
    ++entries;
  }

  AliInfo(Form("==========L3Polarity==========="));
  TObjArray *aliasL3Polarity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[1]);
  if(!aliasL3Polarity) {
    Log(Form("L3Polarity not found!!!"));
  } else {
    AliGRPDCS *dcs6 = new AliGRPDCS(aliasL3Polarity,fStartTime,fEndTime);
    TString sL3Polarity = dcs6->ProcessDCS(1);  
    if (sL3Polarity) {
      Log(Form("<L3Polarity> for run %d: %s",fRun, sL3Polarity.Data()));
    } else {
      Log("L3Polarity not put in TMap!");
    }
    mapDCS->Add(new TObjString("fL3Polarity"),new TObjString(sL3Polarity));
    ++entries;
  }

  AliInfo(Form("==========DipolePolarity==========="));
  TObjArray *aliasDipolePolarity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[2]);
  if(!aliasDipolePolarity) {
    Log(Form("DipolePolarity not found!!!"));
  } else {
    AliGRPDCS *dcs8 = new AliGRPDCS(aliasDipolePolarity,fStartTime,fEndTime);
    TString sDipolePolarity = dcs8->ProcessDCS(1);  
    if (sDipolePolarity) {
      Log(Form("<DipolePolarity> for run %d: %s",fRun, sDipolePolarity.Data()));
    } else {
      Log("DipolePolarity not put in TMap!");
    }
    mapDCS->Add(new TObjString("fDipolePolarity"),new TObjString(sDipolePolarity));
    ++entries;
  }

  AliInfo(Form("==========LHCLuminosity==========="));
  TObjArray *aliasLHCLuminosity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[3]);
  if(!aliasLHCLuminosity) {
    Log(Form("LHCLuminosity not found!!!"));
  } else {
    AliGRPDCS *dcs3 = new AliGRPDCS(aliasLHCLuminosity,fStartTime,fEndTime);
    TString sMeanLHCLuminosity = dcs3->ProcessDCS(5);  
    if (sMeanLHCLuminosity) {
      Log(Form("<LHCLuminosity> for run %d: %s",fRun, sMeanLHCLuminosity.Data()));
    } else {
      Log("LHCLuminosity not put in TMap!");
    }
    mapDCS->Add(new TObjString("fLHCLuminosity"), new TObjString(sMeanLHCLuminosity));
    ++entries;
  }

  AliInfo(Form("==========BeamIntensity==========="));
  TObjArray *aliasBeamIntensity = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[4]);
  if(!aliasBeamIntensity) {
    Log(Form("BeamIntensity not found!!!"));
  } else {
    AliGRPDCS *dcs4 = new AliGRPDCS(aliasBeamIntensity,fStartTime,fEndTime);
    TString sMeanBeamIntensity = dcs4->ProcessDCS(5);  
    if (sMeanBeamIntensity) {
      Log(Form("<BeamIntensity> for run %d: %s",fRun, sMeanBeamIntensity.Data()));
    } else {
      Log("BeamIntensity not put in TMap!");
    }
    mapDCS->Add(new TObjString("fBeamIntensity"),new TObjString(sMeanBeamIntensity));
    ++entries;
  }

  AliInfo(Form("==========L3Current==========="));
  TObjArray *aliasL3Current = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[5]);
  if(!aliasL3Current) {
    Log(Form("L3Current not found!!!"));
  } else {
    AliGRPDCS *dcs5 = new AliGRPDCS(aliasL3Current,fStartTime,fEndTime);
    TString sMeanL3Current = dcs5->ProcessDCS(5);  
    if (sMeanL3Current) {
      Log(Form("<L3Current> for run %d: %s",fRun, sMeanL3Current.Data()));
    } else {
      Log("L3Current not put in TMap!");
    }
    mapDCS->Add(new TObjString("fL3Current"),new TObjString(sMeanL3Current));
    ++entries;
  }


  AliInfo(Form("==========DipoleCurrent==========="));
  TObjArray *aliasDipoleCurrent = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[6]);
  if(!aliasDipoleCurrent) {
    Log(Form("DipoleCurrent not found!!!"));
  }  else {
    AliGRPDCS *dcs7 = new AliGRPDCS(aliasDipoleCurrent,fStartTime,fEndTime);
    TString sMeanDipoleCurrent = dcs7->ProcessDCS(5);  
    if (sMeanDipoleCurrent) {
      Log(Form("<DipoleCurrent> for run %d: %s",fRun, sMeanDipoleCurrent.Data()));
    } else {
      Log("DipoleCurrent not put in TMap!");
    }
    mapDCS->Add(new TObjString("fDipoleCurrent"),new TObjString(sMeanDipoleCurrent));
    ++entries;
  }

  AliInfo(Form("==========CavernTemperature==========="));
  TObjArray *aliasCavernTemperature = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[7]);
  if(!aliasCavernTemperature) {
    Log(Form("CavernTemperature not found!!!"));
  }  else {
    AliGRPDCS *dcs9 = new AliGRPDCS(aliasCavernTemperature,fStartTime,fEndTime);
    TString sMeanCavernTemperature = dcs9->ProcessDCS(5);  
    if (sMeanCavernTemperature) {
      Log(Form("<CavernTemperature> for run %d: %s",fRun, sMeanCavernTemperature.Data()));
    } else {
      Log("CavernTemperature not put in TMap!");
    }
    mapDCS->Add(new TObjString("fCavernTemperature"),new TObjString(sMeanCavernTemperature));
    ++entries;
  }

  AliInfo(Form("==========CavernPressure==========="));
  TObjArray *aliasCavernPressure = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[8]);
  if(!aliasCavernPressure) {
    Log("CavernPressure not found!!!");
  } else {
    AliGRPDCS *dcs10 = new AliGRPDCS(aliasCavernPressure,fStartTime,fEndTime);
    TString sMeanCavernPressure = dcs10->ProcessDCS(5);  
    if (sMeanCavernPressure) {
      Log(Form("<CavernPressure> for run %d: %s",fRun, sMeanCavernPressure.Data()));
    } else {
      Log("CavernPressure not put in TMap!");
    }
    mapDCS->Add(new TObjString("fCavernPressure"),new TObjString(sMeanCavernPressure));
    ++entries;
  }


   // NEEDS TO BE REVISED, CONFIRMED
   AliInfo(Form("==========P2PressureMap==========="));
   AliDCSSensorArray *dcsSensorArray = GetPressureMap(valueMap);
   if( fPressure->NumFits()==0 ) {
     Log("Problem with the pressure sensor values!!!");
   } 
   else {
      AliDCSSensor* sensorP2 = dcsSensorArray->GetSensor(fgkDCSDataPoints[9]);
      if( sensorP2->GetFit() ) {
        Log(Form("<P2Pressure> for run %d: Sensor Fit found",fRun));
        mapDCS->Add( new TObjString("fP2Pressure"), sensorP2 );
        ++entries;
      } 
      else {
        Log(Form("ERROR Sensor Fit for %s not found: ", fgkDCSDataPoints[9] ));
      }
      
      /*
      AliDCSSensor* sensorMeyrin = dcsSensorArray->GetSensor(fgkDCSDataPoints[10]);
      if( sensorMeyrin->GetFit() ) {
        Log(Form("<MeyrinPressure> for run %d: Sensor Fit found",fRun));
        mapDCS->Add( new TObjString("fMeyrinPressure"), sensorMeyrin );
        ++entries;
      } else {
        Log(Form("ERROR Sensor Fit for %s not found: ", fgkDCSDataPoints[10] ));
      }
      */
   }

  return entries;
}

//_______________________________________________________________
AliDCSSensorArray *AliGRPPreprocessor::GetPressureMap(TMap* dcsAliasMap)
{
  // extract DCS pressure maps. Perform fits to save space
  
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    Double_t fitFraction = fPressure->NumFits()/fPressure->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Pressure values extracted, %d fits performed.", fPressure->NumFits()));
    } else { 
      AliInfo("Too few pressure maps fitted!!!");
    }
  } else {
    AliInfo("no atmospheric pressure map extracted!!!");
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
   sqlQuery.Form("SELECT DAQ_time_start, run_type, detectorMask FROM logbook WHERE run = %d", run);
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
   grpData.Add(new TObjString("DAQ_time_start"), new TObjString(row->GetField(0)));
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

   // Receive trigger information
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
