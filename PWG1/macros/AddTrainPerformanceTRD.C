// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   AddTrainPerformanceTRD.C(MC, friends, tasks)
//   tasks : "ALL" or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "MULT"  : TRD single track selection
//     "RES"  : TRD tracking Resolution
//     "CLRES": clusters Resolution
//     "CAL"  : TRD calibration
//     "ALGN" : TRD alignment
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
//     "DET"  : Basic TRD Detector checks
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//
// In compiled mode : 
// Don't forget to load first the libraries
// gSystem->Load("libMemStat.so")
// gSystem->Load("libMemStatGui.so")
// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libANALYSISalice.so")
// gSystem->Load("libPWG1.so");
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#if ! defined (__CINT__) || defined (__MAKECINT__)
//#ifndef __CINT__
#include <Riostream.h>

#include "TStopwatch.h"
#include "TMemStat.h"
#include "TMemStatViewerGUI.h"

#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TError.h"
#include "TChain.h"
#include "TGrid.h"
#include "TAlienCollection.h"
#include "TGridCollection.h"
#include "TGridResult.h"
#include "TGeoGlobalMagField.h"

#include "AliMagF.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

#include "TRD/AliTRDtrackerV1.h"
#include "TRD/AliTRDcalibDB.h"

#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/macros/AddTRDcheckESD.C"
#include "PWG1/TRD/macros/AddTRDinfoGen.C"
#include "PWG1/TRD/macros/AddTRDcheckDET.C"
#include "PWG1/TRD/macros/AddTRDefficiency.C"
#include "PWG1/TRD/macros/AddTRDresolution.C"
#include "PWG1/TRD/macros/AddTRDcheckPID.C"
#endif

#include "../TRD/macros/AliTRDperformanceTrain.h"


Bool_t AddTrainPerformanceTRD(Bool_t mc, Bool_t fr, Char_t *trd="ALL")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTrainPerformanceTRD", "AliAnalysisManager not set!");
    return kFALSE;
  }

  // TRD data containers
  AliAnalysisDataContainer *ci[] = {0x0, 0x0, 0x0};

  // initialize TRD settings
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDtrackerV1::SetNTimeBins(cal->GetNumberOfTimeBins());

  for(Int_t it=0; it<NTRDQATASKS; it++){
    if(gROOT->LoadMacro(Form("$ALICE_ROOT/PWG1/TRD/macros/Add%s.C+", TString(fgkTRDtaskClassName[it])(3,20).Data()))) {
      Error("AddTrainPerformanceTRD()", Form("Error loading %s task.", fgkTRDtaskClassName[it]));
      return kFALSE;
    } 

    switch(it){
    case kCheckESD:
      AddTRDcheckESD(mgr); break;
    case kInfoGen:
      AddTRDinfoGen(mgr, trd, 0x0, ci); break;
    case kCheckDET:
      AddTRDcheckDET(mgr, trd, ci); break;
    case kEfficiency:
      AddTRDefficiency(mgr, trd, ci); break;
    case kResolution:
      AddTRDresolution(mgr, trd, ci); break;
    case kCheckPID:
      AddTRDcheckPID(mgr, trd, ci); break;
    default:
      Warning("AddTrainPerformanceTRD()", Form("No performance task registered at slot %d.", it)); 
    }
  }
  return kTRUE;
}

