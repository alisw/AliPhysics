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
// gSystem->Load("libPWGPP.so");
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#if ! defined (__CINT__) || defined (__MAKECINT__)
//#ifndef __CINT__
#include <Riostream.h>

#include <TStopwatch.h>
#include <TMemStat.h>
//#include <TMemStatViewerGUI.h>

#include <TROOT.h>
#include <TClass.h>
#include <TSystem.h>
#include <TString.h>
#include <TError.h>
#include <TChain.h>
#include <TGrid.h>
#include <TGridCollection.h>
#include <TGridResult.h>
#include <TGeoGlobalMagField.h>

#include <AliMagF.h>
#include <AliTracker.h>
#include <AliLog.h>
#include <AliCDBManager.h>
#include <AliGRPManager.h>
#include <AliGeomManager.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMCEventHandler.h>
#include <AliESDInputHandler.h>

#include <AliTRDtrackerV1.h>
#include <AliTRDcalibDB.h>

#include <AliTRDpwgppHelper.h>
#include "PWGPP/TRD/macros/AddTRDcheckESD.C"
#include "PWGPP/TRD/macros/AddTRDinfoGen.C"
#include "PWGPP/TRD/macros/AddTRDcheckDET.C"
#include "PWGPP/TRD/macros/AddTRDefficiency.C"
#include "PWGPP/TRD/macros/AddTRDresolution.C"
#include "PWGPP/TRD/macros/AddTRDcheckPID.C"
#include "PWGPP/TRD/macros/AddTRDcheckTRK.C"
#include "PWGPP/TRD/macros/AddTRDv0Monitor.C"
#endif

TString opt("");
const Char_t* Translate(Bool_t doCheckESD=kTRUE, Bool_t doCheckDET=kTRUE, Bool_t doEffic=kTRUE, Bool_t doResolution=kTRUE, Bool_t doCheckPID=kTRUE, Bool_t doV0Monitor=kTRUE);
Bool_t AddTrainPerformanceTRD(Char_t *trd="ALL", const Char_t *addMacroPath = "$ALICE_PHYSICS/PWGPP/TRD/macros")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTrainPerformanceTRD", "AliAnalysisManager not set!");
    return kFALSE;
  }

  // TRD data containers
  AliAnalysisDataContainer *ci[AliTRDpwgppHelper::kNOutSlots];
  AliAnalysisDataContainer *ce[6];

  Info("AddTrainPerformanceTRD", Form("Add Macros taken from %s", addMacroPath));
  Info("AddTrainPerformanceTRD", Form("TRD wagons \"%s\"", trd));
  Int_t bitmap = AliTRDpwgppHelper::ParseOptions(trd);
  for(Int_t it=0; it<AliTRDpwgppHelper::kNTRDQATASKS; it++){
    if(gROOT->LoadMacro(Form("%s/Add%s.C", addMacroPath, TString(AliTRDpwgppHelper::TaskClassName(it))(3,20).Data()))) {
      Error("AddTrainPerformanceTRD()", Form("Error loading %s task.", AliTRDpwgppHelper::TaskClassName(it)));
      return kFALSE;
    } 
    if(!AliTRDpwgppHelper::DoTask(it, bitmap)) continue;

    switch(it){
    case AliTRDpwgppHelper::kCheckESD:
      AddTRDcheckESD(mgr); break;
    case AliTRDpwgppHelper::kInfoGen:
      AddTRDinfoGen(mgr, 0, NULL, ci); break;
    case AliTRDpwgppHelper::kCheckDET:
      // map slots
      ce[0]=ci[AliTRDpwgppHelper::kTracksBarrel];
      ce[1]=ci[AliTRDpwgppHelper::kTracksSA];
      ce[2]=ci[AliTRDpwgppHelper::kTracksKink];
      ce[3]=ci[AliTRDpwgppHelper::kEventInfo];
      ce[4]=ci[AliTRDpwgppHelper::kTracklets];
      ce[5]=ci[AliTRDpwgppHelper::kClusters];
      AddTRDcheckDET(mgr, bitmap, ce);
      break;
    case AliTRDpwgppHelper::kEfficiency:
      // map slots
      ce[0]=ci[AliTRDpwgppHelper::kTracksBarrel];
      ce[1]=ci[AliTRDpwgppHelper::kTracksITS];
      ce[2]=ci[AliTRDpwgppHelper::kTracksKink];
      ce[3]=ci[AliTRDpwgppHelper::kEventInfo];
      ce[4]=ci[AliTRDpwgppHelper::kTracklets];
      ce[5]=ci[AliTRDpwgppHelper::kClusters];
      AddTRDefficiency(mgr, bitmap, ce);
      break;
    case AliTRDpwgppHelper::kResolution:
      // map slots
      ce[0]=ci[AliTRDpwgppHelper::kTracksBarrel];
      ce[1]=ci[AliTRDpwgppHelper::kTracksITS];
      ce[2]=ci[AliTRDpwgppHelper::kTracksKink];
      ce[3]=ci[AliTRDpwgppHelper::kEventInfo];
      ce[4]=ci[AliTRDpwgppHelper::kTracklets];
      ce[5]=ci[AliTRDpwgppHelper::kClusters];
      AddTRDresolution(mgr, bitmap, ce); 
      break;
    case AliTRDpwgppHelper::kCheckPID:
      // map slots
      ce[0]=ci[AliTRDpwgppHelper::kTracksBarrel];
      ce[1]=ci[AliTRDpwgppHelper::kEventInfo];
      ce[2]=ci[AliTRDpwgppHelper::kTracklets];
      ce[3]=ci[AliTRDpwgppHelper::kV0List];
      AddTRDcheckPID(mgr, bitmap, ce, &ce[4]);
      break;
    case AliTRDpwgppHelper::kCheckTRK:
      // map slots
      //ce[0]=ci[AliTRDpwgppHelper::kTracksBarrel];
      //ce[1]=ci[AliTRDpwgppHelper::kEventInfo];
      //ce[2]=ci[AliTRDpwgppHelper::kTracklets];
      //ce[3]=ci[AliTRDpwgppHelper::kClusters];
      AddTRDcheckTRK(mgr, 0, ce);
      break;
    case AliTRDpwgppHelper::kV0Monitor:
      // slots already mapped by checkPID
      ce[2] = ce[3];
      ce[3] = ce[4];
      AddTRDv0Monitor(mgr, 0, ce);
      break;
    default:
      Warning("AddTrainPerformanceTRD()", Form("No performance task registered at slot %d.", it)); 
    }
  }
  return kTRUE;
}

const Char_t* Translate(Bool_t doCheckESD, Bool_t doCheckDET, Bool_t doEffic, Bool_t doResolution, Bool_t doCheckPID, Bool_t doCheckV0)
{
  opt.Clear();
  if( doCheckESD==kTRUE &&
      doCheckDET==kTRUE &&
      doEffic==kTRUE &&
      doResolution==kTRUE &&
      doCheckPID==kTRUE &&
      doCheckV0==kTRUE
  ){
    opt="ALL";
  } else {
    Bool_t kINDENT(kFALSE);
    if(doCheckESD){ 
      opt.Append("ESD");
      kINDENT=kTRUE;
    }
    if(doCheckDET){ 
      if(kINDENT) opt.Append(" ");
      opt.Append("DET"); 
      kINDENT = kTRUE;
    }
    if(doEffic){ 
      if(kINDENT) opt.Append(" ");
      opt.Append("EFF");
      kINDENT=kTRUE;
    }
    if(doResolution){ 
      if(kINDENT) opt.Append(" ");
      opt.Append("RES");
      kINDENT=kTRUE;
    }
    if(doCheckPID){ 
      if(kINDENT) opt.Append(" ");
      opt.Append("PID");
      kINDENT=kTRUE;
    }
    if(doCheckV0){ 
      if(kINDENT) opt.Append(" ");
      opt.Append("V0");
      kINDENT=kTRUE;
    }
  }

  return (const Char_t*)opt.Data();
}


