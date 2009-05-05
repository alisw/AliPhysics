///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance to be run on QA train
// 24.04.2009 -  J.Otwinowski@gsi.de
///////////////////////////////////////////////////////////////////////////////
#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <Riostream.h>

#include "TROOT.h"
#include "TClass.h"
//#include "TSystem.h"
#include "TError.h"

#include "AliLog.h"
#include "AliAnalysisManager.h"
//#include "AliAnalysisDataContainer.h"
//#include "AliMCEventHandler.h"
//#include "AliESDInputHandler.h"

#include "PWG1/AliPerformanceTask.h"
#include "PWG1/AliPerformanceObject.h"
#include "PWG1/AliPerformanceEff.h"
#include "PWG1/AliPerformanceDEdx.h"
#include "PWG1/AliPerformanceTPC.h"
#include "PWG1/AliPerformanceDCA.h"
#include "PWG1/AliPerformanceRes.h"
#include "PWG1/AliMCInfoCuts.h"
#include "PWG1/AliRecInfoCuts.h"
#include "PWG1/macros/AddPerformanceTask.h"
#endif

Int_t ParseTPC(Char_t *tpc);

//____________________________________________
void AddPerformanceTask(AliAnalysisManager *mgr=0, Char_t *tpc="ALL")
{
  if(!mgr) { 
    Error("AddPerformanceTask","AliAnalysisManager not set!");
    return;
  }
  // parse options
  Int_t fSteerTPC = ParseTPC(tpc);

  // Add task
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");
  if (!task) return;
  if ( mgr->GetMCtruthEventHandler() ) task->SetUseMCInfo(kTRUE);
  mgr->AddTask(task);

  // Create TPC-ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  }
  // Create TPC-MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  }

  //
  // Create performance objects for TPC and set cuts 
  //

  // TPC at DCA
  TClass ctask; AliPerformanceObject *perf = 0x0;
  for(Int_t icomp=0; icomp<NTPCTASKS; icomp++) {
    if(!(TSTTPCBIT(fSteerTPC, icomp))) continue;
    TString  s(fgkTPCtaskClassName[icomp]);
    if(s.CompareTo("AliPerformanceEff") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceEff(fgkTPCtaskClassName[icomp],fgkTPCtaskClassName[icomp],kTPCMode,fHpt)));
    } else if (s.CompareTo("AliPerformanceRes") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceRes(fgkTPCtaskClassName[icomp],fgkTPCtaskClassName[icomp],kTPCMode,fHpt)));
    } else if (s.CompareTo("AliPerformanceTPC") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceTPC(fgkTPCtaskClassName[icomp],fgkTPCtaskClassName[icomp],kTPCMode,fHpt)));
    } else if (s.CompareTo("AliPerformanceDCA") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceDCA(fgkTPCtaskClassName[icomp],fgkTPCtaskClassName[icomp],kTPCMode,fHpt)));
    } else {
      Warning("No TPC at DCA mode for ",fgkTPCtaskClassName[icomp]);
    }
    perf->SetAliMCInfoCuts(pMCInfoCuts);
    perf->SetAliRecInfoCuts(pRecInfoCuts);

    /*
    new(&ctask) TClass(fgkTPCtaskClassName[icomp]);
    task->AddPerformanceObject((perf = (AliPerformanceObject*)ctask.New()));
    perf->SetAnalysisMode(kTPCMode);
    perf->SetHptGenerator(fHpt);
    */
  }

  // TPC at inner wall
  for(Int_t icomp=0; icomp<NTPCTASKS; icomp++) {
    if(!(TSTTPCBIT(fSteerTPC, icomp))) continue;
    TString  s(fgkTPCtaskClassName[icomp]);
    if(s.CompareTo("AliPerformanceRes") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceRes("AliPerformanceResTPCInner","AliPerformanceResTPCInner",kTPCInnerMode,fHpt)));
    } else if (s.CompareTo("AliPerformanceDEdx") == 0) {
      task->AddPerformanceObject((perf = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInnerMode,fHpt)));
    } else {
      Warning("No TPC at inner wall mode for ",fgkTPCtaskClassName[icomp]);
    } 
    perf->SetAliMCInfoCuts(pMCInfoCuts);
    perf->SetAliRecInfoCuts(pRecInfoCuts);
  }

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  // Create containers for output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(0);
}

//____________________________________________
Int_t ParseTPC(Char_t *tpc)
{
  Printf("---------------------------------------");
  Printf("TPC Performance Task Configuration Options:");
  Printf("---------------------------------------");
  for(int i=0; i<NTPCTASKS+1; i++) {
    Printf("%s",fgkTPCtaskOpt[i]);
  } 
  Printf("%s","HPT");
  Printf("---------------------------------------");
  Printf("Used Options:");
  Printf("---------------------------------------");

  Int_t fSteerTask = 0;
  TObjArray *tasksArray = TString(tpc).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0) {
      Printf("%s","ALL");
      for(Int_t itask = 0; itask < NTPCPERFORMANCE; itask++) SETTPCBIT(fSteerTask, itask);
      continue;
    }
    else if(s.CompareTo("HPT") == 0) {
      fHpt = kTRUE;
      Printf("%s","HPT");
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 0; itask < NTPCTASKS; itask++) {
        if(s.CompareTo(fgkTPCtaskOpt[itask]) != 0) continue;
        SETTPCBIT(fSteerTask, itask); //SETTPCBIT(fSteerTask, 0);
        foundOpt = kTRUE;
        Printf("%s", fgkTPCtaskOpt[itask]);
        break;
      }
      if(!foundOpt) Info("run.C", Form("TPC task %s not implemented (yet).", s.Data()));
    }
  }
  Printf("---------------------------------------");
  return fSteerTask;
}
