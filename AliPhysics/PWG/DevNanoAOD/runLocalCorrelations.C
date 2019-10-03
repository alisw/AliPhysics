#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "AliAnalysisGrid.h"
#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisGrid.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "TRegexp.h"
#include "AliTriggerAnalysis.h"
#include "TChain.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskNanoAODFilter.h"
#include "AliESEHelpers.h"

#endif
void LoadLibs();

class AliAnalysisGrid;



//______________________________________________________________________________
void runLocalCorrelations(
			  const int iMCtruth = 2, 
			  const char * addTaskString = ".x AddTaskNanoAODFilter.C(%d,0)" // 
			  )
{
  LoadLibs();
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
//   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
//   AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(iMCtruth);
//   taskPID->SetUseTPCEtaCorrection(kTRUE); 

  // create task
  cout << "Macro: "<< addTaskString << " " << Form(addTaskString, iMCtruth) << endl;

  AliAnalysisTaskNanoAODFilter * task = (AliAnalysisTaskNanoAODFilter*) gROOT->ProcessLine(Form(addTaskString, iMCtruth));

  // Set Track event and vertex cuts here!
  AliAnalysisNanoAODTrackCuts* trk = new AliAnalysisNanoAODTrackCuts;
  trk->SetBitMask((1 << 4) | (1 << 8)); // hybrid 2010
  trk->SetMaxEta(0.9);
  trk->SetMinPt(0.2);

  AliAnalysisNanoAODEventCuts* evt = new AliAnalysisNanoAODEventCuts;
  evt->SetVertexRange(10);

  task->SetTrkCuts(trk);
  task->SetEvtCuts(evt);
  task->AddSetter(new AliNanoAODSimpleSetter);
  task->SetVarListTrack("pt,theta,phi");
  task->SetVarListHeader("cstCentr,cstMagField");

  task->SelectCollisionCandidates(AliVEvent::kMB);

  // enable debug printouts
  mgr->SetDebugLevel(10);
  //    mgr->SetNSysInfo(100);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // start analysis
  // Always read the same file:
  TChain * chain = new TChain("aodTree");
  chain->Add("./AliAOD.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain,123456789);

}

//______________________________________________________________________________

void LoadLibs() {
  gSystem->Load("libCore");  
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  //  return;
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");

  //  gSystem->Load("libNanoAOD");
  gSystem->Load("libPWGLFspectra");
  gSystem->Load("libPWGDevNanoAOD");

}
