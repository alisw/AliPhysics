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
class AliESETrkCut;
class AliESEEvtCut;

AliESETrkCut * TrkCuts() {

  AliESETrkCut * trk = new AliESETrkCut;

  AliSpectraAODTrackCuts  * trcuts = new AliSpectraAODTrackCuts("TrackCuts");  
  trcuts->SetDCA(100000);
  trcuts->SetTrackBits(1024);
  trcuts->SetPt(15);
  trcuts->SetPtTOFMatching(0.6);   
  trcuts->SetEta(-0.8,0.8);
  trcuts->SetMinTPCcls(70);
  trcuts->PrintCuts();

  trk->SetTrackCuts(trcuts);
  trk->Init();

  return trk;

}

AliESEEvtCut * EvtCuts(Int_t mc) {

  AliESEEvtCut * evt = new AliESEEvtCut;

  AliSpectraAODEventCuts * evcuts = new AliSpectraAODEventCuts("EventCuts");
  evcuts->SetQVectorCut(0,100);
  evcuts->SetCentralityCutMax(100);  
  evcuts->SetCentralityCutMin(0);
  if(mc>0)evcuts->SetIsMC(kTRUE);
  TFile * fCalib = new TFile("./calibV0New.root");
  evcuts->SetCalibFile(fCalib);
  evcuts->SetIsLHC10h(kTRUE);
  evcuts->PrintCuts();
  //  evcuts->SetEventSelectionBit(AliVEvent::kAny);

  evt->SetEventCuts(evcuts);
  evt->Init();

  return evt;

}

//______________________________________________________________________________
void runLocalESE(
		 const int iMCtruth = 2, 
		 const char * addTaskString = ".x AddTaskNanoAODFilter.C(%d,1)" // 
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
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(iMCtruth);
  taskPID->SetUseTPCEtaCorrection(kTRUE); 
  taskPID->SetUserDataRecoPass(2);
  // create task
  cout << "Macro: "<< addTaskString << " " << Form(addTaskString, iMCtruth) << endl;

  AliAnalysisTaskNanoAODFilter * task = (AliAnalysisTaskNanoAODFilter*) gROOT->ProcessLine(Form(addTaskString, iMCtruth));

  // Set Track event and vertex cuts here!
  //  task->SetVarList("pt,theta,phi,cstNSigmaTPCPi,cstNSigmaTPCKa,cstNSigmaTPCPr,cstNSigmaTOFPi,cstNSigmaTOFKa,cstNSigmaTOFPr,cstBayesTPCPi,cstBayesTPCKa,cstBayesTPCPr,cstBayesTOFPi,cstBayesTOFKa,cstBayesTOFPr");
  task->SetVarList("pt,theta,phi,cstNSigmaTPCPi,cstNSigmaTPCKa,cstNSigmaTPCPr,cstNSigmaTOFPi,cstNSigmaTOFKa,cstNSigmaTOFPr");
  task->SetVarListHead("cstCentr,cstQVec");
  AliESETrkCut * trkCuts = TrkCuts();
  AliESEEvtCut * evtCuts = EvtCuts(iMCtruth);
  evtCuts->SetTrackCuts(trkCuts->GetTrackCuts());
  AliAnalysisESESetter * setter  = new AliAnalysisESESetter;
  setter->SetEventCuts(evtCuts->GetEventCuts());

  task->SetTrkCuts(trkCuts);
  task->SetEvtCuts(evtCuts);
  task->SetSetter(setter);

  //task->SelectCollisionCandidates(AliVEvent::kMB);// FIXME
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
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
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
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");

  //  gSystem->Load("libNanoAOD.so");
  gSystem->Load("libPWGLFspectra");
  gSystem->Load("libPWGDevNanoAOD");

}
