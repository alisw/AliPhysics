#include "../AddTaskNanoAODFilter.C"

void filterAOD_Tracks()
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  // Multiplicity selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AddTaskMultSelection();
  
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) AddTaskNanoAODFilter(0, kFALSE);
  AliNanoAODSimpleSetter* setter = new AliNanoAODSimpleSetter;
  task->SetSetter(setter);
  
  // Event selection
  // filter bit
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  AliAnalysisNanoAODEventCuts* evtCuts = new AliAnalysisNanoAODEventCuts;
  evtCuts->SetVertexRange(8);
  evtCuts->SetCutPileUpMV(kTRUE);

  // Track selection
  AliAnalysisNanoAODTrackCuts* trkCuts = new AliAnalysisNanoAODTrackCuts;
  trkCuts->SetBitMask((1 << 8) | (1 << 9)); // hybrid 2011
  trkCuts->SetMaxEta(0.9);
  trkCuts->SetMinPt(1.0);
  
  // Fields to store
  // event level
  // Note: vertices are kept by default
  task->SetVarListHead("OfflineTrigger,MagField,MultSelection.RefMult08");
  // track level
  task->SetVarList("pt,theta,phi");

  task->SetTrkCuts(trkCuts);
  task->SetEvtCuts(evtCuts);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("aodTree");
  chain->Add("AliAOD.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 100);
}
