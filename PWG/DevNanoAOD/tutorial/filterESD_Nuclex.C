#ifndef __CLING__
#include "AliAnalysisManager.h"
#include "AliAnalysisNanoAODCuts.h"
#include "AliAnalysisTaskNanoAODFilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAODHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNanoAODTrack.h"
#include "AliAnalysisTaskNanoAODskimming.h"

#include <TChain.h>
#include <TInterpreter.h>
#endif

void filterESD_Nuclex()
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliESDInputHandler* iH = new AliESDInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  // Physics selection

  AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));


  // Multiplicity selection
  AliAnalysisTaskSE* multSelTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));


  // PID response
  AliAnalysisTaskSE* pidRespTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
  
  AliAnalysisTaskNanoAODskimming* myFilterTask = AliAnalysisTaskNanoAODskimming::AddTask();
  
  ///TODO: put these cuts in a separate AliAnalysisCut class, the AliAnalysisTaskNanoAODskimming should just apply some generic track cuts (maybe a list)
  AliNanoFilterPID* myFilterCuts = new AliNanoFilterPID;
  double nucleiTPCpt[4][2]{{10.6,1.4},{10.,1.8},{1.,10.},{1.,10.}};
  double nucleiTOFpt[4][2]{{10.6,10.},{10.,10.},{1.,10.},{1.,10.}};
  double nucleiTOFsigma[4]{10.,10.,-1.,-1.};
  AliESDtrackCuts* nucleiCuts[4]{AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false),AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false), AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false),AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false)};

  for (int iN{2}; iN < 4; ++iN) {
    myFilterCuts->TriggerOnSpecies(AliPID::EParticleType(AliPID::kDeuteron+iN), nucleiCuts[iN], 5., nucleiTPCpt[iN], nucleiTOFsigma[iN], nucleiTOFpt[iN]);
  }
  myFilterTask->AddEventCut(myFilterCuts);

  // OCDB
  AliAnalysisTaskSE* ocdbTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C(\"raw://\")"));
  
  // V0 finder
  AliAnalysisTaskWeakDecayVertexer* v0Finder = (AliAnalysisTaskWeakDecayVertexer*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
  v0Finder->SetUseImprovedFinding();
  // v0Finder->SetTrackTriggerCuts(nucleiCuts[2], AliPID::kHe3, 5.); ///TODO: fix this, it's buggy.
  
  // ESD filter
  AliAnalysisTaskSE* aodFilterTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kTRUE, kFALSE, kFALSE, kFALSE)"));
  aodFilterTask->SelectCollisionCandidates(AliVEvent::kAny);

  AliAnalysisTaskNanoAODFilter* nanoFilterTask = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/AddTaskNanoAODFilter.C(0, kFALSE)");
  nanoFilterTask->SelectCollisionCandidates(AliVEvent::kAny);  
  nanoFilterTask->AddSetter(new AliNanoAODSimpleSetter);
  
  // Track selection
  AliAnalysisNanoAODTrackCuts* trkCuts = new AliAnalysisNanoAODTrackCuts;
  
  // TODO adapt settings
  trkCuts->SetBitMask(1 << 4); // hybrid 2011
  trkCuts->SetMaxEta(0.9);
  trkCuts->SetMinPt(0.6);
  
  // Fields to store
  // event level
  // Note: vertices are kept by default
  // TODO add other fields
  nanoFilterTask->SetVarListHeader("CentrV0M,OfflineTrigger,MagField");
  // track level
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kPion);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kDeuteron);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kTriton);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kHe3);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kAlpha);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kPion);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kDeuteron);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kTriton);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kHe3);
  nanoFilterTask->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kAlpha);
  // TODO add other fields
  nanoFilterTask->SetVarListTrack("pt,theta,phi");

  nanoFilterTask->SetTrkCuts(trkCuts);

  AliAnalysisNanoAODEventCuts* evtCuts = new AliAnalysisNanoAODEventCuts;

  // V0s
  nanoFilterTask->SaveV0s(kTRUE, new AliAnalysisNanoAODV0Cuts);
  nanoFilterTask->SaveCascades(kTRUE);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 20);
}
