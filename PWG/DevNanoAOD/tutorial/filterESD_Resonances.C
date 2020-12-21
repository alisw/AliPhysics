#ifndef __CLING__
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisNanoAODCuts.h"
#include "AliAnalysisTaskNanoAODFilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNanoAODTrack.h"
#include "AliNanoFilterPID.h"
#include "AliNanoSkimmingPID.h"
#include "AliAnalysisTaskNanoAODskimming.h"
#include "AliAnalysisTaskXi1530.h"

#include <TChain.h>
#include <TInterpreter.h>
#endif

void filterESD_Resonances()
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliESDInputHandler* iH = new AliESDInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.Nano.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  // Physics selection
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  
  // Multiplicity selection
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  
  // PID response
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");

  // OCDB
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C(\"raw://\")"); //for ESD: ->AOD->nanoAOD (2 steps)
  
  // // V0 finder
  AliAnalysisTaskWeakDecayVertexer* v0Finder = (AliAnalysisTaskWeakDecayVertexer*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
  v0Finder->SetUseImprovedFinding();
  v0Finder->SetForceResetV0s(true);
  v0Finder->SetForceResetCascades(true);


  //V0-Related topological selections
  v0Finder->SetV0VertexerDCAFirstToPV(0.05);
  v0Finder->SetV0VertexerDCASecondtoPV(0.05);
  v0Finder->SetV0VertexerDCAV0Daughters(2.0);
  v0Finder->SetV0VertexerCosinePA(0.95);
  v0Finder->SetV0VertexerMinRadius(0.2);
  v0Finder->SetV0VertexerMaxRadius(500);

  // Cascade-Related topological selections
  v0Finder->SetCascVertexerMinV0ImpactParameter(0.01);
  v0Finder->SetCascVertexerV0MassWindow(0.01);
  v0Finder->SetCascVertexerDCABachToPV(0.01);
  v0Finder->SetCascVertexerDCACascadeDaughters(2.0);
  v0Finder->SetCascVertexerCascadeMinRadius(0.2);
  v0Finder->SetCascVertexerCascadeMaxRadius(500);
  v0Finder->SetCascVertexerCascadeCosinePA(0.96);

  // Test1 track selection
  v0Finder->SetExtraCleanup(kFALSE);
  // Test2 pre-selection in dE/dx
  v0Finder->SetPreselectDedx(kFALSE);
  v0Finder->SetPreselectDedxLambda(kFALSE);
  //v0Finder-> SetUseMonteCarloAssociation(kFALSE);

  // ESD filter
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kFALSE, kFALSE, kFALSE, kFALSE)");
  
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/macros/AddTaskNanoAODFilter.C(0, kFALSE)");
  task->AddSetter(new AliNanoAODSimpleSetter);
  task->fUseAliEventCuts = true;
  task->fEventCuts.fCentralityFramework = 1;
  task->fEventCuts.SelectOnlyInelGt0(false);
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  // NOTE filter bit set in AliEventCuts automatically

  // Track selection
  AliAnalysisNanoAODTrackCuts* trkCuts = new AliAnalysisNanoAODTrackCuts;
  trkCuts->SetBitMask((1 << 5)); // ITS+TPC standard cuts
  trkCuts->SetMaxEta(0.8);
  trkCuts->SetMinPt(0.15);
  
  // Fields to store
  // event level
  // Note: vertices are kept by default
  task->SetVarListHeader("OfflineTrigger,MagField,CentrV0M,CentrTRK,RunNumber,T0Spread,NumberOfESDTracks,MultSelection.V0M.Value,MultSelection.SPDTracklets.Value,cstINELgt0");
  task->AddSetter(new AliNanoAODSimpleSetter);
  task->AddSetter(new AliNanoAODINELgt0setter);
  // track level
  task->SetVarListTrack("pt,theta,phi,TPCmomentum,TOFsignal,TPCsignal,integratedLength,DCA,posDCAz,ID,FilterMap,covmat,posx,posy,posz");
  task->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kPion);
  task->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kKaon);
  task->AddPIDField(AliNanoAODTrack::kSigmaTPC, AliPID::kProton);
  task->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kPion);
  task->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kKaon);
  task->AddPIDField(AliNanoAODTrack::kSigmaTOF, AliPID::kProton);
  task->SetTrkCuts(trkCuts);

  // V0s Cuts
  AliAnalysisNanoAODV0Cuts* v0cut = new AliAnalysisNanoAODV0Cuts;
  v0cut->SetTransverseRadius(0.2,500); // 500 is large limit
  v0cut->Setv0EtaMax(1.0);
  v0cut->SetCPAMin(0.95);
  v0cut->SetMaxDCADaughtersToV0Vtx(2.0);
  v0cut->SetMinDCADaughtersToPrimVtx(0.05);
  //v0cut->SetLambdaDaugnSigTPCMax(5); // Disable for the non-PID study

  // V0s
  task->SaveV0s(kTRUE, v0cut);
  task->SaveCascades(kTRUE);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 1000000);
  // gSystem->Exec("mv AnalysisResults.root AnalysisResultsESD.root");

}
