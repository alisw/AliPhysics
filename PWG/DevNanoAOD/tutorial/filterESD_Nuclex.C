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
#include "AliNanoFilterPID.h"
#include "AliNanoSkimmingPID.h"
#include "AliNanoSkimmingV0s.h"
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

  // V0 Hypertriton vertexer
  AliAnalysisTaskHypV0s* hypV0sTask = AliAnalysisTaskHypV0s::AddTask();
  hypV0sTask->SetCVMFSPath("/cvmfs/alice.cern.ch/data/analysis/2020/vAN-20200217/PWGLF/NUCLEX/HypertritonAnalysis/Cuts/splines.root");
  
  AliAnalysisTaskNanoAODskimming* mySkimmingTask = AliAnalysisTaskNanoAODskimming::AddTask();
  AliNanoSkimmingPID* mySkimmingCuts = new AliNanoSkimmingPID;
  double nucleiTPCpt[4][2]{{0.6,1.4},{1.,1.8},{1.,10.},{1.,10.}};
  double nucleiTOFpt[4][2]{{1.4,10.},{1.8,10.},{1.,10.},{1.,10.}};
  double nucleiTOFsigma[4]{10.,10.,-1.,-1.};
  AliESDtrackCuts* nucleiCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false);
  nucleiCuts->SetEtaRange(-0.8,0.8);
  for (int iN{0}; iN < 4; ++iN)
    mySkimmingCuts->fTrackFilter.TriggerOnSpecies(AliPID::EParticleType(AliPID::kDeuteron+iN), nucleiCuts, 1 << 4,  5., nucleiTPCpt[iN], nucleiTOFsigma[iN], nucleiTOFpt[iN]);
  mySkimmingTask->AddEventCut(mySkimmingCuts);
  mySkimmingTask->AddEventCut(new AliNanoSkimmingV0s);

  // OCDB
  AliAnalysisTaskSE* ocdbTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C(\"raw://\")"));
  
  // ESD filter
  AliAnalysisTaskSE* aodFilterTask = reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kTRUE, kFALSE, kFALSE, kFALSE)"));
  aodFilterTask->SelectCollisionCandidates(AliVEvent::kAny);

  AliAnalysisTaskNanoAODFilter* nanoFilterTask = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/macros/AddTaskNanoAODFilter.C(0, kFALSE)");
  nanoFilterTask->SelectCollisionCandidates(AliVEvent::kAny);  
  nanoFilterTask->AddSetter(new AliNanoAODSimpleSetter);
  
  AliNanoFilterPID* myFilterCuts = new AliNanoFilterPID;
  for (int iN{0}; iN < 4; ++iN)
    myFilterCuts->TriggerOnSpecies(AliPID::EParticleType(AliPID::kDeuteron+iN), nullptr, 1 << 4,  5., nucleiTPCpt[iN], nucleiTOFsigma[iN], nucleiTOFpt[iN]);

  nanoFilterTask->SelectCollisionCandidates(AliVEvent::kAny);

  // NOTE no event cuts. Those are in the Skimming task! nanoFilterTask->AddEvtCuts(evtCuts);

  nanoFilterTask->SetTrkCuts(myFilterCuts);
  nanoFilterTask->AddSetter(new AliNanoAODSimpleSetter);

  nanoFilterTask->SetVarListTrack("pt,theta,phi,TPCsignalN,TPCncls,TPCnclsF,TPCNCrossedRows,chi2perNDF,TRDntrackletsPID,TPCmomentum,TOFsignal,TPCsignal,integratedLength,DCA,posDCAz");
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

  nanoFilterTask->SetVarListHeader("OfflineTrigger,MagField,CentrV0M,RunNumber,T0Spread,NumberOfESDTracks");

  nanoFilterTask->SaveV0s(kTRUE);

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
