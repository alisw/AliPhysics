// $Id$

void AddTasksFlavourJet()
{
  const UInt_t uTriggerMask = AliVEvent::kMB;

  const TString sCutDzero = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root";
  const TString sCutDstar = "cutsHF/DStartoKpipiCuts.root"; 

  const TString sInputTrk  = "tracks";
  const TString sUsedTrks  = "PicoTracks";
  const TString sUsedClus  = "";
  const TString sRunPeriod = "-1:-1";

  const Int_t iJetAlgo = 1;
  const Int_t iJetType = 1;

  const Int_t    uBeamType   = 0; 
  const Double_t dJetPtCut   = 0.;
  const Double_t dJetAreaCut = 0.;

  const Int_t    nRadius = 3;
  const Double_t aRadius[] = {  0.2,   0.4,   0.6  };
  const TString  sRadius[] = { "R02", "R04", "R06" };
//=============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "Task manager to have an ESD or AOD input handler.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "This task requires an input event handler");
    return NULL;
  }
//=============================================================================

  UInt_t uAnaType = (((iJetType==0) ||     (iJetType==2)) ? 1 : 0);
  Int_t  iLeading =  ((iJetType==0) ? 3 : ((iJetType==1)  ? 0 : 1));

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bIsMC);
  taskRespPID->SelectCollisionCandidates(uTriggerMask);

  gROOT->LoadMacro("AddTaskSEDmesonsFilterCJ.C");
  AliAnalysisTaskSEDmesonsFilterCJ *taskDmesonsFilter = AddTaskSEDmesonsFilterCJ(sCutDzero,sCutDstar);
  taskDmesonsFilter->SelectCollisionCandidates(uTriggerMask);

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *taskSetupEMCal = AddTaskEmcalSetup();
  taskSetupEMCal->SetGeoPath("$ALICE_ROOT/OADB/EMCAL");
  taskSetupEMCal->SelectCollisionCandidates(uTriggerMask);

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *taskPicoTrack = AddTaskEmcalPicoTrackMaker(sUsedTrks.Data(),sInputTrk.Data(),sRunPeriod.Data());
  taskPicoTrack->SelectCollisionCandidates(uTriggerMask);

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("AddTaskFlavourJetCorrelations.C");

  for (Int_t i=0; i<nRadius; i++) {
    AliEmcalJetTask *taskFJ = AddTaskEmcalJet(sUsedTrks.Data(),sUsedClus.Data(),iJetAlgo,aRadius[i],iJetType);
    taskFJ->SelectCollisionCandidates(uTriggerMask);

    AliAnalysisTaskFlavourJetCorrelations *taskFlavourCJ = AddTaskFlavourJetCorrelations(sUsedTrks,sUsedClus,taskFJ->GetName());
    taskFlavourCJ->SetName(Form("AliAnalysisTaskFlavourJetCorrelations_%s",sRadius[i].Data()));
    taskFlavourCJ->SetForceBeamType(uBeamType);
    taskFlavourCJ->SetAnaType(uAnaType);
    taskFlavourCJ->SetLeadingHadronType(iLeading);
    taskFlavourCJ->SetJetRadius(aRadius[i]);
    taskFlavourCJ->SetJetPtCut(dJetPtCut);
    taskFlavourCJ->SetPercAreaCut(dJetAreaCut);
    taskFlavourCJ->SelectCollisionCandidates(uTriggerMask);
  }

  return;
}
