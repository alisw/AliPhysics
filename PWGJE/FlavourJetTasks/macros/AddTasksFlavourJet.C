// $Id$

void AddTasksFlavourJet(const Int_t iCandType = 1 /*0 = D0, 1=Dstar...*/,
		        const TString sCutFile = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root",
                        const Double_t dJetPtCut   = 1.,
                        const Double_t dJetAreaCut = 0.,
                        const Int_t iAccCut = 1,
                        const TString sRunPeriod = "LHC10b",
                        const Int_t    uBeamType = 0,
                        const UInt_t uTriggerMask = AliVEvent::kMB, /*for jets; the D mesons trigger is defined in the cut object*/
                        const Bool_t bIsMC = kFALSE,
                        TString sText=""/*completes the name of the candidate task lists*/
)
{
  const TString sInputTrk  = "tracks";
  const TString sUsedTrks  = "PicoTracks";
  const TString sUsedClus  = "";

  const Int_t iJetAlgo = 1;
  const Int_t iJetType = 1;

  const Int_t    nRadius = 3;
  const Double_t aRadius[] = {  0.2,   0.4,   0.6  };
  const TString  sRadius[] = { "R02", "R04", "R06" };
//=============================================================================

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "No analysis manager to connect to.");
    return;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "Task manager to have an ESD or AOD input handler.");
    return;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "This task requires an input event handler");
    return;
  }
//=============================================================================

  UInt_t uAnaType = (((iJetType==0) ||     (iJetType==2)) ? 1 : 0);
  Int_t  iLeading =  ((iJetType==0) ? 3 : ((iJetType==1)  ? 0 : 1));

  //D mesons -- PID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bIsMC);

  // -- D meson selection
  gROOT->LoadMacro("AliAnalysisTaskSEDmesonsForJetCorrelations.cxx++g");
  gROOT->LoadMacro("AddTaskDmesonsForJetCorrelations.C");
  AliAnalysisTaskSEDmesonsForJetCorrelations *taskDmesonsFilter = AddTaskDmesonsForJetCorrelations(iCandType,sCutFile,bIsMC,sText);

  // EMCal framework
  // -- Physics selection task
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, uTriggerMask, 5, 5, 10, kTRUE, -1, -1, -1, -1);

  if (!physSelTask) {
    cout << "no physSelTask"; 
    return; 
  }

  // -- 
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *taskSetupEMCal = AddTaskEmcalSetup();
  taskSetupEMCal->SetGeoPath("$ALICE_ROOT/OADB/EMCAL");
  taskSetupEMCal->SelectCollisionCandidates(uTriggerMask);

  // Jet preparation
  gROOT->LoadMacro("/data/Work/jets/testEMCalJetFramework/code/v4/AddTaskJetPreparation.C");
  AddTaskJetPreparation(type,bIsMC,sRunPeriod);

/*gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *taskPicoTrack = AddTaskEmcalPicoTrackMaker(sUsedTrks.Data(),sInputTrk.Data(),sRunPeriod.Data());
  taskPicoTrack->SelectCollisionCandidates(uTriggerMask);*/

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
  gROOT->LoadMacro("AliAnalysisTaskRecoJetCorrelations.cxx++g");
  gROOT->LoadMacro("AddTaskRecoJetCorrelations.C");

  for (Int_t i=0; i<nRadius; i++) {
    AliEmcalJetTask *taskFJ = AddTaskEmcalJet(sUsedTrks.Data(),sUsedClus.Data(),iJetAlgo,aRadius[i],iJetType);
    taskFJ->SelectCollisionCandidates(uTriggerMask);

    AliAnalysisTaskFlavourJetCorrelations *taskDmesonCJ = AddTaskFlavourJetCorrelations(iCandType,
                                                                                        sCutFile,
                                                                                        bIsMC,
                                                                                        taskFJ->GetName(),
                                                                                        Form("JetR%s",sRadius[i].Data()),
                                                                                        iLeading,
                                                                                        aRadius[i],
                                                                                        dJetPtCut,
                                                                                        iAccCut,
                                                                                        dJetAreaCut);

    taskDmesonCJ->SetName(Form("AliAnalysisTaskSEEmcalJetDmesonsCJ_%s",sRadius[i].Data()));
    taskDmesonCJ->SetForceBeamType(uBeamType);
    taskDmesonCJ->SetAnaType(uAnaType);
    taskDmesonCJ->SetLeadingHadronType(iLeading);
//  taskDmesonCJ->SelectCollisionCandidates(uTriggerMask);
  }

  return;
}
