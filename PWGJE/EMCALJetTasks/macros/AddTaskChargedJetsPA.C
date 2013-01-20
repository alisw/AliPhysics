AliAnalysisTaskChargedJetsPA* AddTaskChargedJetsPA(
  Double_t            jetRadius               = 0.4,
  Int_t               trigger                 = AliVEvent::kINT7,
  Bool_t              isMC                    = kFALSE,
  Double_t            randomConeR             = 0.4,
  Double_t            trackBgrdConeR          = 0.6,
  const char*         usedTracks              = "PicoTracks",
  const char*         usedClusters            = "CaloClustersCorr",
  Double_t            trackEtaWindow          = 0.9,
  Double_t            vertexWindow            = 10.0,
  Double_t            vertexMaxR              = 1.0,
  Int_t               minVertexContributors   = 1,
  Double_t            minJetPt                = 5.0, // signal jet min pt
  Double_t            dijetLeadingMinPt       = 10.0,
  Double_t            dijetMaxAngleDev        = 10.0,
  Int_t               numberOfPtHardBins      = 0
)
{
  // #### Detect the demanded trigger with its readable name
  TString triggerName(Form("Trigger_%i", trigger));
  if (trigger == AliVEvent::kAnyINT)
    triggerName = "kAnyINT";
  else if (trigger == AliVEvent::kAny)
    triggerName = "kAny";
  else if(trigger == AliVEvent::kINT7)
    triggerName = "kINT7";
  else if(trigger == AliVEvent::kEMC7)
    triggerName = "kEMC7";
  else if(trigger == AliVEvent::kEMCEJE)
    triggerName = "kEMCEJE";
  else if(trigger == AliVEvent::kEMCEGA)
    triggerName = "kEMCEGA";

  // #### Define manager and data container names
  const char*         kFileName               = "ChargedJetsPA.root"; // hard coded
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskChargedJetsPA", "No analysis manager to connect to.");
    return NULL;
  }
  TString myContName("");
  if(isMC)
    myContName = Form("MCChargedJets_pA_R0%2.0f_%s",jetRadius*100,triggerName.Data());
  else
    myContName = Form("ChargedJets_pA_R0%2.0f_%s",jetRadius*100,triggerName.Data());

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,"",1,jetRadius,1,0.150,0.300);// anti-kt
  AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(usedTracks,"",0,jetRadius,1,0.150,0.300); // kt

  // #### Define analysis task
  AliAnalysisTaskChargedJetsPA *task = NULL;
  contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, kFileName);
  task = new AliAnalysisTaskChargedJetsPA(Form("Analysis_pA_%s_%s", jetFinderTask->GetName(), triggerName.Data()), usedTracks, usedClusters, jetFinderTask->GetName(),jetFinderTaskKT->GetName());

  // #### Task preferences
  task->SetAcceptanceWindows(trackEtaWindow, vertexWindow, vertexMaxR, minVertexContributors, jetRadius, jetRadius);
  task->SetSignalJetMinPt(minJetPt);
  task->SetSignalJetMinArea(0.6*jetRadius*jetRadius*TMath::Pi());
  task->SetDijetLeadingMinPt(dijetLeadingMinPt);
  task->SetDijetMaxAngleDeviation(dijetMaxAngleDev);
  task->SetRandConeRadius(randomConeR);
  task->SetTrackBackgroundConeRadius(trackBgrdConeR);
  task->SelectCollisionCandidates(trigger);
  if(numberOfPtHardBins)
    task->SetNumberOfPtHardBins(numberOfPtHardBins);

  // #### Add analysis task
  manager->AddTask(task);
  manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
  manager->ConnectOutput(task, 1, contHistos);
  return task;
}
