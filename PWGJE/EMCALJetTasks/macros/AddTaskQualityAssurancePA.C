AliAnalysisTaskQualityAssurancePA* AddTaskQualityAssurancePA(
  Double_t            jetRadius               = 0.4,
  Int_t               trigger                 = AliVEvent::kINT7,
  Bool_t              isMC                    = kFALSE,
  const char*         usedTracks              = "PicoTracks",
  const char*         usedClusters            = "CaloClustersCorr",
  Double_t            trackEtaWindow          = 0.9,
  Double_t            vertexWindow            = 10.0,
  Double_t            vertexMaxR              = 1.0,
  Double_t            minJetPt                = 5.0, // signal jet min pt
  Int_t               numberOfPtHardBins      = 0,
  TString             runNumbers              = "195344 195346 195351 195389 195390 195391 195478 195479 195480 195481 195482 195483",
  Bool_t              isEMCalTrain            = kTRUE
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
  else if(trigger == AliVEvent::kMB)
    triggerName = "kMB";
  else if(trigger == AliVEvent::kEMC7)
    triggerName = "kEMC7";
  else if(trigger == AliVEvent::kEMCEJE)
    triggerName = "kEMCEJE";
  else if(trigger == AliVEvent::kEMCEGA)
    triggerName = "kEMCEGA";

  // #### On EMCaltrain automatically determine the run numbers
  if (gSystem->Getenv("ETRAIN_RUNNO") && runNumbers == "")
    runNumbers = gSystem->Getenv("ETRAIN_RUNNO");


  Int_t numberOfAnalyzedRuns = runNumbers.Tokenize(" ")->GetEntries();

  ::Info("AddTaskQualityAssurancePA", Form("Running with event selection %s on N=%d runs (%s)",triggerName.Data(), numberOfAnalyzedRuns, runNumbers.Data()));

  // #### Define manager and data container names

  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskQualityAssurancePA", "No analysis manager to connect to.");
    return NULL;
  }
  TString myContName("");
  if(isMC)
    myContName = Form("QualityAssurancePA_R0%2.0f_%s_MC",jetRadius*100,triggerName.Data());
  else
    myContName = Form("QualityAssurancePA_R0%2.0f_%s",jetRadius*100,triggerName.Data());

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,"",1,jetRadius,1,0.150,0.300);// anti-kt

  // #### Define analysis task
  AliAnalysisTaskQualityAssurancePA *task = NULL;
  contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJetsPA", AliAnalysisManager::GetCommonFileName()));

  task = new AliAnalysisTaskQualityAssurancePA(Form("QualityAssurancePA_%s_%s", jetFinderTask->GetName(), triggerName.Data()), usedTracks, usedClusters, jetFinderTask->GetName());

  if(isEMCalTrain)
    RequestMemory(task,200*1024);

  // #### Task preferences
  task->SetAcceptanceWindows(trackEtaWindow, vertexWindow, vertexMaxR, jetRadius);
  task->SetRunNumbers(runNumbers);
  task->SetSignalJetMinPt(minJetPt);
  task->SetSignalJetMinArea(0.6*jetRadius*jetRadius*TMath::Pi());
  task->SelectCollisionCandidates(trigger);
  if(numberOfPtHardBins)
    task->SetNumberOfPtHardBins(numberOfPtHardBins);

  // #### Add analysis task
  manager->AddTask(task);

  manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
  manager->ConnectOutput(task, 1, contHistos);
  return task;
}
