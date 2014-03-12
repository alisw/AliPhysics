AliAnalysisTaskChargedJetsPA* AddTaskChargedJetsPA(
  Double_t            jetRadius               = 0.4,
  Int_t               trigger                 = AliVEvent::kINT7,
  Bool_t              isMC                    = kFALSE,
  Int_t               ptHardBin               = -1,
  Double_t            randomConeR             = 0.4,
  Double_t            trackBgrdConeR          = 0.6,
  const char*         containerSuffix         = "",
  const char*         usedTracks              = "PicoTracks",
  const char*         centralityType          = "V0A",
  Double_t            trackEtaWindow          = 0.9,
  Double_t            minJetPt                = 0.15, // signal jet min pt
  Double_t            minBackgroundJetPt      = 0.0, // background jet min pt
  Double_t            dijetLeadingMinPt       = 10.0,
  Double_t            dijetMaxAngleDev        = 10.0,
  Int_t               numberOfPtHardBins      = 0,
  const char*         externalMacro           = NULL,
  Bool_t              useVertexCut            = kTRUE,
  Bool_t              usePileUpCut            = kTRUE,
  Bool_t              isEMCalTrain            = kFALSE,
  Bool_t              calculateExternalRho    = kTRUE,
  Bool_t              analyzeDeprecatedBackgrounds = kTRUE,
  Int_t               numberOfCentralityBins  = 20
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

  // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskChargedJetsPA", "No analysis manager to connect to.");
    return NULL;
  }

  TString stringPtHard("");
  TString containerNameSuffix("");

  if (ptHardBin!=-1)
    stringPtHard = Form("_PtHard_%d",ptHardBin);
  if (strcmp(containerSuffix,""))
    containerNameSuffix = Form("_%s", containerSuffix);

  TString myContName("");
  if(isMC)
    myContName = Form("AnalysisR0%2.0f_%s_MC%s%s", jetRadius*100, triggerName.Data(), stringPtHard.Data(), containerNameSuffix.Data());
  else
    myContName = Form("AnalysisR0%2.0f_%s%s%s", jetRadius*100, triggerName.Data(), stringPtHard.Data(), containerNameSuffix.Data());

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,"",1,jetRadius,1,0.150,0.300); // anti-kt
  AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(usedTracks,"",0,jetRadius,1,0.150,0.300); // kt

  if(jetRadius < 0.1)
  {
    jetFinderTask->SetMinJetArea(0.0);
    jetFinderTaskKT->SetMinJetArea(0.0);
    jetFinderTask->SetMinJetPt(0.15);
    jetFinderTaskKT->SetMinJetPt(0.15);
    jetFinderTask->SetGhostArea(0.001);
    jetFinderTaskKT->SetGhostArea(0.001);
  }

  if(minBackgroundJetPt == -1.0)
  {
    if(analyzeDeprecatedBackgrounds)
      minBackgroundJetPt = 0.0;
    else
      minBackgroundJetPt = 0.15;
  }


  jetFinderTaskKT->SetMinJetPt(minBackgroundJetPt);

  // #### Define extern rho task
  if(calculateExternalRho)
  {
    TString myRhoName("ExternalRhoTask");

    AliAnalysisTaskRhoSparse* mgrTask = manager->GetTask(myRhoName.Data());
    if(!mgrTask) // rho task does not yet exist
    {
      AliEmcalJetTask* jetFinderRho = AddTaskEmcalJet(usedTracks,"",1,0.4,1,0.150,0.300); // anti-kt
      AliEmcalJetTask* jetFinderRhoKT = AddTaskEmcalJet(usedTracks,"",0,0.4,1,0.150,0.300); // kt
      jetFinderRhoKT->SetMinJetPt(0);

      gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
      AliAnalysisTaskRhoSparse* rhotask = AddTaskRhoSparse(jetFinderRhoKT->GetName(), jetFinderRho->GetName(), usedTracks, "", myRhoName.Data(), jetRadius,"TPC", 0., 5., 0, 0,0,kFALSE,myRhoName.Data(),kTRUE);
    }    
  }

  // #### Define analysis task
  AliAnalysisTaskChargedJetsPA *task = NULL;
  contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJetsPA", AliAnalysisManager::GetCommonFileName()));
  task = new AliAnalysisTaskChargedJetsPA(Form("AnalysisPA_%s_%s", jetFinderTask->GetName(), triggerName.Data()), usedTracks, jetFinderTask->GetName(),jetFinderTaskKT->GetName());

  // #### Task preferences
  task->SetAcceptanceWindows(trackEtaWindow, jetRadius, jetRadius);
  task->SetAnalyzeQA(kTRUE);
  task->SetAnalyzeBackground(kTRUE);
  task->SetAnalyzeDeprecatedBackgrounds(analyzeDeprecatedBackgrounds);
  task->SetUsePileUpCut(usePileUpCut);
  task->SetUseDefaultVertexCut(useVertexCut);
  task->SetSignalJetMinPt(minJetPt);
  task->SetSignalJetMinArea(0.6*jetRadius*jetRadius*TMath::Pi());
  task->SetDijetLeadingMinPt(dijetLeadingMinPt);
  task->SetDijetMaxAngleDeviation(dijetMaxAngleDev);
  task->SetRandConeRadius(randomConeR);
  task->SetBackgroundJetMinPt(minBackgroundJetPt);
  task->SetTRBackgroundConeRadius(trackBgrdConeR);
  task->SelectCollisionCandidates(trigger);
  task->SetCentralityType(centralityType);
  task->SetNumberOfCentralityBins(numberOfCentralityBins);
  task->SetUsePtHardBin(ptHardBin);
  if(calculateExternalRho)
    task->SetExternalRhoTaskName(myRhoName.Data());

  if(numberOfPtHardBins)
    task->SetNumberOfPtHardBins(numberOfPtHardBins);

  // #### Add analysis task
  manager->AddTask(task);
  manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
  manager->ConnectOutput(task, 1, contHistos);

  if(isEMCalTrain)
    RequestMemory(task,200*1024);

  // #### Do some nasty piggybacking on demand
  if (externalMacro)
    gROOT->Macro(externalMacro);

  return task;
}
