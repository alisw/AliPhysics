AliAnalysisTaskChargedJetsPA* AddTaskChargedJetsPA(
  Double_t            jetRadius               = 0.4,
  Int_t               trigger                 = AliVEvent::kINT7,
  Bool_t              isMC                    = kFALSE,
  Double_t            randomConeR             = 0.4,
  Double_t            trackBgrdConeR          = 0.6,
  const char*         usedTracks              = "PicoTracks",
  const char*         centralityType          = "V0A",
  Double_t            trackEtaWindow          = 0.9,
  Double_t            vertexWindow            = 10.0,
  Double_t            vertexMaxR              = 1.0,
  Double_t            minJetPt                = 5.0, // signal jet min pt
  Double_t            dijetLeadingMinPt       = 10.0,
  Double_t            dijetMaxAngleDev        = 10.0,
  Int_t               numberOfPtHardBins      = 0,
  const char*         fileEtaCorrectionFactors= "alien:///alice/cern.ch/user/r/rhaake/pA/EtaCorrectionFactors.root",
  const char*         externalMacro           = NULL
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
  TString myContName("");
  if(isMC)
    myContName = Form("AnalysisR0%2.0f_%s_MC",jetRadius*100,triggerName.Data());
  else
    myContName = Form("AnalysisR0%2.0f_%s",jetRadius*100,triggerName.Data());

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,"",1,jetRadius,1,0.150,0.300); // anti-kt
  AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(usedTracks,"",0,jetRadius,1,0.150,0.300); // kt

  // #### Load correction factors from alien
  TH1D* corrFactorsKT = NULL;
  TH1D* corrFactorsRC = NULL;
  TH1D* corrFactorsTR = NULL;
  if (fileEtaCorrectionFactors)
  {
    // trying to connect to alien
    if (!TGrid::Connect("alien://"))
      ::Warning("AddTaskChargedJetsPA", "AliEn connection failed!");
    else
    { 
      ::Info("AddTaskChargedJetsPA", "AliEn connection successful!");
      // Copy eta correction file
      Bool_t copied = TFile::Cp(fileEtaCorrectionFactors,"file:EtaCorrectionFactors.root");
      if(copied)
      {
        TFile* tmpFile= new TFile("EtaCorrectionFactors.root","READ");
        corrFactorsKT = static_cast<TH1D*>(tmpFile->Get("EtaCorrectionFactorsKT"));
        corrFactorsRC = static_cast<TH1D*>(tmpFile->Get("EtaCorrectionFactorsRC"));
        corrFactorsTR = static_cast<TH1D*>(tmpFile->Get("EtaCorrectionFactorsTR"));
      }
      else
        ::Warning("AddTaskChargedJetsPA", "AliEn copying failed!");
    }
  }
  // #### Define analysis task
  AliAnalysisTaskChargedJetsPA *task = NULL;
  contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJetsPA", AliAnalysisManager::GetCommonFileName()));
  task = new AliAnalysisTaskChargedJetsPA(Form("AnalysisPA_%s_%s", jetFinderTask->GetName(), triggerName.Data()), usedTracks, jetFinderTask->GetName(),jetFinderTaskKT->GetName());

  // #### Task preferences
  task->SetAcceptanceWindows(trackEtaWindow, vertexWindow, vertexMaxR, jetRadius, jetRadius);
  task->SetSignalJetMinPt(minJetPt);
  task->SetSignalJetMinArea(0.6*jetRadius*jetRadius*TMath::Pi());
  task->SetDijetLeadingMinPt(dijetLeadingMinPt);
  task->SetDijetMaxAngleDeviation(dijetMaxAngleDev);
  task->SetRandConeRadius(randomConeR);
  task->SetTRBackgroundConeRadius(trackBgrdConeR);
  task->SelectCollisionCandidates(trigger);
  task->SetCentralityType(centralityType);
  if(numberOfPtHardBins)
    task->SetNumberOfPtHardBins(numberOfPtHardBins);

  if(corrFactorsKT)
    task->SetKTEtaCorrectionFactors(corrFactorsKT);
  if(corrFactorsRC)
    task->SetRCEtaCorrectionFactors(corrFactorsRC);
  if(corrFactorsTR)
    task->SetTREtaCorrectionFactors(corrFactorsTR);

  // #### Add analysis task
  manager->AddTask(task);
  manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
  manager->ConnectOutput(task, 1, contHistos);

  // #### Do some nasty piggybacking on demand
  if (externalMacro)
    gROOT->LoadMacro(externalMacro);


  return task;
}
