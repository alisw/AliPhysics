AliAnalysisTaskChargedJetsPA* AddTaskChargedJetsPA(
  const char*         containerSuffix         = "",
  Double_t            jetRadius               = 0.4,
  const char*         centralityType          = "V0A",
  Int_t               trigger                 = AliVEvent::kINT7,
  Bool_t              isPA                    = kTRUE,
  Bool_t              isMC                    = kFALSE,
  Bool_t              doJetProfileAnalysis    = kFALSE,
  Bool_t              doTrackcutAnalysis      = kFALSE,
  Bool_t              doJetAnalysis           = kTRUE,
  const char*         usedTracks              = "PicoTracks",
  Int_t               numberOfCentralityBins  = 20,
  Double_t            areaPercentage          = 0.6,
  Double_t            ktJetRadius             = 0.4,
  Double_t            trackBgrdConeR          = 0.6,
  Double_t            minJetTrackPt           = 0.150,
  Double_t            minEta                  = -0.9,
  Double_t            maxEta                  = +0.9,
  Double_t            minJetEta               = -0.5,
  Double_t            maxJetEta               = +0.5,
  Int_t               recombscheme            = 1,
  Bool_t              isEMCalTrain            = kFALSE
)
{
  cout << " ############ MACRO EXECUTION STARTED ############\n";
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

  TString containerNameSuffix("");

  if (strcmp(containerSuffix,""))
    containerNameSuffix = Form("_%s", containerSuffix);

  TString bgrdName("Background");
  TString myContName("");
  TString myContJPName("");
  TString myContTCName("");
  if(isMC)
  {
    bgrdName = Form("BackgroundR0%2.0f_%s_MC%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContName = Form("AnalysisR0%2.0f_%s_MC%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContJPName = Form("JetProfileR0%2.0f_%s_MC%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContTCName = Form("TrackcutsR0%2.0f_%s_MC%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
  }
  else
  {
    bgrdName = Form("BackgroundR0%2.0f_%s%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContName = Form("AnalysisR0%2.0f_%s%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContJPName = Form("JetProfileR0%2.0f_%s%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
    myContTCName = Form("TrackcutsR0%2.0f_%s%s", jetRadius*100, triggerName.Data(), containerNameSuffix.Data());
  }

  if(doJetAnalysis)
  {
    // #### Add necessary jet finder tasks
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,"",1,jetRadius,1,minJetTrackPt,0.300,0.005,recombscheme, "Jet", 0.0); // anti-kt
    AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(usedTracks,"",0,ktJetRadius,1,minJetTrackPt,0.300,0.005,recombscheme, "Jet", 0.0); // kt
    cout << " Jet finder tasks added: " <<  jetFinderTask << " + " <<  jetFinderTaskKT << endl;
    jetFinderTask->SelectCollisionCandidates(AliVEvent::kAny);
    jetFinderTaskKT->SelectCollisionCandidates(AliVEvent::kAny);

    // #### Define external rho task
    AliEmcalJetTask* jetFinderRho = AddTaskEmcalJet(usedTracks,"",1,0.4,1,minJetTrackPt,0.300,0.005,recombscheme, "Jet", 0.0); // anti-kt
    AliEmcalJetTask* jetFinderRhoKT = AddTaskEmcalJet(usedTracks,"",0,0.4,1,minJetTrackPt,0.300,0.005,recombscheme, "Jet", 0.0); // kt
    jetFinderRho->SelectCollisionCandidates(AliVEvent::kAny);
    jetFinderRhoKT->SelectCollisionCandidates(AliVEvent::kAny);
    cout << " Jet finder tasks (used for bgrd) added: " <<  jetFinderRho << " + " <<  jetFinderRhoKT << endl;
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
    AliAnalysisTaskRhoSparse* rhotask = AddTaskRhoSparse(jetFinderRhoKT->GetName(), NULL, usedTracks, "", bgrdName.Data(), 0.4,"TPC", 0., 5., 0, 0,2,kFALSE,bgrdName.Data(),kTRUE);
    rhotask->SelectCollisionCandidates(AliVEvent::kAny);
    cout << " Background task added: " <<  rhotask << endl;
  }

  // #### Define analysis task
  AliAnalysisTaskChargedJetsPA *task = NULL;
  AliAnalysisDataContainer* contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJets", AliAnalysisManager::GetCommonFileName()));
  if(doJetProfileAnalysis)
    AliAnalysisDataContainer* contJetProfile = manager->CreateContainer(myContJPName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJets", AliAnalysisManager::GetCommonFileName()));
  if(doTrackcutAnalysis)
    AliAnalysisDataContainer* contTrackcuts = manager->CreateContainer(myContTCName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJets", AliAnalysisManager::GetCommonFileName()));

  if(doJetAnalysis)
  {
    task = new AliAnalysisTaskChargedJetsPA(Form("AnalysisPA%s_%s_%s", containerNameSuffix.Data(), jetFinderTask->GetName(), triggerName.Data()), usedTracks, jetFinderTask->GetName(),jetFinderTaskKT->GetName(), doJetProfileAnalysis, doTrackcutAnalysis);
    task->SetExternalRhoTaskName(bgrdName.Data());
  }
  else
    task = new AliAnalysisTaskChargedJetsPA(Form("AnalysisPA%s_%s", containerNameSuffix.Data(), triggerName.Data()), usedTracks, "","", doJetProfileAnalysis, doTrackcutAnalysis);

  cout << " Main task created: " <<  task << endl;


  // #### Task preferences
  task->SetIsPA(isPA);
  task->SetAnalyzeJetProfile(doJetProfileAnalysis);
  task->SetAnalyzeTrackcuts(doTrackcutAnalysis);
  task->SetAcceptanceEta(minEta,maxEta);
  task->SetAcceptanceJetEta(minJetEta,maxJetEta);
  task->SetSignalJetRadius(jetRadius);
  task->SetBackgroundJetRadius(jetRadius);
  task->SetSignalJetMinArea(areaPercentage*jetRadius*jetRadius*TMath::Pi());
  task->SetRandConeRadius(jetRadius);
  task->SelectCollisionCandidates(trigger);
  task->SetCentralityType(centralityType);
  task->SetNumberOfCentralityBins(numberOfCentralityBins);
  task->SetDoJetAnalysis(doJetAnalysis);
  // for the case of pp analysis, set special settings
  if(!isPA)
  {
    task->SetCentralityToOne(1);
    task->SetUsePileUpCut(0);
    task->SetUseDefaultVertexCut(kFALSE);
    cout << " Using pp settings: No pileup correction, simple vertex correction, centrality=1" << endl;
  }
  cout << " Settings set." << endl;

  // #### Add analysis task
  manager->AddTask(task);
  cout << " Task added to manager" << endl;
  manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
  cout << " Input connected, common input container: " << manager->GetCommonInputContainer() << endl;
  manager->ConnectOutput(task, 1, contHistos);
  cout << " Output connected, contHistos: " << contHistos << endl;

  if(doJetProfileAnalysis)
  {
    manager->ConnectOutput(task, 2, contJetProfile);
  }
  if(doTrackcutAnalysis && !doJetProfileAnalysis)
  {
    manager->ConnectOutput(task, 2, contTrackcuts);
  }
  else if(doTrackcutAnalysis && doJetProfileAnalysis)
  {
    manager->ConnectOutput(task, 3, contTrackcuts);
  }

  if(isEMCalTrain)
    RequestMemory(task,600*1024);

  cout << " ############ MACRO EXECUTION SUCCESSFUL, will return " << task << " ############\n";
  return task;
}
