// AddTaskSEDmesonsFilterCJ.C

AliAnalysisTaskSEDmesonsFilterCJ *AddTaskSEDmesonsFilterCJ(AliAnalysisTaskSEDmesonsFilterCJ::ECandidateType cand = AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi,
                                                           TString filename = "DStartoKpipiCuts.root",
                                                           Bool_t theMCon = kFALSE,
                                                           Bool_t reco = kTRUE, /*must be true if theMCon is false*/
                                                           TString suffix = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSEDmesonsFilterCJ", "No analysis manager to connect to.");
    return NULL;
  } 

  Bool_t useStdC = kFALSE;
  TFile* filecuts = TFile::Open(filename);
  if (!filecuts || (filecuts && !filecuts->IsOpen())) {
    ::Warning("AddTaskSEDmesonsFilterCJ", "Input file not found: use std cuts");
    useStdC = kTRUE;
  }

  if (!theMCon && !reco) {
    ::Warning("AddTaskSEDmesonsFilterCJ", "'reco' must be kTRUE if 'theMCon' is kFALSE!");
    reco = kTRUE;
  }

  if (theMCon) {
    suffix += "MC";
    if (reco) suffix += "rec";  
  }

  TString candname("DStar"); 
  if (cand==0) candname="D0";

  AliRDHFCuts *analysiscuts = 0x0;
  switch (cand) {
  case 0 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    break;
  case 1 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    analysiscuts->SetName("DStartoKpipiCuts");
    break;
  }
  
  if (!analysiscuts) {
    AliFatal("Specific AliRDHFCuts not found!");
    return;
  }

  TString taskFiltername("DmesonsFilterCJ");
  taskFiltername += candname;
  taskFiltername += suffix;
  if (theMCon) taskFiltername += "MC";
  if (!reco)   taskFiltername += "gen";
  
  AliAnalysisTaskSEDmesonsFilterCJ* task = mgr->GetTask(taskFiltername.Data());
  if (task) {
    ::Info("AddTaskSEDmesonsFilterCJ", Form("Task %s already exist, continue",taskFiltername.Data()));
    return task;
  }
  else {
    ::Info("AddTaskSEDmesonsFilterCJ", "Creating the task");

    // create the task
    task = new AliAnalysisTaskSEDmesonsFilterCJ(taskFiltername.Data(), analysiscuts, cand);
    task->SetMC(theMCon);
    task->SetUseReco(reco);
    mgr->AddTask(task);
  }
  
  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DmesonsForJetCorrelations";
  outputfile += suffix;

  TString nameContainer0 = "histograms";
  TString nameContainer1 = "cuts";
  TString nameContainer2 = "Dcandidates";
  TString nameContainer3 = "DSBcandidates";
  TString nameContainer4 = "DcandidatesAndTracks";
  TString nameContainer5 = "DSBcandidatesAndTracks";
  TString nameContainer6 = "MCDcandidatesAndTracks";
  TString nameContainer7 = "NormalizationCounter";

  nameContainer0 += candname;
  nameContainer1 += candname;
  nameContainer2 += candname;
  nameContainer3 += candname;
  nameContainer4 += candname;
  nameContainer5 += candname;
  nameContainer6 += candname;
  nameContainer7 += candname;  

  nameContainer0 += suffix;
  nameContainer1 += suffix;
  nameContainer2 += suffix;
  nameContainer3 += suffix;
  nameContainer4 += suffix;
  nameContainer5 += suffix;
  nameContainer6 += suffix;
  nameContainer7 += suffix;

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer0, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainer1, AliRDHFCuts::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(nameContainer2, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(nameContainer3, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(nameContainer4, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(nameContainer5, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(nameContainer6, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer(nameContainer7, AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()); 
 
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);
  mgr->ConnectOutput(task,5,coutput5);
  mgr->ConnectOutput(task,6,coutput6);
  mgr->ConnectOutput(task,7,coutput7);
  mgr->ConnectOutput(task,8,coutput8);

  ::Info("AddTaskSEDmesonsFilterCJ", "Input and Output connected to the manager");
  return task;
}

