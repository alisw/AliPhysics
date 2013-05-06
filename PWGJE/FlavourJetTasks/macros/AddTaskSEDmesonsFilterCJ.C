// $Id$

AliAnalysisTaskSEDmesonsFilterCJ *AddTaskSEDmesonsFilterCJ(
  AliAnalysisTaskSEDmesonsForJetCorrelations::ECandidateType cand = AliAnalysisTaskSEDmesonsForJetCorrelations::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  TString suffix = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSEDmesonsFilterCJ", "No analysis manager to connect to.");
    return NULL;
  } 

  Bool_t useStdC = kFALSE;
  TFile* filecuts=TFile::Open(filename);
  if(!filecuts || (filecuts && !filecuts->IsOpen())) {
    cout<<"Input file not found: use std cuts"<<endl;
    useStdC = kTRUE;
  }

  AliRDHFCuts *analysiscuts=0x0;
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
  
  if (!analysiscuts) { // mm let's see if everything is ok
    AliFatal("Specific AliRDHFCuts not found");
    return;
  } 

  printf("CREATE TASK\n"); //CREATE THE TASK

  // create the task
  AliAnalysisTaskSEDmesonsFilterCJ *task = new AliAnalysisTaskSEDmesonsFilterCJ("AnaTaskSEDmesonsFilterCJ",analysiscuts,cand);
  task->SetMC(theMCon); //D meson settings
  mgr->AddTask(task);

  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DmesonsForJetCorrelations";
  outputfile += suffix;

  TString nameContainer0="histograms";
  TString nameContainer1="cuts";
//TString nameContainer2="DcandidatesSel";

  nameContainer0 += suffix;
  nameContainer1 += suffix;
//nameContainer2 += suffix;

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainer1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
//AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(nameContainer2, TList::Class(),AliAnalysisManager::kExchangeContainer, outputfile.Data());
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
//mgr->ConnectOutput(task,3,coutput3);

  Printf("Input and Output connected to the manager");
  return task ;
}

