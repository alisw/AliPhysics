void *AddTaskDFilterAndCorrelations(
  AliAnalysisTaskSEDmesonsFilterCJ::ECandidateType cand = AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  Bool_t reco = kTRUE /*must be true if theMCon is false*/,
  TString suffix = "",
  TString jetArrname = "",
  Bool_t triggerOnLeadingJet = kFALSE,
  Int_t leadingHadType = 0 /*0 = charged, 1 = neutral, 2 = both*/,
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  const char *cutType = "TPC",
  Double_t percjetareacut = 1.)
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
  AliAnalysisTaskSEDmesonsFilterCJ *taskFilter = new AliAnalysisTaskSEDmesonsFilterCJ("AnaTaskSEDmesonsFilterCJ",analysiscuts,cand);
  if(!theMCon) reco=kTRUE;
  taskFilter->SetMC(theMCon); //D meson settings
  taskFilter->SetUseReco(reco);
  taskFilter->SetName("AliAnalysisTaskSEDmesonsFilterCJ");
  mgr->AddTask(taskFilter);
  
    // create the task
    
  AliAnalysisTaskFlavourJetCorrelations *taskCorr = new AliAnalysisTaskFlavourJetCorrelations("AnaTaskFlavourJetCorrelations", 
     analysiscuts, cand);
  taskCorr->SetName("AliAnalysisTaskFlavourJetCorrelations");
  taskCorr->SetJetsName(jetArrname);
  taskCorr->SetMC(theMCon);
  taskCorr->SetUseReco(reco);
  taskCorr->SetTriggerOnLeadingJet(triggerOnLeadingJet);
  taskCorr->SetJetAcceptanceType(cutType);
  taskCorr->SetJetPtCut(jptcut);
  taskCorr->SetPercAreaCut(percjetareacut);
  
  mgr->AddTask(taskCorr);

  if(theMCon) {
     suffix+="MC";
     if(reco) suffix+="rec";  
  }
  
  TString candname="DStar"; 
  if(cand==0)  candname="D0";
  
  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputfileF = outputfile, outputfileC = outputfile;
  outputfileF += ":PWG3_D2H_DmesonsForJetCorrelations";
  outputfileC += ":PWG3_D2H_DEmcalJet";
  
  outputfileF += suffix;
  outputfileC += suffix;

  TString nameContainerF0="histograms";
  TString nameContainerF1="cuts";
  
  TString nameContainerC0="hCor";
  TString nameContainerC1="cutsJ";

  TString nameContainerFC2="Dcandidates";
  TString nameContainerFC3="DSBcandidates";

  nameContainerF0  += candname;
  nameContainerF1  += candname;
  nameContainerC0  += candname;
  nameContainerC1  += candname;
  nameContainerFC2 += candname;
  nameContainerFC3 += candname;
  
  nameContainerF0  += suffix;
  nameContainerF1  += suffix;
  nameContainerC0  += suffix;
  nameContainerC1  += suffix;
  nameContainerFC2 += suffix;
  nameContainerFC3 += suffix;
  

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  AliAnalysisDataContainer *coutputF0 = mgr->CreateContainer(nameContainerF0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfileF.Data());
  
  AliAnalysisDataContainer *coutputF1 = mgr->CreateContainer(nameContainerF1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfileF.Data());
  
  AliAnalysisDataContainer *coutputC0 = mgr->CreateContainer(nameContainerC0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfileC.Data());

  AliAnalysisDataContainer *coutputC1 = mgr->CreateContainer(nameContainerC1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfileC.Data());
  
  AliAnalysisDataContainer *coutputFC2 = mgr->CreateContainer(nameContainerFC2, TClonesArray::Class(),AliAnalysisManager::kExchangeContainer, outputfileF.Data()); //
  
  AliAnalysisDataContainer *coutputFC3 = mgr->CreateContainer(nameContainerFC3, TClonesArray::Class(),AliAnalysisManager::kExchangeContainer, outputfileF.Data()); //
  
  mgr->ConnectInput(taskFilter,0,mgr->GetCommonInputContainer());
  mgr->ConnectInput(taskCorr,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(taskFilter,1,coutputF0);
  mgr->ConnectOutput(taskFilter,2,coutputF1);
  mgr->ConnectOutput(taskFilter,3,coutputFC2);
  mgr->ConnectOutput(taskFilter,4,coutputFC3);
  
  
  mgr->ConnectInput(taskCorr,1,coutputFC2);
  mgr->ConnectInput(taskCorr,2,coutputFC3);
  mgr->ConnectOutput(taskCorr,1,coutputC0);
  mgr->ConnectOutput(taskCorr,2,coutputC1);
  //if(cand==1) mgr->ConnectOutput(task,4,coutput4);

  Printf("Input and Output connected to the manager");
  return; 
}

