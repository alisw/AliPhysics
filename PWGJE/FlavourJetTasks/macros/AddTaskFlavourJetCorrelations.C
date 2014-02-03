AliAnalysisTaskFlavourJetCorrelations *AddTaskFlavourJetCorrelations(
  AliAnalysisTaskFlavourJetCorrelations::ECandidateType cand = AliAnalysisTaskFlavourJetCorrelations::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  Bool_t reco = kTRUE /*must be true if theMCon is false*/,
  TString jetArrname = "",
  TString suffix = "",
  Bool_t triggerOnLeadingJet = kFALSE,
  Int_t leadingHadType = 0 /*0 = charged, 1 = neutral, 2 = both*/,
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  const char *cutType = "TPC",
  Double_t percjetareacut = 1.)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFlavourJetCorrelations::AddTaskFlavourJetCorrelations", "No analysis manager to connect to.");
    return NULL;
  } 

  Bool_t useStdC = kFALSE;
  TFile* filecuts=TFile::Open(filename);
  if (!filecuts || (filecuts && !filecuts->IsOpen())) {
    cout<<"Input file not found: use std cuts"<<endl;
    useStdC = kTRUE;
  }

  AliRDHFCuts *analysiscuts = 0x0;
  switch (cand) {
  case 0 :
    if (useStdC) {
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
  AliAnalysisTaskFlavourJetCorrelations *task = new AliAnalysisTaskFlavourJetCorrelations("AnaTaskFlavourJetCorrelations", 
     analysiscuts, cand);
  task->SetJetsName(jetArrname);
  task->SetMC(theMCon);
  task->SetUseReco(reco);
  task->SetTriggerOnLeadingJet(triggerOnLeadingJet);
  task->SetJetAcceptanceType(cutType);
  task->SetJetPtCut(jptcut);
  task->SetPercAreaCut(percjetareacut);
  
  mgr->AddTask(task);

  if(theMCon) {
     suffix+="MC";
     if(reco) suffix+="rec";  
  }

  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DEmcalJet";
  outputfile += suffix;

  TString candname="DStar"; 
  if(cand==0)  candname="D0";

  TString nameContainer1="hCor";
  TString nameContainer2="cutsJ";

  nameContainer1 += candname;
  nameContainer2 += candname;
  
  nameContainer1 += suffix;
  nameContainer2 += suffix;
  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();

  // ----- output data -----
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainer2, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);


  Printf("Input and Output connected to the manager");
  return task ;
}
