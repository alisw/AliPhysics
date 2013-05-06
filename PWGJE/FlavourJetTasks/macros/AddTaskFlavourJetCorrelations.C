// $Id$

AliAnalysisTaskFlavourJetCorrelations *AddTaskFlavourJetCorrelations(
  AliAnalysisTaskRecoJetCorrelations::ECandidateType cand = AliAnalysisTaskRecoJetCorrelations::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  TString jetArrname = "",
  TString suffix = "",
  Bool_t triggerOnLeadingJet = kFALSE,
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  Int_t acceptance = 1 /*1= 0-2pi, 2=emcal cut*/,
  Double_t areaCut = 0.)
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
  AliAnalysisTaskFlavourJetCorrelations *task = new AliAnalysisTaskFlavourJetCorrelations("AnaTaskFlavourJetCorrelations", analysiscuts, cand, jetArrname);

  //D meson settings
  task->SetMC(theMCon);
  task->SetTriggerOnLeadingJet(triggerOnLeadingJet);

  //jet settings
  task->SetJetRadius(R);
  task->SetJetPtCut(jptcut);

  Float_t etaCov=0.9;
  if (acceptance==2) etaCov=0.7; //EMCal

  Float_t minEta = -etaCov+R;
  Float_t maxEta =  etaCov-R;
  if (acceptance) task->SetJetEtaLimits(minEta, maxEta);

  Float_t minPhi = 0.;
  Float_t maxPhi = 2.*TMath::Pi();
//if (acceptance==2) { /*80-180 degree*/ }
  if (acceptance) task->SetJetPhiLimits(minPhi, maxPhi);

  //Float_t area=0.6*TMath::Pi()*R*R;
  task->SetJetAreaCut(areaCut);
  mgr->AddTask(task);

  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DEmcalJet";
  outputfile += suffix;

  TString nameContainer0="histogramsCorrelations";
  TString nameContainer1="cuts";

  nameContainer0 += suffix;
  nameContainer1 += suffix;

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();

  // ----- output data -----

  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainer1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  Printf("Input and Output connected to the manager");
  return task ;
}
