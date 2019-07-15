AliAnalysisTaskSELbtoLcpi4 *AddTaskLbUPG(TString finname="LcdaughtersCut.root",Int_t ndebug=0, Bool_t lccuts=kFALSE,  Double_t *lccutvalues=NULL, Double_t *ptcutvalues=NULL,const char*  postname="")
 {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAliAnalysisTaskSELbtoLcpi", "No analysis manager to connect to.");
    return 0;
  }

  TFile* filecuts=TFile::Open(finname.Data());
  if(!filecuts->IsOpen()){
   printf("Input file not found: using std cut object\n");
   auto stdcuts=kTRUE;
 }
  AliRDHFCutsLctopKpi* prodcuts=new AliRDHFCutsLctopKpi();
  prodcuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiProdCuts");
  prodcuts->SetName("LctopKpiProdCuts");
  prodcuts->SetMinPtCandidate(-1.);
  prodcuts->SetMaxPtCandidate(10000.);

  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
  analysiscuts->SetName("LctopKpiAnalysisCuts");
  analysiscuts->SetMinPtCandidate(-1.);
  analysiscuts->SetMaxPtCandidate(10000.);

  //analysis task
  TString combinedName="LambdabAnalysis";
  
  combinedName.Form("LambdabAnalysis%s",(char*)postname);
  AliAnalysisTaskSELbtoLcpi4 *task
    =new AliAnalysisTaskSELbtoLcpi4( combinedName,
                                     kTRUE,
                                     analysiscuts,
                                     prodcuts,
                                     ndebug);

  task->SetCutsond0Lcdaughters(lccuts);
  if (lccuts){
  task->ApplyD0CutLcdaughters(lccutvalues);
  }
  task->SetPtConfiguration(ptcutvalues);

  mgr->AddTask(task);

  
  char* outputFileName=(char*)AliAnalysisManager::GetCommonFileName();
  TString container1;
  TString container2;
  TString container3;
  container1.Form("HistosUPG%s", (char*)postname);
  container2.Form("ntupleUPGLb%s",(char*)postname);
  container3.Form("ntupleUPGLc%s",(char*)postname);
  

  AliAnalysisDataContainer *coutput1
     =mgr->CreateContainer(container1,
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", outputFileName, "ITSImproverUpg"));


  AliAnalysisDataContainer *coutput2
     =mgr->CreateContainer(container2,
                           TNtuple::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s",outputFileName, "ITSImproverUpg"));

  AliAnalysisDataContainer *coutput3
     =mgr->CreateContainer(container3,
                           TNtuple::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s",outputFileName, "ITSImproverUpg"));


  mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  return task;
}

