AliAnalysisTaskSELbtoLcpi4 *AddTaskLbUPG(TString finname="LcdaughtersCut.root",Int_t ndebug=0,Int_t applyFixesITS3Analysis=-1,const char*  postname="")
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
     
  task->SetCutsond0Lcdaughters(kFALSE);
  task->SetPtConfiguration(0.,30.,2.,14.,0.,999.,4.);

  //temporary to check effect of fixes separately
  if(applyFixesITS3Analysis == 0){
    task->SetApplyFixesITS3AnalysisBit(kTRUE);
    task->SetApplyFixesITS3AnalysiskAll(kTRUE);
    task->SetApplyFixesITS3AnalysisHijing(kTRUE);
  } else if(applyFixesITS3Analysis == 1){
    task->SetApplyFixesITS3AnalysisBit(kTRUE);
  } else if(applyFixesITS3Analysis == 2){
    task->SetApplyFixesITS3AnalysiskAll(kTRUE);
  } else if(applyFixesITS3Analysis == 3){
    task->SetApplyFixesITS3AnalysisHijing(kTRUE);
  } else {
    task->SetApplyFixesITS3AnalysisBit(kFALSE);
    task->SetApplyFixesITS3AnalysiskAll(kFALSE);
    task->SetApplyFixesITS3AnalysisHijing(kFALSE);
  }
 
  mgr->AddTask(task);

  
  TString outputFileName=(char*)AliAnalysisManager::GetCommonFileName();
  TString container1;
  TString container2;
  TString directory;
  directory=Form(":Lbtask%s",(char*)postname);
  outputFileName+=directory;
  container1.Form("HistosUPG%s", (char*)postname);
  container2.Form("ntupleUPGLb%s",(char*)postname);
  

  AliAnalysisDataContainer *coutput1
     =mgr->CreateContainer(container1,
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputFileName.Data());


  AliAnalysisDataContainer *coutput2
     =mgr->CreateContainer(container2,
                           TNtuple::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputFileName.Data());
 

  coutput2->SetSpecialOutput();
  mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  return task;
}

