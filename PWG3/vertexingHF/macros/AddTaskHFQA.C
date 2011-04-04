AliAnalysisTaskSEHFQA* AddTaskHFQA(AliAnalysisTaskSEHFQA::DecChannel ch,TString filecutsname="D0toKpiCuts.root",Bool_t readMC=kFALSE, Bool_t simplemode=kFALSE){
  //
  // Test macro for the AliAnalysisTaskSE for HF mesons quality assurance
  //Author: C.Bianchin chiara.bianchin@pd.infn.it

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFQA", "No analysis manager to connect to.");
    return NULL;
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(filecutsname.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: using std cut object"<<endl;
    stdcuts=kTRUE;
  }

  Bool_t onoff[4]={kTRUE,kTRUE,kFALSE,kTRUE};

  AliRDHFCuts *analysiscuts=0x0;

  TString filename="",out1name="nEntriesQA",out2name="outputPid",out3name="outputTrack",out4name="cuts",out5name="countersCentrality",out6name="outputCentrCheck",out7name="outputEvSel",inname="input",suffix="",cutsobjname="",centr="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_QA";

  switch (ch){
  case 0:
    cutsobjname="AnalysisCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dplus";
    break;
  case 1:
    cutsobjname="D0toKpiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix="D0";
    break;
  case 2:
    cutsobjname="DStartoKpipiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsDstartoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dstar";
    break;
  case 3:
    cutsobjname="DstoKKpiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
    suffix="Ds";
    break;
  case 4:
    cutsobjname="D0toKpipipiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpipipi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsD0toKpipipi*)filecuts->Get(cutsobjname);
    suffix="D04";
    break;
  case 5:
    cutsobjname="LctopKpiAnalysisCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsLctopKpi();
      analysiscuts->SetStandardCutsPP2010();
    }
    else analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get(cutsobjname);
    suffix="Lc";
    break;
  }

  inname+=suffix;
  out1name+=suffix;
  out2name+=suffix;
  out3name+=suffix;
  out4name=cutsobjname;
  out5name+=suffix;
  out6name+=suffix;
  out7name+=suffix;

  if(!analysiscuts && filecutsname!="none"){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  }

  centr=Form("%.0f%.0f",analysiscuts->GetMinCentrality(),analysiscuts->GetMaxCentrality());
  inname+=centr;
  out1name+=centr;
  out2name+=centr;
  out3name+=centr;
  out4name+=centr;
  out5name+=centr;
  out6name+=centr;
  out7name+=centr;

 
  AliAnalysisTaskSEHFQA* taskQA=new AliAnalysisTaskSEHFQA(Form("QA%s",suffix.Data()),ch,analysiscuts);

  taskQA->SetReadMC(readMC);
  taskQA->SetSimpleMode(simplemode); // set to kTRUE to go faster in PbPb
  taskQA->SetTrackOn(onoff[0]);
  taskQA->SetPIDOn(onoff[1]);
  taskQA->SetCentralityOn(onoff[2]);
  taskQA->SetEvSelectionOn(onoff[3]);
  mgr->AddTask(taskQA);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(taskQA,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(out1name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //events analysed
  mgr->ConnectOutput(taskQA,1,coutput1);

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //PID
  if(onoff[1]) mgr->ConnectOutput(taskQA,2,coutput2);

  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(out3name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of tracks
  if(onoff[0]) mgr->ConnectOutput(taskQA,3,coutput3);

  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(out4name,AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  mgr->ConnectOutput(taskQA,4,coutput4);

  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of centrality
  if(onoff[2]) mgr->ConnectOutput(taskQA,5,coutput5);

  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of centrality
  if(onoff[2]) mgr->ConnectOutput(taskQA,6,coutput6);

  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(out7name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //event selection
  if(onoff[3]) mgr->ConnectOutput(taskQA,7,coutput7);

 return taskQA;
}

