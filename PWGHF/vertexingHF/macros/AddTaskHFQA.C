AliAnalysisTaskSEHFQA* AddTaskHFQA(AliAnalysisTaskSEHFQA::DecChannel ch,TString filecutsname="",Bool_t readMC=kFALSE, Bool_t simplemode=kFALSE, Int_t system=1 /*0=pp, 1=PbPb*/, TString finDirname="",Bool_t trackon=kTRUE,Bool_t pidon=kTRUE,Bool_t centralityon=kTRUE, Bool_t eventselon=kTRUE, Bool_t flowobson=kFALSE){
  //
  // Test macro for the AliAnalysisTaskSE for HF mesons quality assurance
  //Author: C.Bianchin chiara.bianchin@pd.infn.it

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFQA", "No analysis manager to connect to.");
    return NULL;
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filecutsname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(filecutsname.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }
 
  if(system==0) centralityon=kFALSE;

  AliRDHFCuts *analysiscuts=0x0;

  TString filename="",out1name="nEntriesQA",out2name="outputPid",out3name="outputTrack",out4name="cuts",out5name="countersCentrality",out6name="outputCentrCheck",out7name="outputEvSel",out8name="outputFlowObs",inname="input",suffix="",cutsobjname="",centr="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_QA";
  filename += finDirname.Data();

  switch (ch){
  case 0:
    cutsobjname="AnalysisCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dplus";
    break;
  case 1:
    cutsobjname="D0toKpiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix="D0";
    break;
  case 2:
    cutsobjname="DStartoKpipiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDstartoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dstar";
    break;
  case 3:
    cutsobjname="DstoKKpiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
    suffix="Ds";
    break;
  case 4:
    cutsobjname="D0toKpipipiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpipipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsD0toKpipipi*)filecuts->Get(cutsobjname);
    suffix="D04";
    break;
  case 5:
    cutsobjname="LctopKpiAnalysisCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsLctopKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get(cutsobjname);
    suffix="Lc";
    break;
  case 6:
    cutsobjname="LctoV0AnalysisCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsLctoV0bachelor();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsLctoV0*)filecuts->Get(cutsobjname);
    suffix="LcToV0x";
    break;
  }

  inname+=suffix;
  out1name+=suffix;
  out2name+=suffix;
  out3name+=suffix;
  out4name=cutsobjname;
  out4name+=suffix;
  out5name+=suffix;
  out6name+=suffix;
  out7name+=suffix;
  out8name+=suffix;

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
  out8name+=centr;
  inname+= finDirname.Data();
  out1name+= finDirname.Data();
  out2name+= finDirname.Data();
  out3name+= finDirname.Data();
  out4name+= finDirname.Data();
  out5name+= finDirname.Data();
  out6name+= finDirname.Data();
  out7name+= finDirname.Data();
  out8name+= finDirname.Data();

 
  AliAnalysisTaskSEHFQA* taskQA=new AliAnalysisTaskSEHFQA(Form("QA%s",suffix.Data()),ch,analysiscuts);

  taskQA->SetReadMC(readMC);
  taskQA->SetSimpleMode(simplemode); // set to kTRUE to go faster in PbPb
  taskQA->SetTrackOn(trackon);
  taskQA->SetPIDOn(pidon);
  taskQA->SetCentralityOn(centralityon);
  taskQA->SetEvSelectionOn(eventselon);
  taskQA->SetFlowObsOn(flowobson);
  mgr->AddTask(taskQA);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(taskQA,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(out1name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //events analysed
  mgr->ConnectOutput(taskQA,1,coutput1);

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //PID
  if(pidon) mgr->ConnectOutput(taskQA,2,coutput2);

  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(out3name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of tracks
  if(trackon) mgr->ConnectOutput(taskQA,3,coutput3);

  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(out4name,AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  mgr->ConnectOutput(taskQA,4,coutput4);

  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of centrality
  if(centralityon) mgr->ConnectOutput(taskQA,5,coutput5);

  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of centrality
  if(centralityon) mgr->ConnectOutput(taskQA,6,coutput6);

  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(out7name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //event selection
  if(eventselon) mgr->ConnectOutput(taskQA,7,coutput7);

  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer(out8name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //flow observables
  if(flowobson) mgr->ConnectOutput(taskQA,8,coutput8);

 return taskQA;
}

