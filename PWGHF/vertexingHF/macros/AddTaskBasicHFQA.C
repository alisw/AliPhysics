AliAnalysisTaskSEHFQA* AddTaskBasicHFQA(AliAnalysisTaskSEHFQA::DecChannel ch, 
					Bool_t readMC, 
					TString suffix=""){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBasicHFQA", "No analysis manager to connect to.");
    return NULL;
  }

  AliRDHFCuts *analysiscuts=0x0;
  switch (ch){
  case 0:
    analysiscuts = new AliRDHFCutsDplustoKpipi();
    suffix="Dplus";
    break;
  case 1:
    analysiscuts = new AliRDHFCutsD0toKpi();
    suffix="D0";
    break;
  case 2:
    analysiscuts = new AliRDHFCutsDstoKKpi();
    suffix="Ds";
    break;
  case 4:
    analysiscuts = new AliRDHFCutsD0toKpipipi();
    suffix="D04";
    break;
  case 5:
    analysiscuts = new AliRDHFCutsLctopKpi();
    suffix="Lc";
    break;
  case 6:
    analysiscuts = new AliRDHFCutsLctoV0();
     suffix="LcToV0x";
    break;
  }
  
  analysiscuts->SetTriggerClass("");
  analysiscuts->SetTriggerMask(AliVEvent::kAnyINT);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff);

  AliAODPidHF* pid4hf=new AliAODPidHF();
  Double_t plim[2]={0.6,0.8};
  Double_t nsigma[5]={2.,1.,2.,3.,0.};
  pid4hf->SetPLimit(plim,2);
  pid4hf->SetAsym(kTRUE);
  pid4hf->SetSigma(nsigma);
  pid4hf->SetMatch(1);
  pid4hf->SetTPC(1);
  pid4hf->SetTOF(1);
  pid4hf->SetITS(0);
  pid4hf->SetTRD(0);
  pid4hf->SetCompat(kTRUE);
  analysiscuts->SetPidHF(pid4hf);

  AliAnalysisTaskSEHFQA* taskQA=new AliAnalysisTaskSEHFQA(Form("QA%s",suffix.Data()),ch,analysiscuts);
  taskQA->SetReadMC(readMC);
  taskQA->SetSimpleMode(kTRUE);
  taskQA->SetTrackOn(kTRUE);
  taskQA->SetPIDOn(kTRUE);
  taskQA->SetCentralityOn(kFALSE);
  taskQA->SetEvSelectionOn(kTRUE);
  taskQA->SetFillDistributionsForTrackEffChecks(kFALSE);
  mgr->AddTask(taskQA);

  TString filename="",out1name="nEntriesQA",out2name="outputPid",out3name="outputTrack",out4name="cuts",out5name="countersCentrality",out6name="outputCentrCheck",out7name="outputEvSel",out8name="outputFlowObs",inname="input",cutsobjname="",centr="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_QA";
  inname+=suffix;
  out1name+=suffix;
  out2name+=suffix;
  out3name+=suffix;
  out4name+=suffix;
  out5name+=suffix;
  out6name+=suffix;
  out7name+=suffix;
  out8name+=suffix;


  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(taskQA,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //events analysed
  mgr->ConnectOutput(taskQA,1,coutput1);

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //PID
  mgr->ConnectOutput(taskQA,2,coutput2);

  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(out3name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //quality of tracks
  mgr->ConnectOutput(taskQA,3,coutput3);

  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(out4name,AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  mgr->ConnectOutput(taskQA,4,coutput4);

  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(out7name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //event selection
  mgr->ConnectOutput(taskQA,7,coutput7);

 return taskQA;
}
