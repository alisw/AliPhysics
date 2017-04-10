AliAnalysisTaskCombinHF *AddTaskCombinHF(Int_t meson = 0,
                                         Bool_t readMC = kTRUE,
                                         TString containerStr = "",
                                         TString cutObjFile = "",
                                         TString cutObjNam = "",
                                         Int_t filterMask = 1,
                                         Double_t ptcut = 0.1,
                                         Double_t etacut = 0.9,
                                         Int_t pidStrategy = 0,
                                         Int_t casePID = 0,
                                         Double_t bayesThresKaon = 0.4,
                                         Double_t bayesThresPion = 0.4,
                                         Double_t minMass = 1.6,
                                         Double_t maxMass = 2.15,
                                         Double_t maxPt = 20.,
                                         Double_t ptBinWidth = 0.5)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCombinHF", "No analysis manager to connect to.");
  }
  
  //Analysis Task
  
  
  AliRDHFCuts* analysiscuts=0x0;
  AliAODPidHF* pid=0x0;
  if(!cutObjFile.IsNull()){
    TFile *f=TFile::Open(cutObjFile.Data(),"READ");
    if(f){
      analysiscuts=(AliRDHFCuts*)f->Get(cutObjNam.Data());
      AliESDtrackCuts *esdc=analysiscuts->GetTrackCuts();
      pid=analysiscuts->GetPidHF();
    }
  }
  else {
    if(meson==0) analysiscuts=new AliRDHFCutsD0toKpi();
    else if(meson==1) analysiscuts=new AliRDHFCutsDplustoKpipi();
    else if(meson==3) analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetStandardCutsPP2010();
    pid=new AliAODPidHF();
    pid->SetMatch(5);
    pid->SetTPCnSigmaRangeForPions(-3.,3.);
    pid->SetTPCnSigmaRangeForKaons(-2.,3.);
    pid->SetTPCnSigmaRangeForProtons(-3.,3.);
    pid->SetTOFnSigmaRangeForPions(-3.,3.);
    pid->SetTOFnSigmaRangeForKaons(-2.,2.);
    pid->SetTOFnSigmaRangeForProtons(-3.,3.);
    
  }
  if(!analysiscuts){
    Printf("Wrong file or cut object name set");
    return 0x0;
  }
  
  AliAnalysisTaskCombinHF *dTask = new AliAnalysisTaskCombinHF(meson,analysiscuts);
  if(!cutObjFile.IsNull()){
    dTask->SetKaonTrackCuts(esdc);
    dTask->SetPionTrackCuts(esdc);
    dTask->SetTrackCuts(esdc);
  }
  dTask->SetReadMC(readMC);
  dTask->SetDebugLevel(0);
  
  dTask->SetFilterMask(filterMask);
  dTask->SetPtAccCut(ptcut);
  dTask->SetEtaAccCut(etacut);
  
  // mass and pt range for histograms
  dTask->SetMassWindow(minMass, maxMass);
  dTask->SetMaxPt(maxPt);
  dTask->SetPtBinWidth(ptBinWidth);
  
  // PID settings
  dTask->SetPIDHF(pid);
  dTask->SetPIDstrategy(pidStrategy);
  dTask->SetPIDselCaseZero(casePID);
  dTask->SetBayesThres(bayesThresKaon, bayesThresPion);
  
  
  mgr->AddTask(dTask);
  
  
  // Create containers for input/output
  
  TString mesname="Dzero";
  if(meson==1) mesname="Dplus";
  else if(meson==3) mesname="Ds";
  TString inname = Form("cinput%s%s",mesname.Data(),containerStr.Data());
  TString outname = Form("coutput%s%s",mesname.Data(),containerStr.Data());
  TString normname = Form("coutput%sNorm%s",mesname.Data(),containerStr.Data());
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG3_D2H_InvMass%sLowPt%s",mesname.Data(),containerStr.Data());
  
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               outputfile.Data());
  
  mgr->ConnectInput(dTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dTask,1,coutput);
  
  mgr->ConnectOutput(dTask,2,coutputNorm);  
  
  return dTask;
}
