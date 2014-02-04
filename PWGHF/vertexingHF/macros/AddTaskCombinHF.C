AliAnalysisTaskCombinHF *AddTaskCombinHF(Int_t meson=0, Bool_t readMC=kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCombinHF", "No analysis manager to connect to.");
  }

  //Analysis Task

  AliRDHFCuts* analysiscuts;
  if(meson==0) analysiscuts=new AliRDHFCutsD0toKpi();
  else analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetStandardCutsPP2010();
  
  AliAnalysisTaskCombinHF *dTask = new AliAnalysisTaskCombinHF(meson,analysiscuts);
  dTask->SetReadMC(readMC);
  dTask->SetDebugLevel(0);
  AliAODPidHF* pid=new AliAODPidHF();
  pid->SetMatch(5);
  pid->SetTPCnSigmaRangeForPions(-3.,3.);
  pid->SetTPCnSigmaRangeForKaons(-2.,3.);
  pid->SetTPCnSigmaRangeForProtons(-3.,3.);
  pid->SetTOFnSigmaRangeForPions(-3.,3.);
  pid->SetTOFnSigmaRangeForKaons(-2.,2.);
  pid->SetTOFnSigmaRangeForProtons(-3.,3.);
  dTask->SetPIDHF(pid);
  mgr->AddTask(dTask);
  
  // Create containers for input/output 

  TString mesname="Dzero";
  if(meson==1) mesname="Dplus";
  TString inname = Form("cinput%s",mesname.Data());
  TString outname = Form("coutput%s",mesname.Data());
  TString normname = Form("coutput%sNorm",mesname.Data());

  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),
							  AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG3_D2H_InvMass%sLowPt",mesname.Data());
  
  
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
