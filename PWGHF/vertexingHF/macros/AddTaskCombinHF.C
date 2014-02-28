AliAnalysisTaskCombinHF *AddTaskCombinHF(Int_t meson=0, TString containerStr="",Bool_t readMC=kTRUE, TString cutObjFile="",TString cutObjNam="")
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
      pid=analysiscuts->GetPidHF();
    }
  }
  else {
    if(meson==0) analysiscuts=new AliRDHFCutsD0toKpi();
    else analysiscuts=new AliRDHFCutsDplustoKpipi();
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
  dTask->SetReadMC(readMC);
  dTask->SetDebugLevel(0);
  dTask->SetPIDHF(pid);


  mgr->AddTask(dTask);
  
  // Create containers for input/output 

  TString mesname="Dzero";
  if(meson==1) mesname="Dplus";
  TString inname = Form("cinput%s",mesname.Data());
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
