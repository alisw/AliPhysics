AliAnalysisTaskDmesonMCPerform *AddTaskDmesonMCPerform(TString suffix="", 
						       Int_t centOpt=AliRDHFCuts::kCentOff,
						       TString dpluscutfilename="")
{
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMCPerform", "No analysis manager to connect to.");
  }

  AliRDHFCutsDplustoKpipi* analysiscutsdp=0x0;
  if(dpluscutfilename!="") {
    TFile* fileDpcuts=TFile::Open(dpluscutfilename.Data());
    if(!fileDpcuts ||(fileDpcuts&& !fileDpcuts->IsOpen())){
      AliFatal("Input file not found : check your cut object");
    }else{
      analysiscutsdp = (AliRDHFCutsDplustoKpipi*)fileDpcuts->Get("AnalysisCuts");
    }
  }
  
  AliAnalysisTaskDmesonMCPerform *task = new AliAnalysisTaskDmesonMCPerform();
  //  task->SetAODMismatchProtection(-1);
  task->SetUseCentrality(centOpt);
  if(analysiscutsdp) task->SetDplusAnalysisCuts(analysiscutsdp);
  mgr->AddTask(task);
  

  // Create containers for input/output

  TString outname = "coutputDperf";
  outname += suffix.Data();

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGHF_D2H_MCPerform";


  AliAnalysisDataContainer *coutputDmc = mgr->CreateContainer(outname,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputfile.Data());
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutputDmc);

  return task;
}
