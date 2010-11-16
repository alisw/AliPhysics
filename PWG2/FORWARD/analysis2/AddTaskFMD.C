/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
AliForwardMultiplicity* 
AddTaskFMD(Int_t nCutBins=1, Float_t correctionCut=0.1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMD", "No analysis manager to connect to.");
    return NULL;
  }   

  AliForwardMultiplicity* task = new AliForwardMultiplicity("FMD");
  task->GetHistCollector().SetNCutBins(nCutBins);
  task->GetHistCollector().SetCorrectionCut(correctionCut);
  mgr->AddTask(task);
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  AliMCEventHandler* mcHandler = 
    dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  Info("AddTaskFMD", "MC handler %p", mcHandler);
  if(mcHandler) {
    pars->SetRealData(kFALSE);
    pars->SetProcessPrimary(kTRUE);
    pars->SetProcessHits(kFALSE);
  }
  else {
    pars->SetRealData(kTRUE);
    pars->SetProcessPrimary(kFALSE);
    pars->SetProcessHits(kFALSE);
  }
  pars->Init();
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);


  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  return task;
}
