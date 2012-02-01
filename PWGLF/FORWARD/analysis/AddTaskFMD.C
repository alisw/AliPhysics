// This is the macro to include the FMD analysis in the train.  It
// depends on two libraries: libFMDanalysis.so in AliRoot and
// libPhysics.so in ROOT. It has been tested to work with the example
// scripts in the ANALYSIS webpages.  
//
// Author: Hans Hjersing Dalsgaard <hans.dalsgaard@cern.ch>


AliFMDAnalysisTaskSE* AddTaskFMD(Float_t cmsGeV=900.F, 
				 const char* col="p-p", 
				 Float_t bkG=5.F) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskESDFilter", "No analysis manager to connect to.");
    return NULL;
  }   
  
  
  // --- Generate parameter manager ----------------------------------
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->SetEnergy(cmsGeV);
  pars->SetCollisionSystem(col);
  pars->SetMagField(bkG);
  pars->SetProcessPrimary(kFALSE);
  pars->SetProcessHits(kFALSE);
  pars->SetRealData(kTRUE);

  // --- Check if we have an MC handler ------------------------------
  AliMCEventHandler* eventHandler = 
    dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()
				      ->GetMCtruthEventHandler());
  if(eventHandler) {
    pars->SetRealData(kFALSE);
    pars->SetProcessPrimary(kTRUE);
    pars->SetProcessHits(kFALSE);
  }
  pars->Init();
  pars->PrintStatus();
  
  // --- Make the task -----------------------------------------------
  AliFMDAnalysisTaskSE *taskfmd = new AliFMDAnalysisTaskSE("TaskFMD");
  mgr->AddTask(taskfmd);

  // --- Connect input/output ----------------------------------------
  TString outputfile = Form("%s:%s", 
			    AliAnalysisManager::GetCommonFileName(),
			    pars->GetDndetaAnalysisName());

  AliAnalysisDataContainer *cout_fmd = 
    mgr->CreateContainer("BackgroundCorrected", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(taskfmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskfmd, 1, cout_fmd);

  return taskfmd;
}
//
// EOF
//

