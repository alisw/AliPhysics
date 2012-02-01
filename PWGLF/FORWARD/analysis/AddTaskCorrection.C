// This is the macro to include the FMD analysis in the train.  It
// depends on two libraries: libFMDanalysis.so in AliRoot and
// libPhysics.so in ROOT. It has been tested to work with the example
// scripts in the ANALYSIS webpages.  
//
// Author: Hans Hjersing Dalsgaard <hans.dalsgaard@cern.ch>


AliAnalysisTask* AddTaskCorrection(Float_t blow =0,
				   Float_t bhigh =100,				   
				   Float_t cmsGeV=900.F, 
				   const char* col="p-p", 
				   Float_t bkG=5.F) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskESDFilter", "No analysis manager to connect to.");
    return NULL;
  }   
  
  
  // --- Generate parameter manager ----------------------------------
  /* AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->SetEnergy(cmsGeV);
  pars->SetCollisionSystem(col);
  pars->SetMagField(bkG);
  pars->SetProcessPrimary(kFALSE);
  pars->SetProcessHits(kFALSE);
  pars->SetRealData(kTRUE);
  */
  // --- Check if we have an MC handler ------------------------------
  //AliMCEventHandler* eventHandler = 
  //  dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()
  //				      ->GetMCtruthEventHandler());
  //if(eventHandler) {
  //  pars->SetRealData(kFALSE);
  //  pars->SetProcessPrimary(kTRUE);
  //  pars->SetProcessHits(kFALSE);
  // }
  //pars->Init();
  //pars->PrintStatus();
  const Char_t* outFile = AliAnalysisManager::GetCommonFileName();
  // --- Make the task -----------------------------------------------
  AliFMDAnalysisTaskGenerateCorrection *taskfmd = new AliFMDAnalysisTaskGenerateCorrection("TaskFMD");
  taskfmd->SetBLow(blow);
  taskfmd->SetBHigh(bhigh);
  mgr->AddTask(taskfmd);
  
  // --- Connect input/output ----------------------------------------
  AliAnalysisDataContainer *cout_fmd1 = mgr->CreateContainer("Hits", 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, outFile);
  AliAnalysisDataContainer *cout_fmd2 = mgr->CreateContainer("Primaries", 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, outFile);
  
  AliAnalysisDataContainer *cout_fmd3 = mgr->CreateContainer("vertex", 
							     TH1F::Class(),
							     AliAnalysisManager::kOutputContainer, outFile);
  AliAnalysisDataContainer *cout_fmd4 = mgr->CreateContainer("Correction", 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, outFile);
  
  
  // Dummy AOD output container for jet analysis (no client yet)
  //AliAnalysisDataContainer *c_aodfmd = mgr->CreateContainer("cAODfmd", 
  //						    TTree::Class(),
  //						    AliAnalysisManager::kExchangeContainer);
  // Connect to data containers
  // mgr->ConnectInput  (taskfmd,     0, cin_esd  );
  // mgr->ConnectOutput (fmdana,     0, c_aodfmd );
   mgr->ConnectOutput (taskfmd,     1, cout_fmd1 );
  mgr->ConnectOutput (taskfmd,     2, cout_fmd2 );
  mgr->ConnectOutput (taskfmd,     3, cout_fmd3 );
  mgr->ConnectOutput (taskfmd,     4, cout_fmd4 );

  
  TString outputfile = Form("%s:%s", 
			    AliAnalysisManager::GetCommonFileName(),
			    "Background");

  //  AliAnalysisDataContainer *cout_fmd = 
  // mgr->CreateContainer("BackgroundCorrected", TList::Class(), 
  //		 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(taskfmd, 0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput(taskfmd, 1, cout_fmd);

  return taskfmd;
}
//
// EOF
//

