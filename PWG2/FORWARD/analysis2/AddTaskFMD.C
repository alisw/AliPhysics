/** 
 * @defgroup pwg2_forward_scripts Scripts used in the analysis
 *
 * @ingroup pwg2_forward
 */
/**
 * @file 
 * @ingroup pwg2_forward_scripts
 * 
 */
/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_scripts
 */
AliAnalysisTask*
AddTaskFMD(Bool_t mc, UShort_t sys=0, UShort_t sNN=0, Short_t field=0)
{
  gSystem->Load("libPWG2forward2");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMD", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliForwardMultiplicityBase* task = 0;
  if (mc) task = new AliForwardMCMultiplicityTask("FMD");
  else    task = new AliForwardMultiplicityTask("FMD");
  mgr->AddTask(task);
  
  // --- Do a local initialisation with assumed values ---------------
  if (sys > 0 && sNN > 0) 
    AliForwardCorrectionManager::Instance().Init(sys,sNN,field);

  // --- Configure the task ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2", 
			   gROOT->GetMacroPath()));
  const char* config = gSystem->Which(gROOT->GetMacroPath(),
				      "ForwardAODConfig.C");
  if (!config) 
    Warning("AddTaskFMD", "ForwardAODConfig.C not found in %s",
	    gROOT->GetMacroPath());
  else {
    Info("AddTaskFMD", 
	 "Loading configuration of AliForwardMultiplicityTask from %s",
	 config);
    gROOT->Macro(Form("%s((AliForwardMultiplicityBase*)%p)", config, task));
    delete config;
  }
  
  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  // outputfile += ":PWG2forwardDnDeta"; 
  // Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  return task;
}
