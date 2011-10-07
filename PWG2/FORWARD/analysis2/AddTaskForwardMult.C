/** 
 * @defgroup pwg2_forward_scripts Scripts used in the analysis
 *
 * @ingroup pwg2_forward
 */
/** 
 * @defgroup pwg2_forward_scripts_tasks Scripts to add tasks to manager 
 * @ingroup pwg2_forward_scripts
 */
/**
 * @file   AddTaskForwardMult.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:54 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_aod
 */
AliAnalysisTask*
AddTaskForwardMult(Bool_t mc, UShort_t sys=0, UShort_t sNN=0, Short_t field=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWG2forward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardMult", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliForwardMultiplicityBase* task = 0;
  
  if (mc)
    task = new AliForwardMCMultiplicityTask("FMD");
  else    
    task = new AliForwardMultiplicityTask("FMD");
  mgr->AddTask(task);
  
  // --- Do a local initialisation with assumed values ---------------
  if (sys > 0 && sNN > 0) 
    AliForwardCorrectionManager::Instance().Init(sys,sNN,field,mc);

  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_ROOT)/PWG2/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_ROOT)/PWG2/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  const char* config = gSystem->Which(gROOT->GetMacroPath(),
				      "ForwardAODConfig.C");
  if (!config) 
    Warning("AddTaskForwardMult", "ForwardAODConfig.C not found in %s",
	    gROOT->GetMacroPath());
  else {
    Info("AddTaskForwardMult", 
	 "Loading configuration of '%s' from %s",
	 task->ClassName(), config);
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
