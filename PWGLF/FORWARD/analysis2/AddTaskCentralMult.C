/**
 * @file   AddTaskCentralMult.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:25 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */

/**
 * This is the macro to include the Central multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_aod
 */
AliAnalysisTask* 
AddTaskCentralMult(Bool_t mc=false, 
		   UShort_t sys=0, UShort_t sNN=0, Short_t field=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWG2forward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCentralMult", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliCentralMultiplicityTask* task = 0;
  if (!mc) task = new AliCentralMultiplicityTask("Central");
  else     task = new AliCentralMCMultiplicityTask("Central");
  if(sys>0 && sNN > 0)
    task->GetManager().Init(sys, sNN, field);
  mgr->AddTask(task);

  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_ROOT)/PWG2/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_ROOT)/PWG2/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  const char* config = gSystem->Which(gROOT->GetMacroPath(),
				      "CentralAODConfig.C");
  if (!config) 
    Warning("AddTaskCentralMult", "CentralAODConfig.C not found in %s",
	    gROOT->GetMacroPath());
  else {
    Info("AddTaskCentralMult", 
	 "Loading configuration of '%s' from %s",
	 task->ClassName(), config);
    gROOT->Macro(Form("%s((AliCentralMultiplicityTask*)%p)", config, task));
    delete config;
  }

  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Central", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);
  
  return task;
}
//
// EOF
//
