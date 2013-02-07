/**
 * @file   AddTaskCentralMult.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:25 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */

/**
 * This is the macro to include the Central multiplicity in a train.  
 * 
 * @param mc       If true, assume MC input 
 * @param sys      Pre-set collision system
 * @param sNN      Pre-set collition energy
 * @param field    Pre-set magnetic field
 * @param config   Configuration file to use 
 * @param corrs    Corrections to use 
 * 
 * @return Newly created task 
 *
 * @ingroup pwglf_forward_aod
 */
AliAnalysisTask* 
AddTaskCentralMult(Bool_t      mc=false, 
		   UShort_t    sys=0, 
		   UShort_t    sNN=0, 
		   Short_t     field=0, 
		   const char* config="CentralAODConfig.C", 
		   const char* corrs=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCentralMult", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task -----------------------------------------------
  AliCentralMultiplicityTask* task = 0;
  if (!mc) task = new AliCentralMultiplicityTask("Central");
  else     task = new AliCentralMCMultiplicityTask("Central");
  task->Configure(config);
  mgr->AddTask(task);

  // --- Set optional corrections path -------------------------------
  AliCentralMultiplicityTask::Manager& cm = task->GetManager();
  if (corrs && corrs[0] != '\0') { 
    cm->SetAcceptancePath(Form("%s/CentralAcceptance", corrs));
    cm->SetSecMapPath(Form("%s/CentralSecMap", corrs));
  }

  // --- Prime the corrections ---------------------------------------
  if(sys>0 && sNN > 0) {
    cm.Init(sys, sNN, field);
    if (!cm.HasSecondaryCorrection()) 
      Fatal("AddTaskCentralMult", "No secondary correction defined!");
    if (!cm.HasAcceptanceCorrection()) 
      Fatal("AddTaskCentralMult", "No acceptance correction defined!");
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
