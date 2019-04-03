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
 * This generates a histogram on the output AOD 
 * 
 * @param mc       If true, assume MC input
 * @param runNo    Pre-set run number
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
		   ULong_t     runNo=0,
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

  // --- Set optional corrections path -------------------------------
  AliCentralCorrectionManager& cm = 
    AliCentralCorrectionManager::Instance();
  if (corrs && corrs[0] != '\0') cm.SetPrefix(corrs); 

  // --- Prime the corrections ---------------------------------------
  if(sys>0 && sNN > 0) {
    cm.Init(runNo, sys, sNN, field);
  }

  // --- Make the output container and connect it --------------------
  task->Connect(0,0);
  
  return task;
}
//
// EOF
//
