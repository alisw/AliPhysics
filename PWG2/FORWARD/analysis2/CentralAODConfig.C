/**
 * @file   CentralAODConfig.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 22:10:18 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
/**
 * Configuration script for central multiplicity task.  
 *
 * You can copy this to your working directory or to some other
 * directory up-front in your ROOT macro path, and edit it to suit your
 * needs.
 * 
 * @ingroup pwg2_forward_aod
 */
void
CentralAODConfig(AliCentralMultiplicityTask* task)
{

  // --- Set options on task -----------------------------------------
  // Whether to do correction for secondaries
  task->SetUseSecondary(true);
  // Whether to do correction for acceptance
  task->SetUseAcceptance(true);
  // task->GetInspector().SetDebug(4);

  // task->GetManager().SetSecMapPath(".");
}

//
// EOF
//
