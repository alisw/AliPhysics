/**
 * @file   CentralAODConfig.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 22:10:18 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * Configuration script for central multiplicity task.  
 *
 * You can copy this to your working directory or to some other
 * directory up-front in your ROOT macro path, and edit it to suit your
 * needs.
 * 
 * @param task  Task to configure 
 *
 * @ingroup pwglf_forward_aod
 */
void
CentralAODConfig(AliCentralMultiplicityTask* task)
{

  // --- Set options on task -----------------------------------------
  // Whether to do correction for secondaries
  task->SetUseSecondary(true);
  // Whether to do correction for acceptance
  task->SetUseAcceptance(true);

  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetInspector().SetMaxVzErr(0.2);
  // Least number of constributors to 2nd pile-up vertex
  task->GetInspector().SetMinPileupContributors(3);
  // Least distance from primary to 2nd pile-up vertex (cm)
  task->GetInspector().SetMinPileupDistance(.8);
  // V0-AND triggered events flagged as NSD 
  task->GetInspector().SetUseV0AndForNSD(false);
  // Use primary vertex selection from 1st physics WG
  // task->GetInspector().SetUseFirstPhysicsVtx(true);
  // Use satellite collisions
  // task->GetInspector().SetUseDisplacedVertices(true);
  // task->GetInspector().SetDebug(4);

  // task->GetManager().SetSecMapPath(".");
}

//
// EOF
//
