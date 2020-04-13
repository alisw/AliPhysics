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
  task->SetUseSecondary(false/*true*/);
  // Whether to do correction for acceptance - deprecated
  // The tasks stores the per-event phi acceptance in the overflow bin 
  // and the eta coverage in the underflow bin.  Use these to correct the 
  // data for acceptance. 
  task->SetUseAcceptance(false);
  // Whether to make diagnostics or not - off by default
  // task->SetMakeDiagnostics(true);

  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Least number of constributors to 2nd pile-up vertex
  task->GetEventInspector().SetMinPileupContributors(3);
  // Least distance from primary to 2nd pile-up vertex (cm)
  task->GetEventInspector().SetMinPileupDistance(.8);
  // V0-AND triggered events flagged as NSD 
  // task->GetEventInspector().SetUseV0AndForNSD(false);
  // Use primary vertex selection from 1st physics WG
  // task->GetEventInspector().SetUseFirstPhysicsVtx(true);
  // Use satellite collisions
  // task->GetEventInspector().SetUseDisplacedVertices(true);
  // task->GetEventInspector().SetDebug(4);
}

//
// EOF
//
