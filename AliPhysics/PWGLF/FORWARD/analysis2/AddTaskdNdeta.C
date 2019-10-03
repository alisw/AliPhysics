/**
 * @file   AddTaskdNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * Create a @f$ dN/d\eta@f$ analysis task 
 * 
 * 
 * @param config Configuration script  
 * @param which  Which type of task to add (Forward, Central, or MCTruth)
 * 
 * @return Newly created and configured task
 *
 * @ingroup pwglf_forward_dndeta
 */
AliAnalysisTask*
AddTaskdNdeta(const char* which    = "Forward",
	      const char* config   = "dNdetaConfig.C")
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Analysis manager --------------------------------------------
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  // --- Check that we have an AOD input handler ---------------------
  UShort_t aodInput = 0;
  if (!(aodInput = AliForwardUtil::CheckForAOD())) 
    Fatal("","Cannot proceed without and AOD handler");
  if (aodInput == 2 &&
      !AliForwardUtil::CheckForTask("AliForwardMultiplicityBase")) 
    Fatal("","The relevant task wasn't added to the train");


  // --- Make our object ---------------------------------------------
  AliBasedNdetaTask*    task = 0;
  TString w(which);
  if (w.EqualTo("Forward",TString::kIgnoreCase))
    task = new AliForwarddNdetaTask("Forward");
  else if (w.EqualTo("Central",TString::kIgnoreCase))
    task = new AliCentraldNdetaTask("Central");
  else if (w.EqualTo("MCTruth", TString::kIgnoreCase))
    task = new AliMCTruthdNdetaTask("MCTruth");
  else
    Fatal("AddTaskdNdeta.C", "Unknown dN/deta task: %s", which);

  // Set-up task using a script 
  task->Configure(config);
  
  // Connect to manager 
  task->Connect(0,0);

  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
