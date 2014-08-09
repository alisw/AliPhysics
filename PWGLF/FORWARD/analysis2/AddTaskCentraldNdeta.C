/**
 * @file   AddTaskCentraldNdeta.C
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
 * Create the central @f$ dN/d\eta@f$ analysis task 
 * 
 * @param trig      Trigger to use 
 * @param vzMin     Smallest @f$ v_z@f$
 * @param vzMax     Biggest @f$ v_z@f$
 * @param useCent   Whether to use the centrality or not
 * @param scheme    Normalisation scheme
 * @param cutEdges  Whether to cut edges when rebinning 
 * @param mcanalysisfilename Take final MC corrections from this - if present
 * @param trigEff   Trigger efficiency 
 * @param trigEff0  Trigger efficiency for 0-bin
 * @param corrEmpty Correct for empty bins 
 * 
 * @return Newly created and configured task
 *
 * @ingroup pwglf_forward_dndeta
 */
AliAnalysisTask*
AddTaskCentraldNdeta(const char* config    = "dNdetaConfig.C",
		     const char* trig      = "INEL", 
		     Double_t    vzMin     = -10, 
		     Double_t    vzMax     = +10, 
		     const char* cent      = "",
		     const char* scheme    = 0,
		     Double_t    trigEff   = 1, 
		     Double_t    trigEff0  = 1,
		     Bool_t      satVtx    = false)
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
      !AliForwardUtil::CheckForTask("AliCentralMultiplicityTask")) 
    Fatal("","The relevant task wasn't added to the train");

  // --- Make our object ---------------------------------------------
  AliCentraldNdetaTask* task = new AliCentraldNdetaTask("Central");
  
  // Set the vertex range to use 
  task->SetIpZRange(vzMin, vzMax);
  // Set the trigger mask to use (INEL,INEL>0,NSD)
  task->SetTriggerMask(trig);
  task->SetTriggerEff(trigEff); // 0.997535);
  task->SetTriggerEff0(trigEff0); 

  // Bit mask of 
  // 
  //    kNone               Normalise to accepted events 
  //    kEventLevel         Normalise to all events in selected range 
  //    kBackground         Also correct for background triggers 
  //    kTriggerEfficiency  Correct for trigger efficiency 
  // 
  // kNone and kEventLevel are mutually exclusive.  If kEventLevel is
  // not specified, then kNone is assumed.  kBackground only makes
  // sense with kEventLevel. Furthermore, there
  // are some constants that encode the common cases
  //     
  //    kFull    = kEventLevel |  kBackground | kTriggerEfficiency
  // 
  // Default is kFull
  task->SetNormalizationScheme(AliBasedNdetaTask::kFull);
  if (scheme) task->SetNormalizationScheme(scheme);
  // Set the centrality bins to use.  These are mutually exclusive.
  // Note, that a bin specified as a-b, covers the interval from a,
  // inclusive to b exclusive.  An upper bound of 100 is treated
  // especially, and the upper bound is inclusive in that case .
  if (cent) {
    if (task->SetCentralityMethod(cent)) {
      Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
      task->SetCentralityAxis(10, bins);
    }
  }

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
