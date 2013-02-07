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
AddTaskCentraldNdeta(const char* trig      = "INEL", 
		     Double_t    vzMin     = -10, 
		     Double_t    vzMax     = +10, 
		     Bool_t      useCent   = false,
		     const char* scheme    = 0,
		     Bool_t      cutEdges  = false,
		     Double_t    trigEff   = 1, 
		     Double_t    trigEff0  = 1,
		     Bool_t      corrEmpty = true,
		     const char* mcanalysisfilename = "none")
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
  task->SetMCFinalCorrFilename(mcanalysisfilename);
  
  // Set the vertex range to use 
  task->SetVertexRange(vzMin, vzMax);
  // Set the trigger mask to use (INEL,INEL>0,NSD)
  task->SetTriggerMask(trig);
  task->SetTriggerEff(trigEff); // 0.997535);
  task->SetTriggerEff0(trigEff0); 
  // Whether to cut edges when re-binning 
  task->SetCutEdges(cutEdges);
  // Whether to correct for empty bins when projecting 
  task->SetCorrEmpty(corrEmpty);
  // Whether to use TH2::ProjectionX 
  task->SetUseROOTProjectX(false);
  // Bit mask of 
  // 
  //    kNone               Normalise to accepted events 
  //    kEventLevel         Normalise to all events in selected range 
  //    kBackground         Also correct for background triggers 
  //    kTriggerEfficiency  Correct for trigger efficiency 
  //    kShape              Correct shape 
  // 
  // kNone and kEventLevel are mutually exclusive.  If kEventLevel is
  // not specified, then kNone is assumed.  kBackground only makes
  // sense with kEventLevel. Furthermore, there
  // are some constants that encode the common cases
  //     
  //    kFull    = kEventLevel |  kBackground | kShape | kTriggerEfficiency
  // 
  // Default is kFull
  task->SetNormalizationScheme(AliBasedNdetaTask::kFull);
  if (scheme) task->SetNormalizationScheme(scheme);
  // Set the centrality bins to use.  These are mutually exclusive.
  // Note, that a bin specified as a-b, covers the interval from a,
  // inclusive to b exclusive.  An upper bound of 100 is treated
  // especially, and the upper bound is inclusive in that case .
  if (useCent) {
    Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    task->SetCentralityAxis(11, bins);
  }
  mgr->AddTask(task);

  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("CentralSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("CentralResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  
  // --- connect input/output ----------------------------------------
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);
  
  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
