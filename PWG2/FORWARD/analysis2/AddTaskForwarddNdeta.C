/**
 * @file   AddTaskForwarddNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
/** 
 * Create the Forward @f$ dN/d\eta@f$ analysis task 
 * 
 * @param trig      Trigger to use 
 * @param vzMin     Smallest @f$ v_z@f$
 * @param vzMax     Biggest @f$ v_z@f$
 * @param useCent   Whether to use the centrality or not
 * @param scheme    Normalisation scheme
 * @param cutEdges  Whether to cut edges when rebinning 
 * 
 * @return Newly created and configured task
 *
 * @ingroup pwg2_forward_dndeta
 */
AliAnalysisTask*
AddTaskForwarddNdeta(const char* trig     = "INEL", 
		     Double_t    vzMin    = -10, 
		     Double_t    vzMax    = +10, 
		     Bool_t      useCent  = false,
		     const char* scheme   = 0,
		     Bool_t      cutEdges = false)
{
  // --- Analysis manager --------------------------------------------
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  // --- Fix up physics selection to give proper A,C, and E triggers -
  AliInputEventHandler* ih =
    static_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());

  // Ignore trigger class when selecting events.  This mean that we
  // get offline+(A,C,E) events too
  ps->SetSkipTriggerClassSelection(true);
  
  // --- Make our object ---------------------------------------------
  AliForwarddNdetaTask* task = new AliForwarddNdetaTask("Forward");
  // Set the vertex range to use 
  task->SetVertexRange(vzMin, vzMax);
  // Set the trigger mask to use (INEL,INEL>0,NSD)
  task->SetTriggerMask(trig);
  // Whether to cut edges when re-binning 
  task->SetCutEdges(cutEdges);
  // Bit mask of 
  // 
  //    kNone           Normalise to accepted events 
  //    kEventLevel     Normalise to all events in selected range 
  //    kAltEventLevel  Normalise to all events in selected range 
  //    kBackground     Also correct for background triggers 
  //    kShape          Correct shape 
  // 
  // kNone, kEventLevel, and kAltEventLevel are mutually exclusive.
  // If neither kEventLevel, nor kAltEventLevel is specified, then
  // kNone is assumed.  kBackground (when implemented) only makes
  // sense with kEventLevel and kAltEventLevel.  Furthermore, there
  // are some constants that encode the common cases
  //     
  //    kFull    = kEventLevel |  kBackground | kShape 
  //    kAltFull = kAltEventLevel |  kBackground | kShape 
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
    mgr->CreateContainer("ForwardSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
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
