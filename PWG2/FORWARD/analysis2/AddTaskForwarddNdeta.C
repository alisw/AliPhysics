/**
 * @file   AddTaskForwarddNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * 
 */
AliAnalysisTask*
AddTaskForwarddNdeta(const char* trig     = "INEL", 
		     Double_t    vzMin    = -10, 
		     Double_t    vzMax    = +10, 
		     Bool_t      useCent  = false,
		     Bool_t      cutEdges = false)
{
  // analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwarddNdetaTask* task = new AliForwarddNdetaTask("Forward");
  task->SetVertexRange(vzMin, vzMax);
  task->SetTriggerMask(trig);
  task->SetCutEdges(cutEdges);
  task->SetUseShapeCorrection(false);
  if (useCent) {
    Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    task->SetCentralityAxis(11, bins);
  }
  mgr->AddTask(task);

  // create containers for input/output
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("ForwardSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  
  // connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
