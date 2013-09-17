/**
 * @file   AddTaskForwarddNdeta.C
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
 * Create the Forward @f$ dN/d\eta@f$ analysis task 
 * 
 * @param trig      Trigger to use 
 * @param vzMin     Smallest @f$ v_z@f$
 * @param vzMax     Biggest @f$ v_z@f$
 * @param maxN      Maximum Nch 
 * @param usePhiAcc Use stored phi acceptance 
 * @param useAsymm  Make asymmetric bins 
 *
 * @return Newly created and configured task
 *
 * @ingroup pwglf_forward_dndeta
 */
AliAnalysisTask*
AddTaskForwardMultDists(const char* trig      = "V0AND", 
			Double_t    vzMin     = -4, 
			Double_t    vzMax     = +4,
			UShort_t    maxN      = 150,
			UShort_t    nDiv      = 1,
			Bool_t      usePhiAcc = true,
			Bool_t      useAsymm  = false)
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
  AliForwardMultDists* task = new AliForwardMultDists("Forward");
  
  // Set the vertex range to use 
  task->SetIpZRange(vzMin, vzMax);
  // Set the trigger mask to use (INEL,INEL>0,NSD)
  task->SetTriggerMask(trig);
  // Set the maximum number of charged particles
  task->SetMaxN(maxN);
  // Set number of divisions per particle number 
  task->SetNDivisions(nDiv);
  // Set whether to use stored phi acceptance 
  task->SetUsePhiAcc(usePhiAcc);
  // add the task 
  mgr->AddTask(task);

  // Add bins 
  Double_t limits[] = { 1, 1.5, 2, 3, 0 };
  Double_t* plimit  = limits;
  while (*plimit) {
    Double_t eta = *plimit;
    task->AddBin(-eta, +eta);
    if (useAsymm) { 
      task->AddBin(-eta, 0);
      task->AddBin(0,    +eta);
    }
    plimit++;
  }

  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("ForwardMultSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardMultResults", TList::Class(), 
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
