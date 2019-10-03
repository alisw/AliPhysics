/**
 * @file   AddTaskForwardMultDists.C
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
 * Create the Forward @f$ dN/d\eta@f$ analysis task.
 * Christian's sketch of a task 
 * 
 * @param trig      Trigger to use 
 * @param vzMin     Smallest @f$ v_z@f$
 * @param vzMax     Biggest @f$ v_z@f$
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
  // Set whether to use stored phi acceptance 
  task->SetUsePhiAcc(usePhiAcc);
  // add the task 
  mgr->AddTask(task);

  // Variable size axis objects
  // for |eta|<0.5 from CMD
  AliForwardMultDists::BinSpec b05(-0.5, 0.5, -0.5);
  b05.Push(21, 1);
  b05.Push(1, 3); 
  b05.Push(1, 5);
  b05.Push(3, 5); // <-- Extra 

  // for |eta|<0.5 from ALICE
  AliForwardMultDists::BinSpec a05(-0.5, 0.5, -0.5);
  a05.Push(21, 1);
  a05.Push(3, 2);

  // For |eta|<1 from CMS 
  AliForwardMultDists::BinSpec b10(-1.0, 1.0, -0.5);
  b10.Push(35, 1);
  b10.Push(1, 3);
  b10.Push(4, 5);
  b10.Push(2, 5); // <-- Extra 

  // For |eta|<1 from ALICE 
  AliForwardMultDists::BinSpec a10(-1.0, 1.0, -0.5);
  a10.Push(41, 1);
  a10.Push(1,  2);

  // For |eta|<1.3 from ALICE 
  AliForwardMultDists::BinSpec a13(-1.3, 1.3, -0.5);
  a13.Push(41, 1);
  a13.Push(7,  2);

  // For |eta|<1.5 from CMS 
  AliForwardMultDists::BinSpec b15(-1.5, 1.5, -0.5);
  b15.Push(46, 1);
  b15.Push(1,  2);
  b15.Push(1,  3); 
  b15.Push(4,  5);
  b15.Push(4,  5); // <-- Extra

  // for |eta|<2.0 from CMS 
  AliForwardMultDists::BinSpec b20(-2.0, 2.0, -0.5);
  b20.Push(62, 1);
  b20.Push(2,  2);
  b20.Push(1,  4);
  b20.Push(3,  10);
  b20.Push(2,  10); // <-- Extra

  // for |eta|<2 from CMS 
  AliForwardMultDists::BinSpec b24(-2.4, 2.4, -0.5);
  b24.Push(57, 1);
  b24.Push(1,  2);
  b24.Push(1,  3); 
  b24.Push(3,  10);
  b24.Push(4,  10); // <-- Extra
  
  // Add bins 
  AliForwardMultDists::BinSpec*  bs[] = { &b05, &b10, &b15, &b20, &b24, 0 };
  AliForwardMultDists::BinSpec** pb   = bs;
  while (*pb) {
    AliForwardMultDists::BinSpec* b = *pb;
    task->AddBin(*b);
    if (useAsymm) { 
      task->AddBin(b->fEtaMin, 0,          b->Axis());
      task->AddBin(0,          b->fEtaMax, b->Axis());
    }
    pb++;
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
