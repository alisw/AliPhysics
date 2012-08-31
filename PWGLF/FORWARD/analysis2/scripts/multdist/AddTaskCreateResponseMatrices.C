/**
 * @file   AddTaskCentraldNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * 
 */
AliAnalysisTask*
AddTaskCreateResponseMatrices(const char* trig="INEL", Double_t vzMin=-10, Double_t vzMax=10)
{
  // analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardCreateResponseMatrices* task = new AliForwardCreateResponseMatrices("ResponseMatrices");
  //AliForwarddNdetaTask* task = new AliForwarddNdetaTask("Forward");
  
  task->SetVertexRange(vzMin, vzMax);
  task->SetTriggerMask(trig);


  //task->SetTriggerMask(AliAODForwardMult::kMCNSD); // trig);
  
  mgr->AddTask(task);
  //Add Full eta-ranges
  task->AddBin(-3.4,5.1);
  
  //Add Symmetric eta bins.
  Double_t limits[] = { 3.4, 3.0, 2.5, 2.4, 2.0, 1.5, 1.4, 1.0, 0.5, 0. };
  Double_t* limit = limits;
  while ((*limit) > 0.1) { 
    task->AddBin(-(*limit), +(*limit));
    limit++;
  }
  
  //Add 0-<eta> ranges
  task->AddBin(0,5.1);
  task->AddBin(0,5.0);
  task->AddBin(0,4.5);
  task->AddBin(0,4.0);
  task->AddBin(0,3.5);

 limit = limits;
 while ((*limit) > 0.1) { 
   task->AddBin(0, +(*limit));
   limit++;
 }

 limit = limits;
 while ((*limit) > 0.1) { 
   task->AddBin(-(*limit),0);
   limit++;
 }


 //Add 0.5 eta intervals
 for (Double_t l = -3; l < 5; l += 0.5){ 
   task->AddBin(l, l+.5);
 }
 
 //Add 0.20 eta intervals
 for (Double_t l = -3; l < 5; l += 0.2){ 
   task->AddBin(l, l+.2);
 }
  
   
  // create containers for input/output
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("ResponseMatrices", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			AliAnalysisManager::GetCommonFileName());
 /*
  AliAnalysisDataContainer *output = 
  mgr->CreateContainer("CentralResults", TList::Class(), 
  AliAnalysisManager::kParamContainer, 
  AliAnalysisManager::GetCommonFileName());
 */
  // connect input/output
 mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
 mgr->ConnectOutput(task, 1, sums);
 // mgr->ConnectOutput(task, 2, output);
 
 return task;
}

  
//________________________________________________________________________
//
// EOF
// 
