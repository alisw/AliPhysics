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
AddTaskMultiplicity(const char* trig="INEL", Double_t vzMin=-10, Double_t vzMax=10, Int_t lowCent, Int_t highCent, Int_t nBins)
{
  // analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardMultiplicityDistribution* task = new AliForwardMultiplicityDistribution("Mult");
  //AliForwarddNdetaTask* task = new AliForwarddNdetaTask("Forward");

  task->SetVertexRange(vzMin, vzMax);
  task->SetTriggerMask(trig);
  task->SetCentrality(lowCent, highCent);
  task->SetNBins(nBins);
  mgr->AddTask(task);
  /*
 if (useCent) {
   Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
   task->SetCentralityAxis(11, bins);
   task->InitCentBins();
 }
  */
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
 /*
 task->AddBin(-3.0,-2.5);
 task->AddBin(-2.5,-2.0);
 task->AddBin(-2.0,-1.5);
 task->AddBin(-1.5,-1.0);
 task->AddBin(-1.0,-0.5);
 
 task->AddBin(0.5,1.0);
 task->AddBin(1.0,1.5);
 task->AddBin(1.5,2.0);
 task->AddBin(2.0,2.5);
 task->AddBin(2.5,3.0);
 task->AddBin(3.0,3.5);
 task->AddBin(3.5,4.0);
 task->AddBin(4.0,4.5);
 task->AddBin(4.5,5.0);
 */

 
 
 // create containers for input/output
 AliAnalysisDataContainer *sums = 
   mgr->CreateContainer("CentralSums", TList::Class(), 
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
