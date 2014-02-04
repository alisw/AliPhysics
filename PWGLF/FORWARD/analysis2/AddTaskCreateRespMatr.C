/**
 * @file   AddTaskCreateRespMatr.C
 * @author Valentina Zaccolo
 * @date   Fri Jan 11 15:34:26 2013
 * 
 * @brief Script to add a task to create response matrices
 * 
 * 
 */
AliAnalysisTask*
AddTaskCreateRespMatr(const char* trig="V0AND",
                      Double_t vzMin=-4, 
                      Double_t vzMax=4)
{
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardCreateResponseMatrices* task = 
    new AliForwardCreateResponseMatrices("ResponseMatrices");  
  // Set the Vertex Range to Use
  task->SetIpZRange(vzMin, vzMax);
  // Set the Trigger Mask to Use (INEL, NSD, ...)
  task->SetTriggerMask(trig);
  // Set the Number of Bins
  //  task->SetNBins(nBins);

  //Add Full eta-ranges
  task->AddBin(-3.4,5.1);
  
  //Add Symmetric eta bins.
  Double_t limits[] = { 3.4, 3.0, 2.5, 2.4, 2.0, 1.5, 1.4, 1.0, 0.5, 0. };
  Double_t* limit = limits;
  while ((*limit) > 0.1) { 
    task->AddBin(-(*limit), +(*limit));
    // task->AddBin(0,+(*limit));
    // task->AddBin(0,-(*limit));
    limit++;
  }
  // task->AddBin(0,5.0);
  // task->AddBin(0,4.5);
  // task->AddBin(0,4.0); 
  // task->AddBin(0,3.5); 

  // Add 0.5 eta intervals
  // for (Double_t l = -3; l < 5; l += 0.5) task->AddBin(l, l+.5);
 
  // Add 0.20 eta intervals
  // for (Double_t l = -3; l < 5; l += 0.2) task->AddBin(l, l+.2);

  task->Connect();
  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
