/**
 * @file   AddTaskCreateRespMatr.C
 * @author Valentina Zaccolo
 * @date   Fri Jan 11 15:34:26 2013
 * 
 * @brief Script to add a task to create response matrices
 * 
 * 
 */
AliAnalysisTaskSE*
AddTaskCreateRespMatr(const char* trig="V0AND",
                      Double_t vzMin=-4, 
                      Double_t vzMax=4)
{
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardCreateResponseMatrices* task = 
    new AliForwardCreateResponseMatrices("ResponseMatrices");  
  task->SetIpZRange(vzMin, vzMax);
  task->SetTriggerMask(trig);
  task->DefaultBins();
  task->Connect();
  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
