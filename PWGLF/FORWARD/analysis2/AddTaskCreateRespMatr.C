/**
 * @file   AddTaskCreateRespMatr.C
 * @author Valentina Zaccolo
 * @date   Fri Jan 11 15:34:26 2013
 * 
 * @brief Script to add a task to create response matrices
 * 
 * This is using Valentina's modified  class 
 *
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
/** 
 * Add task to create response matrix 
 * 
 * @param trig     Trigger class 
 * @param vzMin    Least Z-coordiante of interaction point
 * @param vzMax    Largest Z-coordiante of interaction point
 * 
 * @return Pointer to task 
 *
 * @ingroup pwglf_forward_scripts_tasks_vz
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
