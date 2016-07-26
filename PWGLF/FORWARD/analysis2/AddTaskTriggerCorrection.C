/**
 * @file   AddTaskTriggerCorrection.C
 * @author Valentina Zaccolo
 * @date   Mon Feb 3 13:56:26 2014
 * 
 * @brief Script to add a task to create trigger bias correction
 * 
 * 
 */
AliAnalysisTaskSE*
AddTaskTriggerCorrection(const char* trig="INEL",
			 Double_t vzMin=-4,
			 Double_t vzMax=4)
{

  // Make our object.  2nd argumenent is absolute max Eta 3rd argument
  // is absolute max Vz
  AliForwardTriggerBiasCorrection* task =
    new AliForwardTriggerBiasCorrection("TriggerCorrection");
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

