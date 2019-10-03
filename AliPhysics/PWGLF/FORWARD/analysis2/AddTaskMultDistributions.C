/**
 * @defgroup pwglf_forward_scripts_tasks_vz Valentina's code
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/*
 * @file   AddTaskMultDistributions.C
 * @author Valentina Zaccolo
 * @date   Thu Nov 22 11:29:26 2012
 * 
 * @brief Script to add a multiplicity task
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
/** 
 * Script to add a multiplicity task
 * 
 * @param trig     Trigger class 
 * @param vzMin    Least Z-coordiante of interaction point
 * @param vzMax    Largest Z-coordiante of interaction point
 * @param lowCent  Not used yet
 * @param highCent Not used yet 
 * @param nBins    Number of bins 
 * 
 * @return Pointer to task 
 *
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
AliAnalysisTaskSE*
AddTaskMultDistributions(const char* trig = "INEL",
			 Double_t    vzMin = -4,
			 Double_t    vzMax = 4,
			 Int_t       lowCent = 0,
			 Int_t       highCent = 0,
			 Int_t       nBins = 600)
{
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardMultiplicityDistribution* task =
    new AliForwardMultiplicityDistribution("Mult");
  task->SetIpZRange(vzMin, vzMax);
  task->SetTriggerMask(trig);
  task->SetCentralityAxis(lowCent, highCent);
  task->DefaultBins();
  task->Connect();
  return task;
}


//________________________________________________________________________
//
// EOF
// 

