/**
 * @file   AddTaskMultDists.C
 * @author Valentina Zaccolo
 * @date   Thu Nov 22 11:29:26 2012
 * 
 * @brief Script to add a multiplicity task
 * 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
/** 
 * Function to add task to train 
 * 
 * @param trig      Trigger to use 
 * @param vzMin     Least z-coordinate of the interaction point
 * @param vzMax     Largest z-coordinate of the interaction point
 * @param lowCent   Least centrality to consider 
 * @param highCent  Largest centrality to consider 
 * @param nBins     Number of bins to use 
 * 
 * @return Newly allocated task, or null 
 *
 * @deprecated Use AddTaskMultDistribution.C instead 
 *
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
AliAnalysisTask*
AddTaskMultDists(const char* trig     = "V0AND",
		 Double_t    vzMin    = -4,
		 Double_t    vzMax    = 4,
		 Int_t 	     lowCent  = 0,
		 Int_t 	     highCent = 0, 
		 Int_t 	     nBins    = 400)
{
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  AliForwardMultiplicityDistribution* task = 
    new AliForwardMultiplicityDistribution("Mult");
  // Set the Vertex Range to Use
  task->SetIpZRange(vzMin, vzMax);
  // Set the Trigger Mask to Use (INEL, NSD, ...)
  task->SetTriggerMask(trig);
  // Set the Centrality limits
  task->SetCentrality(lowCent, highCent);
  // Set the Number of Bins
  task->SetNBins(nBins);

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
