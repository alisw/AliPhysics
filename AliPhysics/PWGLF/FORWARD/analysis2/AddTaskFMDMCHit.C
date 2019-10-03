/**
 * @file   AddTaskFMDMCHit.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Mar 23 12:14:03 2011
 * 
 * @brief  Investigate energy loss in simulations 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * This is the macro to include the FMD energy fitter in a train.  
 * 
 * @param useTuple  If true, create NTuple of hits
 * @param debug     Debug level
 *
 * @return Newly created task 
 *
 * @ingroup pwglf_forward_eloss
 */
AliAnalysisTask*
AddTaskFMDMCHit(Bool_t useTuple=false, 
		Int_t  debug=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliFMDDigit",             "FMDbase");
  gROOT->LoadClass("AliFMDHit",               "FMDsim");
  gROOT->LoadClass("AliAODForwardMult",       "PWGLFforward2");
  gROOT->LoadClass("AliFMDMCHitEnergyFitter", "PWGLFforwardhit");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMDELoss", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliFMDMCHitEnergyFitterTask* task = 
    new AliFMDMCHitEnergyFitterTask("ForwardHitELoss", 
				    useTuple);
  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(100, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(15); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(500);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kTRUE);
  // Set whether to make the correction object 
  // task->GetEnergyFitter().SetDoMakeObject(kTRUE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.05);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(1);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data - 10K seems the absolute minimum
  task->GetEnergyFitter().SetMinEntries(3000 /*10000*/);
  // If set, only collect statistics for MB.  This is to prevent a
  // bias when looping over data where the MB trigger is downscaled.
  // task->SetOnlyMB(onlyMB);
  // Debug
  task->SetDebug(debug);

  // --- Set limits on fits the energy -------------------------------
  // DO NOT CHANGE THESE UNLESS YOU KNOW WHAT YOU ARE DOING
  // Maximum relative error on parameters 
  // AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  // AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  // AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;
  
  // --- Make the output container and connect it --------------------
  task->Connect(0,0);
  if (useTuple) { 
    AliAnalysisDataContainer* tuple = 
      mgr->CreateContainer("tuple", TNtuple::Class(), 
			   AliAnalysisManager::kOutputContainer,
			   "forward_tuple.root"
			   /*AliAnalysisManager::GetCommonFileName()*/);
    mgr->ConnectOutput(task, 3, tuple);
  }

  Printf("Returning task %p", task);
  return task;
}
//
// EOF
//
