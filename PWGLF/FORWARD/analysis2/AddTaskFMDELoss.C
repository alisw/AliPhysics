/**
 * @file   AddTaskFMDELoss.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:14:03 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * @defgroup pwglf_forward_eloss Energy Loss Fits
 *
 * Fitting the energy loss @f$\Delta/\Delta_{mip}@f$ spectra 
 *
 * @ingroup pwglf_forward_topical
 */

/**
 * This is the macro to include the FMD energy fitter in a train.  
 * 
 * @param mc        Assume MC input 
 * @param useCent   Use centrality information 
 * @param onlyMB    Only collect statistics for MB (INEL) events
 * @param debug     Debug level
 * @param residuals If set, also do residuals 
 *
 * @return Newly created task 
 *
 * @ingroup pwglf_forward_eloss
 */
AliAnalysisTask*
AddTaskFMDELoss(Bool_t        mc, 
		Bool_t        useCent,
		Bool_t        onlyMB=false, 
		Int_t         debug=0,
		const Char_t* residuals="")
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMDELoss", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // --- Make the task and add it to the manager ---------------------
  AliFMDEnergyFitterTask* task = new AliFMDEnergyFitterTask("ForwardELoss");
  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(15); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(500);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kTRUE);
  // Set whether to make the correction object 
  task->GetEnergyFitter().SetDoMakeObject(kTRUE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data - 10K seems the absolute minimum
  task->GetEnergyFitter().SetMinEntries(10000);
  // If set, only collect statistics for MB.  This is to prevent a
  // bias when looping over data where the MB trigger is downscaled.
  task->SetOnlyMB(onlyMB);
  // Debug
  task->SetDebug(debug);

  TString resi(residuals);
  resi.ToUpper();
  AliFMDEnergyFitter::EResidualMethod rm = AliFMDEnergyFitter::kNoResiduals;
  if (!resi.IsNull() && !resi.BeginsWith("no")) {
    if (resi.BeginsWith("square")) 
      rm = AliFMDEnergyFitter::kResidualSquareDifference;
    else if (resi.BeginsWith("scale")) 
      rm = AliFMDEnergyFitter::kResidualScaledDifference;
    else // Anything else gives plain difference and errors in errors
      rm = AliFMDEnergyFitter::kResidualDifference;
  }
  task->GetEnergyFitter().SetStoreResiduals(rm);

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

  Printf("Returning task %p", task);
  return task;
}
//
// EOF
//
