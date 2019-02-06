/**
 * @file   AddTaskForwardQA.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:14:03 2011
 * 
 * @brief  Include the Forward QA task in a train.  
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * @defgroup pwglf_forward_qa Quality Assurance
 * 
 * Code to deal with Quality Assurance 
 *
 * @ingroup pwglf_forward_topical
 */
/**
 * This is the macro to include the Forward QA task in a train.  
 * 
 * @param mc       Monte-carlo input 
 * @param useCent  Use centrality 
 *
 * @return newly constructured task object 
 *
 * @ingroup pwglf_forward_eloss
 */
 AliForwardQATask*
AddTaskForwardQA(Bool_t mc=false, Bool_t useCent=false)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Make the task and add it to the manager ---------------------
  AliForwardQATask* task = new AliForwardQATask("forwardQA");

  // --- Cuts --------------------------------------------------------
  // Note, the absolute lowest signal to consider - ever -
  // irrespective of MC or real data, is 0.15.  Signals below this
  // will contain remenants of the pedestal (yes, the width of the
  // pedestal is small, but there are many _many_ channels with only
  // pedestal value in them, so the absolute number of high-value
  // pedestal signals is large - despite the low probablity).
  AliFMDMultCuts cSharingLow(AliFMDMultCuts::kFixed,0.15);
  AliFMDMultCuts cSharingHigh(AliFMDMultCuts::kLandauSigmaWidth,1);
  AliFMDMultCuts cDensity(AliFMDMultCuts::kLandauSigmaWidth,1);

  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Disable use of code from 1st Physics
  task->GetEventInspector().SetUseFirstPhysicsVtx(kFALSE);

  // --- Set parameters on energy loss fitter ------------------------
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
  task->GetEnergyFitter().SetDoMakeObject(kFALSE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.45);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(10000);
  // Set reqularization cut 
  task->GetEnergyFitter().SetRegularizationCut(1e8);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(4);

  // --- Sharing filter ----------------------------------------------
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetUseAngleCorrectedSignals(true);
  // Disable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetZeroSharedHitsBelowThreshold(false);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetHCuts(cSharingHigh);
  // Lower cuts from object
  task->GetSharingFilter().SetLCuts(cSharingLow);
  // Whether to use simple merging algorithm
  task->GetSharingFilter().SetUseSimpleSharing(true);
  // Whether to allow for 3 strip hits - deprecated
  task->GetSharingFilter().SetAllow3Strips(false);

  // --- Density calculator ------------------------------------------
  // Least acceptable quality of ELoss fits
  task->GetDensityCalculator()
    .SetMinQuality(AliFMDCorrELossFit::kDefaultQuality);
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(3);
  // Wet whether to use poisson statistics to estimate N_ch
  task->GetDensityCalculator().SetUsePoisson(true);
  // How to lump for Poisson calculator - 64 strips, 4 regions
  task->GetDensityCalculator().SetLumping(32,4);
  // Set whether or not to include sigma in cut
  task->GetDensityCalculator().SetCuts(cDensity);
  // Set the maximum ratio of outlier bins to the total number of bins
  // task->GetDensityCalculator().SetMaxOutliers(.10);
  task->GetDensityCalculator().SetMaxOutliers(1.0);//Disable filter
  // Set whether or not to use the phi acceptance
  //   AliFMDDensityCalculator::kPhiNoCorrect
  //   AliFMDDensityCalculator::kPhiCorrectNch
  //   AliFMDDensityCalculator::kPhiCorrectELoss
  task->GetDensityCalculator()
    .SetUsePhiAcceptance(AliFMDDensityCalculator::kPhiCorrectNch);

  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  // AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  // AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  // AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;

  // --- Debug -------------------------------------------------------
  // Set the overall debug level (1: some output, 3: a lot of output)
  task->SetDebug(1);
    
  // --- Make the output container and connect it --------------------
  task->Connect(AliAnalysisManager::GetCommonFileName(), "trending.root");

  return task;
}
//
// EOF
// 
