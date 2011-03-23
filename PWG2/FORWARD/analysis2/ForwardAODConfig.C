/**
 * @file   ForwardAODConfig.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 13:56:02 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_scripts_tasks
 * 
 */
/**
 * Configuration script for forward multiplicity task.  
 *
 * You can copy this to your working directory or to some other
 * directory up-front in your ROOT macro path, and edit it to suit your
 * needs.
 * 
 * @ingroup pwg2_forward_aod
 */
void
ForwardAODConfig(AliForwardMultiplicityBase* task)
{
  if (!task) return;

  Info("ForwardAODConfig", "Setting up task %s (%p)", task->GetName(), task);

  // --- General parameters ------------------------------------------
  // Whether to enable low flux specific code 
  task->SetEnableLowFlux(kFALSE);

  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  
  // --- Energy Loss Fitter ------------------------------------------
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(10); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(300);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kFALSE);
  // Set whether to make the correction object 
  task->GetEnergyFitter().SetDoMakeObject(kFALSE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(1000);

  // --- Sharing filter ----------------------------------------------
  // Set the low cut used for sharing - overrides settings in eloss fits
  task->GetSharingFilter().SetLowCut(0.3);
  // Set the number of xi's (width of landau peak) to stop at 
  task->GetSharingFilter().SetNXi(1);

  // --- Density calculator 
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(3);
  // Wet whether to use poisson statistics to estimate N_ch
  task->GetDensityCalculator().SetUsePoisson(false);
  // Set the lower multiplicity cut.  Overrides setting in energy loss fits.
  task->GetDensityCalculator().SetMultCut(0.3); //was 0.3

  // --- Corrector ---------------------------------------------------
  // Whether to use the secondary map correction
  task->GetCorrections().SetUseSecondaryMap(true);
  // Whether to use the vertex bias correction
  task->GetCorrections().SetUseVertexBias(false);
  // Whether to use the vertex bias correction
  task->GetCorrections().SetUseAcceptance(true);
  // Whether to use the merging efficiency correction 
  task->GetCorrections().SetUseMergingEfficiency(false);

  // --- Histogram Collector -----------------------------------------
  // Set the number of extra bins (beyond the secondary map border) 
  task->GetHistCollector().SetNCutBins(2);
  // Set the correction cut, that is, when bins in the secondary map 
  // is smaller than this, they are considered empty 
  task->GetHistCollector().SetCorrectionCut(0.5);
  // How to calculate the value of overlapping bins. 
  // Possible values are 
  //    kStraightMean 
  //    kStraightMeanNoZero 
  //    kWeightedMean 
  //    kLeastError 
  task->GetHistCollector().SetMergeMethod(AliFMDHistCollector::kStraightMean);
  // How to find the fiducial area of the secondary maps 
  // Possible values are 
  //   kByCut    Only bins larger that cut are trusted 
  //   kDistance Only bins that are more than half the size of it neighbors
  task->GetHistCollector().SetFiducialMethod(AliFMDHistCollector::kByCut);

  // --- Debug -------------------------------------------------------
  // Set the overall debug level (1: some output, 3: a lot of output)
  task->SetDebug(0);
  // Set the debug level of a single algorithm 
  // task->GetEventInspector().SetDebug(4);

  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;
}
//
// EOF
//
