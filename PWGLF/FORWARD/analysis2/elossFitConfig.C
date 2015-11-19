/**
 * @file   elossFitConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Aug 14 15:25:13 2014
 * 
 * @brief  Configure Energy loss fitter 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 * 
 */
/** 
 * Configure ther energy loss fitter task 
 * 
 * @param task Task to configure 
 *
 * @ingroup pwglf_forward_eloss
 */
void
elossFitConfig(AliFMDEnergyFitterTask* task)
{
  if (!task) return;

  Info("elossFitConfig", "Setting up task %s (%p)", task->GetName(), task);

  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // How to tag events as pile-up.  Bit pattern of 
  //
  //   - 0x1:      SPD multi-vertex 
  //   - 0x2:      Track multi-vertex 
  //   - 0x4:      Out-of-bunch
  // 
  task->GetEventInspector().SetPileupFlags(0x7);

  // --- ESD Fixer ---------------------------------------------------
  // IF the noise correction is bigger than this, flag strip as dead 
  task->GetESDFixer().SetMaxNoiseCorrection(0.04);
  // Dead region in FMD2i
  task->GetESDFixer().AddDeadRegion(2, 'I', 16, 17, 256, 511);  

  // --- Energy loss fitter ------------------------------------------
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
  task->GetEnergyFitter().SetLowCut(0.45);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data - 10K seems the absolute minimum
  task->GetEnergyFitter().SetMinEntries(10000);
  // Set reqularization cut 
  task->GetEnergyFitter().SetRegularizationCut(1e8);
  // Check if we're to store the residuals.  This can be one of
  // AliFMDEnergyFitter::EResidualMethod:
  //   
  // - AliFMDEnergyFitter::kNoResiduals - no residuals calculated
  // - AliFMDEnergyFitter::kResidualSquareDifference
  // - AliFMDEnergyFitter::kResidualScaledDifference
  // - AliFMDEnergyFitter::kResidualDifference
  // 
  AliFMDEnergyFitter::EResidualMethod rm = AliFMDEnergyFitter::kNoResiduals;
  task->GetEnergyFitter().SetStoreResiduals(rm);

  // --- Set limits on fits the energy -------------------------------
  // DO NOT CHANGE THESE UNLESS YOU KNOW WHAT YOU ARE DOING
  // Maximum relative error on parameters 
  // AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  // AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  // AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;
}
// 
// EOF
// 
