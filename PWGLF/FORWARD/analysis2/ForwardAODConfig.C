/**
 * @file   ForwardAODConfig.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 13:56:02 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_scripts_tasks
 * 
 */
/**
 * Configuration script for forward multiplicity task.  
 *
 * You can copy this to your working directory or to some other
 * directory up-front in your ROOT macro path, and edit it to suit your
 * needs.
 * 
 * @param task  Task to configure 
 *
 * @ingroup pwglf_forward_aod
 */
void
ForwardAODConfig(AliForwardMultiplicityBase* task)
{
  if (!task) return;

  Info("ForwardAODConfig", "Setting up task %s (%p)", task->GetName(), task);

  // --- General parameters ------------------------------------------
  // Whether to enable low flux specific code 
  task->SetEnableLowFlux(kFALSE);

  // Would like to use dynamic cast but CINT interprets that as a 
  // static cast - sigh!
  Bool_t mc = false;
  if (task->IsA() == AliForwardMCMultiplicityTask::Class()) 
    mc = true;
  
#if 0
  if (mc) {
    AliForwardMCMultiplicityTask* mcTask = 
      static_cast<AliForwardMCMultiplicityTask*>(task);
    mcTask->SetOnlyPrimary(true);
  }
#endif

  Double_t nXi = mc ? 2 : 2;   //HHD test
  Bool_t   includeSigma = false; //true;

  AliFMDMultCuts cSharing;
  //c.SetNXi(mc ? 1 : 1);
  //c.SetIncludeSigma(true);
  //c.SetMPVFraction(0.5);
  Double_t factor = 1.;
  //if(mc) factor = 1.15;
  cSharing.SetMultCuts(0.3, 0.3, 0.3, 0.3, 0.3);
 
  AliFMDMultCuts cDensity;
  //c2.SetNXi(mc ? 1 : 1);
  //  c2.SetIncludeSigma(false);
  //c2.SetMPVFraction(0.5);
  //Double_t factor = 1.2;
  cDensity.SetMultCuts(0.3, 0.3, 0.3, 0.3, 0.3);
  
  
  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Least number of constributors to 2nd pile-up vertex
  task->GetEventInspector().SetMinPileupContributors(3);
  // Least distance from primary to 2nd pile-up vertex (cm)
  task->GetEventInspector().SetMinPileupDistance(.8);
  // V0-AND triggered events flagged as NSD 
  task->GetEventInspector().SetUseV0AndForNSD(true);
  // Use primary vertex selection from 1st physics WG
  // task->GetEventInspector().SetUseFirstPhysicsVtx(true);

  // --- Sharing filter ----------------------------------------------
  // Set the low cut used for sharing - overrides settings in eloss fits
  //  Float_t factor = 1.;
  //if(mc) factor = 1.2;
  //task->GetSharingFilter().SetLowCut(0.3*factor);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetUseAngleCorrectedSignals(true);
  // Disable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetZeroSharedHitsBelowThreshold(false);
  // Whether to use simple merging algorithm
  task->GetSharingFilter().SetUseSimpleSharing(false);
  // Whether to allow for 3 strip hits 
  task->GetSharingFilter().SetAllow3Strips(true);
  // Do not cut fixed/hard on multiplicity 
  task->GetSharingFilter().GetHCuts().SetMultCuts(-1);
  // Set the number of xi's (width of landau peak) to stop at 
  task->GetSharingFilter().GetHCuts().SetNXi(nXi);
  // Set whether or not to include sigma in cut
  task->GetSharingFilter().GetHCuts().SetIncludeSigma(includeSigma);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetLCuts(cSharing);
  
   
  // --- Density calculator ------------------------------------------
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(10);
  // Wet whether to use poisson statistics to estimate N_ch
  task->GetDensityCalculator().SetUsePoisson(true);
  // Set to use the running average in Poisson 
  // task->GetDensityCalculator().SetUseRunningAverage(false);
  // Set whether or not to include sigma in cut
  task->GetDensityCalculator().SetCuts(cDensity);
  // Set lumping (nEta,nPhi)
  task->GetDensityCalculator().SetLumping(32,4);
  // Set whether or not to use the phi acceptance
  //   AliFMDDensityCalculator::kPhiNoCorrect
  //   AliFMDDensityCalculator::kPhiCorrectNch
  //   AliFMDDensityCalculator::kPhiCorrectELoss
  task->GetDensityCalculator()
    .SetUsePhiAcceptance(AliFMDDensityCalculator::kPhiCorrectNch);

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
  // task->SetDebug(0);
  // Set the debug level of a single algorithm 
  // task->GetSharingFilter().SetDebug(3);

  // --- Eventplane Finder -------------------------------------------
  task->GetEventPlaneFinder().SetUsePhiWeights(false);

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
