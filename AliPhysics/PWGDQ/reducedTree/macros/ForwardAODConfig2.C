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
ForwardAODConfig2(AliForwardMultiplicityBase* task)
{
  if (!task) return;

  Info("ForwardAODConfig", "Setting up task %s (%p)", task->GetName(), task);

  // --- General parameters ------------------------------------------
  // Whether to enable low flux specific code 
  task->SetEnableLowFlux(kFALSE);

  // --- Check for MC ------------------------------------------------
  // Would like to use dynamic cast but CINT interprets that as a 
  // static cast - sigh!
  Bool_t         mc = (task->IsA()==AliForwardMCMultiplicityTask::Class());

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
  
  // --- Event inspector ---------------------------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Least number of constributors to 2nd pile-up vertex -was 3
  task->GetEventInspector().SetMinPileupContributors(1000);
  // Least distance from primary to 2nd pile-up vertex (cm)
  task->GetEventInspector().SetMinPileupDistance(50);
  // V0-AND triggered events flagged as NSD 
  task->GetEventInspector().SetUseV0AndForNSD(false);
  // Set the kind of vertex to look for.  Can be one of 
  //  
  //   - kNormal:    SPD vertex 
  //   - kpA2012:    Selection tuned for 2012 pA data 
  //   - kpA2013:    Selection tuned for 2013 pA data 
  //   - kPWGUD:     Selection used by 'first physics' 
  //   - kDisplaced: Satellite collisions, with kNormal fall-back 
  // 
  task->GetEventInspector().SetVertexMethod(AliFMDEventInspector::kNormal);
  // Which centrality estimator to use 
  task->GetEventInspector().SetCentralityMethod("V0M");
  // How to tag events as pile-up.  Bit pattern of 
  //
  //   - 0x1:      SPD multi-vertex 
  //   - 0x2:      Track multi-vertex 
  //   - 0x4:      Out-of-bunch
  //   - 0x8:      SPD multi-vertex in mult bins 
  // 
  task->GetEventInspector().SetPileupFlags(0xf);

  // --- Event classifier --------------------------------------------
  // Enable/Disable centrality estimation from AliPPVsMultUtils
  task->GetMultEventClassifier().SetUseCentrality(true);
  
  // --- ESD fixer ---------------------------------------------------
  // Sets the noise factor that was used during reconstruction.  If
  // this is set to 4 or more, then this correction will be disabled.
  task->GetESDFixer().SetRecoNoiseFactor(1);
  // IF the noise correction is bigger than this, flag strip as dead 
  task->GetESDFixer().SetMaxNoiseCorrection(0.05);
  // Sets whether to recalculate eta 
  task->GetESDFixer().SetRecalculateEta(false);
  // If true, consider AliESDFMD::kInvalidMult as a zero signal.  This
  // has the unfortunate side effect, that we cannot use the
  // on-the-fly calculation of the phi acceptance.  
  // 
  // *IMPORTANT*
  // 
  // Before revision 43711 of AliFMDReconstructor, all strips with no
  // signal where set to kInvalidMult.  From revision 43711 (Release
  // 4-19-Rev-09) empty strips that are not marked as bad have a
  // signal of 0 (zero).  That means, that for any reconstruction done
  // with releases prior to 4-19-Rev-09 this flag _must_ be defined as
  // true. 
  // 
  // The unfortunate side effect mentioned above is especially cruel
  // in this case, since we would benefit a lot from this setting for
  // that data.  However, one can add dead strips here using
  // AliFMDSharingFilter::AddDeadStrip or
  // AliFMDSharingFilter::AddDeadRegion to remedy the situation, since
  // strips added explicitly here are always ignored.  In the future,
  // the acceptance maker script should generate the list
  // automaticallu.
  //
  // LHC10c-900Gev is effected up-to and including pass3 
  // LHC10c-7TeV is effected up-to and including pass2
  // LHC10c-CPass0 should be OK, but has limited statistics 
  // LHC10c_11a_FMD should be OK, but has few runs  
  task->GetESDFixer().SetInvalidIsEmpty(false);
  // Dead region in FMD2i
  task->GetESDFixer().AddDeadRegion(2, 'I', 16, 17, 256, 511);  
  // One can add extra dead strips from a script like 
  // 
  //   void deadstrips(AliFMDSharingFilter* filter)
  //   {
  //     filter->AddDead(...);
  //     // ... and so on 
  //   }
  //
  // and then do here 
  // 
  // task->GetESDFixer().AddDead("deadstrips.C");

  // --- Sharing filter ----------------------------------------------
  // If the following is set to true, then the merging of shared
  // signals is disabled completely
  // task->GetSharingFilter().SetMergingDisabled(false);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetUseAngleCorrectedSignals(true);
  // Ignore the ESD information when angle correcting.
  // 
  // *IMPORTANT* 
  // 
  // This is to counter a known issue with AliESDFMD with ClassDef 
  // version < 4, where the angle correction flag is incorrectly set.
  // A fix is coming to AliESDFMD to handle it directly in the class. 
  // Only set the flag below to true if you know it to be necessary for
  // your data set.
  task->GetSharingFilter().SetIgnoreESDWhenAngleCorrecting(false);
  // Disable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetZeroSharedHitsBelowThreshold(false);
  // Whether to use simple merging algorithm
  task->GetSharingFilter().SetUseSimpleSharing(true);
  // Whether to allow for 3 strip hits - deprecated
  task->GetSharingFilter().SetAllow3Strips(false);
  // Set upper sharing cut 
  task->GetSharingFilter().SetHCuts(cSharingHigh);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetLCuts(cSharingLow);
   
  // --- Density calculator ------------------------------------------
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(10);
  // Wet whether to use poisson statistics to estimate N_ch
  task->GetDensityCalculator().SetUsePoisson(true);
  // Set whether or not to include sigma in cut
  task->GetDensityCalculator().SetCuts(cDensity);
  // Set lumping (nEta,nPhi)
  task->GetDensityCalculator().SetLumping(32,4);
  // Recalculate phi taking (x,y) offset of IP into account 
  task->GetDensityCalculator().SetRecalculatePhi(true);
  // Least acceptable quality of ELoss fits
  task->GetDensityCalculator()
    .SetMinQuality(AliFMDCorrELossFit::kDefaultQuality);
  // Set the maximum ratio of outlier bins to the total number of bins
  // task->GetDensityCalculator().SetMaxOutliers(.10);
  task->GetDensityCalculator().SetMaxOutliers(1.0);//Disable filter
  // Set the maximum relative diviation between N_ch from Eloss and Poisson
  task->GetDensityCalculator().SetOutlierCut(0.5);
  // Set whether or not to use the phi acceptance
  //   AliFMDDensityCalculator::kPhiNoCorrect
  //   AliFMDDensityCalculator::kPhiCorrectNch
  //   AliFMDDensityCalculator::kPhiCorrectELoss
  task->GetDensityCalculator()
    .SetUsePhiAcceptance(AliFMDDensityCalculator::kPhiCorrectNch);

  // --- Corrector ---------------------------------------------------
  // Whether to use the secondary map correction.  By default we turn
  // off secondary correction for normal data, and on for simulated
  // data.
  task->GetCorrections().SetUseSecondaryMap(mc);
  // Whether to use the vertex bias correction (deprecated)
  task->GetCorrections().SetUseVertexBias(false);
  // Whether to use the acceptance correction from dead-strips (deprecated)
  task->GetCorrections().SetUseAcceptance(false);
  // Whether to use the merging efficiency correction  (deprecated)
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
  //    kSum
  //    kPreferInner
  //    kPreferOuter
  task->GetHistCollector().SetMergeMethod(AliFMDHistCollector::kStraightMean);
  // How to find the fiducial area of the secondary maps 
  // Possible values are 
  //   kByCut    Only bins larger that cut are trusted 
  //   kDistance Only bins that are more than half the size of it neighbors
  task->GetHistCollector().SetFiducialMethod(AliFMDHistCollector::kByCut);
  // Additional diagnostics output - off by default
  // 
  // If this option is enabled, then the summed per-vertex, per-ring
  // d2N/detadphi histograms will be stored in the output, as well as
  // copies of the secondary maps
  task->GetHistCollector().SetMakeBGHitMaps(false);
  //
  // If this option is enabled, then a 3D histogram will be made for
  // each ring, summing dN/deta for each centrality bin.
  task->GetHistCollector().SetMakeCentralitySums(false);

  // --- Eventplane Finder -------------------------------------------
  task->GetEventPlaneFinder().SetUsePhiWeights(false);

  // --- Ring AOD output ---------------------------------------------
  // If set to true, then 5 additional branches will be created on the
  // output AOD - one for each FMD ring.  The branches each contain a
  // TH2D object of the (primary) charged particle multiplicity per
  // (eta,phi)-bin in that event 
  task->SetStorePerRing(true);

  // --- Set limits on fits the energy -------------------------------
  // DO NOT CHANGE THESE UNLESS YOU KNOW WHAT YOU ARE DOING
  // Maximum relative error on parameters 
  // AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  // AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  // AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 10;

  // --- Debug -------------------------------------------------------
  // Set the overall debug level (1: some output, 3: a lot of output)
  // task->SetDebug(0);
  // Set the debug level of a single algorithm 
  // task->GetSharingFilter().SetDebug(3);
}
//
// EOF
//
