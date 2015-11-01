/**
 * @file   dNdetaConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu May 15 20:33:03 2014
 * 
 * @brief  Configuration file for Forward dN/deta task 
 * 
 * 
 */
void dNdetaConfig(AliBasedNdetaTask* task)
{
  // - Whether to correct for empty bins when projecting 
  // task->SetCorrEmpty(true);
  task->SetCorrEmpty(false); 
  // task->IsA()->InheritsFrom(AliCentraldNdetaTask::Class()));
  // - Whether to use TH2::ProjectionX 
  task->SetUseROOTProjectX(false);

  // --- Special for Hans' analysis ----------------------------------
  // - Set the filename of the corresponding MC analysis
  // const char* mcanalysisfilename = 
  //    "/home/hehi/alex/work/dispVtxDNdeta/mcCorrectionPos.root"
  // task->SetMCFinalCorrFilename(mcanalysisfilename);

  // --- Other things we may overwrite -------------------------------
  // - Set the vertex range to use 
  // task->SetIpZRange(vzMin, vzMax);

  // - Set the trigger mask to use (INEL,INEL>0,NSD)
  // task->SetTriggerMask(trig);

  // - set the rejection mask (PILEUP OUTLIER)
  // task->SetFilterMask(filter);

  // - Set the trigger efficiency 
  // task->SetTriggerEff(trigEff); // 0.997535);
  // task->SetTriggerEff0(trigEff0); 

  // - Set how to normalize. Bit mask of
  // 
  //    kNone           Normalise to accepted events 
  //    kEventLevel     Normalise to all events in selected range 
  //    kAltEventLevel  Normalise to all events in selected range 
  //    kBackground     Also correct for background triggers 
  // 
  //   kNone, kEventLevel, and kAltEventLevel are mutually exclusive.
  //   If neither kEventLevel, nor kAltEventLevel is specified, then
  //   kNone is assumed.  kBackground (when implemented) only makes
  //   sense with kEventLevel and kAltEventLevel.  Furthermore, there
  //   are some constants that encode the common cases
  //     
  //    kFull    = kEventLevel |  kBackground  
  //    kAltFull = kAltEventLevel |  kBackground  
  // 
  //   Default is kFull
  // task->SetNormalizationScheme(scheme);

  // - Whether to use util for pile-up filtering 
  // task->SetUseUtilPileup(false);
  // 
  // If the above is set to true, then one can do 
  // AliAnalysisUtils& au = task->GetAnalysisUtils();
  // 
  // - Use 'multi-vertex' pile-up selection.  If true use track
  //   vertices rather than tracklet vertices
  // au.SetUseMVPlpSelection(false);
  // - Use 'out-of-bunch' pile-up selection
  // au.SetUseOutOfBunchPileUp(false);
  // - If using track pile-up vertex, check for pile-up from other BC
  // au.SetCheckPlpFromDifferentBCMV(false);
  // - Least number of contributors to track pile-up vertex
  // au.SetMinPlpContribMV(5);
  // - Largest chi^2/nu of track pile-up vertex
  // au.SetMaxPlpChi2MV(5);
  // - Least weighted distance between primary and track pile-up
  //   vertex (cm)
  // au.SetMinWDistMV(15.);
  // - Wether to use an adaptive algorithm for the tracklet pile-up
  //   flagging. If this is enabled, the parameters MinPlpContribSPD,
  //   MinPlpZdistSPD, nSigmaPlpZdistSPD, nSigmaPlpDiamXYSPD,
  //   nSigmaPlpDiamZSPD are not used.  Instead the parameters are set
  //   according the the number of tracklets (nTrkL):
  // 
  //   - nTrkL < 20 MinPlpContribSPD=3, 
  //   - nTrkL < 50 MinPlpContribSPD=4
  //   - else       MinPlpContribSPD=5
  // 
  //   and MinPlpZdistSPD=0.8, nSigmaZdistSPD=3, nSigmaDiamXYSPD=2,
  //   nSigmaDiamZSPD=5.
  // au.SetUseSPDCutInMultBins(false);
  // - Least number of contributors to tracklet pile-up vertex
  // au.SetMinPlpContribSPD(5);
  // - Least distance along Z between primary and tracklet pile-up
  //   vertex (cm)
  // au.SetMinPlpZdistSPD(0.8);
  // - Least distance along Z between primary and tracklet pile-up
  //   vertex in terms of the errors on the vertices
  // au.SetnSigmaPlpZdistSPD(3);
  // - Only consider tracklet vertices within this number of times the
  //   sigma in the XY-plane of the interaction diamond.
  // au.SetnSigmaPlpDiamXYSPD(2);
  // - Only consider tracklet vertices within this number of times the
  //   sigma along Z of the interaction diamond.
  // au.SetnSigmaPlpDiamZ(5);

  // - Set the centrality estimator to use 
  // task->SetCentralityMethod(cent);

  // - Set the centrality bins to use.  These are mutually exclusive.
  //   Note, that a bin specified as a-b, covers the interval from a,
  //   inclusive to b exclusive.  An upper bound of 100 is treated
  //   especially, and the upper bound is inclusive in that case .
  // Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  // task->SetCentralityAxis(10, bins);
  // task->SetCentralityAxis("pbpb");
  // task->SetCentralityAxis("ppb");
  // task->SetCentralityAxis("0-5-10-20-30-40-50-60-80-90");

  // - Set Re-weighting function based on IPz
  //
  // TF1* f = new TF1("reweight",
  //                  "TMath::Gaus(x,[0],[1],true)/TMath::Gaus(x,[2],[3],true",
  //                  -10,10,4);
  // f->SetParNames("#mu_{emp},#sigma_{emp},#mu_{this},#sigma_{this}");
  // f->SetParameters(0.592,6.836,muIpz,sigmaIpz);
  // f->SetParErrors(0.023,0.029,eMuIpz,eSigmaIpz);
  // task->SetIpzReweight(f);
  
  // - Set satellite vertex flag
  // task->SetSatelliteVertices(satVtx);
}
// EOF

