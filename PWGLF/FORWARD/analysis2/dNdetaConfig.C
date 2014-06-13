/**
 * @file   ForwarddNdetaConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu May 15 20:33:03 2014
 * 
 * @brief  Configuration file for Forward dN/deta task 
 * 
 * 
 */
void dNdetaConfig(AliBasedNdetaTask* task)
{
  // - Whether to cut edges when re-binning 
  task->SetCutEdges(false);
  // - Whether to correct for empty bins when projecting 
  // task->SetCorrEmpty(true);
  task->SetCorrEmpty(task->IsA()->InheritsFrom(AliCentraldNdetaTask::Class()));
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

  // - Set the trigger efficiency 
  // task->SetTriggerEff(trigEff); // 0.997535);
  // task->SetTriggerEff0(trigEff0); 

  // - Set how to normalize. Bit mask of
  // 
  //    kNone           Normalise to accepted events 
  //    kEventLevel     Normalise to all events in selected range 
  //    kAltEventLevel  Normalise to all events in selected range 
  //    kBackground     Also correct for background triggers 
  //    kShape          Correct shape 
  // 
  //   kNone, kEventLevel, and kAltEventLevel are mutually exclusive.
  //   If neither kEventLevel, nor kAltEventLevel is specified, then
  //   kNone is assumed.  kBackground (when implemented) only makes
  //   sense with kEventLevel and kAltEventLevel.  Furthermore, there
  //   are some constants that encode the common cases
  //     
  //    kFull    = kEventLevel |  kBackground | kShape 
  //    kAltFull = kAltEventLevel |  kBackground | kShape 
  // 
  //   Default is kFull
  // task->SetNormalizationScheme(scheme);

  // - Set the centrality estimator to use 
  // task->SetCentralityMethod(cent);

  // - Set the centrality bins to use.  These are mutually exclusive.
  //   Note, that a bin specified as a-b, covers the interval from a,
  //   inclusive to b exclusive.  An upper bound of 100 is treated
  //   especially, and the upper bound is inclusive in that case .
  // Short_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  // task->SetCentralityAxis(10, bins);
  
  // - Set satellite vertex flag
  // task->SetSatelliteVertices(satVtx);
}
// EOF

