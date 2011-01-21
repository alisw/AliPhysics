//
// This macro is defined in order to have a unique point
// where the standard cuts are configured, in order to be sure
// that all common parts of the cuts will be defined coherently
//
AliESDtrackCuts QualityCutsITS()
{
  // create output variable
  AliESDtrackCuts cuts;
  
  // general acceptance/pt cuts
  cuts.SetPtRange ( 0.15, 1.0e+10);
  cuts.SetEtaRange(-0.8 , 0.8);
  
  // DCA cuts
  cuts.SetMaxDCAToVertexXYPtDep("0.02289+0.03136/pt^1.3");
  cuts.SetMaxDCAToVertexZ(2.0);
  cuts.SetDCAToVertex2D(kFALSE);
  cuts.SetRequireSigmaToVertex(kFALSE);
  
  // ITS related cuts for TPC+ITS tracks
  cuts.SetRequireITSStandAlone(kTRUE);
  cuts.SetRequireITSPureStandAlone(kFALSE);
  cuts.SetRequireITSRefit(kTRUE); 
  cuts.SetMinNClustersITS(4);
  cuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cuts.SetMaxChi2PerClusterITS(3.0);
  
  // finished
  return cuts;
}
