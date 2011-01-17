//
// This macro is defined in order to have a unique point
// where the standard cuts are configured, in order to be sure
// that all common parts of the cuts will be defined coherently
//
AliESDtrackCuts QualityCutsTPC()
{
  // create output variable
  AliESDtrackCuts cuts;
  
  // general acceptance/pt cuts
  cuts.SetPtRange ( 0.15, 1.0e+10);
  cuts.SetEtaRange(-0.8 , 0.8);
  
  // DCA cuts
  cuts.SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  cuts.SetMaxDCAToVertexZ(2.0);
  cuts.SetDCAToVertex2D(kFALSE);
  cuts.SetRequireSigmaToVertex(kFALSE);
  
  // TPC related cuts for TPC+ITS tracks
  cuts.SetMinNClustersTPC(70);
  cuts.SetMaxChi2PerClusterTPC(4);
  cuts.SetAcceptKinkDaughters(kFALSE);
  cuts.SetRequireTPCRefit(kTRUE);
  
  // ITS related cuts for TPC+ITS tracks
  cuts.SetRequireITSRefit(kTRUE);
  cuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  
  // finished
  return cuts;
}
