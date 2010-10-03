//__________________________________________________//
AliBalance *GetBalanceFunctionObject(const char* analysisLevel = "ESD") {
  //Function to setup the AliProtonAnalysis object and return it
  AliBalance *gBalance = new AliBalance();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetAnalysisType(1);
  gBalance->SetNumberOfBins(18);
  gBalance->SetInterval(-0.9,0.9);

  return gBalance;
}

//__________________________________________________//
AliESDtrackCuts *GetTrackCutsObject() {
  //Function to setup the AliESDtrackCuts object and return it
  AliESDtrackCuts *cuts = new AliESDtrackCuts("bfTrackCuts","bfTrackCuts");
  cuts->SetMinNClustersTPC(80);
  cuts->SetMinNClustersITS(2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  cuts->SetRequireITSRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(0.3);
  
  cuts->SetPtRange(0.3,1.5);
  cuts->SetEtaRange(-0.8,0.8);

  cuts->SaveHistograms("trackCuts");

  return cuts;
}
