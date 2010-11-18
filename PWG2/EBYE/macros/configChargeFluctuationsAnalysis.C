//__________________________________________________//
AliESDtrackCuts *GetTrackCutsObject() {
  //Function to setup the AliESDtrackCuts object and return it
  AliESDtrackCuts *cuts = new AliESDtrackCuts("trackCuts","trackCuts");
  cuts->SetMinNClustersTPC(80);
  cuts->SetMinNClustersITS(2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  cuts->SetRequireITSRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  cuts->SetMaxDCAToVertexZ(0.3);
  
  cuts->SetPtRange(0.3,1.5);
  cuts->SetEtaRange(-0.7,0.7);

  cuts->SaveHistograms("trackCuts");

  return cuts;
}
