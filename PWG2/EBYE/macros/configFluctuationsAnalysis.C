//__________________________________________________//
AliESDtrackCuts *GetTrackCutsObject() {
  //Function to setup the AliESDtrackCuts object and return it
  AliESDtrackCuts *cuts = new AliESDtrackCuts("ebyeTrackCuts","ebyeTrackCuts");
  cuts->SetMinNClustersTPC(80);
  //cuts->SetMinNClustersITS(2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  //cuts->SetRequireITSRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXY(3.0);
  cuts->SetMaxDCAToVertexZ(3.0);
  
  cuts->SetPtRange(0.3,1.5);
  cuts->SetEtaRange(-0.8,0.8);

  cuts->SaveHistograms("trackCuts");

  return cuts;
}
