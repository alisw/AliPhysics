/* $Id$ */

// this macro creates the track and event cuts used in this analysis

AliESDtrackCuts* CreateTrackCuts(Bool_t hists = kTRUE)
{
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  if (hists)
    esdTrackCuts->DefineHistograms(1);

  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  //esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,1e99);
  //Printf("WARNING: no cut on 1/pt");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);

  esdTrackCuts->SetMinNsigmaToVertex(3);
  esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
  esdTrackCuts->SetAcceptKingDaughters(kFALSE);

  return esdTrackCuts;
}
