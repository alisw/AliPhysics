/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track cuts used in this analysis

AliESDtrackCuts* CreateTrackCuts(AliPWG0Helper::AnalysisMode analysisMode,  Bool_t hists = kTRUE,  Float_t ptMin = 0,  Float_t etacut =1e10)
{
  AliESDtrackCuts* esdTrackCuts = 0;
  
  // see https://twiki.cern.ch/twiki/bin/view/ALICE/SelectionOfPrimaryTracksForPp2009DataAnalysis
  
  if (analysisMode & AliPWG0Helper::kTPC)
  {
    esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

    TString tag("TPC-only tracking");

    esdTrackCuts->SetMaxDCAToVertexZ(3.2);
    esdTrackCuts->SetMaxDCAToVertexXY(2.4);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
  
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }

  if (analysisMode & AliPWG0Helper::kTPCITS)
  {
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kTRUE);
    esdTrackCuts->SetPtRange(ptMin);  // adding pt cut
    esdTrackCuts->SetEtaRange(-etacut,etacut);  
    TString tag("Global tracking");
  }
  if ( analysisMode & AliPWG0Helper::kTPCSPD) {

    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kFALSE);
    TString tag("Global tracking+tracklets");

  }

  if (hists)
    esdTrackCuts->DefineHistograms(1);

  // cuts for data without field
  if (!(analysisMode & AliPWG0Helper::kFieldOn))
  {
    tag += " without field";
  }

  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}

