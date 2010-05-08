/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track cuts used in this analysis

AliESDtrackCuts* CreateTrackCuts(AliPWG0Helper::AnalysisMode analysisMode, Bool_t hists = kTRUE)
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
  
    TString tag("Global tracking");
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

