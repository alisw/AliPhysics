/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track and event cuts used in this analysis

AliESDtrackCuts* CreateTrackCuts(AliPWG0Helper::AnalysisMode analysisMode, Bool_t fieldOn = kTRUE, Bool_t hists = kTRUE)
{
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  if (hists)
    esdTrackCuts->DefineHistograms(1);

  // default cuts for ITS+TPC
  Double_t cov1 = 2;
  Double_t cov2 = 2;
  Double_t cov3 = 0.5;
  Double_t cov4 = 0.5;
  Double_t cov5 = 2;
  Double_t nSigma = 4;

  Bool_t tpcRefit = kTRUE;
  Bool_t sigmaToVertex = kTRUE;

  TString tag("Global tracking");

  // TPC-only cuts
  if (analysisMode == AliPWG0Helper::kTPC)
  {
    cov1 = 9;
    cov2 = 9;
    cov3 = 1e10;
    cov4 = 1e10;
    cov5 = 1e10;
    
    tpcRefit = kFALSE;
    sigmaToVertex = kFALSE;
    
    tag = "TPC-only tracking";
  }

  // cuts for data without field
  if (!fieldOn)
  {
    cov5 = 1e10;
    tag += " without field";
  }

  esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);

  esdTrackCuts->SetRequireSigmaToVertex(sigmaToVertex);

  if (sigmaToVertex) {
    esdTrackCuts->SetMaxNsigmaToVertex(nSigma);
  }
  else{
    
    esdTrackCuts->SetMaxDCAToVertexZ(3.2);
    esdTrackCuts->SetMaxDCAToVertexXY(2.4);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
  }

  esdTrackCuts->SetRequireTPCRefit(tpcRefit);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);

  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);

  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}
