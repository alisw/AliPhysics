/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track and event cuts used in this analysis

// last modified: 2011-03-28 
// m.l.knichel@gsi.de
// E.PerezLezama@gsi.de
// j.gronefeld@gsi.de
// added cut modes 200,201: replacing TPCNcluster cut


AliESDtrackCuts* CreatedNdPtTrackCuts(Int_t cutMode=1, Bool_t isMC = kFALSE, Double_t scaleChi2 = 1, Bool_t fieldOn = kTRUE, Bool_t hists = kTRUE)
{
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  if (hists)
    esdTrackCuts->DefineHistograms(1);

  Double_t cov1, cov2, cov3, cov4, cov5;
  Double_t nSigma;
  Double_t maxDCAtoVertex, maxDCAtoVertexXY, maxDCAtoVertexZ;
  Double_t minNClustersTPC;
  Double_t maxChi2PerClusterTPC;
  Double_t minPt, maxPt;

  // default cuts for ITS+TPC
  if (cutMode == 0)
  {
    cov1 = 2;
    cov2 = 2;
    cov3 = 0.5;
    cov4 = 0.5;
    cov5 = 2;
    nSigma = 3;
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
    esdTrackCuts->SetMinNsigmaToVertex(nSigma);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag("Global tracking");
  }

  // TPC-only cuts (vertex n sigma cut)
  if (cutMode == 1)
  {
    // beta cuts (still under investigation)
    //cov1 = 4;
    //cov2 = 4;
    cov1 = 2;
    cov2 = 2;
    cov3 = 0.5;
    cov4 = 0.5;
    cov5 = 2;
    nSigma = 4;
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
    esdTrackCuts->SetMinNsigmaToVertex(nSigma);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only cuts (vertex maxDCAtoVertex cut)
  if (cutMode == 2)
  {
    // beta cuts (still under investigation)
    maxDCAtoVertex = 3.0; // cm
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no vertex cuts
  if (cutMode == 3)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no cuts at all
  if (cutMode == 4)
  {

    // beta cuts (still under investigation)
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no chi2
  if (cutMode == 5)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no kink removal
  if (cutMode == 6)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no minNClustersTPC
  if (cutMode == 7)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }
  // TPC-only no kink removal no minNClustersTPC
  if (cutMode == 8)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertex = 3.0; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no minNClustersTPC no maxChi2PerClusterTPC
  if (cutMode == 9)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertex = 3.0; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (loose cuts, absolute DCA cut)
  if (cutMode == 10)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertex = 2.8; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }


  // TPC-only (loose cuts, no DCA cut)
  if (cutMode == 11)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts, no DCA cut)
  if (cutMode == 12)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.2;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (tight cuts, no DCA cut)
  if (cutMode == 13)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.3;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (loose cuts, no pt cut)
  if (cutMode == 14)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts, no pt cut)
  if (cutMode == 15)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (tight cuts, no pt cuts)
  if (cutMode == 16)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }
  // TPC-only (loose cuts)
  if (cutMode == 17)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.6; // cm
    maxDCAtoVertexZ  = 2.1; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts)
  if (cutMode == 18)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 1.8; // cm
    minPt=0.2;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

 // TPC-only (tight cuts)
  if (cutMode == 19)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 1.8; // cm
    minPt=0.3;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts included)
  if (cutMode == 20)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts excluded)
  if (cutMode == 21)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts excluded, no chi2, no DCA)
  if (cutMode == 22)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  // TPC-only
  if (cutMode == 23)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    TString tag = "TPC-only tracking";
  }

  // TPC-only tight cuts
  if (cutMode == 230)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 0.3; // cm
    maxDCAtoVertexZ  = 0.3; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    TString tag = "TPC-only tracking";
  }


  // TPC (no pt cut, no eta cut)
  if (cutMode == 24)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);

    TString tag = "TPC tracking";
  }

  // TPC-only (no pt cut, no eta cut) updated 2011
  if (cutMode == 201)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    Float_t minNCrossedRowsTPC = 120;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);

    esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);

    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking (2011)";
  }

  // TPC multiplicity cuts (test 2013)
  if (cutMode == 203)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    Float_t minNCrossedRowsTPC = 80;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);

    esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);

    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC Multiplicity Cuts (2013)";
  }


  //
  // systematic errors DCA cut studies
  //
  // TPC-only
  if (cutMode == 25)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 2.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 26)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.6; // cm
    maxDCAtoVertexZ  = 2.4; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  //
  // systematic errors cut studies
  //
  // TPC-only
  if (cutMode == 27)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.8; // cm
    maxDCAtoVertexZ  = 2.6; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 28)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.0; // cm
    maxDCAtoVertexZ  = 2.8; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 29)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.2; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 30)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 31)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.6; // cm
    maxDCAtoVertexZ  = 3.4; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }


  if (cutMode == 32)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.8; // cm
    maxDCAtoVertexZ  = 3.6; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 33)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.8; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 34)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.2; // cm
    maxDCAtoVertexZ  = 4.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 35)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.4; // cm
    maxDCAtoVertexZ  = 4.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

//
// cut stability systematics
//

  if (cutMode == 36)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

 if (cutMode == 37)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 90;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 38)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 39)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 5.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 40)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 2.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 41)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.4; // cm
    maxDCAtoVertexZ  = 4.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }

  if (cutMode == 42)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    //esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    TString tag = "TPC-only tracking";
  }
  // test
  if (cutMode == 43)
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    //minPt=0.15;
    //maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    //esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    //esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    //esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    TString tag = "TPC-only tracking";
  }

  // TPC-only + pt cut + eta cut
  if (cutMode == 45)
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    //minPt=0.15;
    //maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    TString tag = "TPC-only tracking";
  }

  // TPC-tracks + SPD point + ITS refit
  if (cutMode == 50)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    //Double_t maxEtaInAcc=0.8;
    Double_t maxdcaxyITSTPC=0.2;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //esdTrackCuts->SetEtaRange(-maxEtaInAcc,maxEtaInAcc);

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster";
  }

  // TPC-tracks + SPD point + ITS refit
  if (cutMode == 60)
  {
    Int_t    minclsITS=4;
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcaxyITSTPC=0.2;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMinNClustersITS(minclsITS);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);

    TString tag = "Global tracking: TPC refit + ITS refit + >3 ITS clusters + >=1 SPD cluster";
  }

  /*
  // TPC-tracks + SPD point + ITS refit + DCAr(pt)
  if (cutMode == 70)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcaxyITSTPC=1.e9;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt)";
  }
  */

  // TPC-tracks + SPD point + ITS refit + DCAr(pt)
  if (cutMode == 70)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt)";
  }

  // TPC+ITS combine tracking + DCAr(pt) + DCAz(pt)
  if (cutMode == 71)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.01+0.011/pt^0.72)
    esdTrackCuts->SetMaxDCAToVertexZPtDep("0.07+0.077/pt^0.72");

    TString tag = "TPC+ITS combine tracking + DCAr(pt) + DCAz(pt)";
  }

  // TPC+ITS combine tracking + DCAr(pt) (2010)
  if (cutMode == 72)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    TString tag = "TPC+ITS combine tracking + DCAr(pt) (2010)";
  }

  // TPC+ITS combine tracking + DCAr(pt) (2011)
  if (cutMode == 200)
  {
    //Int_t    minclsTPC=70;
    Float_t minNCrossedRowsTPC = 120;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);

    //esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    TString tag = "TPC+ITS combine tracking + DCAr(pt) (2011)";
  }

// TPC+ITS combine tracking + DCAr(pt) (2011)
  if (cutMode == 222)
  {
    //Int_t    minclsTPC=70;
    Float_t minNCrossedRowsTPC = 120;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);

    //esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMaxChi2PerClusterITS(36.);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    // tpcc cut
    esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);

    TString tag = "TPC+ITS combine tracking + DCAr(pt) + Chi2TPCcc + Chi2ITS";
  }

// Default track cuts (2015) for 5TeV data
  if (cutMode == 223)
  {

    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);


    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMaxChi2PerClusterITS(36.);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    // tpcc cut
    if(isMC) esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
    else     esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(scaleChi2*36.);

    // Geometrical-Length Cut
    esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7); 
    
    TString tag = "Default track cuts (2015) for 5TeV data";
  }



  // TPC-tracks + SPD point + ITS refit + DCAr(pt) 4-sigma
  if (cutMode == 75)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 4*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.02+0.024/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) 4-sigma";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) 10-sigma
  if (cutMode == 80)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 10*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.05+0.06/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) 10 sigma";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 60 TPCclust
  if (cutMode == 85)
  {
    Int_t    minclsTPC=60;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 80 clusters
  if (cutMode == 90)
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 80 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + TPCchi2=3.5
  if (cutMode == 95)
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=3.5;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + TPCchi2 3.5";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + TPCchi2=4.5
  if (cutMode == 100)
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=4.5;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + TPCchi2 4.5";
  }

  // TPC-tracks
  if (cutMode == 110)
  {

    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e9; // cm
    maxDCAtoVertexZ  = 1.e9; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);

    TString tag = "TPC-tracks loose criteria";
  }


  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 50 TPCclust
  if (cutMode == 120)
  {
    Int_t    minclsTPC=50;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 70 TPCclust + accept kink daughters
  if (cutMode == 130)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 30 TPCclust + accept kink daughters
  if (cutMode == 140)
  {
    Int_t    minclsTPC=30;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // Adam Kisiel track selectiion
  if (cutMode == 150)
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=0.25;
    Double_t maxdcaxyITSTPC=0.2;

    //
    // TPC
    //
    //esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    //esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    //esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    TString tag = "Adam Kisiel track selection";
  }

  // TPC+ITS refit + SPD any
  // for cut studies
  if (cutMode == 151)
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    TString tag = "TPC+ITS refit required - for cut studies";
  }

  // TPC refit
  // for cut studies
  if (cutMode == 152)
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    //esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    TString tag = "TPC refit required - for cut studies";
  }

  // TPC
  // for cut studies
  if (cutMode == 153)
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    //esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    TString tag = "TPC stand alone - for cut studies";
  }

  // TPC+ITS refit
  // for cut studies
  if (cutMode == 154)
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    TString tag = "TPC+ITS refit and KinkRejection required - for cut studies";
  }

  // TPC+ITS refit  + TPC DCA rough cuts
  // for cut studies
  if (cutMode == 155)
  {
    //
    // TPC
    //
    maxDCAtoVertexXY = 5.0; // cm
    maxDCAtoVertexZ  = 5.0; // cm

    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);

    // ITS
    esdTrackCuts->SetRequireITSRefit(kTRUE);
  }

  // Only TPC refit and KinksRemoval required
  if (cutMode == 156)
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    TString tag = "TPC refit + Kink rejection required - for cut studies";
  }


  // TPC+ITS combine tracking + DCAr(pt) (2011)
  if ((cutMode >= 2000) && (cutMode <= 2100))
  {
    //Int_t    minclsTPC=70;
    Float_t minNCrossedRowsTPC = 120;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;
    Double_t maxdaczTPC=3.0;
    Double_t maxdcaxyTPC=3.0;

    //
    // TPC
    //
    if (cutMode >= 2001) { esdTrackCuts->SetRequireTPCRefit(kTRUE); }


    if (cutMode >= 2002) { esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC); }
    if (cutMode >= 2003) { esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC); }
    if (cutMode >= 2004) { esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl); }
    if (cutMode >= 2005) { esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster); }
    if (cutMode >= 2006) { esdTrackCuts->SetMaxDCAToVertexZ(maxdaczTPC); }
    if (cutMode >= 2007) { esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyTPC); }
    //
    // ITS
    //
    if (cutMode >= 2008) { esdTrackCuts->SetRequireITSRefit(kTRUE); }
    if (cutMode >= 2009) { esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); }
    if (cutMode >= 2010) { esdTrackCuts->SetMaxChi2PerClusterITS(36.); }
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    if (cutMode >= 2011) { esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC); }

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    if (cutMode >= 2012) { esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); }
    if (cutMode >= 2013) { esdTrackCuts->SetAcceptKinkDaughters(kFALSE); }

    if (cutMode >= 2014) { esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); }

    TString tag = "for cut/efficiency studies)";
  }

  if ((cutMode >= 3000) && (cutMode <= 3100))
  {
    //Int_t    minclsTPC=70;
    Float_t minNCrossedRowsTPC = 120;
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
    Float_t maxFractionSharedTPCCluster = 0.4;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;
    Double_t maxdaczTPC=3.0;
    Double_t maxdcaxyTPC=3.0;

    //
    // TPC
    //
    if (cutMode >= 3001) { esdTrackCuts->SetRequireTPCRefit(kTRUE); }

    if (cutMode >= 3002) { esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl); }
    if (cutMode >= 3003) { esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC); }
    if (cutMode >= 3004) { esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC); }
    if (cutMode >= 3005) { esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster); }
    if (cutMode >= 3006) { esdTrackCuts->SetMaxDCAToVertexZ(maxdaczTPC); }
    if (cutMode >= 3007) { esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyTPC); }
    //
    // ITS
    //
    if (cutMode >= 3008) { esdTrackCuts->SetRequireITSRefit(kTRUE); }
    if (cutMode >= 3009) { esdTrackCuts->SetMaxChi2PerClusterITS(36.); }
    if (cutMode >= 3010) { esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); }
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    if (cutMode >= 3011) { esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC); }

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    if (cutMode >= 3012) { esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); }
    if (cutMode >= 3013) { esdTrackCuts->SetAcceptKinkDaughters(kFALSE); }

    if (cutMode >= 3014) { esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); }

    TString tag = "for cut/efficiency studies (version 3)";
  }


   if ((cutMode >= 5000) && (cutMode <= 5400))
  {
    //According to 223 for study of systematic uncertanties. Change to fit the default cuts for 5TeV analysis 
    //Just like the 4000 but now with intiger increasing numbers. 
    //Easier to use. 
    //
    // TPC
    //
    TString tag = "Study of systematic uncertanties + Geomtrical cut";
    
    esdTrackCuts->SetRequireTPCRefit(kTRUE); 
    //esdTrackCuts->SetMinNCrossedRowsTPC(120); 
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); 
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4); 
    esdTrackCuts->SetMaxDCAToVertexXY(3.0); 
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE); 
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); 
    esdTrackCuts->SetMaxChi2PerClusterITS(36.);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(2.0); 
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); 
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE); 
    // tpcc cut
    if(isMC) esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
    else     esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(scaleChi2*36.);
    
    // Geometrical-Length Cut
    esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7); 
    
    //
    // Swich Low/High for study of systematics
    //
    if(cutMode==5001){esdTrackCuts->SetMaxChi2PerClusterITS(25.);}
    if(cutMode==5002){esdTrackCuts->SetMaxChi2PerClusterITS(49.);}
    if(cutMode==5003){esdTrackCuts->SetMaxChi2PerClusterTPC(3); }
    if(cutMode==5004){esdTrackCuts->SetMaxChi2PerClusterTPC(5); }
    if(cutMode==5005){esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
    if(cutMode==5006){esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
    if(cutMode==5007){esdTrackCuts->SetMaxFractionSharedTPCClusters(0.2);}
    if(cutMode==5008){esdTrackCuts->SetMaxFractionSharedTPCClusters(1.0);}
    if(cutMode==5009){
      if(isMC) esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
      else     esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(scaleChi2*25.);
    } 
    if(cutMode==5010){
      if(isMC) esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);
      else     esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(scaleChi2*49.);
    }
    if(cutMode==5011){esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");}
    if(cutMode==5012){esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");}
    if(cutMode==5013){esdTrackCuts->SetMaxDCAToVertexZ(1.0); }
    if(cutMode==5014){esdTrackCuts->SetMaxDCAToVertexZ(5.0); }
    if(cutMode==5015){esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff); }
    if(cutMode==5016){esdTrackCuts->SetCutGeoNcrNcl(3,120,1.5,0.85,0.7);}	
    if(cutMode==5017){esdTrackCuts->SetCutGeoNcrNcl(3,140,1.5,0.85,0.7);}	

    if(cutMode==5018){esdTrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);}	//Make a varaition of cut on the width of the dead zone
    if(cutMode==5019){esdTrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);}	// Make a varaition of cut on the width of the dead zone
    
    if(cutMode==5020){esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.80,0.65);}  //Make a variation of cut Nc,Ncl  THE EFFECT IS NEGLIGIBLE
    }

 




  // cuts for data without field
  if (!fieldOn)
  {
    cov5 = 1e10;
    tag += " without field";
  }


//Matching Systematic Uncertainty
  if ((cutMode >= 2100) && (cutMode <= 2199))
  {
    
    //Float_t minNCrossedRowsTPC = 120; 
    Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8; 
    Float_t maxFractionSharedTPCCluster = 0.4;    
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm

    esdTrackCuts->SetRequireTPCRefit(kTRUE);

    //esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);

    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    
    if (cutMode==2100){      
      esdTrackCuts->SetMinNCrossedRowsTPC(120);
      TString tag = "Calculate matching efficiency: TPC only with Crossed Rows";
    }
    
    if (cutMode==2101){
      esdTrackCuts->SetMinNCrossedRowsTPC(120);
      esdTrackCuts->SetRequireITSRefit(kTRUE); 
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); 
      TString tag = "Calculate matching efficiency: TPC + ITS with Crossed Rows";
    }
    
    if (cutMode==2102){
      esdTrackCuts->SetMinNCrossedRowsTPC(120);
      esdTrackCuts->SetRequireITSRefit(kTRUE); 
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff); 
      TString tag = "Calculate matching efficiency: TPC + ITS without SPC hit with Crossed Rows";
    }
    
    if (cutMode==2103){
      esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
      TString tag = "Calculate matching efficiency: Include geometric length cut. TPC only";
    }
    
    if (cutMode==2104){
      esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      TString tag = "Calculate matching efficiency: Include geometric length cut. TPC + ITS";
    }

    if (cutMode==2105){
      esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
      TString tag = "Calculate matching efficiency: Include geometric length cut. TPC + ITS without SPD hit";
    }

    if (cutMode==2106){
      esdTrackCuts->SetCutGeoNcrNcl(0,130,1.5,0.85,0.7);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =0). TPC only";
    }
    
    if (cutMode==2107){
      esdTrackCuts->SetCutGeoNcrNcl(0,130,1.5,0.85,0.7);
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =0). TPC + ITS";
    }

    if (cutMode==2108){
      esdTrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =2). TPC only";
    }
    
    if (cutMode==2109){
      esdTrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =2). TPC + ITS";
    }
   
    if (cutMode==2110){
      esdTrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =4). TPC only";
    }
    
    if (cutMode==2111){
      esdTrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);
      esdTrackCuts->SetRequireITSRefit(kTRUE);
      esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      TString tag = "Calculate matching efficiency: Include geometric length cut(DeadZone =4). TPC + ITS";
    }
  
  
  }


  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}
