/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track and event cuts used in this analysis

AliESDtrackCuts* CreateTrackCuts(Int_t cutMode=1, Bool_t fieldOn = kTRUE, Bool_t hists = kTRUE)
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
    esdTrackCuts->SetMaxDCAToVertex(maxDCAtoVertex);    
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
    esdTrackCuts->SetMaxDCAToVertex(maxDCAtoVertex);    
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
    esdTrackCuts->SetMaxDCAToVertex(maxDCAtoVertex);    
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
    esdTrackCuts->SetMaxDCAToVertex(maxDCAtoVertex);    
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

  // TPC-only + pt cut + eta cut 
  if (cutMode == 23) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    //minPt=0.15;
    //maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    TString tag = "TPC-only tracking";
  }

  // TPC-only (no pt cut, no eta cut)
  if (cutMode == 24) 
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


  // cuts for data without field
  if (!fieldOn)
  {
    cov5 = 1e10;
    tag += " without field";
  }

  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}
