/* $Id: CreateCuts.C,v 1.5 2008/01/11 08:28:52 jgrosseo Exp $ */

// this macro creates the track and event cuts used in this analysis

// last modified: 2011-03-28 
// m.l.knichel@gsi.de
// added cut modes 200,201: replacing TPCNcluster cut


AliESDtrackCuts* CreatedNdPtTrackCuts(Int_t cutMode=1, Bool_t fieldOn = kTRUE, Bool_t hists = kTRUE)
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

    TString tag = "TPC+ITS refit required - for cut studies";
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
