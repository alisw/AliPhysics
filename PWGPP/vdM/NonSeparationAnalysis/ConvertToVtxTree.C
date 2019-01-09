// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>
#include <TArrayI.h>
#include <TMath.h>
// #include <T.h>
// #include <T.h>
// #include <T.h>
// #include <T.h>
// #include <T.h>
#include "AliESDVertex.h"
#include "AliAnalysisTaskVdM.h"

void ConvertToVtxTree() {
  UInt_t run, timestamp, bx;
  UShort_t triggered, ntrkletsS, ntrksTRKnc;
  Bool_t isV0and, isSPI7, isA, isC, isT0, isV0M;
  Float_t xTRKnc, yTRKnc, zTRKnc, covMtxXX, covMtxXY, covMtxYY, covMtxXZ, covMtxYZ, covMtxZZ, chi2;
  Float_t timeAD[2], timeV0[2];
  UShort_t decV0online[2], decADonline[2];
  UShort_t decV0offline[2], decADoffline[2];

  //  TFile *fSave = TFile::Open("root/5533/test.root", "RECREATE");
  TFile *fSave = TFile::Open("test.root", "RECREATE");

  fSave->mkdir("Vertex_Performance");
  fSave->cd("Vertex_Performance");
  TTree *tVtx = new TTree("cOutputVtxESD", "cOutputVtxESD");
  tVtx->Branch("run",        &run);
  tVtx->Branch("timestamp",  &timestamp);
  tVtx->Branch("bx",         &bx);
  tVtx->Branch("triggered",  &triggered);
  tVtx->Branch("ntrkletsS",  &ntrkletsS);
  tVtx->Branch("xTRKnc",     &xTRKnc);
  tVtx->Branch("yTRKnc",     &yTRKnc);
  tVtx->Branch("zTRKnc",     &zTRKnc);
  tVtx->Branch("ntrksTRKnc", &ntrksTRKnc);
  tVtx->Branch("isV0and",    &isV0and);
  tVtx->Branch("isSPI7",     &isSPI7);
  tVtx->Branch("isA",        &isA);
  tVtx->Branch("isC",        &isC);
  tVtx->Branch("isT0",       &isT0);
  tVtx->Branch("isV0M",      &isV0M);
  tVtx->Branch("covMtxXX",   &covMtxXX);
  tVtx->Branch("covMtxXY",   &covMtxXY);
  tVtx->Branch("covMtxYY",   &covMtxYY);
  tVtx->Branch("covMtxXZ",   &covMtxXZ);
  tVtx->Branch("covMtxYZ",   &covMtxYZ);
  tVtx->Branch("covMtxZZ",   &covMtxZZ);
  tVtx->Branch("chi2",       &chi2);
  tVtx->Branch("timeAD",     &timeAD, "C/F:A");
  tVtx->Branch("timeV0",     &timeV0, "C/F:A");
  tVtx->Branch("decADonline",   &decADonline,  "C/S:A");
  tVtx->Branch("decADoffline",  &decADoffline, "C/S:A");
  tVtx->Branch("decV0online",   &decV0online,  "C/S:A");
  tVtx->Branch("decV0offline",  &decV0offline, "C/S:A");

  TFile::Open("282027.root");
  gDirectory->cd("results.root");
  TTree *TE = (TTree*)gDirectory->Get("TE");
#if 0
  TE->SetEstimate(-1);
  const Int_t n = TE->Draw("fEventInfo.fTimeStamp", "", "GOFF");
  TArrayI idx(n);
  //  TMath::Sort(TE->GetEntries(), TE->GetV1(), idx.GetArray(), kFALSE);
  TMath::Sort(n, TE->GetV1(), idx.GetArray(), kFALSE);
#endif
  AliESDVertex *esdVtx = NULL;
  TE->SetBranchAddress("VertexTracks", &esdVtx);

  AliAnalysisTaskVdM::TreeData *td = NULL;
  TE->SetBranchAddress("AliAnalysisTaskVdM::TreeData", &td);

  Double_t covMatrix[6] = { 0 };

  for (Long64_t i=0, n=TE->GetEntries(); i<n; ++i) {
    //    TE->GetEntry(idx[i]);
    TE->GetEntry(i);
    run = 0;
    timestamp = td->fEventInfo.fTimeStamp;
    bx        = td->fEventInfo.fBCID;
    triggered = 0;
    ntrkletsS = 0;
    isV0and   = Bool_t(td->fEventInfo.fClassMask & (1ULL<<0));
    isT0      = Bool_t(td->fEventInfo.fClassMask & (1ULL<<1));
    isSPI7    = 0;
    isA       = 0;
    isC       = 0;
    isV0M     = 0;
    ntrksTRKnc  = esdVtx->GetNContributors();
    xTRKnc   = esdVtx->GetX();
    yTRKnc   = esdVtx->GetY();
    zTRKnc   = esdVtx->GetZ();
    esdVtx->GetCovMatrix(covMatrix);
    covMtxXX = covMatrix[0];
    covMtxXY = covMatrix[1];
    covMtxYY = covMatrix[2];
    covMtxXZ = covMatrix[3];
    covMtxYZ = covMatrix[4];
    covMtxZZ = covMatrix[5];
    chi2     = esdVtx->GetChi2();
    for (Int_t j=0; j<2; ++j) {
      timeAD[j] = td->fADInfo.fTime[j];
      timeV0[j] = td->fV0Info.fTime[j];
      decADonline[j] = td->fADInfo.fDecisionOnline[j];
      decV0online[j] = td->fV0Info.fDecisionOnline[j];
      decADoffline[j] = td->fADInfo.fDecisionOffline[j];
      decV0offline[j] = td->fV0Info.fDecisionOffline[j];
    }
    if ((i%100000)==0)
      Printf("%9lld/%9lld %f %f %f",
             i, n,
             esdVtx->GetX(),
             esdVtx->GetY(),
             esdVtx->GetZ());
    tVtx->Fill();
  }
  fSave->Write();
}
