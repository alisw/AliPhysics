// -*- C++ -*-
// $Id$

#include <TCut.h>
#include <TDirectory.h>
#include <TEventList.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TSystem.h>

#include <Math/Functor.h>
#include <TMinuit.h>
#include <TMinuitMinimizer.h>
#include <TRobustEstimator.h>
#include <TStopwatch.h>

#include "AliLog.h"
#include "AliLuminousRegionFit.h"

ClassImp(AliLuminousRegionFit);

Bool_t AliLuminousRegionFit::CutOutliers(TTree *t, Bool_t doCut, Double_t sigmaThreshold,
					 TCut &sel, TVectorD &mu, TMatrixDSym &cov)
{
  if (doCut) {
    Double_t       sigmaCut[3]  = { 0.00, 0.00, 0.00 };
    const Double_t sigmaCut0[3] = { 0.01, 0.01, 6.00 };
    for (Int_t i=0; i<3; ++i) {
      sigmaCut[i] = 0.5*sigmaThreshold*TMath::Sqrt(cov(i,i));
      sigmaCut[i] = TMath::Max(1.0*sigmaThreshold*sigmaCut0[i], sigmaCut[i]);
      sigmaCut[i] = TMath::Min(1.2*sigmaThreshold*sigmaCut0[i], sigmaCut[i]);
    }
    sel = TCut(TString::Format("((xTRKnc%+.8f)/%f)**2 + ((yTRKnc%+.8f)/%f)**2 + ((zTRKnc%+.8f)/%f)**2 < 1",
			       -mu[0], sigmaCut[0],
			       -mu[1], sigmaCut[1],
			       -mu[2], sigmaCut[2]));
  }

  Int_t m = t->Draw("xTRKnc:yTRKnc:zTRKnc", sel, "GOFF");
  if (m < 10)
    return kFALSE;

  // update mu
  mu[0] = TMath::Mean(m, t->GetV1());
  mu[1] = TMath::Mean(m, t->GetV2());
  mu[2] = TMath::Mean(m, t->GetV3());

  // update sigma
  m = t->Draw(TString::Format("(xTRKnc%+.8f)**2:(yTRKnc%+.8f)**2:(zTRKnc%+.8f)**2", -mu[0], -mu[1], -mu[2]), sel, "GOFF");
  cov(0,0) = TMath::Mean(m, t->GetV1());
  cov(1,1) = TMath::Mean(m, t->GetV2());
  cov(2,2) = TMath::Mean(m, t->GetV3());

  // compute covariances
  m = t->Draw(TString::Format("(xTRKnc%+.8f)*(yTRKnc%+.8f):(xTRKnc%+.8f)*(zTRKnc%+.8f):(yTRKnc%+.8f)*(zTRKnc%+.8f)",
			      -mu[0], -mu[1],
			      -mu[0], -mu[2],
			      -mu[1], -mu[2]), sel, "GOFF");
  cov(0,1) = cov(1,0) = TMath::Mean(m, t->GetV1());
  cov(0,2) = cov(2,0) = TMath::Mean(m, t->GetV2());
  cov(1,2) = cov(2,1) = TMath::Mean(m, t->GetV3());

  return kTRUE;
}

Double_t AliLuminousRegionFit::MinuitFunction(const Double_t *par)
{
  // see ATLAS-CONF-2010-027
  TMatrixD sigma(2,2);

  Double_t result = 0.0;
  for (Int_t i=0; i<fN; ++i) {
    sigma(0,0) = par[3]*par[3];
    sigma(1,1) = par[4]*par[4];
    sigma(0,1) = par[3]*par[4]*par[6];
    sigma(1,0) = par[3]*par[4]*par[6];
//     if (!hasVtxCovariance) {
//       const Double_t sigmaV = 500*TMath::Power(vtx_ntrk[i], -1.4/2.)*1e-6*100; // microns -> cm
//       sigma(0,0) += par[9]*par[9]*sigmaV*sigmaV;
//       sigma(1,1) += par[9]*par[9]*sigmaV*sigmaV;
//     } else {
    sigma(0,0) += par[9]*par[9]*fCov[0][i]; //vtx_cxx[i];
    sigma(1,1) += par[9]*par[9]*fCov[2][i]; //vtx_cyy[i];
    sigma(0,1) += par[9]*par[9]*fCov[1][i]; //vtx_cxy[i];
    sigma(1,0) += par[9]*par[9]*fCov[1][i]; //vtx_cxy[i];
    // if (fCov[0][i] > 2e-4 ||
    //     fCov[2][i] > 2e-4)
    //   continue;
    //     }
    Double_t det(0);
    sigma.InvertFast(&det);
    Double_t sum = TMath::Log(TMath::Power(TMath::TwoPi(), 1.5));
    sum += 0.5*TMath::Log(det);
    sum +=     TMath::Log(TMath::Abs(par[5]));
    const Double_t x = fX[0][i] - par[0] - par[7]*(fX[2][i] - par[2]);
    const Double_t y = fX[1][i] - par[1] - par[8]*(fX[2][i] - par[2]);
    const Double_t z = fX[2][i] - par[2];

    sum += 0.5*(x*x*sigma(0,0) + y*y*sigma(1,1) + x*y*sigma(0,1) + x*y*sigma(1,0));
    sum += 0.5*z*z/(par[5]*par[5]);
    result += sum;
  }
  return result;
}

void AliLuminousRegionFit::ComputeMoments(TTree *t,
					  const TCut& sel, const TVectorD &mu, const TMatrixDSym &cov,
					  TVectorD &x, TMatrixDSym &cx, Double_t &llRatio)
{
  fN = t->Draw("xTRKnc:yTRKnc:zTRKnc:covMtxXX:covMtxXY:covMtxYY:chi2", sel, "PARA GOFF");
  for (Int_t i=0; i<3; ++i) {
    fX[i]   = t->GetVal(  i);
    fCov[i] = t->GetVal(3+i);
  }
  fChi2 = t->GetVal(6);
  ROOT::Math::Functor fcn(this, &AliLuminousRegionFit::MinuitFunction, 10);
  TMinuitMinimizer m("minimize", 10);
  m.UseStaticMinuit(kTRUE);
  m.SetFunction(fcn);

  m.SetLimitedVariable( 0, "#mu_x", mu[0], 0.01*TMath::Abs(mu[0]),  -1,   1);
  m.SetLimitedVariable( 1, "#mu_y", mu[1], 0.01*TMath::Abs(mu[1]),  -1,   1);
  m.SetLimitedVariable( 2, "#mu_z", mu[2], 0.01*TMath::Abs(mu[2]), -10,  10);

  const Double_t sigma[3] = {
    TMath::Sqrt(cov(0,0)),
    TMath::Sqrt(cov(1,1)),
    TMath::Sqrt(cov(2,2))
  };
  m.SetLimitedVariable( 3, "#sigma_x", sigma[0], 0.01*sigma[0], 0,   1);
  m.SetLimitedVariable( 4, "#sigma_y", sigma[1], 0.01*sigma[1], 0,   1);
  m.SetLimitedVariable( 5, "#sigma_z", sigma[2], 0.01*sigma[2], 0, 200);

  m.SetLimitedVariable( 6, "C_xy", cov(0,1)/(sigma[0]*sigma[1]), .01, -1, 1);
  m.SetLimitedVariable( 7, "s_x",  0.0, 0.0001, -1.0, 1.0);
  m.SetLimitedVariable( 8, "s_y",  0.0, 0.0001, -1.0, 1.0);
  m.SetLimitedVariable( 9, "k",    1.0, 0.01,    0.5, 5.0);

  m.Minimize();
  m.PrintResults();

  // log likelihood ratio test
  const Double_t llNormal = 0.5*(1.0 + TMath::Log(TMath::TwoPi()));
  llRatio = 2*(m.MinValue() + 3*fN*llNormal)/(fN-m.NDim());

  const Double_t *p = m.X();
  for (Int_t i=0; i<10; ++i)
    x[i] = p[i];

  gMinuit->mnemat(cx.GetMatrixArray(), 10);
}

TGraph* AliLuminousRegionFit::ConfGraph(TGraph *g, const char *name, Int_t color, Int_t marker, Int_t lineWidth)
{
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetLineColor(color);
  g->SetLineWidth(lineWidth);
  g->SetName(name);
  fListSave->Add(g);
  return g;
}
TGraphErrors* AliLuminousRegionFit::ConfGraphErr(TGraphErrors *g, const char *name, Int_t color, Int_t marker)
{
  g->SetMarkerStyle(kFullDotLarge);
  ConfGraph(g, name, color, marker);
  return g;
}

TTree* AliLuminousRegionFit::SkimTTree(TTree *t, const TArrayI *idxOutliers)
{
  const Int_t m = idxOutliers->GetSize();
  TArrayI idx(m);
  TMath::Sort(m, idxOutliers->GetArray(), idx.GetArray(), kFALSE);

  TTree *tNew = t->CloneTree(0);
  tNew->SetDirectory(NULL);
  for (Int_t i=0,j=0, n=t->GetEntries(); i<n; ++i) {
    const Bool_t skip = (m && idxOutliers->At(idx[j]) == i);
    if (skip) {
      j += (j != m-1);
      continue;
    }
    t->GetEntry(i);
    tNew->Fill();
  }
  delete t;

  return tNew;
}

TTree* AliLuminousRegionFit::RemoveOutliersRobust(TTree *t, TVectorD &mu, TMatrixDSym &cov)
{
  const Int_t m = t->GetEntries();
  if (m < 10) {
    delete t; t = NULL;
    return t;
  }

  const Int_t hh = TMath::Nint(0.9*m);
  TRobustEstimator robEst(m, 3, hh);
  Float_t xyz[3] = { 0,0,0 };
  t->SetBranchAddress("xTRKnc", &xyz[0]);
  t->SetBranchAddress("yTRKnc", &xyz[1]);
  t->SetBranchAddress("zTRKnc", &xyz[2]);
  for (Int_t j=0; j<m; ++j) {
    t->GetEntry(j);
    Double_t data[3] = {xyz[0], xyz[1], xyz[2]};
    robEst.AddRow(data);
  }
  t->ResetBranchAddresses();

  robEst.Evaluate();

  t = SkimTTree(t, robEst.GetOuliers()); // typo in ROOT (Ouliers instead of Outliers)
  if (t->GetEntries() < 10) {
    delete t; t = NULL;
    return t;
  }

  mu.ResizeTo(3);
  cov.ResizeTo(3,3);
  robEst.GetMean(mu);
  robEst.GetCovariance(cov);
  return t;
}

TTree* AliLuminousRegionFit::RemoveOutliersOld(TTree *t, TCut &sel, TVectorD &mu, TMatrixDSym &cov)
{
  Bool_t success = kTRUE;
  for (Int_t j=0; j<4 && success; ++j)
    success = CutOutliers(t, j!=0, 7.0, sel, mu, cov);

  if (!success) {
    delete t;
    t = NULL;
  }
  return t;
}

Bool_t AliLuminousRegionFit::DoFit(TString  scanName,
				   Int_t    scanType,
				   Double_t offset,
                                   const TCut& cut,
                                   Int_t    bcSel)
{
  const char* graphNamesMoment[9] = {  "gX", "gY", "gZ",  "gSX", "gSY", "gSZ",  "gCXY", "gCXZ", "gCYZ" };
  const char* graphNamesModel[10] = { "geX","geY","geZ", "geSX","geSY","geSZ", "geCXY", "gesX","gesY","gek" };

  if (fListSave) delete fListSave;
  fListSave = new TList;
  fListSave->SetOwner(kTRUE);

  TGraph *gMoments[9];
  for (Int_t i=0; i<9; ++i)
    gMoments[i] = ConfGraph(new TGraph, graphNamesMoment[i]);

  TGraphErrors *gModel[10];
  for (Int_t i=0; i<10; ++i)
    gModel[i] = ConfGraphErr(new TGraphErrors, graphNamesModel[i]);

  // fTE->SetBranchStatus("*",         1);
  // fTE->SetBranchStatus("timestamp", 1);
  // fTE->SetBranchStatus("*TRK*",     1);
  // fTE->SetBranchStatus("isV0and",   1);
  // fTE->SetBranchStatus("isV0M",     1);
  // fTE->SetBranchStatus("covMtx*",   1);
  // fTE->SetBranchStatus("bx",        1);
  // fTE->SetBranchStatus("chi2",   1);

  UInt_t   timeStamp(0);
  Short_t  ntrksTRKnc(0);
  Bool_t   isV0and(0);
  Bool_t   isV0M(0);
  UInt_t   bcid(0);
  Float_t  chi2(0);
  fTE->SetBranchAddress("timestamp",  &timeStamp);
  fTE->SetBranchAddress("ntrksTRKnc", &ntrksTRKnc);
  fTE->SetBranchAddress("isV0and",    &isV0and);
  fTE->SetBranchAddress("isV0M",      &isV0M);
  fTE->SetBranchAddress("bx",         &bcid);
  fTE->SetBranchAddress("chi2",       &chi2);

  const Int_t n = fTSep->Draw("timeStart:timeEnd:sep", "", "GOFF");
  Printf("n=%d / %lld", n, fTSep->GetEntries());
  const Double_t *timeStart = fTSep->GetV1();
  const Double_t *timeEnd   = fTSep->GetV2();
  const Double_t *sep       = fTSep->GetV3();
  TTree **tTemp = new TTree*[n];
  for (Int_t i=0; i<n; ++i) {
    AliInfoF("CloneTree for sep=%+.3f mm (%.0f - %.0f) (%.0f)", sep[i], timeStart[i], timeEnd[i], timeEnd[i]-timeStart[i]);
    tTemp[i] = fTE->CloneTree(0);
    tTemp[i]->SetDirectory(0);
  }

  const Double_t ttMin = TMath::MinElement(n, timeStart);
  const Double_t ttMax = TMath::MaxElement(n, timeEnd);

  fTE->LoadBaskets(2000*1000*1000);

  TCut cutTime(TString::Format("timestamp>=%.0f && timestamp<=%.0f", ttMin, ttMax));
  Printf("START >>elist");
  if (bcSel >= 0) {
    TCut cutBCID(TString::Format("bx==%d", bcSel));
    cutTime *= cutBCID;
  }
  fTE->Draw(">>elist", cutTime*cut);
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  Printf("STOP >>elist %p %d/%lld", elist, elist->GetN(), fTE->GetEntries());
  elist->Print();

  for (Long64_t j=0, m=elist->GetN(); j<m; ++j) {
    fTE->GetEntry(elist->GetEntry(j));
    if ((j%(100*1000)) == 0)
      Printf("%7lld/%lld %d", j, m, timeStamp);
    for (Long64_t i=0; i<n; ++i) {
      if (timeStamp > timeStart[i] && timeStamp < timeEnd[i])
	tTemp[i]->Fill();
    }
  }

  TVectorD    beamSep(2); // beam separation in X and in Y
  TVectorD    mu(3);      // mean
  TMatrixDSym cov(3);     // covariance matrix
  TVectorD    x(10);      // model parameters
  TMatrixDSym cx(10);     // model covariance matrix
  Double_t    llRatio(0); // LL ratio (goodness of fit)

  TTree *tBeamSpot = new TTree;
  tBeamSpot->SetName("TBeamSpot");
  tBeamSpot->Branch("scanType", &scanType);
  tBeamSpot->Branch("beamSep",   beamSep.GetMatrixArray(), "X/D:Y");
  tBeamSpot->Branch("mu",             mu.GetMatrixArray(), "X/D:Y:Z");
  tBeamSpot->Branch("cov",           cov.GetMatrixArray(), "XX/D:XY:XZ:YX:YY:YZ:YX:YY:YZ");
  tBeamSpot->Branch("modelPar",        x.GetMatrixArray(), "muX/D:muY/D:muZ:sigmaX:sigmaY:sigmaZ:rhoXY:sX:sY:k");
  tBeamSpot->Branch("modelCov",       cx.GetMatrixArray(), "cov[100]/D");
  tBeamSpot->Branch("llRatio",  &llRatio);
  fListSave->Add(tBeamSpot);

  for (Int_t i=0; i<n; ++i) {
    tTemp[i]->SetEstimate(-1);

    TCut sel("1");
//     tTemp[i] = RemoveOutliersRobust(tTemp[i], mu, cov);
    // old version of outlier removal:
#if 1
    tTemp[i] = RemoveOutliersOld(tTemp[i], sel, mu, cov);
    if (NULL == tTemp[i])
      continue;
#endif
    const Int_t idx[3][2] = { {0,1}, {0,2}, {1,2} };
    for (Int_t j=0; j<3; ++j) {
      gMoments[  j]->SetPoint(gMoments[  j]->GetN(), sep[i],
			      mu[j]);
      gMoments[3+j]->SetPoint(gMoments[3+j]->GetN(), sep[i],
			      TMath::Sqrt(cov(j,j)));
      gMoments[6+j]->SetPoint(gMoments[6+j]->GetN(), sep[i],
			      cov(idx[j][0], idx[j][1])/TMath::Sqrt(cov(idx[j][0],idx[j][0])*cov(idx[j][1],idx[j][1])));
    }

    ComputeMoments(tTemp[i], sel, mu, cov, x, cx, llRatio);

    for (Int_t j=0; j<10; ++j) {
      gModel[j]->SetPoint(gModel[j]->GetN(),    sep[i], x[j]);
      gModel[j]->SetPointError(gModel[j]->GetN()-1, 0, TMath::Sqrt(cx(j,j)));
    }

    // conversions from mm -> cm and from um -> cm
    beamSep[0] = (scanType == 0 ? 0.1*sep[i] : 1e-4*offset);
    beamSep[1] = (scanType == 1 ? 0.1*sep[i] : 1e-4*offset);

    tBeamSpot->Fill();
  }

  tBeamSpot->ResetBranchAddresses();

  for (Int_t i=0; i<n; ++i)
    delete tTemp[i];
  delete[] tTemp;

  gSystem->Exec(TString::Format("mkdir -p root/%d", fFillNumber));
  TString saveFileName = TString::Format("root/%d/lumiRegion_", fFillNumber) + scanName;
  saveFileName += (bcSel < 0
		   ? TString(".root")
		   : TString::Format("_bcid%d.root", bcSel));

  TFile *fSave = TFile::Open(saveFileName, "RECREATE");
  fListSave->Write();
  fSave->Write();
  fSave->Close();
  SafeDelete(fListSave);

  fTE->ResetBranchAddresses();
  return kTRUE;
}
