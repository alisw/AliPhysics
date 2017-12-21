// -*- C++ -*-

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TParameter.h>
#include <TEventList.h>
#include <TEntryList.h>

#include "AliLog.h"
#include "AliDoubleGaussianBeamProfile.h"
#include "AliNonseparationModelFit.h"

ClassImp(AliNonseparationModelFit);

AliNonseparationModelFit::AliNonseparationModelFit()
  : TObject()
  , fMoments(6)
  , fRates(6)
  , fPar(29)
  , fNRates(0)
  , fFitToRates(kFALSE)
  , fFixCrossingAngles(!kTRUE)
  , fScaleRateError(1)
  , fFitToOffsetScans(kTRUE)
  , fMuOffsetsX(3)
  , fMuOffsetsY(3)
  , fChi2Moments(6)
  , fChi2Rates(6)
  , fNDFMoments(6)
  , fNDFRates(6)
  , fMinimizer("minimize", 29)
  , fFcn(this, &AliNonseparationModelFit::MinuitFunction, 29)
{
  fMoments.SetOwner(kTRUE);
  fRates.SetOwner(kTRUE);

  fPar         = 0.0;
  fMuOffsetsX  = 0.0;
  fMuOffsetsY  = 0.0;
  fChi2Moments = 0.0;
  fChi2Rates   = 0.0;
  fNDFMoments  = 0.0;
  fNDFRates    = 0.0;

  ROOT::Math::MinimizerOptions opt = fMinimizer.Options();
  fMinimizer.SetPrintLevel(5);
  fMinimizer.SetOptions(opt);
  fMinimizer.UseStaticMinuit(kTRUE);
  fMinimizer.SetFunction(fFcn);
}

AliNonseparationModelFit::~AliNonseparationModelFit() {
  //
}

void AliNonseparationModelFit::Add(Int_t idx, TTree *tMoments, const TCut& cut, TGraph *gRate) {
  const TString eventListName = TString::Format("elist_%d", idx);
  tMoments->Draw(">>"+eventListName, cut, "GOFF");
  TEventList *elist = NULL;
  gDirectory->GetObject(eventListName, elist);
  elist->Print();
  tMoments->SetEventList(elist);
  fMoments.AddAt(tMoments, idx);

  fRates.AddAt(gRate, idx);
  fNRates += (gRate != 0);
}

Double_t AliNonseparationModelFit::GetNDFMoments() const {
  Double_t sum=0;
  for (Int_t i=0; i<6; ++i)
    sum += fNDFMoments[i];
  return sum;
}
Double_t AliNonseparationModelFit::GetNDFRates() const {
  Double_t sum=0;
  for (Int_t i=0; i<6; ++i)
    sum +=  fNDFRates[i];
  return sum;
}
Double_t AliNonseparationModelFit::GetChi2Moments() const {
  Double_t sum=0;
  for (Int_t i=0; i<6; ++i) {
    sum += fChi2Moments[i];
    AliInfoF("chi2Moments[%d]/ndf = %.1f (%f/%.0f)", i, fNDFMoments[i] ? fChi2Moments[i]/fNDFMoments[i] : 0.0, fChi2Moments[i], fNDFMoments[i]);
  }
  return sum;
}
Double_t AliNonseparationModelFit::GetChi2Rates() const {
  Double_t sum=0;
  for (Int_t i=0; i<6; ++i) {
    sum += fChi2Rates[i];
    AliInfoF("chi2Rate[%d]/ndf = %.1f", i, fNDFRates[i] ? fChi2Rates[i]/fNDFRates[i] : 0.0);
  }
  return sum;
}


const char* AliNonseparationModelFit::GetParName(Int_t idx) {
  if (idx >= 29)
    AliFatalClassF("%d >= 29", idx);

  if (idx < 20)
    return AliDoubleGaussianBeamProfile::GetParName(idx);

  static const char* parNames[] = {
    "#DeltaX",
    "#DeltaY",
    "#DeltaZ",
    "#DeltaX^{O}",
    "#DeltaY^{O}",
    "#DeltaZ^{O}",
    "sR",
    "LengthScaleX",
    "LengthScaleY"
  };
  return parNames[idx-20];
}

Int_t FindPoint(Double_t x, const TGraph *g, Double_t dx) {
  const Double_t *xx = g->GetX();
  for (Int_t i=0, n=g->GetN(); i<n; ++i) {
    if (TMath::Abs(xx[i]-x) < dx)
      return i;
  }
  return -1;
}

Double_t AliNonseparationModelFit::MinuitFunction(const Double_t *par) {
  Int_t       scanType(0);
  TVectorD    sep(2);
  TVectorD    mom(10);
  TMatrixDSym cov(10);

  TVectorD profile(8);
  TVectorD modelDGPar(20);

  for (Int_t i=0; i<20; ++i)
    modelDGPar[i] = par[i];

  Double_t chi2 = 0.0;

  TTree        *t = NULL;
  TGraphErrors *g = NULL;
  for (Int_t k=0; k<6; ++k) { // loop over scans
    fChi2Moments[k] = fChi2Rates[k] = fNDFMoments[k] = fNDFRates[k] = 0.0;

    t = (TTree*)fMoments.At(k);
    if (!t)
      continue;
    TEntryList *elist = t->GetEntryList();

    t->SetBranchAddress("scanType", &scanType);
    t->SetBranchAddress("beamSep",  sep.GetMatrixArray());
    t->SetBranchAddress("modelPar", mom.GetMatrixArray());
    t->SetBranchAddress("modelCov", cov.GetMatrixArray());

    for (Long64_t i=0, m=elist->GetN(); i<m; ++i) {
      t->GetEntry(elist->GetEntry(i));
      //Printf("%p %p %p", sep.GetMatrixArray(), mom.GetMatrixArray(), cov.GetMatrixArray());
#if 0
      AliDoubleGaussianBeamProfile::Eval(par[24]*sep(0) - fMuOffsetsX[k/2],
					 par[25]*sep(1) - fMuOffsetsY[k/2],
					 modelDGPar, profile);
      TVectorD profileP(8);
      AliDoubleGaussianBeamProfile::Eval(par[24]*sep(0) - fMuOffsetsX[k/2],
					 par[25]*sep(1) - fMuOffsetsY[k/2],
					 modelDGPar, profileP, 1e-3, kTRUE);
      Printf("TEST: ");
      profile.Print();
      profileP.Print();
#else
      AliDoubleGaussianBeamProfile::Eval(par[24+3]*sep(0) - fMuOffsetsX[k/2],
					 par[25+3]*sep(1) - fMuOffsetsY[k/2],
					 modelDGPar, profile, 1e-5, !kTRUE);
#endif
      // compute chi2 including correlations
      const TMatrixD covData_Inv(TMatrixD::kInverted, cov);
      Double_t xx = 0;
      const Double_t cx[7] = { 0.0005,  0.0005, 0.1, // 0 1 2
			       0.00001, 0.00001, 0.1, // 3 4 5
			       0.01};
      for (Int_t j1=0; j1<7; ++j1) {
	for (Int_t j2=0; j2<7; ++j2) {
	  if (j2 != j1) continue;
	  xx +=
	    (mom(j1) - (j1<3 ? par[20+j1+3*(k>=4)] : 0.0) - profile(1+j1)) *
	    // (mom(j2) - (j2<3 ? par[20+j2] : 0.0) - profile(1+j2)) / (j1>=0 && j1<=6 ? cx[j1]*cx[j2] : cov(j1,j2));
            (mom(j2) - (j2<3 ? par[20+j2+3*(k>=4)] : 0.0) - profile(1+j2)) * covData_Inv(j1,j2);
	}
      }
      fChi2Moments[k] += xx;
      //Printf("k=%d sep=%f,%f k=%f %f", k, sep(0), sep(1), mom[9], xx);
      fNDFMoments[k] += 7;

      // chi2 from rates
      g = (TGraphErrors*)fRates.At(k);
      if (FitToRates() && NULL != g) {
	const Int_t l = FindPoint(10*sep(scanType), g, 0.002);
	if (l < 0)
	  continue;
	const Double_t *g_x  = g->GetX();
	const Double_t *g_y  = g->GetY();
	const Double_t *g_ey = g->GetEY();
	if (g_ey[l] == 0)
	  continue;

	fChi2Rates[k] += TMath::Power((g_y[l] - par[23+3]*profile(0)) / g_ey[l], 2) * fScaleRateError;
	fNDFRates[k]  += 1;
      }
    } // next separation

    chi2 += fChi2Moments[k] + fChi2Rates[k];

    t->ResetBranchAddresses();
  } // next scan

  return chi2;
}
void AliNonseparationModelFit::SetVar(Int_t idx, Double_t val, Double_t step, Double_t min, Double_t max) {
  fMinimizer.SetLimitedVariable(idx, GetParName(idx), val, step, min, max);
}

void AliNonseparationModelFit::DoFit(const TString &saveFileName) {

  for (Int_t i=0; i<6; ++i) {
    TGraph *g = (TGraph*)fRates.At(i);
    if (NULL == g)
      continue;
    TF1 *fg = (TF1*)g->FindObject("gaus");
    if ((i%2)==0)
      fMuOffsetsX[i/2] = 0.1*fg->GetParameter(1); // conversion mm -> cm
    else
      fMuOffsetsY[i/2] = 0.1*fg->GetParameter(1); // conversion mm -> cm
  }

  fFitToOffsetScans = !kTRUE;

  // (1) fit without rates
  SetFitToRates(kFALSE);
  fMinimizer.FixVariable(3);
  fMinimizer.FixVariable(7);
  fMinimizer.FixVariable(12);
  fMinimizer.FixVariable(16);

  if (fMoments.At(4) == NULL) { // do not use DeltaXYZ^O without offset scans
    fMinimizer.FixVariable(23);
    fMinimizer.FixVariable(24);
    fMinimizer.FixVariable(25);
  }

  //X crossing angles
  fMinimizer.FixVariable(18);
  fMinimizer.FixVariable(19);

  fMinimizer.FixVariable(23+3);
  fMinimizer.FixVariable(24+3);
  fMinimizer.FixVariable(25+3);

  fMinimizer.Minimize();
  fMinimizer.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
   	 GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));

  fMinimizer.ReleaseVariable(3);
  fMinimizer.ReleaseVariable(7);
  fMinimizer.ReleaseVariable(12);
  fMinimizer.ReleaseVariable(16);

  fMinimizer.ReleaseVariable(18);
  fMinimizer.ReleaseVariable(19);

  fMinimizer.Clear();
  fMinimizer.Minimize();
  fMinimizer.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	 GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));

  fFitToOffsetScans = kTRUE;

  //X
  // fMinimizer.FixVariable(3);
  // fMinimizer.FixVariable(7);
  // fMinimizer.FixVariable(12);
  // fMinimizer.FixVariable(16);

  fMinimizer.Clear();
  fMinimizer.Minimize();
  fMinimizer.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	 GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));


  if (fNRates) {
    // (2a) fit with rates leaving only the scale free
    SetFitToRates(kTRUE);
    for (Int_t i=0; i<29; ++i)
      fMinimizer.FixVariable(i);

    fMinimizer.ReleaseVariable(23);

    fMinimizer.Minimize();
    fMinimizer.PrintResults();

    Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	   GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));

    Printf("chi2Rates = %.0f NDFRates=%.0f chi2/ndf = %.1f", GetChi2Rates(), GetNDFRates(), GetChi2Rates()/GetNDFRates());

    // (2b) scale the errors of the rates to have <chi2> = 6
    Double_t meanChi2Rates = GetChi2Rates()/GetNDFRates();
    SetScaleRateError(6.0/meanChi2Rates);
    if (fScaleRateError > 1.0)
      fScaleRateError = 1.0;

    // (2c) do a full fit
    for (Int_t i=0; i<23; ++i)
      if (!fFixCrossingAngles)
	fMinimizer.ReleaseVariable(i);
      else
	if (i!=19 && i!=18) fMinimizer.ReleaseVariable(i);

    fMinimizer.Clear();
    fMinimizer.Minimize();
    fMinimizer.PrintResults();

    Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	   GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));

    if (fScaleRateError != 1.0) {

      // (2d) adjust the rate errors a 2nd time
      meanChi2Rates = GetChi2Rates()/GetNDFRates()/fScaleRateError;
      SetScaleRateError(6.0/meanChi2Rates);

      if (fScaleRateError > 1.0) {
	fScaleRateError = 1.0;
      } else {
	fMinimizer.Clear();
	fMinimizer.Minimize();
	fMinimizer.PrintResults();

	Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	       GetChi2(), GetNDF()-fMinimizer.NFree(), GetChi2()/(GetNDF()-fMinimizer.NFree()));
      }
    }
  }
  printf("   ");
  for (Int_t j=0; j<24; ++j)
    printf("%5d ", 1+j);
  printf("\n");

  for (Int_t j=0; j<24; ++j) {
    printf("%2d ", 1+j);
    for (Int_t k=0; k<=j; ++k) {
      printf("%5.2f ", fMinimizer.Correlation(j,k));
    }
    printf("\n");
  }
  Printf("chi2Rates = %.0f NDFRates=%.0f chi2/ndf = %.1f", GetChi2Rates(), GetNDFRates(), GetChi2Rates()/GetNDFRates());

  TParameter<Double_t>            pNDF("ndf",            GetNDF()-fMinimizer.NFree());
  TParameter<Double_t>           pChi2("chi2",           GetChi2());
  TParameter<Double_t> pScaleRateError("scaleRateError", fScaleRateError);

  Printf("saveFileName.Data=%s", saveFileName.Data());
  TFile *fSave = TFile::Open(saveFileName, "RECREATE");
  gMinuit->Write("m");
  pNDF.Write();
  pChi2.Write();
  pScaleRateError.Write();
  fMuOffsetsX.Write("muOffsetsX");
  fMuOffsetsY.Write("muOffsetsY");
  for (Int_t i=0; i<6; ++i) {
    if (fRates.At(i))
      fRates.At(i)->Write(Form("gRateScan%c_%d", (i%2) ? 'Y' : 'X', i/2));
    if (fMoments.At(i)) {
      TTree *tt = (TTree*)fMoments.At(i);
      tt->GetEntryList()->Print();
      tt->GetEntryList()->SetDirectory(fSave);
      fMoments.At(i)->Write(Form("fMoments%c_%d", (i%2) ? 'Y' : 'X', i/2));
    }
  }
  fSave->Write();
  gDirectory->ls();
  fSave->Close();
}
