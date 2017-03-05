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

#include <Math/Functor.h>
#include <TMinuit.h>
#include <TFitter.h>
#include <TMinuitMinimizer.h>

#include "AliLog.h"
#include "AliDoubleGaussianBeamProfile.h"
#include "AliNonseparationModelFit.h"

ClassImp(AliNonseparationModelFit);

AliNonseparationModelFit::AliNonseparationModelFit()
  : TObject()
  , fMoments(new TObjArray(6))
  , fRates(new TObjArray(6))
  , fPar(26)
  , fNRates(0)
  , fFitToRates(kFALSE)
  , fScaleRateError(1)
  , fFitToOffsetScans(kTRUE)
  , fMuOffsetsX(3)
  , fMuOffsetsY(3)
{
  fMoments->SetOwner(kTRUE);
  fRates->SetOwner(kTRUE);

  for (Int_t i=0; i<6; ++i)
    fChi2Moments[i] = fChi2Rates[i] = fNDFMoments[i] = fNDFRates[i] = 0.0;
}

AliNonseparationModelFit::~AliNonseparationModelFit() {
  if (fMoments)
    delete fMoments;
   fMoments = NULL;

  if (fRates)
    delete fRates;
   fRates = NULL;
}

void AliNonseparationModelFit::Add(Int_t idx, TTree *tMoments, TGraph *gRate) {
  fMoments->AddAt(tMoments, idx);
  fRates->AddAt(gRate, idx);
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
  if (idx >= 26)
    AliFatalClassF("%d >= 26", idx);

  if (idx < 20)
    return AliDoubleGaussianBeamProfile::GetParName(idx);

  static const char* parNames[] = {
    "#DeltaX",
    "#DeltaY",
    "#DeltaZ",
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

Bool_t AliNonseparationModelFit::IsGoodLumiRegionFit(Int_t k, Int_t scanType, const TVectorD& sep,
						     const TVectorD& mom, const TMatrixDSym& cov) {
  const Double_t maxSep = (k>=4 ? 0.04 : 0.06);
  if (k<4 && (TMath::Abs(sep(0)) > maxSep || TMath::Abs(sep(1)) > maxSep))
    return kFALSE;
  if (k==4 && TMath::Abs(sep(0)) > maxSep)
    return kFALSE;
  if (k==5 && TMath::Abs(sep(1)) > maxSep)
    return kFALSE;
  if (mom(9) < 1 || mom(9) > 1.6)
    return kFALSE;
  return kTRUE;
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

    t = (TTree*)fMoments->At(k);
    if (NULL == t)
      continue;
    t->SetBranchAddress("scanType", &scanType);
    t->SetBranchAddress("beamSep",  sep.GetMatrixArray());
    t->SetBranchAddress("modelPar", mom.GetMatrixArray());
    t->SetBranchAddress("modelCov", cov.GetMatrixArray());

    for (Int_t i=0, n=t->GetEntries(); i<n; ++i) { // loop over separations
      t->GetEntry(i);
      if (!IsGoodLumiRegionFit(k, scanType, sep, mom, cov))
	continue;
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
      AliDoubleGaussianBeamProfile::Eval(par[24]*sep(0) - fMuOffsetsX[k/2],
					 par[25]*sep(1) - fMuOffsetsY[k/2],
					 modelDGPar, profile, 5e-3, !kTRUE);
#endif
      // compute chi2 including correlations
      const TMatrixD covData_Inv(TMatrixD::kInverted, cov);
      for (Int_t j1=0; j1<7; ++j1) {
	for (Int_t j2=0; j2<7; ++j2) {
	  fChi2Moments[k] += 
	    (mom(j1) - (j1<3 ? par[20+j1] : 0.0) - profile(1+j1)) *
	    (mom(j2) - (j2<3 ? par[20+j2] : 0.0) - profile(1+j2)) * covData_Inv(j1,j2);
	}
      }
      fNDFMoments[k] += 7;

      // chi2 from rates
      g = (TGraphErrors*)fRates->At(k);
      if (FitToRates() && NULL != g) {
	const Int_t l = FindPoint(10*sep(scanType), g, 0.002);
	if (l < 0)
	  continue;
	const Double_t *g_x  = g->GetX();
	const Double_t *g_y  = g->GetY();
	const Double_t *g_ey = g->GetEY();
	if (g_ey[l] == 0)
	  continue;
	
	fChi2Rates[k] += TMath::Power((g_y[l] - par[23]*profile(0)) / g_ey[l], 2) * fScaleRateError;
	fNDFRates[k]  += 1;
      }
    } // next separation

    chi2 += fChi2Moments[k] + fChi2Rates[k];

    t->ResetBranchAddresses();
  } // next scan

  return chi2;
}

void AliNonseparationModelFit::DoFit(TVectorD& par,
				     const TString &saveFileName) {

  ROOT::Math::Functor fcn(this, &AliNonseparationModelFit::MinuitFunction, 26);
  TMinuitMinimizer m("minimize", 26); 
  m.UseStaticMinuit(kTRUE);
  m.SetFunction(fcn);

  m.SetLimitedVariable( 0, GetParName( 0),  par[ 0], 0.001, 0.008,  0.015);
  m.SetLimitedVariable( 1, GetParName( 1),  par[ 1], 0.001, 0.008,  0.015);
  m.SetLimitedVariable( 2, GetParName( 2),  par[ 2], 0.10, 4.0, 10.0);
  m.SetLimitedVariable( 3, GetParName( 3),  par[ 3], 0.01, -1, 1);

  m.SetLimitedVariable( 4, GetParName( 4),  par[ 4], 0.01, 1.0, 2.0);
  m.SetLimitedVariable( 5, GetParName( 5),  par[ 5], 0.01, 1.0, 2.0);
  m.SetLimitedVariable( 6, GetParName( 6),  par[ 6], 0.01, 0.5, 2.0);
  m.SetLimitedVariable( 7, GetParName( 7),  par[ 7], 0.01, -1, 1);

  m.SetLimitedVariable( 8, GetParName( 8),  par[ 8], 0.01, 0.5, 0.95);

  m.SetLimitedVariable( 9, GetParName( 9),  par[ 9], 0.001, 0.008,  0.015);
  m.SetLimitedVariable(10, GetParName(10),  par[10], 0.001, 0.008,  0.015);
  m.SetLimitedVariable(11, GetParName(11),  par[11], 0.10, 4.0, 10.0);
  m.SetLimitedVariable(12, GetParName(12),  par[12], 0.01, -1, 1);

  m.SetLimitedVariable(13, GetParName(13),  par[13], 0.01, 1.0, 2.0);
  m.SetLimitedVariable(14, GetParName(14),  par[14], 0.01, 1.0, 2.0);
  m.SetLimitedVariable(15, GetParName(15),  par[15], 0.01, 0.5, 2.0);
  m.SetLimitedVariable(16, GetParName(16),  par[16], 0.01, -1, 1);

  m.SetLimitedVariable(17, GetParName(17),  par[17], 0.01, 0.5, 0.95);

  m.SetLimitedVariable(18, GetParName(18),  par[18], 1e-7, -0.1, +0.1);
  m.SetLimitedVariable(19, GetParName(19),  par[19], 1e-6, -0.1, +0.1);

  m.SetLimitedVariable(20, GetParName(20),  par[20], 1e-3, par[20]-0.09, par[20]+0.09);
  m.SetLimitedVariable(21, GetParName(21),  par[21], 1e-3, par[21]-0.09, par[21]+0.09);
  m.SetLimitedVariable(22, GetParName(22),  par[22], 0.1, -2, 2);

  m.SetLimitedVariable(23, GetParName(23),  par[23], 1e-4, 0.0, 0.001);

  m.SetLimitedVariable(24, GetParName(24),  par[24], 0.001, 0.95, 1.05);
  m.SetLimitedVariable(25, GetParName(25),  par[25], 0.001, 0.95, 1.05);

  for (Int_t i=0; i<6; ++i) {
    TGraph *g = (TGraph*)fRates->At(i);
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
  m.FixVariable(3);
  m.FixVariable(7);
  m.FixVariable(12);
  m.FixVariable(16);

  m.FixVariable(18);
  m.FixVariable(19);

  m.FixVariable(23);
  m.FixVariable(24);
  m.FixVariable(25);

  m.Minimize();
  m.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	 GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));

  m.ReleaseVariable(3);
  m.ReleaseVariable(7);
  m.ReleaseVariable(12);
  m.ReleaseVariable(16);
  m.Clear();
  m.Minimize();
  m.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	 GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));


  fFitToOffsetScans = kTRUE;

  m.FixVariable(3);
  m.FixVariable(7);
  m.FixVariable(12);
  m.FixVariable(16);

  m.Clear();
  m.Minimize();
  m.PrintResults();

  Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	 GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));


  if (fNRates) {
    // (2a) fit with rates leaving only the scale free
    SetFitToRates(kTRUE);
    for (Int_t i=0; i<26; ++i)
      m.FixVariable(i);

    m.ReleaseVariable(23);

    m.Minimize();
    m.PrintResults();

    Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	   GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));

    Printf("chi2Rates = %.0f NDFRates=%.0f chi2/ndf = %.1f", GetChi2Rates(), GetNDFRates(), GetChi2Rates()/GetNDFRates());

    // (2b) scale the errors of the rates to have <chi2> = 6 
    Double_t meanChi2Rates = GetChi2Rates()/GetNDFRates();
    SetScaleRateError(6.0/meanChi2Rates);
    if (fScaleRateError > 1.0)
      fScaleRateError = 1.0;  

    // (2c) do a full fit
    for (Int_t i=0; i<23; ++i)
      if (i!=19 && i!=18) m.ReleaseVariable(i);
    
    m.Clear();
    m.Minimize();
    m.PrintResults();

    Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	   GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));

    if (fScaleRateError != 1.0) {

      // (2d) adjust the rate errors a 2nd time
      meanChi2Rates = GetChi2Rates()/GetNDFRates()/fScaleRateError;
      SetScaleRateError(6.0/meanChi2Rates );

      if (fScaleRateError > 1.0) {
	fScaleRateError = 1.0;
      } else {    
	m.Clear();
	m.Minimize();
	m.PrintResults();
	
	Printf("chi2=%.0f NDF=%.0f chi2/NDF=%.1f",
	       GetChi2(), GetNDF()-m.NFree(), GetChi2()/(GetNDF()-m.NFree()));
      }
    }
  }
  Printf("chi2Rates = %.0f NDFRates=%.0f chi2/ndf = %.1f", GetChi2Rates(), GetNDFRates(), GetChi2Rates()/GetNDFRates());

  TParameter<Double_t>            pNDF("ndf",            GetNDF()-m.NFree());
  TParameter<Double_t>           pChi2("chi2",           GetChi2());
  TParameter<Double_t> pScaleRateError("scaleRateError", fScaleRateError);

  Printf("%s", saveFileName.Data());
  TFile *fSave = TFile::Open(saveFileName, "RECREATE");
  gMinuit->Write("m");
  pNDF.Write();
  pChi2.Write();
  pScaleRateError.Write();
  fMuOffsetsX.Write("muOffsetsX");
  fMuOffsetsY.Write("muOffsetsY");
  for (Int_t i=0; i<6; ++i) {
    if (fRates->At(i))
      fRates->At(i)->Write(Form("gRateScan%c_%d", (i%2) ? 'Y' : 'X', i/2));
    if (fMoments->At(i))
      fMoments->At(i)->Write(Form("fMoments%c_%d", (i%2) ? 'Y' : 'X', i/2));
  }

  fSave->Write();
  fSave->Close();
}
// this comment serves no purpose
