/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                Dielectron SignalFunc                                  //
//                                                                       //
//                                                                       //
/*
Dielectron signal extraction class using functions as input.

A function to describe the signal as well as one to describe the background
has to be deployed by the user. Alternatively on of the default implementaions
can be used.

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>
#include <TList.h>
#include <TFitResult.h>
//#include <../hist/hist/src/TF1Helper.h>

#include <AliLog.h>

#include "AliDielectronSignalFunc.h"

ClassImp(AliDielectronSignalFunc)

TH1F* AliDielectronSignalFunc::fgHistSimPM=0x0;

AliDielectronSignalFunc::AliDielectronSignalFunc() :
AliDielectronSignalBase(),
fFuncSignal(0x0),
fFuncBackground(0x0),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE),
fPolDeg(0),
fChi2Dof(0.0)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalFunc::AliDielectronSignalFunc(const char* name, const char* title) :
AliDielectronSignalBase(name, title),
fFuncSignal(0x0),
fFuncBackground(0x0),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE),
fPolDeg(0),
fChi2Dof(0.0)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronSignalFunc::~AliDielectronSignalFunc()
{
  //
  // Default Destructor
  //
  if(fFuncSignal) delete fFuncSignal;
  if(fFuncBackground) delete fFuncBackground;
  if(fFuncSigBack) delete fFuncSigBack;
}


//______________________________________________
void AliDielectronSignalFunc::Process(TObjArray * const arrhist)
{
  //
  // Fit the invariant mass histograms and retrieve the signal and background
  //
  switch(fMethod) {
  case kFittedMC :
	ProcessFitIKF(arrhist);
	break;
		  
  case kFitted :
    ProcessFit(arrhist);
    break;
    
  case kLikeSign :
    ProcessLS(arrhist);
    break;
    
  case kEventMixing :
    ProcessEM(arrhist);
    break;
    
  default :
    AliError("Background substraction method not known!");
  }
}

//______________________________________________
void AliDielectronSignalFunc::ProcessFitIKF(TObjArray * const arrhist) {
  //
  // Fit the +- invariant mass distribution only
  // 
  //
  
  const Double_t bigNumber = 100000.;
  Double_t chi2ndfm1 = bigNumber;
  Double_t ratiom1   = bigNumber;
  Double_t chi2ndf   = bigNumber;
  Int_t nDP =0;
  
  Int_t maxPolDeg = 8;
      
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
    if(fRebin>1) fHistDataPM->Rebin(fRebin);
  
  fgHistSimPM = (TH1F*)(arrhist->At(3))->Clone("histPMsim");  // +- mc shape
  if (!fgHistSimPM) {
	AliFatal("No mc peak shape found at idx 3.");
	return;
  }
  if(fRebin>1) fgHistSimPM->Rebin(fRebin);
  
  // try out the polynomial degrees
  for (Int_t iPD=0; iPD<=maxPolDeg; iPD++) {
    TH1F *hf1 = (TH1F *) fHistDataPM->Clone(Form("hf1_PD%d",iPD));
	
	FitOneMinv(hf1, fgHistSimPM, iPD);
	if (fChi2Dof > 0) chi2ndf = fChi2Dof;
    AliInfo(Form("nDP: %d, iPD: %d, chi2ndf: %f", nDP, iPD, chi2ndf));
	
    ratiom1 = TMath::Abs(fChi2Dof - 1);
    if (chi2ndfm1 > ratiom1) { // search for the closest to 1.
      chi2ndfm1 = ratiom1;
      nDP       = iPD;
    }
  }
    
  
  // fit again with the best polynomial degree
  TH1F *h2 = (TH1F *) fHistDataPM->Clone(Form("h2_PD%d",nDP));
  
  FitOneMinv(h2, fgHistSimPM, nDP);
  AliInfo(Form("Best Fit: PD %d, chi^2/ndf %.3f, S/B %.3f",nDP,fChi2Dof,fValues(3)));
    
}
//______________________________________________
void AliDielectronSignalFunc::FitOneMinv(TH1F *hMinv, TH1F *hSim, Int_t pod) {
  //
  // main function to fit an inv mass spectrum
  //
  
  TObjArray *arrResults = new TObjArray;
  arrResults->SetOwner();
  arrResults->AddAt(hMinv,0);
  
  // Degree of polynomial
  fPolDeg = pod;
    
  // inclusion and exclusion areas (values)
  const Double_t kJPMass    = 3.096916;
  // inclusion and exclusion areas (bin numbers)
  const Int_t binIntLo = hMinv->FindBin(fIntMin);
  const Int_t binIntHi = hMinv->FindBin(fIntMax);
  // for error calculation
  Double_t intAreaEdgeLo = hMinv->GetBinLowEdge(binIntLo);
  Double_t intAreaEdgeHi = hMinv->GetBinLowEdge(binIntHi)+hMinv->GetBinWidth(binIntHi);
  Double_t norm = (binIntHi-binIntLo)/(fIntMax - fIntMin); 
  
  TH1F  *hBfit = (TH1F *) hMinv->Clone(); // for bg only fit (excluding peak region)
  TH1F  *hSigF = (TH1F *) hMinv->Clone(); // signal with subtracted bg
  TH1F  *hBgrd = (TH1F *) hMinv->Clone(); // bg histogram
  
  hBfit->Reset();
  hSigF->Reset();
  hBgrd->Reset();
  
  // extract start parameter for the MC signal fit
  Double_t bgvalAv = (hMinv->Integral(1,hMinv->GetNbinsX()+1) - hMinv->Integral(binIntLo,binIntHi)) / (hMinv->GetNbinsX()+1 - (binIntHi-binIntLo));
  Double_t pkval     = hMinv->GetBinContent(hMinv->FindBin(kJPMass)) - bgvalAv;
  Double_t heightMC  = hSim->GetBinContent(hSim->FindBin(kJPMass));
  Double_t peakScale = (heightMC > 0. ? pkval/heightMC : 0.0);
  
  Int_t nBgnd = 2 + fPolDeg; // degree + c1st oefficient + higher coefficients
  Int_t nMinv = nBgnd + 1;  // bgrd + peakscale
  
  // Create the spectra without peak region
  for (Int_t iBin = 0; iBin <= hMinv->GetNbinsX(); iBin++) {
    if ((iBin < binIntLo) || (iBin > binIntHi)) {
      hBfit->SetBinContent(iBin,hMinv->GetBinContent(iBin));
      hBfit->SetBinError(iBin,hMinv->GetBinError(iBin));
	}
  }
    
  
  // =======
  //   1.
  // =======
  // Do the fit to the background spectrum
  TF1 *fBo = new TF1("bgrd_fit",BgndFun,fFitMin,fFitMax,nBgnd);
  for (Int_t iPar=0; iPar<nBgnd; iPar++) {
    if (iPar == 0) fBo->FixParameter(0, fPolDeg);
	if (iPar == 1) fBo->SetParameter(iPar, bgvalAv);
    if (iPar >= 2) fBo->SetParameter(iPar, 0.);
  }
  hBfit->Fit(fBo,"0qR");
  //hBfit->SetNameTitle("bgrd_fit");
  arrResults->AddAt(fBo,1);
  
  
  // =======
  //   2.
  // =======
  // Fit the whole spectrum with peak and background
  TF1 *fSB = new TF1("bgrd_peak_fit",MinvFun,fFitMin,fFitMax,nMinv);
  fSB->FixParameter(0, fPolDeg);
  fSB->SetParameter(1, peakScale);
  // copy the polynomial parameters
  for (Int_t iPar=0; iPar<=fPolDeg; iPar++) 
	fSB->SetParameter(2+iPar, fBo->GetParameter(iPar+1));
  hMinv->Fit(fSB,"0qR");
  arrResults->AddAt(fSB,2);
  
  
  // =======
  //   3.
  // =======
  // Create the background function
  TF1 *fB = new TF1("bgrdOnly_fkt",BgndFun,fFitMin,fFitMax,nBgnd);
  fB->FixParameter(0,fPolDeg);
  for (Int_t iDeg=0; iDeg<=fPolDeg; iDeg++) {
	fB->SetParameter(1+iDeg,fSB->GetParameter(2+iDeg));
	fB->SetParError(1+iDeg,fSB->GetParError(2+iDeg));
  }
  // create background histogram from background function
  hBgrd->Eval(fB);
  hBgrd->Fit(fB,"0qR");
  // calculate the integral and integral error from fit function
  Double_t intc = fB->Integral(intAreaEdgeLo, intAreaEdgeHi) * norm;
  Double_t inte = fB->IntegralError(intAreaEdgeLo, intAreaEdgeHi) * norm;
  arrResults->AddAt(fB,3);

  // Fill the background spectrum erros. Use the error from the fit function for the background fB
  for (Int_t iBin = 0; iBin <= hBgrd->GetNbinsX(); iBin++) {
    Double_t x = hBgrd->GetBinCenter(iBin);
    if ((x >= fFitMin) && (x <= fFitMax)) {
	  Double_t binte = inte / TMath::Sqrt((binIntHi-binIntLo)+1);
	  hBgrd->SetBinError(iBin,binte); 
    }
  }
  arrResults->AddAt(hBgrd,4);
  
  // =======
  //   4.
  // =======
  // Subtract the background  
  hSigF->Add(hMinv,hBgrd,1.0,-1.0);
  for (Int_t iBin = 0; iBin <= hSigF->GetNbinsX(); iBin++) {
    Double_t x = hSigF->GetBinCenter(iBin);
    if ((x < fFitMin) || (x > fFitMax)) {
      hSigF->SetBinContent(iBin,0.0);
      hSigF->SetBinError(iBin,0.0);
    }
  }
  hSigF->SetNameTitle("peak_only","");
  arrResults->AddAt(hSigF,5);
  
  // =======
  //   5.
  // =======
  // Fit the background-subtracted spectrum
  TF1 *fS  = new TF1("peakOnly_fit",PeakFunCB,fFitMin,fFitMax,5);
  fS->SetParameters(-.05,1,kJPMass,.003,700);
  fS->SetParNames("alpha","n","meanx","sigma","N");
  hSigF->Fit(fS,"0qR");
  arrResults->AddAt(fS,6);
  
  
  // connect data members 
  fFuncSignal     = (TF1*) arrResults->At(6)->Clone();
  fFuncBackground = (TF1*) arrResults->At(3)->Clone();
  fFuncSigBack    = (TF1*) arrResults->At(2)->Clone();
  fHistSignal     = (TH1F*)arrResults->At(5)->Clone();
  fHistBackground = (TH1F*)arrResults->At(4)->Clone();
  
  // fit results
  Double_t chi2 = fSB->GetChisquare();
  Int_t    ndf  = fSB->GetNDF();
  fChi2Dof = chi2/ndf;
  
  // signal + signal error
  fValues(0) = hSigF->IntegralAndError(binIntLo, binIntHi, fErrors(0));
  fValues(1) = intc; // background
  fErrors(1) = inte; // background error
  // S/B (2) and significance (3)
  SetSignificanceAndSOB();
  fValues(4) = fS->GetParameter(2); // mass
  fErrors(4) = fS->GetParError(2);  // mass error
  fValues(5) = fS->GetParameter(3); // mass wdth
  fErrors(5) = fS->GetParError(3);  // mass wdth error
  
    
  delete arrResults;
  
}
//______________________________________________________________________________
Double_t AliDielectronSignalFunc::BgndFun(const Double_t *x, const Double_t *par) {
  // parameters
  // [0]:   degree of polynomial
  // [1]:   constant polynomial coefficient
  // [2]..: higher polynomial coefficients
  
  Int_t    deg   = ((Int_t) par[0]);
  
  Double_t f     = 0.0;
  Double_t yy    = 1.0;
  Double_t xx    = x[0];

  for (Int_t i = 0; i <= deg; i++) {
    f  += par[i+1] * yy;
    yy *= xx;
  } 
  
    
  return f;
}
//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakFun(const Double_t *x, const Double_t *par) {
  // Fit MC signal shape
  // parameters
  // [0]:   scale for simpeak
  
  Double_t xx  = x[0];
  
  TH1F *hPeak = fgHistSimPM;
  if (!hPeak) {
    printf("F-AliDielectronSignalFunc::PeakFun: No histogram for peak fit defined!\n");
  }
  
  Int_t idx = hPeak->FindBin(xx);
  if ((idx <= 1) ||
      (idx >= hPeak->GetNbinsX())) {
    return 0.0;
  }
  
  return (par[0] * hPeak->GetBinContent(idx));
  
}

//______________________________________________________________________________
Double_t AliDielectronSignalFunc::MinvFun(const Double_t *x, const Double_t *par) {
  // parameters
  // [0]:   degree of polynomial             -> [0]   for BgndFun
  // [1]:   scale for simpeak                -> [0]   for PeakFun
  // [2]:   constant polynomial coefficient  -> [1]   for BgndFun
  // [3]..: higher polynomial coefficients   -> [2].. for BgndFun
  
  Int_t    deg = ((Int_t) par[0]);
  Double_t parPK[25], parBG[25];
  
  parBG[0] = par[0]; // degree of polynomial
  
  parPK[0] = par[1]; // MC minv scale
  for (Int_t i = 0; i <= deg; i++) parBG[i+1] = par[i+2]; // polynomial coefficients
  
  Double_t peak = PeakFun(x,parPK);
  Double_t bgnd = BgndFun(x,parBG);
  Double_t f    = peak + bgnd;
  
  return f;
}

//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakFunCB(const Double_t *x, const Double_t *par) {
  // Crystal Ball function fit
  
  Double_t alpha = par[0];
  Double_t     n = par[1];
  Double_t meanx = par[2];
  Double_t sigma = par[3];
  Double_t    nn = par[4];
   
  Double_t a = TMath::Power((n/TMath::Abs(alpha)), n) * TMath::Exp(-.5*alpha*alpha);
  Double_t b = n/TMath::Abs(alpha) - TMath::Abs(alpha);
  
  Double_t arg = (x[0] - meanx)/sigma;
  Double_t fitval = 0;
  
  if (arg > -1.*alpha) {
    fitval = nn * TMath::Exp(-.5*arg*arg);
  } else {
    fitval = nn * a * TMath::Power((b-arg), (-1*n));
  }
  
  return fitval;
}


//______________________________________________
void AliDielectronSignalFunc::ProcessFit(TObjArray * const arrhist) {
  //
  // Fit the +- invariant mass distribution only
  // Here we assume that the combined fit function is a sum of the signal and background functions
  //    and that the signal function is always the first term of this sum
  //
  
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataPM->Sumw2();
  if(fRebin>1)
    fHistDataPM->Rebin(fRebin);
  
  fHistSignal = new TH1F("HistSignal", "Fit substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground = new TH1F("HistBackground", "Fit contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  // the starting parameters of the fit function and their limits can be tuned
  // by the user in its macro
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  TFitResultPtr pmFitPtr = fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  //TFitResult *pmFitResult = pmFitPtr.Get(); // used only with TF1Helper
  fFuncSignal->SetParameters(fFuncSigBack->GetParameters());
  fFuncBackground->SetParameters(fFuncSigBack->GetParameters()+fFuncSignal->GetNpar());
  
  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    Double_t m = fHistDataPM->GetBinCenter(iBin);
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t epm = fHistDataPM->GetBinError(iBin);
    Double_t bknd = fFuncBackground->Eval(m);
    Double_t ebknd = 0;
    for(Int_t iPar=fFuncSignal->GetNpar(); iPar<fFuncSigBack->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSigBack->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        ebknd += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    Double_t signal = pm-bknd;
    Double_t error = TMath::Sqrt(epm*epm+ebknd);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, error);
    fHistBackground->SetBinContent(iBin, bknd);
    fHistBackground->SetBinError(iBin, TMath::Sqrt(ebknd));
  }
  
  if(fUseIntegral) {
    // signal
    fValues(0) = fFuncSignal->Integral(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
    fErrors(0) = 0;
    for(Int_t iPar=0; iPar<fFuncSignal->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSignal->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,fFuncSignal),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,fFuncSignal),0,0,0);
        fErrors(0) += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Integral(fIntMin,fIntMax)*gradientJpar.Integral(fIntMin,fIntMax)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    // background
    fValues(1) = fFuncBackground->Integral(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
    fErrors(1) = 0;
    for(Int_t iPar=fFuncSignal->GetNpar(); iPar<fFuncSigBack->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSigBack->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        fErrors(1) += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Integral(fIntMin, fIntMax)*gradientJpar.Integral(fIntMin, fIntMax)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
  }
  else {
    // signal
    fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
                                               fHistSignal->FindBin(fIntMax), fErrors(0));
    // background
    fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
      fHistBackground->FindBin(fIntMax),
      fErrors(1));
  }
  // S/B and significance
  SetSignificanceAndSOB();
  fValues(4) = fFuncSigBack->GetParameter(fParMass);
  fErrors(4) = fFuncSigBack->GetParError(fParMass);
  fValues(5) = fFuncSigBack->GetParameter(fParMassWidth);
  fErrors(5) = fFuncSigBack->GetParError(fParMassWidth);
  
  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalFunc::ProcessLS(TObjArray * const arrhist) {
  //
  // Substract background using the like-sign spectrum
  //
  fHistDataPP = (TH1F*)(arrhist->At(0))->Clone("histPP");  // ++    SE
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataMM = (TH1F*)(arrhist->At(2))->Clone("histMM");  // --    SE
  if (fRebin>1) {
    fHistDataPP->Rebin(fRebin);
    fHistDataPM->Rebin(fRebin);
    fHistDataMM->Rebin(fRebin);
  }
  fHistDataPP->Sumw2();
  fHistDataPM->Sumw2();
  fHistDataMM->Sumw2();
  
  fHistSignal = new TH1F("HistSignal", "Like-Sign substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground = new TH1F("HistBackground", "Like-sign contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  // fit the +- mass distribution
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  // declare the variables where the like-sign fit results will be stored
//   TFitResult *ppFitResult = 0x0;
//   TFitResult *mmFitResult = 0x0;
  // fit the like sign background
  TF1 *funcClonePP = (TF1*)fFuncBackground->Clone("funcClonePP");
  TF1 *funcCloneMM = (TF1*)fFuncBackground->Clone("funcCloneMM");
  fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
//   TFitResultPtr ppFitPtr = fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
//   ppFitResult = ppFitPtr.Get();
  fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
//   TFitResultPtr mmFitPtr = fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
//   mmFitResult = mmFitPtr.Get();
  
  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    Double_t m = fHistDataPM->GetBinCenter(iBin);
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t pp = funcClonePP->Eval(m);
    Double_t mm = funcCloneMM->Eval(m);
    Double_t epm = fHistDataPM->GetBinError(iBin);
    Double_t epp = 0;
    for(Int_t iPar=0; iPar<funcClonePP->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<funcClonePP->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,funcClonePP),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,funcClonePP),0,0,0);
        epp += ppFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    Double_t emm = 0;
    for(Int_t iPar=0; iPar<funcCloneMM->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<funcCloneMM->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,funcCloneMM),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,funcCloneMM),0,0,0);
        emm += mmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    
    Double_t signal = pm-2.0*TMath::Sqrt(pp*mm);
    Double_t background = 2.0*TMath::Sqrt(pp*mm);
    // error propagation on the signal calculation above
    Double_t esignal = TMath::Sqrt(epm*epm+(mm/pp)*epp+(pp/mm)*emm);
    Double_t ebackground = TMath::Sqrt((mm/pp)*epp+(pp/mm)*emm);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, esignal);
    fHistBackground->SetBinContent(iBin, background);
    fHistBackground->SetBinError(iBin, ebackground);
  }
  
  // signal
  fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
                                             fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
                                                 fHistBackground->FindBin(fIntMax),
                                                 fErrors(1));
  // S/B and significance
  SetSignificanceAndSOB();
  fValues(4) = fFuncSigBack->GetParameter(fParMass);
  fErrors(4) = fFuncSigBack->GetParError(fParMass);
  fValues(5) = fFuncSigBack->GetParameter(fParMassWidth);
  fErrors(5) = fFuncSigBack->GetParError(fParMassWidth);
  
  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalFunc::ProcessEM(TObjArray * const arrhist) {
  //
  // Substract background with the event mixing technique
  //
  arrhist->GetEntries();   // just to avoid the unused parameter warning
  AliError("Event mixing for background substraction method not implemented!");
}

//______________________________________________
void AliDielectronSignalFunc::SetFunctions(TF1 * const combined, TF1 * const sig, TF1 * const back,
                                           Int_t parM, Int_t parMres)
{
  //
  // Set the signal, background functions and combined fit function
  // Note: The process method assumes that the first n parameters in the
  //       combined fit function correspond to the n parameters of the signal function
  //       and the n+1 to n+m parameters to the m parameters of the background function!!!
  
  if (!sig||!back||!combined) {
    AliError("Both, signal and background function need to be set!");
    return;
  }
  fFuncSignal=sig;
  fFuncBackground=back;
  fFuncSigBack=combined;
  fParMass=parM;
  fParMassWidth=parMres;
}

//______________________________________________
void AliDielectronSignalFunc::SetDefaults(Int_t type)
{
  //
  // Setup some default functions:
  // type = 0: gaus signal + linear background in 2.5 - 4 GeV inv. mass
  // type = 1: gaus signal + exponential background in 2.5 - 4 GeV inv. mass
  // type = 2: half gaussian, half exponential signal function
  // type = 3: Crystal-Ball function
  // type = 4: Crystal-Ball signal + exponential background
  //
  
  if (type==0){
    fFuncSignal=new TF1("DieleSignal","gaus",2.5,4);
    fFuncBackground=new TF1("DieleBackground","pol1",2.5,4);
    fFuncSigBack=new TF1("DieleCombined","gaus+pol1(3)",2.5,4);
    
    fFuncSigBack->SetParameters(1,3.1,.05,2.5,1);
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
  }
  else if (type==1){
    fFuncSignal=new TF1("DieleSignal","gaus",2.5,4);
    fFuncBackground=new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])",2.5,4);
    fFuncSigBack=new TF1("DieleCombined","gaus+[3]*exp(-(x-[4])/[5])",2.5,4);
    
    fFuncSigBack->SetParameters(1,3.1,.05,1,2.5,1);
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
  }
  else if (type==2){
    // half gaussian, half exponential signal function
    // exponential background
    fFuncSignal = new TF1("DieleSignal","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",2.5,4);
    fFuncBackground = new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])+[3]",2.5,4);
    fFuncSigBack = new TF1("DieleCombined","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))+[4]*exp(-(x-[5])/[6])+[7]",2.5,4);
    fFuncSigBack->SetParameters(1.,3.1,.05,.1,1,2.5,1,0);
    
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
    fFuncSigBack->FixParameter(6,2.5);
    fFuncSigBack->FixParameter(7,0);
  }
}


//______________________________________________
void AliDielectronSignalFunc::Draw(const Option_t* option)
{
  //
  // Draw the fitted function
  //
  
  TString drawOpt(option);
  drawOpt.ToLower();
  
  Bool_t optStat=drawOpt.Contains("stat");
  
  fFuncSigBack->SetNpx(200);
  fFuncSigBack->SetRange(fIntMin,fIntMax);
  fFuncBackground->SetNpx(200);
  fFuncBackground->SetRange(fIntMin,fIntMax);
  
  TGraph *grSig=new TGraph(fFuncSigBack);
  grSig->SetFillColor(kGreen);
  grSig->SetFillStyle(3001);
  
  TGraph *grBack=new TGraph(fFuncBackground);
  grBack->SetFillColor(kRed);
  grBack->SetFillStyle(3001);
  
  grSig->SetPoint(0,grBack->GetX()[0],grBack->GetY()[0]);
  grSig->SetPoint(grSig->GetN()-1,grBack->GetX()[grBack->GetN()-1],grBack->GetY()[grBack->GetN()-1]);
  
  grBack->SetPoint(0,grBack->GetX()[0],0.);
  grBack->SetPoint(grBack->GetN()-1,grBack->GetX()[grBack->GetN()-1],0.);
  
  fFuncSigBack->SetRange(fFitMin,fFitMax);
  fFuncBackground->SetRange(fFitMin,fFitMax);
  
  if (!drawOpt.Contains("same")){
    if (fHistDataPM){
      fHistDataPM->Draw();
      grSig->Draw("f");
    } else {
      grSig->Draw("af");
    }
  } else {
    grSig->Draw("f");
  }
  if(fMethod==kFitted || fMethod==kFittedMC) grBack->Draw("f");
  fFuncSigBack->Draw("same");
  fFuncSigBack->SetLineWidth(2);
  if(fMethod==kLikeSign) {
    fHistDataPP->SetLineWidth(2);
    fHistDataPP->SetLineColor(6);
    fHistDataPP->Draw("same");
    fHistDataMM->SetLineWidth(2);
    fHistDataMM->SetLineColor(8);
    fHistDataMM->Draw("same");
  }
  
  if(fMethod==kFitted || fMethod==kFittedMC)
    fFuncBackground->Draw("same");
  
  if (optStat) DrawStats();
  
}
