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

  Class used for extracting the signal from an invariant mass spectrum.
  It implements the AliDielectronSignalBase and -Ext classes and it uses user provided
  functions to fit the spectrum with a combined signa+background fit.
  Used invariant mass spectra are provided via an array of histograms. There are serveral method
  to estimate the background and to extract the raw yield from the background subtracted spectra.

  Example usage:

  AliDielectronSignalFunc *sig = new AliDielectronSignalFunc();


  1) invariant mass input spectra

  1.1) Assuming a AliDielectronCF container as data format (check class for more details)
  AliDielectronCFdraw *cf = new AliDielectronCFdraw("path/to/the/output/file.root");
  TObjArray *arrHists = cf->CollectMinvProj(cf->FindStep("Config"));

  1.2) Assuming a AliDielectronHF grid as data format (check class for more details)
  AliDielectronHFhelper *hf = new AliDielectronHFhelper("path/to/the/output/file.root", "ConfigName");
  TObjArray *arrHists = hf->CollectHistos(AliDielectronVarManager::kM);

  1.3) Assuming a single histograms
  TObjArray *histoArray = new TObjArray();
  arrHists->Add(signalPP);            // add the spectrum histograms to the array
  arrHists->Add(signalPM);            // the order is important !!!
  arrHists->Add(signalMM);


  2) background estimation

  2.1) set the method for the background estimation (methods can be found in AliDielectronSignalBase)
  sig->SetMethod(AliDielectronSignalBase::kFitted);
  2.2) rebin the spectras if needed
  //  sig->SetRebin(2);
  2.3) add any background function you like
  TF1 *fB = new TF1("fitBgrd","pol3",minFit,maxFit);


  3) configure the signal extraction

  3.1) chose one of the signal functions (MCshape, CrystalBall, Gauss)
  TF1 *fS = new TF1("fitSign",AliDielectronSignalFunc::PeakFunCB,minFit,maxFit,5); // has 5 parameters
  //  TF1 *fS = new TF1("fitSign",AliDielectronSignalFunc::PeakFunGaus,minFit,maxFit,3); // has 3 parameters
  //  sig->SetMCSignalShape(hMCsign);
  //  TF1 *fS = new TF1("fitSign",AliDielectronSignalFunc::PeakFunMC,minFit,maxFit,1); // requires a MC shape
  3.2) set the method for the signal extraction (methods can be found in AliDielectronSignalBase)
  depending on the method serveral inputs are needed (e.g. MC shape, PDG code of the particle of interest)
  //  sig->SetParticleOfInterest(443); //default is jpsi
  //  sig->SetMCSignalShape(signalMC);
  //  sig->SetIntegralRange(minInt, maxInt);
  sig->SetExtractionMethod(AliDielectronSignal::BinCounting); // this is the default


  4) combined fit of bgrd+signal

  4.1) combine the two functions
  sig->CombineFunc(fS,fB);
  4.2) apply fitting ranges and the fit options
  sig->SetFitRange(minFit, maxFit);
  sig->SetFitOption("NR");


  5) start the processing

  sig->Process(arrHists);
  sig->Print(""); // print values and errors extracted


  6) access the spectra and values created

  6.1) standard spectra as provided filled in AliDielectronSignalExt
  TH1F *hsign = (TH1F*) sig->GetUnlikeSignHistogram();  // same as the input (rebinned)
  TH1F *hbgrd = (TH1F*) sig->GetBackgroundHistogram();  // filled histogram with fitBgrd
  TH1F *hextr = (TH1F*) sig->GetSignalHistogram();      // after backgound extraction (rebinned)
  TObject *oPeak = (TObject*) sig->GetPeakShape();      // can be a TF1 or TH1 depending on the method
  6.2) flow spectra
  TF1 *fFitSign  = sig->GetCombinedFunction();                // combined fit function
  TF1 *fFitExtr  = sig->GetSignalFunction();                  // signal function
  TF1 *fFitBgrd  = sig->GetBackgroundFunction();              // background function
  6.3) access the extracted values and errors
  sig->GetValues();     or GetErrors();                 // yield extraction

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>
#include <TList.h>
#include <TFitResult.h>
//#include <../hist/hist/src/TF1Helper.h>

#include <AliLog.h>

#include "AliDielectronSignalFunc.h"
#include "AliDielectron.h"

ClassImp(AliDielectronSignalFunc)

TF1*  AliDielectronSignalFunc::fFuncSignal=0x0;
TF1*  AliDielectronSignalFunc::fFuncBackground=0x0;
Int_t AliDielectronSignalFunc::fNparPeak=0;
Int_t AliDielectronSignalFunc::fNparBgnd=0;


AliDielectronSignalFunc::AliDielectronSignalFunc() :
AliDielectronSignalExt(),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE),
fDof(0),
fChi2Dof(0.0),
fHistCombinatorialBackground(0x0)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalFunc::AliDielectronSignalFunc(const char* name, const char* title) :
AliDielectronSignalExt(name, title),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE),
fDof(0),
fChi2Dof(0.0),
fHistCombinatorialBackground(0x0)
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
  if(fFuncSigBack) delete fFuncSigBack;
}


//______________________________________________
void AliDielectronSignalFunc::Process(TObjArray * const arrhist)
{
  //
  // Fit the invariant mass histograms and retrieve the signal and background
  //
  switch(fMethod) {
  case kCombinatorialPlusFit:
    ProcessCombinatorialPlusFit(arrhist);
    break;
    
  case kFitted :
    ProcessFit(arrhist);
    break;

  case kLikeSignFit :
    ProcessFitLS(arrhist);
    break;

  case kEventMixingFit :
    ProcessFitEM(arrhist);
    break;

  case kLikeSign :
  case kLikeSignArithm :
  case kLikeSignRcorr:
  case kLikeSignArithmRcorr:
    ProcessLS(arrhist);    // process like-sign subtraction method
    break;

  case kEventMixing :
    ProcessEM(arrhist);    // process event mixing method
    break;

  case kRotation:
    ProcessRotation(arrhist);
    break;

  default :
    AliError("Background substraction method not known!");
  }
}
//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakFunMC(const Double_t *x, const Double_t *par) {
  // Fit MC signal shape
  // parameters
  // [0]:   scale for simpeak
  
  Double_t xx  = x[0];
  
  TH1F *hPeak = fHistSimPM;
  if (!hPeak) {
    printf("E-AliDielectronSignalFunc::PeakFun: No histogram for peak fit defined!\n");
    return 0.0;
  }
  
  Int_t idx = hPeak->FindBin(xx);
  if ((idx <= 1) ||
      (idx >= hPeak->GetNbinsX())) {
    return 0.0;
  }
  
  return (par[0] * hPeak->GetBinContent(idx));
  
}

//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakFunCB(const Double_t *x, const Double_t *par) {
  // Crystal Ball fit function
  
  Double_t alpha = par[0];
  Double_t     n = par[1];
  Double_t meanx = par[2];
  Double_t sigma = par[3];
  Double_t    nn = par[4];
   
  Double_t a = TMath::Power((n/TMath::Abs(alpha)), n) * TMath::Exp(-.5*alpha*alpha);
  Double_t b = n/TMath::Abs(alpha) - TMath::Abs(alpha);
  
  Double_t arg = (x[0] - meanx)/sigma;
  Double_t fitval = 0.;
  
  if (arg > -1.*alpha) {
    fitval = nn * TMath::Exp(-.5*arg*arg);
  } else {
    fitval = nn * a * TMath::Power((b-arg), (-1.*n));
  }
  
  return fitval;
}

//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakFunGaus(const Double_t *x, const Double_t *par) {
  // Gaussian fit function
  //printf("fNparBgrd %d \n",fNparBgnd);
  Double_t     n = par[0];
  Double_t  mean = par[1];
  Double_t sigma = par[2];
  Double_t    xx = x[0];

  return ( n*TMath::Exp(-0.5*TMath::Power((xx-mean)/sigma,2)) );
}

//______________________________________________
void AliDielectronSignalFunc::ProcessFit(TObjArray * const arrhist) {
  //
  // Fit the +- invariant mass distribution only
  // Here we assume that the combined fit function is a sum of the signal and background functions
  //    and that the signal function is always the first term of this sum
  //
  
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataPM->SetBinErrorOption(TH1::kPoisson);
  fHistDataPM->Sumw2(kFALSE);
  if(fRebin>1)
    fHistDataPM->Rebin(fRebin);
  
//   for(Int_t ibin=1; ibin<=fHistDataPM->GetXaxis()->GetNbins(); ibin++) {
//     if(fHistDataPM->GetBinError(ibin)<1e-30 ) fHistDataPM->SetBinError(ibin, fgkErrorZero);
//   }
  fHistSignal = new TH1F("HistSignal", "Fit substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground = new TH1F("HistBackground", "Fit contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  // the starting parameters of the fit function and their limits can be tuned
  // by the user in its macro
//   fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  TFitResultPtr pmFitPtr = fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  TFitResult *pmFitResult = pmFitPtr.Get(); // used only with TF1Helper
  if(pmFitResult){
    fDof =  pmFitResult->Ndf();
    if(fDof) fChi2Dof = pmFitResult->Chi2() / fDof;
  }
  fFuncSignal->SetParameters(fFuncSigBack->GetParameters());
  fFuncBackground->SetParameters(fFuncSigBack->GetParameters()+fFuncSignal->GetNpar());
  // fill the background spectrum

  fHistBackground->Eval(fFuncBackground);
  // set the error for the background histogram
  fHistBackground->Fit(fFuncBackground,"0qR","",fFitMin,fFitMax);
  Double_t inte  = fFuncBackground->IntegralError(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
  Double_t binte = inte / TMath::Sqrt((fHistDataPM->FindBin(fIntMax)-fHistDataPM->FindBin(fIntMin))+1);
  for(Int_t iBin=fHistDataPM->FindBin(fIntMin); iBin<=fHistDataPM->FindBin(fIntMax); iBin++) {
    fHistBackground->SetBinError(iBin, binte);
  }

  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    //    Double_t m = fHistDataPM->GetBinCenter(iBin);
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t epm = fHistDataPM->GetBinError(iBin);
    Double_t bknd = fHistBackground->GetBinContent(iBin);
    Double_t ebknd = fHistBackground->GetBinError(iBin);
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
    Double_t error = TMath::Sqrt(epm*epm+ebknd*ebknd);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, error);
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
    fErrors(1) = fFuncBackground->IntegralError(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
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
  
  fHistBackground->GetListOfFunctions()->Add(fFuncBackground);
  fHistDataPM->GetListOfFunctions()->Add(fFuncSigBack);
  fHistSignal->GetListOfFunctions()->Add(fFuncSignal);

}


//______________________________________________
void AliDielectronSignalFunc::ProcessCombinatorialPlusFit(TObjArray * const arrhist) {
  //
  // Describe the backgorund by a combinatorial background historgram, multiplied with a fit function
  // e.g. 1 + a * exp( -m / m_0 )
  // the fit function describes the correlated bg.

  // Example usage in: PWGDQ/dielectron/macrosJPSI/multiplicity13TeV/example_hybrid.C



  
  fHistDataPM = (TH1*)(arrhist->At(AliDielectron::kEv1PM))->Clone("histPM");  // +-    SE
  fHistCombinatorialBackground = (TH1*)(arrhist->At(AliDielectron::kEv1PMRot))->Clone("histCombinatorial");
  
  
  fHistDataPM->Sumw2(kFALSE);
  fHistDataPM->SetBinErrorOption(TH1::kPoisson);
  
  
  
  
  if(fHistCombinatorialBackground->GetDefaultSumw2()) fHistCombinatorialBackground->Sumw2();
  fHistDataPM->SetDirectory(0);
  fHistCombinatorialBackground->SetDirectory(0);
  

  // rebin the histograms
  if (fRebin>1) {
    fHistDataPM->Rebin(fRebin);
    fHistCombinatorialBackground->Rebin(fRebin);
  }
  
  
  fHistBackground = new TH1F("HistBackground", "Fit contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());

  
  // Step 1 : scale TR to OS
  // this will be the histogram containing the combinatorial background
  
  if (fScaleMax>fScaleMin && fScaleMax2>fScaleMin2) fScaleFactor=ScaleHistograms(fHistDataPM,fHistCombinatorialBackground,fScaleMin,fScaleMax,fScaleMin2,fScaleMax2);
  else if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistCombinatorialBackground,fScaleMin,fScaleMax);
  else if (fScaleMin>0.){
    fScaleFactor=fScaleMin;
    fHistCombinatorialBackground->Scale(fScaleFactor);
  }

  // Step 3 : Make Ratio of pm histo and combinatorial background to determine the correlated bg
  
  TGraphAsymmErrors* hRatioPMtoCombinatorial = new TGraphAsymmErrors();
  
  Int_t  excludeFrom  = fHistDataPM->GetXaxis()->FindBin( 2.2 );
  Double_t excludeUntil = fHistDataPM->GetXaxis()->FindBin( 3.2 );
  Double_t ipoint =0;
  
  for( int ibin = 1; ibin < fHistDataPM->GetXaxis()->GetNbins(); ++ibin ){
    Bool_t exclude =  (ibin > excludeFrom && ibin < excludeUntil );
    
    if(!exclude){
      Double_t x = fHistDataPM->GetXaxis()->GetBinCenter(ibin);
      Double_t exl = fHistDataPM->GetXaxis()->GetBinWidth(ibin) / 2.;
      Double_t exh = exl;
      
      Double_t denominator = fHistCombinatorialBackground->GetBinContent(ibin) > 0. ? fHistCombinatorialBackground->GetBinContent(ibin) : 1.;
      
      Double_t y = fHistDataPM->GetBinContent(ibin) / denominator ;
      Double_t eyl = fHistDataPM->GetBinErrorLow(ibin)/ denominator;
      Double_t eyh = fHistDataPM->GetBinErrorUp(ibin) / denominator;

      
      hRatioPMtoCombinatorial->SetPoint( ipoint, x, y ) ;
      hRatioPMtoCombinatorial->SetPointError( ipoint, exl, exh, eyl, eyh) ; 
      ipoint++;
    }
    
  }
  // Step 4 : Fit correlated bg
  
  fFuncBackground->SetParameters(fFuncSigBack->GetParameters()+fFuncSignal->GetNpar());
  TFitResultPtr pmFitPtr = hRatioPMtoCombinatorial->Fit(fFuncBackground, fFitOpt.Data(), "", fFitMin, fFitMax);
  TFitResult *pmFitResult = pmFitPtr.Get(); // used only with TF1Helper
  if(pmFitResult){
    fDof =  pmFitResult->Ndf();
    if(fDof) fChi2Dof = pmFitResult->Chi2() / fDof;
  }
  
  
  TH1D hCorrelatedBackground ("HistCorrelatedBg", "Fit substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  hCorrelatedBackground.Sumw2();
  
  hCorrelatedBackground.Eval(fFuncBackground);
  
  
  for( int ibin = 1; ibin < fHistDataPM->GetXaxis()->GetNbins(); ++ibin ) {
    Double_t inte  = fFuncBackground->IntegralError(  fHistDataPM->GetBinLowEdge(ibin), fHistDataPM->GetBinLowEdge(ibin+1) ) / fHistDataPM->GetBinWidth(1);
    hCorrelatedBackground.SetBinError(ibin, inte);
  }

  
  // Step 5 : bg = combinatorial bg * fit of correlated bg
  
  fHistBackground = (TH1*) fHistCombinatorialBackground->Clone("histBackground");
  fHistBackground->Multiply( &hCorrelatedBackground );
  fHistSignal = new TH1F("HistSignal", "Fit substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t epm = (fHistDataPM->GetBinError(iBin) < 1e-30 ? fgkErrorZero : fHistDataPM->GetBinError(iBin));
    Double_t bknd = fHistBackground->GetBinContent(iBin);
    Double_t ebknd = fHistBackground->GetBinError(iBin);

    Double_t signal = pm-bknd;
    Double_t error = TMath::Sqrt(epm*epm+ebknd*ebknd);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, error);
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
void AliDielectronSignalFunc::ProcessFitLS(TObjArray * const arrhist) {
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
void AliDielectronSignalFunc::ProcessFitEM(TObjArray * const arrhist) {
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


//______________________________________________________________________________
void AliDielectronSignalFunc::CombineFunc(TF1 * const peak, TF1 * const bgnd) {
  //
  // combine the bgnd and the peak function
  //

  if (!peak||!bgnd) {
    AliError("Both, signal and background function need to be set!");
    return;
  }
  fFuncSignal=peak;
  fFuncBackground=bgnd;

  fNparPeak     = fFuncSignal->GetNpar();
  fNparBgnd     = fFuncBackground->GetNpar();

  fFuncSigBack = new TF1("BgndPeak",AliDielectronSignalFunc::PeakBgndFun, fFitMin,fFitMax, fNparPeak+fNparBgnd);
  return;
}

//______________________________________________________________________________
Double_t AliDielectronSignalFunc::PeakBgndFun(const Double_t *x, const Double_t *par) {
  //
  // merge peak and bgnd functions
  //
  for(int i = 0 ; i < fNparPeak; ++i ) fFuncSignal->SetParameter(i, par[i]);
  for(int i = fNparPeak ; i < fNparPeak+fNparBgnd; ++i ) fFuncSignal->SetParameter(i, par[i]);
  return (fFuncSignal->Eval(x[0]) + fFuncBackground->Eval(x[0]) );
}

