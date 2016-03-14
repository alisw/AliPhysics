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
//                Dielectron SignalBase                                  //
//                                                                       //
//                                                                       //
/*
Base class for signal extraction from a histogram or an array of histograms
The histogram is assumed to be an inv. mass spectrum,
the array of histograms is assumed to be an array with inv. mass histograms
resulting from single and mixed events, as defined in AliDielectron.cxx

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <TVectorT.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TMath.h>
#include <TH1.h>
#include <TDatabasePDG.h>

#include "AliDielectronSignalFunc.h"
#include "AliDielectronSignalBase.h"

ClassImp(AliDielectronSignalBase)

const char* AliDielectronSignalBase::fgkValueNames[6] = {
  "Signal","Background","Significance","Signal/Background","Mass","MassWidth"};
const Double_t AliDielectronSignalBase::fgkErrorZero = 0.5 * TMath::ChisquareQuantile(0.6827,2);
  
AliDielectronSignalBase::AliDielectronSignalBase() :
  TNamed(),
  fHistSignal(0),
  fHistBackground(0),
  fHistDataPM(0),
  fHistDataPP(0),
  fHistDataMM(0),
  fHistDataME(0),
  fHistRfactor(0),
  fPeakShapeObj(0x0),
  fHistSimPM(0x0),
  fValues(6),
  fErrors(6),
  fIntMin(0),
  fIntMax(0),
  fFitMin(0),
  fFitMax(0),
  fRebin(1),
  fMethod(kLikeSign),
  fScaleMin(0.),
  fScaleMax(0.),
  fScaleMin2(0.),
  fScaleMax2(0.),
  fScaleFactor(1.),
  fMixingCorr(kFALSE),
  fPeakMethod(kBinCounting),
  fProcessed(kFALSE),
  fPOIpdg(443)

{
  //
  // Default Constructor
  //
}

//______________________________________________
AliDielectronSignalBase::AliDielectronSignalBase(const char* name, const char* title) :
  TNamed(name, title),
  fHistSignal(0),
  fHistBackground(0),
  fHistDataPM(0),
  fHistDataPP(0),
  fHistDataMM(0),
  fHistDataME(0),
  fHistRfactor(0),
  fPeakShapeObj(0x0),
  fHistSimPM(0x0),
  fValues(6),
  fErrors(6),
  fIntMin(0),
  fIntMax(0),
  fFitMin(0),
  fFitMax(0),
  fRebin(1),
  fMethod(kLikeSign),
  fScaleMin(0.),
  fScaleMax(0.),
  fScaleMin2(0.),
  fScaleMax2(0.),
  fScaleFactor(1.),
  fMixingCorr(kFALSE),
  fPeakMethod(kBinCounting),
  fProcessed(kFALSE),
  fPOIpdg(443)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronSignalBase::AliDielectronSignalBase(const char* name, const char* title, bool enummaps) :
  TNamed(name, title),
  fHistSignal(0),
  fHistBackground(0),
  fHistDataPM(0),
  fHistDataPP(0),
  fHistDataMM(0),
  fHistDataME(0),
  fHistRfactor(0),
  fPeakShapeObj(0x0),
  fHistSimPM(0x0),
  fValues(6),
  fErrors(6),
  fIntMin(0),
  fIntMax(0),
  fFitMin(0),
  fFitMax(0),
  fRebin(1),
  fMethod(kLikeSign),
  fScaleMin(0.),
  fScaleMax(0.),
  fScaleMin2(0.),
  fScaleMax2(0.),
  fScaleFactor(1.),
  fMixingCorr(kFALSE),
  fPeakMethod(kBinCounting),
  fProcessed(kFALSE),
  fPOIpdg(443),
  MapBackgroundMethod(),
  MapSignalExtractionMethod()
{
  if(enummaps){
    AliDielectronSignalBase::SetBackgroundEnumMap();
    AliDielectronSignalBase::SetSignalExtractionEnumMap();
  }
  //
  // Named + set enum maps Constructor
  //
}
//______________________________________________
AliDielectronSignalBase::~AliDielectronSignalBase()
{
  //
  // Default Destructor
  //
  if(fHistSignal) delete fHistSignal;
  if(fHistBackground) delete fHistBackground;
  if (fHistDataPP) delete fHistDataPP;
  if (fHistDataPM) delete fHistDataPM;
  if (fHistDataMM) delete fHistDataMM;
  if (fHistDataME) delete fHistDataME;
  if (fHistRfactor)delete fHistRfactor;  
  if (fPeakShapeObj) delete fPeakShapeObj;
  if (fHistSimPM) delete fHistSimPM;
}

//______________________________________________
TPaveText* AliDielectronSignalBase::DrawStats(Double_t x1/*=0.*/, Double_t y1/*=0.*/, Double_t x2/*=0.*/, Double_t y2/*=0.*/)
{
  //
  // Draw extracted values in a TPaveText
  // with the corners x1,y2,x2,y2
  //
  if (TMath::Abs(x1)<1e-20&&TMath::Abs(x2)<1e-20){
    x1=.6;
    x2=.9;
    y1=.7;
    y2=.9;
  }
  TPaveText *t=new TPaveText(x1,y1,x2,y2,"brNDC");
  t->SetFillColor(kWhite);
  t->SetBorderSize(1);
  t->SetTextAlign(12);
  t->AddText(Form("Range  : %.2f - %.2f GeV/c^{2}", fIntMin, fIntMax));
  t->AddText(Form("Signal : %.1f #pm %.1f", fValues(0), fErrors(0)));
  t->AddText(Form("Backgnd: %.1f #pm %.1f", fValues(1), fErrors(1)));
  t->AddText(Form("Signif.: %.2f #pm %.2f", fValues(2), fErrors(2)));
  t->AddText(Form("S/B    : %.2f #pm %.2f", fValues(3), fErrors(3)));
  if(fValues(4)>0)
    t->AddText(Form("Mass: %.2f #pm %.2f GeV/c^{2}", fValues(4), fErrors(4)));
  if(fValues(5)>0)
    t->AddText(Form("Mass res.: %.1f #pm %.1f MeV/c^{2}", 1000*fValues(5), 1000*fErrors(5)));
  t->Draw();

  return t;
}

//______________________________________________
void AliDielectronSignalBase::Print(Option_t */*option*/) const
{
  //
  // Print the statistics
  //
  printf("Signal : %.5g #pm %.5g\n",fValues(0), fErrors(0));
  printf("Backgnd: %.5g #pm %.5g\n",fValues(1), fErrors(1));
  printf("Signif.: %.5g #pm %.5g\n",fValues(2), fErrors(2));
  printf("SoB    : %.5g #pm %.5g\n",fValues(3), fErrors(3));
  if(fValues(4)>0)
    printf("Mass: %.5g #pm %.5g\n", fValues(4), fErrors(4));
  if(fValues(5)>0)
    printf("Mass res.: %.5g #pm %.5g\n", fValues(5), fErrors(5));
}

//______________________________________________
Double_t AliDielectronSignalBase::ScaleHistograms(TH1* histRaw, TH1* histBackground, Double_t intMin, Double_t intMax)
{
  //
  // scale histBackground to match the integral of histRaw in the interval intMin, intMax
  //

  //protect using over and underflow bins in normalisation calculation
  if (intMin<histRaw->GetXaxis()->GetXmin()) intMin=histRaw->GetXaxis()->GetXmin();
  if (intMin<histBackground->GetXaxis()->GetXmin()) intMin=histBackground->GetXaxis()->GetXmin();
  
  if (intMax>histRaw->GetXaxis()->GetXmax())
    intMax=histRaw->GetXaxis()->GetXmax()-histRaw->GetBinWidth(histRaw->GetNbinsX())/2.;
  if (intMax>histBackground->GetXaxis()->GetXmax())
    intMax=histBackground->GetXaxis()->GetXmax()-histBackground->GetBinWidth(histBackground->GetNbinsX())/2.;
  
  Double_t intRaw  = histRaw->Integral(histRaw->FindBin(intMin),histRaw->FindBin(intMax));
  Double_t intBack = histBackground->Integral(histBackground->FindBin(intMin),histBackground->FindBin(intMax));
  Double_t scaleFactor=intBack>0?intRaw/intBack:0.;
  if (intBack>0){
    histBackground->Sumw2();
    histBackground->Scale(scaleFactor);
  }

  return scaleFactor;
}
//______________________________________________
Double_t AliDielectronSignalBase::ScaleHistograms(TH1* histRaw, TH1* histBackground, Double_t intMin, Double_t intMax, Double_t intMin2, Double_t intMax2)
{
  //
  // scale histBackground to match the integral of histRaw in the interval intMin, intMax and intMin2, intMax2
  //

  if(intMin2==intMax2) return (ScaleHistograms(histRaw, histBackground, intMin, intMax));

  //protect using over and underflow bins in normalisation calculation
  if (intMin<histRaw->GetXaxis()->GetXmin()) intMin=histRaw->GetXaxis()->GetXmin();
  if (intMin<histBackground->GetXaxis()->GetXmin()) intMin=histBackground->GetXaxis()->GetXmin();
  
  if (intMax2>histRaw->GetXaxis()->GetXmax())
    intMax2=histRaw->GetXaxis()->GetXmax()-histRaw->GetBinWidth(histRaw->GetNbinsX())/2.;
  if (intMax2>histBackground->GetXaxis()->GetXmax())
    intMax2=histBackground->GetXaxis()->GetXmax()-histBackground->GetBinWidth(histBackground->GetNbinsX())/2.;
  
  Double_t intRaw  = histRaw->Integral(histRaw->FindBin(intMin),histRaw->FindBin(intMax));
  Double_t intBack = histBackground->Integral(histBackground->FindBin(intMin),histBackground->FindBin(intMax));
  intRaw  += histRaw->Integral(histRaw->FindBin(intMin2),histRaw->FindBin(intMax2));
  intBack += histBackground->Integral(histBackground->FindBin(intMin2),histBackground->FindBin(intMax2));

  Double_t scaleFactor=intBack>0?intRaw/intBack:0.;
  if (intBack>0){
    histBackground->Sumw2();
    histBackground->Scale(scaleFactor);
  }

  return scaleFactor;
}

TObject* AliDielectronSignalBase::DescribePeakShape(ESignalExtractionMethod method, Bool_t replaceValErr,  TH1F *mcShape) {
  //
  // Describe the extracted peak by the selected method and overwrite signal etc if needed
  //
  AliDielectronSignalFunc* SignalFuncObj = new AliDielectronSignalFunc();

  fPeakMethod=method;
  Double_t data=0.;
  Double_t mc=0.;
  Double_t massPOI=TDatabasePDG::Instance()->GetParticle(fPOIpdg)->Mass();
  Double_t nPOI     = fHistSignal->GetBinContent(fHistSignal->FindBin(massPOI));
  Double_t binWidth = fHistSignal->GetBinWidth(  fHistSignal->FindBin(massPOI));
  TF1 *fit=0x0;
  Int_t parMass =-1;
  Int_t parSigma=-1;

  // do the scaling/fitting
  switch(fPeakMethod) {
  case kBinCounting: /*nothing needs to be done*/ 
    break;

  case kMCScaledMax:
    if(!mcShape) { printf(" ERROR: No MC histogram passed. Returning. \n"); return 0x0; }
    data = fHistSignal->GetBinContent(fHistSignal->FindBin(massPOI));
    mc   = mcShape->GetBinContent(fHistSignal->FindBin(massPOI));
    mcShape->Scale(data / mc );
    break;

  case kMCScaledInt:
    if(!mcShape) { printf(" ERROR: No MC histogram passed. Returning. \n"); return 0x0; }
    if(mcShape->GetBinWidth(1)!=fHistSignal->GetBinWidth(1)) 
      printf(" WARNING: MC and signal histogram have different bin widths. \n");
    data = fHistSignal->Integral(fHistSignal->FindBin(fIntMin),fHistSignal->FindBin(fIntMax));
    mc   = mcShape->Integral(mcShape->FindBin(fIntMin),mcShape->FindBin(fIntMax));
    mcShape->Scale(data / mc );
    break;

  case kMCFitted:
    if(!mcShape && !fHistSimPM) { printf(" ERROR: No MC histogram passed or set. Returning. \n"); return 0x0; }
    if(!fHistSimPM) fHistSimPM=mcShape;
    fit = new TF1("fitMC",SignalFuncObj,&AliDielectronSignalFunc::PeakFunMC,fFitMin,fFitMax,1);
    fit->SetParNames("N");
    fHistSignal->Fit(fit,"RNI0");
    break;

  case kCrystalBall:
    fit = new TF1("fitCB",SignalFuncObj,&AliDielectronSignalFunc::PeakFunCB,fFitMin,fFitMax,5);
    fit->SetParNames("alpha","n","meanx","sigma","N");
    //  fit->SetParameters(-.2,5.,gMjpsi,.06,20);
    //  fit->SetParameters(1.,3.6,gMjpsi,.08,700);
    fit->SetParameters(0.4, 4.0, massPOI, 0.025, 1.3*nPOI);
    fit->SetParLimits(0, 0.0,           1.           );
    fit->SetParLimits(1, 0.01,          10.          );
    fit->SetParLimits(2, massPOI-0.02,  massPOI+0.02 );
    fit->SetParLimits(3, 0.001,          0.2         );
    fit->SetParLimits(4, 0.2*nPOI,      2.0*nPOI     );
    parMass=2;
    parSigma=3;
    fHistSignal->Fit(fit,"RNI0");
    break;

  case kGaus:
    fit = new TF1("fitGaus",SignalFuncObj,&AliDielectronSignalFunc::PeakFunGaus,fFitMin,fFitMax,3);
    //fit = new TF1("fitGaus","gaus",fFitMin,fFitMax);
    fit->SetParNames("N","meanx","sigma");
    fit->SetParameters(1.3*nPOI, massPOI, 0.025);
    fit->SetParLimits(0, 0.2*nPOI,      2.0*nPOI     );
    fit->SetParLimits(1, massPOI-0.02, massPOI+0.02);
    fit->SetParLimits(2, 0.001,         1.           );
    parMass=1;
    parSigma=2;
    fHistSignal->Fit(fit,"RNI0");
    break;

  }

  // overwrite values and errors if requested
  if(replaceValErr) {
    switch(fPeakMethod) {
    case kBinCounting:
      fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin), fHistSignal->FindBin(fIntMax), fErrors(0));
      SetSignificanceAndSOB();
      break;
    case kMCScaledMax:
    case kMCScaledInt:
      fValues(0) = mcShape->IntegralAndError(mcShape->FindBin(fIntMin), mcShape->FindBin(fIntMax), fErrors(0));
      SetSignificanceAndSOB();
      break;

    case kMCFitted:
    case kCrystalBall:
    case kGaus:
      fValues(0) = fit->Integral(fIntMin, fIntMax)/binWidth;
      fErrors(0) = fit->IntegralError(fIntMin, fIntMax)/binWidth;
      SetSignificanceAndSOB();
      break;
    }
    // set mass position and width
    if(parMass>=0) {
      fValues(4) = fit->GetParameter(parMass);
      fErrors(4) = fit->GetParError(parMass);
    }
    if(parSigma>=0) {
      fValues(5) = fit->GetParameter(parSigma);
      fErrors(5) = fit->GetParError(parSigma);
    }
    else {
      // calculate FWHM
      SetFWHM();
    }
  }

  // set the peak method obj
  switch(fPeakMethod) {
  case kBinCounting:
    fPeakShapeObj=(TH1F*)fHistSignal->Clone("BinCount");
    break;
  case kMCScaledMax:
  case kMCScaledInt:
    fPeakShapeObj=mcShape;
    break;
  case kMCFitted:
  case kCrystalBall:
  case kGaus:
    fPeakShapeObj=fit;
    break;
  }

  return fPeakShapeObj;

}
