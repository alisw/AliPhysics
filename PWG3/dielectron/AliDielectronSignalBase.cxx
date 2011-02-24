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

/* $Id$ */

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
#include <TH1F.h>
#include "AliDielectronSignalBase.h"

ClassImp(AliDielectronSignalBase)

AliDielectronSignalBase::AliDielectronSignalBase() :
  TNamed(),
  fHistSignal(0),
  fHistBackground(0),
  fHistDataPM(0),
  fHistDataPP(0),
  fHistDataMM(0),
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
  fScaleFactor(1.),
  fProcessed(kFALSE)
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
  fScaleFactor(1.),
  fProcessed(kFALSE)
{
  //
  // Named Constructor
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
  if(fValues(4)>0) {
    t->AddText(Form("Mass: %.2f #pm %.2f GeV/c^{2}", fValues(4), fErrors(4)));
    t->AddText(Form("Mass res.: %.1f #pm %.1f MeV/c^{2}", 1000*fValues(5), 1000*fErrors(5)));
  }
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
  if(fValues(4)>0){
    printf("Mass: %.5g #pm %.5g\n", fValues(4), fErrors(4));
    printf("Mass res.: %.5g #pm %.5g\n", fValues(5), fErrors(5));
  }
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
