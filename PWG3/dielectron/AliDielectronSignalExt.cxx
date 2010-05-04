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
//                                                                       //
//                      Dielectron SignalExt                             //
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
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>

#include <AliLog.h>

#include "AliDielectronSignalExt.h"

ClassImp(AliDielectronSignalExt)

AliDielectronSignalExt::AliDielectronSignalExt() :
  AliDielectronSignalBase(),
  fSignPM(0x0),
  fSignPP(0x0),
  fSignMM(0x0),
  fBackground(0x0),
  fSignal(0x0),
  fMethod(1),
  fRebin(1),
  fBins(0),
  fDrawMin(0.),
  fDrawMax(0.),
  fFitMin(2.5),
  fFitMax(4)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalExt::AliDielectronSignalExt(const char* name, const char* title) :
  AliDielectronSignalBase(name, title),
  fSignPM(0x0),
  fSignPP(0x0),
  fSignMM(0x0),
  fBackground(0x0),
  fSignal(0x0),
  fMethod(1),
  fRebin(1),
  fBins(0),
  fDrawMin(0.),
  fDrawMax(0.),
  fFitMin(2.5),
  fFitMax(4)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronSignalExt::~AliDielectronSignalExt()
{
  //
  // Default Destructor
  //
  if (fSignPM)     delete fSignPM;
  if (fSignPP)     delete fSignPP;
  if (fSignMM)     delete fSignMM;
  if (fBackground) delete fBackground;
  if (fSignal)     delete fSignal;
}

//______________________________________________
void AliDielectronSignalExt::SetHistograms(TH1F* const unlike, TH1F* const backg, TH1F* const signal)
{
  //
  // set histograms 
  //
  fSignPM = (TH1F*)unlike->Clone("fSignPM");
  fBackground = backg;
  fSignal = signal;
}

//______________________________________________
void AliDielectronSignalExt::Process(TObjArray* const arrhist)
{
  // 
  // signal subtraction. support like-sign subtraction and event mixing method
  //
  switch ( fMethod ){
    case 1 :
      ProcessLS(arrhist);    // process like-sign subtraction method
      break;

    case 2 : 
      ProcessEM(arrhist);    // process event mixing method
      break;

    default :
      AliWarning("Subtraction method not supported. Please check SetMethod() function.");
  }
}

//______________________________________________
void AliDielectronSignalExt::ProcessLS(TObjArray* const arrhist)
{
  //
  // signal subtraction 
  //
  fSignPP = (TH1F*)arrhist->At(0);  // like sign   : plus-plus
  fSignPM = (TH1F*)arrhist->At(1);  // unlike-sign : plus-minus
  fSignMM = (TH1F*)arrhist->At(2);  // like sign   : minus-minus
  fSignPP->Sumw2();
  fSignPM->Sumw2();
  fSignMM->Sumw2();
  
  if ( fRebin>1 ){ Rebin(fRebin); }       // rebinning of histogram

  fBins = fSignPM->GetNbinsX();            // number of bins
  Double_t minX  = fSignPM->GetBinLowEdge(1);       // minimum X value in axis
  Double_t maxX  = fSignPM->GetBinLowEdge(fBins+1); // maximum X value in axis
  
  AliDebug(10,Form("histogram #bin = %d , min = %f , max = %f\n",fBins,minX,maxX));
  TH1F* hBackground = new TH1F("hBackground","Like-sign background",fBins,minX,maxX); 
  TH1F* hSubtracted = new TH1F("hSubtracted","Background subtracted",fBins,minX,maxX);

  // fill out background and subtracted histogram
  for ( Int_t ibin=1; ibin<fBins+1; ibin++ ){
    Double_t mass  = ibin*(maxX-minX)/(Double_t)fBins;
    Double_t backgr = 2.*sqrt( fSignPP->GetBinContent(ibin) * fSignMM->GetBinContent(ibin) );
    Double_t signal = fSignPM->GetBinContent(ibin) - backgr;

    hBackground->Fill(mass, backgr);
    hSubtracted->Fill(mass, signal);
  }
  SetHistograms(fSignPM, hBackground, hSubtracted);


  Double_t signal=0, signal_err=0;
  Double_t background=0, background_err=0;
  
  signal     = fSignal->IntegralAndError(fSignal->FindBin(GetIntegralMin()),
                                         fSignal->FindBin(GetIntegralMax()), signal_err);
  background = fBackground->IntegralAndError(fBackground->FindBin(GetIntegralMin()),
                                             fBackground->FindBin(GetIntegralMax()), background_err);

  //reset result arrays
  Reset();
  //set values
  SetSignal(signal,signal_err);
  SetBackground(background,background_err);
  SetSignificanceAndSOB();

  // cleanup
  //delete hBackground;
  //delete hSubtracted;
}

//______________________________________________
void AliDielectronSignalExt::ProcessEM(TObjArray* const arrhist)
{
  //
  // event mixing. not implemented yet.
  //
  printf("event mixing method is not implemented yet. Like-sign method will be used.\n");
  ProcessLS(arrhist);
}

//______________________________________________
void AliDielectronSignalExt::Rebin(Int_t rebin)
{
  // 
  // rebinning of histograms
  //
  fSignPM->Rebin(rebin);
  fSignPP->Rebin(rebin);
  fSignMM->Rebin(rebin);
}

//______________________________________________
void AliDielectronSignalExt::Draw(const Option_t* option)
{
  //
  // Draw the fitted function
  //
  TString drawOpt(option); 
  drawOpt.ToLower();   

  TCanvas *cSub = new TCanvas("cSub","signal, background subtracted",1400,1000);
  cSub->Divide(2,2);
  cSub->Draw();

  Double_t minX  = fSignPM->GetBinLowEdge(1);       // minimum X value in axis
  Double_t maxX  = fSignPM->GetBinLowEdge(fBins+1); // maximum X value in axis
  if ( TMath::Abs(fDrawMin)<1.e-30 ) fDrawMin = minX;
  if ( TMath::Abs(fDrawMax)<1.e-30 ) fDrawMax = maxX;
  
  cSub->cd(4);
  fSignPM->GetXaxis()->SetRangeUser(fDrawMin, fDrawMax);
  fSignPM->SetLineColor(1);
  fSignPM->SetLineWidth(2);
  fSignPM->SetMarkerStyle(6);
  fSignPM->DrawCopy("P");

  fBackground->SetLineColor(4);
  fBackground->SetMarkerColor(4);
  fBackground->SetMarkerStyle(6);
  fBackground->DrawCopy("Psame");

  fSignal->SetMarkerStyle(20);
  fSignal->SetMarkerColor(2);
  fSignal->DrawCopy("Psame");

  cSub->cd(1);
  fSignal->DrawCopy("P");

  cSub->cd(2);
  fSignPM->DrawCopy("P");

  cSub->cd(3);
  fSignPP->DrawCopy("P");
}

