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
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TLine.h>

#include <AliLog.h>

#include "AliDielectronSignalExt.h"

ClassImp(AliDielectronSignalExt)

AliDielectronSignalExt::AliDielectronSignalExt() :
  AliDielectronSignalBase()
{
  //
  // Default Constructor
  //
}

//______________________________________________
AliDielectronSignalExt::AliDielectronSignalExt(const char* name, const char* title) :
  AliDielectronSignalBase(name, title)
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
}

//______________________________________________
void AliDielectronSignalExt::Process(TObjArray* const arrhist)
{
  // 
  // signal subtraction. support like-sign subtraction and event mixing method
  //
  switch ( fMethod ){
    case kLikeSign :
    case kLikeSignArithm :
      ProcessLS(arrhist);    // process like-sign subtraction method
      break;

    case kEventMixing : 
      ProcessEM(arrhist);    // process event mixing method
      break;

  case kRotation:
      ProcessRotation(arrhist);
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
  fHistDataPP = (TH1*)(arrhist->At(0))->Clone("histPP");  // ++    SE
  fHistDataPM = (TH1*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataMM = (TH1*)(arrhist->At(2))->Clone("histMM");  // --    SE
  fHistDataPP->Sumw2();
  fHistDataPM->Sumw2();
  fHistDataMM->Sumw2();
  fHistDataPP->SetDirectory(0);
  fHistDataPM->SetDirectory(0);
  fHistDataMM->SetDirectory(0);
  
  // rebin the histograms
  if (fRebin>1) { 
    fHistDataPP->Rebin(fRebin);
    fHistDataPM->Rebin(fRebin);
    fHistDataMM->Rebin(fRebin);
  }       

  fHistSignal = new TH1D("HistSignal", "Like-Sign substracted signal",
			 fHistDataPM->GetXaxis()->GetNbins(),
			 fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistSignal->SetDirectory(0);
  fHistBackground = new TH1D("HistBackground", "Like-sign contribution",
			     fHistDataPM->GetXaxis()->GetNbins(),
			     fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground->SetDirectory(0);
  
  // fill out background and subtracted histogram
  for(Int_t ibin=1; ibin<=fHistDataPM->GetXaxis()->GetNbins(); ibin++) {
    Float_t pm = fHistDataPM->GetBinContent(ibin);
    Float_t pp = fHistDataPP->GetBinContent(ibin);
    Float_t mm = fHistDataMM->GetBinContent(ibin);
    Float_t epm = fHistDataPM->GetBinError(ibin);

    Float_t background = 2*TMath::Sqrt(pp*mm);
    Float_t ebackground = TMath::Sqrt(mm+pp);
    if (fMethod==kLikeSignArithm){
      //Arithmetic mean instead of geometric
      background=(pp+mm);
      ebackground=TMath::Sqrt(pp+mm);
      if (TMath::Abs(ebackground)<1e-30) ebackground=1;
    }
//     Float_t signal = pm - background;
//     Float_t error = TMath::Sqrt(epm*epm+mm+pp);

    fHistSignal->SetBinContent(ibin, pm);
    fHistSignal->SetBinError(ibin, epm);
    fHistBackground->SetBinContent(ibin, background);
    fHistBackground->SetBinError(ibin, ebackground);
  }
  //scale histograms to match integral between fScaleMin and fScaleMax
  if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax);

  //subract background
  fHistSignal->Add(fHistBackground,-1);
  
  // signal
  fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
	  			             fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
						 fHistBackground->FindBin(fIntMax), 
						 fErrors(1));
  // S/B and significance
  SetSignificanceAndSOB();

  fProcessed = kTRUE;
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
void AliDielectronSignalExt::ProcessRotation(TObjArray* const arrhist)
{
  //
  // signal subtraction
  //
  fHistDataPM = (TH1*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  if (!fHistDataPM){
    AliError("Unlike sign histogram not available. Cannot extract the signal.");
    return;
  }
  fHistDataPM->Sumw2();

  fHistBackground=(TH1*)arrhist->At(10)->Clone("histRotation");
  if (!fHistBackground){
    AliError("Histgram from rotation not available. Cannot extract the signal.");
    delete fHistDataPM;
    fHistDataPM=0x0;
    return;
  }

  //scale histograms to match integral between fScaleMin and fScaleMax
  if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax);
  
  fHistSignal=(TH1*)fHistDataPM->Clone("histSignal");
  fHistSignal->Add(fHistBackground,-1.);

    // signal
  fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
                                             fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
                                                 fHistBackground->FindBin(fIntMax),
                                                 fErrors(1));
  // S/B and significance
  SetSignificanceAndSOB();
  
  fProcessed = kTRUE;
  
}

//______________________________________________
void AliDielectronSignalExt::Draw(const Option_t* option)
{
  //
  // Draw the fitted function
  //
  TString drawOpt(option); 
  drawOpt.ToLower();   

  Float_t minY = 0.001;
  Float_t maxY = 1.2*fHistDataPM->GetMaximum();
  Float_t minX = 1.001*fHistDataPM->GetXaxis()->GetXmin();
  Float_t maxX = 0.999*fHistDataPM->GetXaxis()->GetXmax();
  Int_t binSize = Int_t(1000*fHistDataPM->GetBinWidth(1));   // in MeV
  Float_t minMinY = fHistSignal->GetMinimum();

  TCanvas *cSub = new TCanvas(Form("%s", fName.Data()),Form("%s", fTitle.Data()),1400,1000);
  cSub->SetLeftMargin(0.15);
  cSub->SetRightMargin(0.0);
  cSub->SetTopMargin(0.002);
  cSub->SetBottomMargin(0.0);
  cSub->Divide(2,2,0.,0.);
  cSub->Draw();

  TVirtualPad* pad = cSub->cd(1);
  pad->SetLeftMargin(0.15);
  pad->SetRightMargin(0.0);
  pad->SetTopMargin(0.005);
  pad->SetBottomMargin(0.0);
  TH2F *range1=new TH2F("range1","",10,minX,maxX,10,minY,maxY);
  range1->SetStats(kFALSE);
  range1->GetYaxis()->SetTitle(Form("entries [counts per %d MeV bin]", binSize));
  range1->GetYaxis()->CenterTitle();
  range1->GetYaxis()->SetLabelSize(0.05);
  range1->GetYaxis()->SetTitleSize(0.06);
  range1->GetYaxis()->SetTitleOffset(0.8);
  range1->Draw();
  fHistDataPM->SetLineColor(1);
  fHistDataPM->SetLineWidth(2);
  //  fHistDataPM->SetMarkerStyle(21);
  fHistDataPM->Draw("Psame");
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.05);
  latex->DrawLatex(0.2, 0.95, "Background un-substracted");
  TLine line;
  line.SetLineWidth(1);
  line.SetLineStyle(2);
  line.DrawLine(fIntMin, minY, fIntMin, maxY);
  line.DrawLine(fIntMax, minY, fIntMax, maxY);

  pad = cSub->cd(2);
  pad->SetLeftMargin(0.);
  pad->SetRightMargin(0.005);
  pad->SetTopMargin(0.005);
  pad->SetBottomMargin(0.0);
  TH2F *range2=new TH2F("range2","",10,minX,maxX,10,minY,maxY);
  range2->SetStats(kFALSE);
  range2->Draw();
  fHistBackground->SetLineColor(4);
  fHistBackground->SetLineWidth(2);
  //  fHistBackground->SetMarkerColor(4);
  //  fHistBackground->SetMarkerStyle(6);
  fHistBackground->Draw("Psame");
  latex->DrawLatex(0.05, 0.95, "Like-sign background");
  line.DrawLine(fIntMin, minY, fIntMin, maxY);
  line.DrawLine(fIntMax, minY, fIntMax, maxY);
  TLegend *legend = new TLegend(0.65, 0.70, 0.98, 0.98);
  legend->SetFillColor(0);
  legend->SetMargin(0.15);
  legend->AddEntry(fHistDataPM, "N_{+-}", "l");
  legend->AddEntry(fHistDataPP, "N_{++}", "l");
  legend->AddEntry(fHistDataMM, "N_{--}", "l");
  legend->AddEntry(fHistSignal, "N_{+-} - 2 #sqrt{N_{++} #times N_{--}}", "l");
  legend->AddEntry(fHistBackground, "2 #sqrt{N_{++} #times N_{--}}", "l");
  legend->Draw();

  
  pad = cSub->cd(3);
  pad->SetLeftMargin(0.15);
  pad->SetRightMargin(0.0);
  pad->SetTopMargin(0.0);
  pad->SetBottomMargin(0.15);
  TH2F *range3=new TH2F("range3","",10,minX,maxX,10,minMinY,maxY);
  range3->SetStats(kFALSE);
  range3->GetYaxis()->SetTitle(Form("entries [counts per %d MeV bin]", binSize));
  range3->GetYaxis()->CenterTitle();
  range3->GetYaxis()->SetLabelSize(0.05);
  range3->GetYaxis()->SetTitleSize(0.06);
  range3->GetYaxis()->SetTitleOffset(0.8);
  range3->GetXaxis()->SetTitle("inv. mass [GeV/c^{2}]");
  range3->GetXaxis()->CenterTitle();
  range3->GetXaxis()->SetLabelSize(0.05);
  range3->GetXaxis()->SetTitleSize(0.06);
  range3->GetXaxis()->SetTitleOffset(1.0);
  range3->Draw();
  fHistDataPM->Draw("Psame");
  fHistDataPP->SetLineWidth(2);
  fHistDataPP->SetLineColor(6);
  fHistDataMM->SetLineWidth(2);
  fHistDataMM->SetLineColor(8);
  fHistDataPP->Draw("Psame");
  fHistDataMM->Draw("Psame");
  line.DrawLine(minX, 0.,maxX, 0.);
  line.DrawLine(fIntMin, minMinY, fIntMin, maxY);
  line.DrawLine(fIntMax, minMinY, fIntMax, maxY);

  pad = cSub->cd(4);
  pad->SetLeftMargin(0.0);
  pad->SetRightMargin(0.005);
  pad->SetTopMargin(0.0);
  pad->SetBottomMargin(0.15);
  TH2F *range4=new TH2F("range4","",10,minX,maxX,10,minMinY,maxY);
  range4->SetStats(kFALSE);
  range4->GetXaxis()->SetTitle("inv. mass [GeV/c^{2}]");
  range4->GetXaxis()->CenterTitle();
  range4->GetXaxis()->SetLabelSize(0.05);
  range4->GetXaxis()->SetTitleSize(0.06);
  range4->GetXaxis()->SetTitleOffset(1.0);
  range4->Draw();
  fHistSignal->SetLineWidth(2);
  fHistSignal->SetLineColor(2);
  fHistSignal->Draw("Psame");
  latex->DrawLatex(0.05, 0.95, "Like-sign background substracted");
  if(fProcessed) DrawStats(0.05, 0.6, 0.5, 0.9);
  line.DrawLine(minX, 0.,maxX, 0.);
  line.DrawLine(fIntMin, minMinY, fIntMin, maxY);
  line.DrawLine(fIntMax, minMinY, fIntMax, maxY);

  cSub->SaveAs(Form("%s_summary.png", fName.Data()));
}

