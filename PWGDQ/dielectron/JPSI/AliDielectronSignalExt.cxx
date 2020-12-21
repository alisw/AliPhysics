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

  Class used for extracting the signal from an invariant mass spectrum.
  Used invariant mass spectra are provided via an array of histograms. There are serveral method
  to estimate the background and to extract the raw yield from the background subtracted spectra.

  Example usage:

  AliDielectronSignalExt *sig = new AliDielectronSignalExt();


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
  sig->SetMethod(AliDielectronSignalBase::kEventMixing);
  2.2) rebin the spectras if needed
  //  sig->SetRebin(2);
  2.3) normalize the backgound spectum to the odd-sign spectrum in the desired range(s)
  sig->SetScaleRawToBackground(minScale, maxScale);
  //  sig->SetScaleRawToBackground(minScale, maxScale, minScale2, maxScale2);


  3) configure the signal extraction

  3.1) set the method for the signal extraction (methods can be found in AliDielectronSignalBase)
  depending on the method serveral inputs are needed (e.g. MC shape, PDG code of the particle of interest)
  //  sig->SetParticleOfInterest(443); //default is jpsi
  //  sig->SetMCSignalShape(signalMC);
  sig->SetIntegralRange(minInt, maxInt);  // range for bin counting
  sig->SetExtractionMethod(AliDielectronSignal::BinCounting); // this is the default


  4) start the processing

  sig->Process(arrHists);
  sig->Print(""); // print values and errors extracted


  5) access the spectra and values created

  5.1) standard spectras
  TH1F *hsign = (TH1F*) sig->GetUnlikeSignHistogram();  // same as the input (rebinned)
  TH1F *hbgrd = (TH1F*) sig->GetBackgroundHistogram();  // scaled input      (rebinned)
  TH1F *hextr = (TH1F*) sig->GetSignalHistogram();      // after backgound extraction (rebinned)
  TObject *oPeak = (TObject*) sig->GetPeakShape();      // can be a TF1 or TH1 depending on the extraction method
  TH1F *hrfac = (TH1F*) sig->GetRfactorHistogram();     // if like-sign correction was activated, o.w. 0x0
  5.2) access the extracted values and errors
  sig->GetValues();     or GetErrors();                 // yield extraction

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

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
#include "AliDielectron.h"

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
AliDielectronSignalExt::AliDielectronSignalExt(const char* name, const char* title, bool enummaps) :
  AliDielectronSignalBase(name, title, enummaps)
{
  //
  // Named + set enum maps Constructor
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
      AliWarning("Subtraction method not supported. Please check SetMethod() function.");
  }
}

//______________________________________________
void AliDielectronSignalExt::ProcessLS(TObjArray* const arrhist)
{
  //
  // signal subtraction 
  //

  fHistDataPP = (TH1*)(arrhist->At(AliDielectron::kEv1PP))->Clone("histPP");  // ++    SE
  fHistDataPM = (TH1*)(arrhist->At(AliDielectron::kEv1PM))->Clone("histPM");  // +-    SE
  fHistDataMM = (TH1*)(arrhist->At(AliDielectron::kEv1MM))->Clone("histMM");  // --    SE
  if(fHistDataPP->GetDefaultSumw2()) fHistDataPP->Sumw2();
  if(fHistDataPM->GetDefaultSumw2()) fHistDataPM->Sumw2();
  if(fHistDataMM->GetDefaultSumw2()) fHistDataMM->Sumw2();
  fHistDataPP->SetDirectory(0);
  fHistDataPM->SetDirectory(0);
  fHistDataMM->SetDirectory(0);
  
  // rebin the histograms
  if (fRebin>1) { 
    fHistDataPP->Rebin(fRebin);
    fHistDataPM->Rebin(fRebin);
    fHistDataMM->Rebin(fRebin);
  }       

  fHistRfactor = new TH1D("HistRfactor", "Rfactor;;N^{mix}_{+-}/2#sqrt{N^{mix}_{++} N^{mix}_{--}}",
                          fHistDataPM->GetXaxis()->GetNbins(),
                          fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistRfactor->Sumw2();
  fHistRfactor->SetDirectory(0);
  
  fHistSignal = new TH1D("HistSignal", "Like-Sign substracted signal",
			 fHistDataPM->GetXaxis()->GetNbins(),
			 fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistSignal->SetDirectory(0);
  fHistBackground = new TH1D("HistBackground", "Like-sign contribution",
			     fHistDataPM->GetXaxis()->GetNbins(),
			     fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground->SetDirectory(0);
  
  // fill out background histogram
  for(Int_t ibin=1; ibin<=fHistDataPM->GetXaxis()->GetNbins(); ibin++) {
    if(fHistDataPM->GetBinError(ibin)<1e-30 ) fHistDataPM->SetBinError(ibin, fgkErrorZero);
    Float_t pp = fHistDataPP->GetBinContent(ibin);
    Float_t mm = fHistDataMM->GetBinContent(ibin);

    Float_t background = 2*TMath::Sqrt(pp*mm);
    Float_t ebackground = TMath::Sqrt(mm+pp);
    if (fMethod==kLikeSignArithm || fMethod==kLikeSignArithmRcorr ){
      //Arithmetic mean instead of geometric
      background=(pp+mm);
      ebackground=TMath::Sqrt(pp+mm);
      if (TMath::Abs(ebackground)<1e-30) ebackground=fgkErrorZero;
    }

    fHistBackground->SetBinContent(ibin, background);
    fHistBackground->SetBinError(ibin, ebackground);
  }

  //correct LS spectrum bin-by-bin with R factor obtained in mixed events
  if(fMixingCorr || fMethod==kLikeSignRcorr || fMethod==kLikeSignArithmRcorr) {
    
    TH1* histMixPP = (TH1*)(arrhist->At(AliDielectron::kEv1PEv2P))->Clone("mixPP");  // ++    ME
    TH1* histMixMM = (TH1*)(arrhist->At(AliDielectron::kEv1MEv2M))->Clone("mixMM");  // --    ME

    TH1* histMixPM = 0x0;
    if (arrhist->At(AliDielectron::kEv1MEv2P)){
      histMixPM   = (TH1*)(arrhist->At(AliDielectron::kEv1MEv2P))->Clone("mixPM");  // -+    ME
    }

    if (arrhist->At(AliDielectron::kEv1PEv2M)){
      TH1 *htmp=(TH1*)(arrhist->At(AliDielectron::kEv1PEv2M));
      if (!histMixPM) fHistDataME   = (TH1*)htmp->Clone("mixPM");                   // +-    ME
      else histMixPM->Add(htmp);
    }
    
    if (!histMixPM){
      AliError("For R-factor correction the mixed event histograms are requires. No +- histogram found");
      return;
    }
    histMixPM->Sumw2();
    
    // rebin the histograms
    if (fRebin>1) { 
      histMixPP->Rebin(fRebin);
      histMixMM->Rebin(fRebin);
      histMixPM->Rebin(fRebin);
    }       
    
    // fill out rfactor histogram
    for(Int_t ibin=1; ibin<=histMixPM->GetXaxis()->GetNbins(); ibin++) {
      Float_t pp  = histMixPP->GetBinContent(ibin);
      Float_t ppe = histMixPP->GetBinError(ibin);
      Float_t mm  = histMixMM->GetBinContent(ibin);
      Float_t mme = histMixMM->GetBinError(ibin);
      Float_t pm  = histMixPM->GetBinContent(ibin);
      Float_t pme = histMixPM->GetBinError(ibin);
      
      Float_t background = 2*TMath::Sqrt(pp*mm);
      // do not use it since ME could be weighted err!=sqrt(entries)
      //      Float_t ebackground = TMath::Sqrt(mm+pp); 
      Float_t ebackground = TMath::Sqrt(mm/pp*ppe*ppe + pp/mm*mme*mme);
      if (fMethod==kLikeSignArithm){
        //Arithmetic mean instead of geometric
        background=(pp+mm);
        ebackground=TMath::Sqrt(ppe*ppe+mme*mme);
        if (TMath::Abs(ebackground)<1e-30) ebackground=fgkErrorZero;
      }
      
      Float_t rcon = 1.0;
      Float_t rerr = 0.0;
      if(background>0.0) {
        rcon = pm/background;
        rerr = TMath::Sqrt((1./background)*(1./background) * pme*pme +
                           (pm/(background*background))*(pm/(background*background)) * ebackground*ebackground);
      }
      fHistRfactor->SetBinContent(ibin, rcon);
      fHistRfactor->SetBinError(ibin, rerr);
    }
    
    fHistBackground->Multiply(fHistRfactor);
    
    if (histMixPP) delete histMixPP;
    if (histMixMM) delete histMixMM;
    if (histMixPM) delete histMixPM;
  }
  

  //scale histograms to match integral between fScaleMin and fScaleMax
  // or if fScaleMax <  fScaleMin use fScaleMin as scale factor
  if (fScaleMax>fScaleMin && fScaleMax2>fScaleMin2) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax,fScaleMin2,fScaleMax2);
  else if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax);
  else if (fScaleMin>0.){
    fScaleFactor=fScaleMin;
    fHistBackground->Scale(fScaleFactor);
  }

  //subract background
  fHistSignal->Add(fHistDataPM);
  fHistSignal->Add(fHistBackground,-1);
  
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
						 fHistBackground->FindBin(fIntMax), 
						 fErrors(1));

  // signal depending on peak description method
  DescribePeakShape(fPeakMethod, kTRUE, fHistSimPM); 
  //printf("%f  %f\n",fValues(0),fValues(1));
  // S/B and significance
  //  SetSignificanceAndSOB();

  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalExt::ProcessEM(TObjArray* const arrhist)
{
  //
  // event mixing of +- and -+
  //

  if (!arrhist->At(AliDielectron::kEv1PM) || !(arrhist->At(AliDielectron::kEv1MEv2P) || arrhist->At(AliDielectron::kEv1PEv2M)) ){
    AliError("Either OS or mixed histogram missing");
    return;
  }

  delete fHistDataPM; fHistDataPM=0x0;
  delete fHistDataME; fHistDataME=0x0;
  delete fHistBackground; fHistBackground=0x0;
  
  fHistDataPM = (TH1*)(arrhist->At(AliDielectron::kEv1PM))->Clone("histPM");  // +-    SE
  fHistDataPM->Sumw2();
  fHistDataPM->SetDirectory(0x0);

  if (arrhist->At(AliDielectron::kEv1MEv2P)){
    fHistDataME   = (TH1*)(arrhist->At(AliDielectron::kEv1MEv2P))->Clone("histMPME");  // -+    ME
  }
  
  if (arrhist->At(AliDielectron::kEv1PEv2M)){
    TH1 *htmp=(TH1*)(arrhist->At(AliDielectron::kEv1PEv2M));
    if (!fHistDataME) fHistDataME   = (TH1*)htmp->Clone("histMPME");  // -+    ME
    else fHistDataME->Add(htmp);
  }

  fHistBackground = (TH1*)fHistDataME->Clone("ME_Background");
  fHistBackground->SetDirectory(0x0);
  fHistBackground->Sumw2();

  // rebin the histograms
  if (fRebin>1) {
    fHistDataPM->Rebin(fRebin);
    fHistDataME->Rebin(fRebin);
    fHistBackground->Rebin(fRebin);
  }
  for(Int_t ibin=1; ibin<=fHistDataPM->GetXaxis()->GetNbins(); ibin++) {
    if(fHistDataPM->GetBinError(ibin)<1e-30 ) fHistDataPM->SetBinError(ibin, fgkErrorZero);
  }

  //scale histograms to match integral between fScaleMin and fScaleMax
  // or if fScaleMax <  fScaleMin use fScaleMin as scale factor
  if (fScaleMax>fScaleMin && fScaleMax2>fScaleMin2) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax,fScaleMin2,fScaleMax2);
  else if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax);
  else if (fScaleMin>0.){
    fScaleFactor=fScaleMin;
    fHistBackground->Scale(fScaleFactor);
  }
  fHistSignal=(TH1*)fHistDataPM->Clone("Signal");
  fHistSignal->Sumw2();
  //  printf(" err: %f %f \n",fHistSignal->GetBinError(75),TMath::Sqrt(fHistSignal->GetBinContent(75)));
  fHistSignal->Add(fHistBackground,-1.);
  //  printf(" err: %f %f \n",fHistSignal->GetBinError(75),TMath::Sqrt(fHistSignal->GetBinContent(75)));
//     // signal
//   fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
//                                              fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
                                                 fHistBackground->FindBin(fIntMax),
                                                 fErrors(1));

  // signal depending on peak description method
  DescribePeakShape(fPeakMethod, kTRUE, fHistSimPM);

  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalExt::ProcessRotation(TObjArray* const arrhist)
{
  //
  // signal subtraction
  //

  if (!arrhist->At(AliDielectron::kEv1PM) || !arrhist->At(AliDielectron::kEv1PMRot) ){
    AliError("Either OS or rotation histogram missing");
    return;
  }
  
  fHistDataPM = (TH1*)(arrhist->At(AliDielectron::kEv1PM))->Clone("histPM");  // +-    SE
  fHistDataPM->Sumw2();
  fHistDataPM->SetDirectory(0x0);

  fHistBackground = (TH1*)(arrhist->At(AliDielectron::kEv1PMRot))->Clone("histRotation");
  fHistBackground->Sumw2();
  fHistBackground->SetDirectory(0x0);

  // rebin the histograms
  if (fRebin>1) {
    fHistDataPM->Rebin(fRebin);
    fHistBackground->Rebin(fRebin);
  }
  for(Int_t ibin=1; ibin<=fHistDataPM->GetXaxis()->GetNbins(); ibin++) {
    if(fHistDataPM->GetBinError(ibin)<1e-30 ) fHistDataPM->SetBinError(ibin, fgkErrorZero);
  } 

  //scale histograms to match integral between fScaleMin and fScaleMax
  // or if fScaleMax <  fScaleMin use fScaleMin as scale factor
  if (fScaleMax>fScaleMin && fScaleMax2>fScaleMin2) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax,fScaleMin2,fScaleMax2);
  else if (fScaleMax>fScaleMin) fScaleFactor=ScaleHistograms(fHistDataPM,fHistBackground,fScaleMin,fScaleMax);
  else if (fScaleMin>0.){
    fScaleFactor=fScaleMin;
    fHistBackground->Scale(fScaleFactor);
  }

  fHistSignal=(TH1*)fHistDataPM->Clone("histSignal");
  fHistSignal->Add(fHistBackground,-1.);
  fHistSignal->SetDirectory(0x0);

  //     // signal
  //   fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
  //                                              fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
                                                 fHistBackground->FindBin(fIntMax),
                                                 fErrors(1));
  // signal depending on peak description method
  DescribePeakShape(fPeakMethod, kTRUE, fHistSimPM);
  
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

