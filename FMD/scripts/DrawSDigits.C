//____________________________________________________________________
//
// $Id: DrawSDigits.C 20907 2007-09-25 08:44:03Z cholm $
//
// Script that contains a class to draw eloss from hits, versus ADC
// counts from digits, using the AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <AliFMDSDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <AliFMDGeometry.h>
#include <iostream>
#include <TStyle.h>
#include <TLatex.h>
#include <TArrayF.h>
#include <AliLog.h>

/** @class DrawSDigits
    @brief Draw hit energy loss versus digit ADC
    @code 
    Root> .L Compile.C
    Root> Compile("DrawSDigits.C")
    Root> DrawSDigits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawSDigits : public AliFMDInput
{
private:
  TH1D* fAdc; // Histogram 
  TProfile2D* fPrimRatio[5]; // Histograms
public:
  void Idx2Det(UShort_t idx, UShort_t& d, Char_t& r) const
  {
    d = 0;
    r = '\0';
    switch (idx) { 
    case 0: d = 1; r = 'I'; break;
    case 1: d = 2; r = 'I'; break;
    case 2: d = 2; r = 'O'; break;
    case 3: d = 3; r = 'I'; break;
    case 4: d = 3; r = 'O'; break;
    }
  }
  Short_t Det2Idx(UShort_t d, Char_t r) const
  { 
    Short_t idx = -1;
    switch (d) { 
    case 1: idx = 0; break;
    case 2: idx = 1; break;
    case 3: idx = 3; break;
    }
    return (idx + ((r == 'O' || r == 'o') ? 1 : 0));
  }
  
    
  //__________________________________________________________________
  DrawSDigits(Int_t m=1100, Double_t amin=-0.5, Double_t amax=1023.5) 
  { 
    AddLoad(kSDigits);
    AddLoad(kGeometry);
    fAdc = new TH1D("adc", "ADC", m, amin, amax);
    fAdc->SetXTitle("ADC value");
    fAdc->Sumw2();

    Int_t   eN   = 50;
    Float_t eMin = -4;
    Float_t eMax = 6;
    Int_t   pN   = 40;
    Float_t pMin = -.5;
    Float_t pMax = 359.5;
    for (UShort_t i = 0; i < 5; i++) { 
      UShort_t d;
      Char_t   r;
      Idx2Det(i, d, r);
      fPrimRatio[i] = new TProfile2D(Form("primRatio_fmd%d%c", d, r), 
				     Form("Primary/Total - FMD%d%c", d, r), 
				     eN, eMin, eMax, pN, pMin, pMax);
      fPrimRatio[i]->SetXTitle("#eta");
      fPrimRatio[i]->SetYTitle("#varphi");
      fPrimRatio[i]->SetZTitle("M_{ch,primary}/M_{ch,total}");
      fPrimRatio[i]->GetXaxis()->SetTitleFont(132);
      fPrimRatio[i]->GetYaxis()->SetTitleFont(132);
      fPrimRatio[i]->GetZaxis()->SetTitleFont(132);
      fPrimRatio[i]->GetXaxis()->SetLabelFont(132);
      fPrimRatio[i]->GetYaxis()->SetLabelFont(132);
      fPrimRatio[i]->GetZaxis()->SetLabelFont(132);
      fPrimRatio[i]->GetXaxis()->SetTitleSize(.06);
      fPrimRatio[i]->GetYaxis()->SetTitleSize(.06);
      fPrimRatio[i]->GetZaxis()->SetTitleSize(.06);
      fPrimRatio[i]->GetXaxis()->SetLabelSize(.06);
      fPrimRatio[i]->GetYaxis()->SetLabelSize(.06);
      fPrimRatio[i]->GetZaxis()->SetLabelSize(.06);
      fPrimRatio[i]->GetXaxis()->SetNdivisions(10);
      fPrimRatio[i]->GetYaxis()->SetNdivisions(10);
      fPrimRatio[i]->GetZaxis()->SetNdivisions(10);
    }
  }
  Bool_t Init() 
  {
    Bool_t ret = AliFMDInput::Init();
    AliFMDGeometry::Instance()->Init();
    AliFMDGeometry::Instance()->InitTransformations();
    return ret;
  }
    
  //__________________________________________________________________
  Bool_t ProcessSDigit(AliFMDSDigit* digit)
  {
    if (!digit) return kTRUE;
    fAdc->Fill(digit->Counts());
    if (digit->NParticles() == 0) return kTRUE;
    

    Short_t primIdx = Det2Idx(digit->Detector(), digit->Ring());
    if (primIdx < 0) return kTRUE;
    
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    Double_t x, y, z;
    geom->Detector2XYZ(digit->Detector(), 
		       digit->Ring(), 
		       digit->Sector(), 
		       digit->Strip(), 
		       x, y, z);
    // We should get the primary vertex and do 
    // z     += fCurrentVertex;
    Double_t phi   =  TMath::ATan2(y, x) * 180. / TMath::Pi();
    Double_t r     =  TMath::Sqrt(y * y + x * x);
    Double_t theta =  TMath::ATan2(r, z);
    Double_t eta   = -TMath::Log(TMath::Tan(theta / 2));
    Double_t ratio = digit->NPrimaries() / digit->NParticles();
    if (phi < 0) phi += 360;
    fPrimRatio[primIdx]->Fill(eta, phi, ratio);
    
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    
    TCanvas* c1 = new TCanvas("fmdSdigitSpectra", "FMD SDigit spectra");
    fAdc->SetStats(kFALSE);
    if (fAdc->GetEntries() != 0) 
      fAdc->Scale(1. / fAdc->GetEntries());
    fAdc->Draw();
    c1->Modified();
    c1->Update();
    c1->cd();

    TCanvas* c2 = new TCanvas("fmdSdigitPrims", "FMD SDigit # prim/# ntotal",
			      800, 800);
    c2->SetTopMargin(0);
    c2->SetRightMargin(0);
    c2->Divide(2, 3, 0, 0);
    for (UShort_t i = 0; i < 5; i++) { 
      Int_t        np = i+1+(i == 0 ? 0 : 1);
      TVirtualPad* p   = c2->cd(np);
      p->SetFillColor(0);
      p->SetTopMargin((np <= 2) ? 0.05 : 0);
      p->SetLeftMargin((np % 2 == 1) ? 0.20 : 0);
      p->SetRightMargin((np % 2 == 1) ? 0 : 0.20);
      p->SetBottomMargin((np >= 5) ? 0.15 : 0);
      fPrimRatio[i]->SetStats(kFALSE);
      fPrimRatio[i]->Draw(Form("COL%c", (np % 2 == 1) ? ' ' : 'Z'));
      UShort_t d;
      Char_t   r;
      Idx2Det(i, d, r);
      TLatex* text = new TLatex(0, 180, Form("FMD%d%c", d, r));
      text->SetTextFont(132);
      text->SetTextAlign(22);
      text->SetTextSize(0.08);
      text->Draw();
    }
    c2->Modified();
    c2->Update();
    c2->cd();

    return kTRUE;
  }

  ClassDef(DrawSDigits,0);
};

//____________________________________________________________________
//
// EOF
//
