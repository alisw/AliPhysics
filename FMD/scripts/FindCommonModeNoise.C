//____________________________________________________________________
//
// $Id: DrawDigits.C 30718 2009-01-22 16:07:40Z cholm $
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
#include <TH1F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TArrayF.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TPad.h>
#include <TLegend.h>

#include <AliCDBManager.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <AliLog.h>

#include <iostream>

/** @class DrawDigits
    @brief Draw pedestal, noise, and both common mode noise corrected
    @code 
    Root> .L Compile.C
    Root> Compile("FindCommonModeNoise.C")
    Root> FindCommonModeNoise fcmn("file");
    Root> fcmn.Run();
    @endcode

    If null is passed instead of a file, then a simulation of common
    mode noise is run instead.

    @ingroup FMD_script
 */
class FindCommonModeNoise : public AliFMDInput
{
private:
  TString fCalibDir;
  Int_t   fCount;
  TFile*  fOut;
  Float_t fShift;
  Float_t fCmn;                // Common mode noise component
  TH1F*   fPed;                // 
  TH1F*   fNoise;              //
  TH1F*   fCenteredPed;        //
  TH1F*   fCorrectedNoise;     //
  TH1F*   fSummedAdc;          //
  TH1F*   fSummedAdc2;         //
  TH1F*   fSummedCenteredAdc;  //
  TH1F*   fSummedCenteredAdc2; //
  TH1F*   fAdc;                //

public:
  //__________________________________________________________________
  FindCommonModeNoise(const char*  file, 
		      const char*  calibDir=0) 
    : AliFMDInput(), 
      fCalibDir(""), 
      fOut(0)
  { 
    fCmn = 2;
    if (calibDir) fCalibDir = calibDir;
    if (file) { 
      AddLoad(kRaw);
      SetRawFile(file);
    }
    else 
      AddLoad(kUser);
    
    TString outName(Form("histo_%s", gSystem->BaseName(fRawFile.Data())));
    if (!outName.EndsWith(".root")) outName.Append(".root");
    fOut  = TFile::Open(outName.Data(),"RECREATE");

    fPed                = MakeHist("ped_dist","Pedestal Distribution");
    fNoise              = MakeHist("noise_dist",
				   "Noise Distribution (Not CMN corrected)");
    fCenteredPed        = MakeHist("cmn_ped_dist",
				   "Pedestal - Average Pedestal Distribution");
    fCorrectedNoise     = MakeHist("noise_corrected_dist",
				   "Noise Distribution (CMN corrected)");
    fSummedAdc          = MakeHist("sum_adc_dist","Sum ADC Distribution");
    fSummedAdc2         = MakeHist("sum_adc2_dist",
				   "Sum ADC Squared Distribution");
    fSummedCenteredAdc  = MakeHist("sum_cmn_adc_dist",
				   "Sum CMN Adjusted ADC Distribution");
    fSummedCenteredAdc2 = MakeHist("sum_cmn_adc2_dist",
				   "Sum CMN Adjusted ADC Squared Distribution");
    fAdc                = MakeHist("adc_dist","Current ADC Distribution");
  }
  //__________________________________________________________________
  Int_t NEvents() const 
  { 
    Int_t nEv = (fReader ? fReader->GetNumberOfEvents() : -1);
    if (nEv > 0) return TMath::Min(nEv, 1000);
    return 1000;
  }
  //__________________________________________________________________
  TH1F* MakeHist(const char* name, const char* title)
  {
    TH1F* h = new TH1F(name,title,51200,-0.5,51199.5);
    h->SetXTitle("Strip number");
    h->SetFillColor(kRed);
    h->SetFillStyle(3001);
    h->SetDirectory(0);
    h->SetStats(kFALSE);
    return h;
  }
  //__________________________________________________________________
  Bool_t Init()
  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetRun(0);
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

    AliFMDParameters* pars = AliFMDParameters::Instance();
    if (fCalibDir.IsNull())
      pars->Init(kFALSE);
    else 
      pars->Init(fCalibDir.Data(), kFALSE);

    fCount = 0;

    return AliFMDInput::Init();
  }

  //__________________________________________________________________
  Bool_t Begin(Int_t e)
  {
    fCount++;
    fAdc->Reset();
    // Common mode noise generation
    fShift = gRandom->Gaus(0, fCmn);
    return AliFMDInput::Begin(e);
  }

  //__________________________________________________________________
  Int_t BinCenter(UShort_t d, Char_t r, UShort_t s, UShort_t t) const
  {
    if (d == 1 && !(r == 'I' || r == 'i')) return -1;

    Int_t idx = -1;
    switch (d) { 
    case 1: idx = 0; break;
    case 2: idx = 10240; break;
    case 3: idx = 3 * 10240; break;
    default: 
      return -1;
    }
    
    UShort_t q    = (r == 'I' || r == 'i' ? 0 : 1);
    // UShort_t nSec = (q == 0 ?  20 :  40);
    UShort_t nStr = (q == 0 ? 512 : 256);

    idx += (q == 0 ? 0 : 10240);
    idx += s * nStr + t;

    if (idx >= 51200) return -1;
    return idx;
  }
  
  //__________________________________________________________________
  Int_t BinCenter(UShort_t chip, UShort_t strip) const 
  {
    return chip * 128 + strip;
  }
  //__________________________________________________________________
  Int_t BinNumber(UShort_t chip, UShort_t strip) const 
  {
    return 1 + BinCenter(chip, strip);
  }
  //__________________________________________________________________
  Float_t GetSignal(UShort_t d, Char_t, UShort_t s, UShort_t t)
  {
    // Calculate signal value.  Only the first 256 strips of FMD are
    // filled. 
    if (d > 1 || s > 0 || t > 256) return 0;
    Float_t c = 0.005;
    Float_t x = (t % 128);
    Float_t p = c * (64 * 64 - 128 * x + x * x) + 100;
    Float_t n = (p-100)/40 + 1;

    // In case we simulate a common mode noise component. 
    if (t < 128) 
      return fShift + gRandom->Gaus(p, n);
   
    // In case all noise is uncorrelated. 
    return gRandom->Gaus(p, TMath::Sqrt(n*n + fCmn*fCmn));
  }    
  //__________________________________________________________________
  Bool_t ProcessRawDigit(AliFMDDigit* digit)
  {
    // From data 
    if (!digit) return kTRUE;

    UShort_t d             =  digit->Detector();
    Char_t   r             =  digit->Ring();
    UShort_t s             =  digit->Sector();
    UShort_t t             =  digit->Strip();
    return ProcessOne(d, r, s, t, digit->Counts());
  }
  //__________________________________________________________________
  Bool_t ProcessUser(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
  {
    // From our simulator in GetValue
    return ProcessOne(d, r, s, t, v);
  }
  //__________________________________________________________________
  Bool_t ProcessOne(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
  {
    Int_t    i             =  BinCenter(d, r, s, t);
    if (i < 0) { 
      std::cout << "Invalid index " << i << " returned for FMD" 
		<< d << r << '[' << s << ',' << t << ']' << std::endl;
      return kFALSE;
    }

    fAdc->Fill(i, v);

    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t End()
  {
    // At the end of each event, fill the summed histograms. 
    for (UShort_t chip = 0; chip < 400; chip++) { 
      Double_t chipAv = 0;
      for (UShort_t strip = 0; strip < 128; strip++) 
	chipAv += fAdc->GetBinContent(1+BinCenter(chip, strip));
      chipAv /= 128;
      
      for (UShort_t strip = 0; strip < 128; strip++) {
	Int_t    i    = BinCenter(chip, strip);
	Double_t adc  = fAdc->GetBinContent(i+1);
	Double_t cAdc = chipAv - adc;
	fSummedAdc->Fill(i, adc);
	fSummedAdc2->Fill(i, adc * adc);
	fSummedCenteredAdc->Fill(i, cAdc);
	fSummedCenteredAdc2->Fill(i, cAdc*cAdc);
      }
    }

    return AliFMDInput::End();
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    if (fEventCount == 0) { 
      std::cerr << "No events!" << std::endl;
      return kFALSE;
    }

    for (UInt_t bin = 1; bin <= 51200; bin++) { 
      if (bin % (51200 / 10) == 0) 
	std::cout << "Bin # " << bin << std::endl;
      Double_t sumAdc  = fSummedAdc->GetBinContent(bin);
      Double_t sumAdc2 = fSummedAdc2->GetBinContent(bin);
      Double_t avAdc   = sumAdc / fEventCount;
      Double_t avAdc2  = sumAdc2 / fEventCount;
      Double_t noise2  = avAdc2 - avAdc * avAdc;
    
      fPed->SetBinContent(bin, avAdc);
      fNoise->SetBinContent(bin, noise2 >= 0 ? TMath::Sqrt(noise2) : 0);
      
      Double_t sumCAdc  = fSummedCenteredAdc->GetBinContent(bin);
      Double_t sumCAdc2 = fSummedCenteredAdc2->GetBinContent(bin);
      Double_t avCAdc   = sumCAdc / fEventCount;
      Double_t avCAdc2  = sumCAdc2 / fEventCount;
      Double_t cNoise2  = avCAdc2 - avCAdc * avCAdc;
      fCenteredPed->SetBinContent(bin, avCAdc);
      fCorrectedNoise->SetBinContent(bin, TMath::Sqrt(cNoise2));
    }
    
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
  
    TCanvas* c = new TCanvas("c", "C", 800, 600);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(0);
    c->SetTopMargin(0);
    c->SetBottomMargin(0);
    c->Divide(2, 1);
    
    TPad* p1 = static_cast<TPad*>(c->cd(1));
    p1->SetTopMargin(0.05);
    p1->SetRightMargin(0.05);
    if (!fReader) fNoise->GetXaxis()->SetRangeUser(-0.5, 255.5);
    fNoise->SetMinimum(0);
    fNoise->SetMaximum(1.5*fNoise->GetMaximum());
    fNoise->DrawCopy();
    fCorrectedNoise->SetMinimum(0);
    fCorrectedNoise->SetFillColor(kBlue);
    fCorrectedNoise->DrawCopy("same");
    TLegend* l1 = new TLegend(0.2, 0.7, 0.945, 0.945);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->SetBorderSize(0);
    l1->AddEntry(fNoise, "Noise", "f");
    l1->AddEntry(fCorrectedNoise, "Corrected noise", "f");
    l1->Draw();
      

    TPad* p2 = static_cast<TPad*>(c->cd(2));
    p2->SetTopMargin(0.05);
    p2->SetRightMargin(0.05);
    if (!fReader) fPed->GetXaxis()->SetRangeUser(-0.5, 255.5);
    fPed->SetMinimum(fCenteredPed->GetMinimum());
    fPed->SetMaximum(1.5*fPed->GetMaximum());
    fPed->DrawCopy();
    fCenteredPed->SetFillColor(kBlue);
    fCenteredPed->DrawCopy("same");
    TLegend* l2 = new TLegend(0.2, 0.7, 0.945, 0.945);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->AddEntry(fPed, "Pedestal", "f");
    l2->AddEntry(fCenteredPed, "Pedestal - #bar{Pedestal}", "f");
    l2->Draw();

    c->Modified();
    c->Update();
    c->cd();
    c->Print("cmn.png");

    if (fOut) { 
      std::cout << "Flusing to disk ... " << std::flush;
      fOut->cd();
      fPed->Write();                //
      fNoise->Write();              //
      fCenteredPed->Write();        //
      fCorrectedNoise->Write();     //
      fSummedAdc->Write();          //
      fSummedAdc2->Write();         //
      fSummedCenteredAdc->Write();  //
      fSummedCenteredAdc2->Write(); //
      fAdc->Write();                //
      fOut->Write();
      fOut->Close();
      std::cout << "done" << std::endl;
    }
    return kTRUE;
  }

  ClassDef(FindCommonModeNoise,0);
};

//____________________________________________________________________
//
// EOF
//
