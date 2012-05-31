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
#include <TH1D.h>
#include <AliCDBManager.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <AliFMDParameters.h>
#include <AliFMDCalibStripRange.h>
#include <AliFMDCalibSampleRate.h>
#include <AliFMDCalibGain.h>
#include <AliFMDCalibPedestal.h>
#include <AliFMDRawReader.h>
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TArrayF.h>
#include <AliLog.h>
#include <TSystem.h>
#include <TFile.h>

/** @class DrawDigits
    @brief Draw hit energy loss versus digit ADC
    @code 
    Root> .L Compile.C
    Root> Compile("DrawCalibRaw.C")
    Root> DrawCalibRaw cr;
    Root> cr.SetRawFile("file");
    Root> cr.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawCalibRaw : public AliFMDInput
{
private:
  Double_t   fFactor;
  TString    fCalibDir;
  TH1D*      fELoss; // Histogram 
  TObjArray* fDets;
  TFile*     fOut;
  Bool_t     fHasData;
  Int_t      fGotNEvents;

public:
  //__________________________________________________________________
  DrawCalibRaw(const       char*  file, 
	       const char* calibDir    = 0,
	       Double_t    noiseFactor = 5, 
	       Bool_t      save        = kTRUE,
	       Int_t       m           = 420, 
	       Double_t    mmin        = -0.5, 
	       Double_t    mmax        = 20.5) 
    : AliFMDInput(), 
      fFactor(noiseFactor),
      fCalibDir(""), 
      fELoss(0), 
      fDets(0), 
      fOut(0), 
      fHasData(kFALSE),
      fGotNEvents(0)
  { 
    if (calibDir) fCalibDir = calibDir;
    AddLoad(kRawCalib);
    fELoss = new TH1D("eLoss", "Scaled Energy loss", m, mmin, mmax);
    fELoss->SetXTitle("#Delta E/#Delta E_{mip}");
    fELoss->Sumw2();
    fELoss->SetDirectory(0);
    
    SetRawFile(file);

    if (!save) return;
    fDets = new TObjArray(3);

  }
  //__________________________________________________________________
  Bool_t CheckFile(const char* prefix, int number, TString& f)
  {
    f = (Form("%s%d.csv", prefix, number));
    f = gSystem->Which(fCalibDir.Data(), f.Data());
    return !f.IsNull();
  }

  //__________________________________________________________________
  Bool_t Init()
  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetRun(0);
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

    AliFMDCalibStripRange* range = 0;
    AliFMDCalibSampleRate* rate  = 0;
    AliFMDCalibPedestal*   peds  = 0;
    AliFMDCalibGain*       gains = 0;
    for (Int_t i = 1; i <= 3; i++) { 
      TString f;
      if (CheckFile("conditions", i, f)) {
	Info("Init", "Reading conditions for FMD%d from %s", i, f.Data());
	if (!range) range = new AliFMDCalibStripRange;
	if (!rate)  rate  = new AliFMDCalibSampleRate;
	std::ifstream in(f.Data());
	range->ReadFromFile(in);
	rate->ReadFromFile(in);
      }
      if (CheckFile("peds", i, f)) {
	Info("Init", "Reading pedestals for FMD%d from %s", i, f.Data());
	if (!peds) peds = new AliFMDCalibPedestal;
	std::ifstream in(f.Data());
	peds->ReadFromFile(in);
      }
      if (CheckFile("gains", i, f)) {
	Info("Init", "Reading gains for FMD%d from %s", i, f.Data());
	if (!gains) gains = new AliFMDCalibGain;
	std::ifstream in(f.Data());
	gains->ReadFromFile(in);
      }
    }

    Int_t mask = (AliFMDParameters::kDeadMap|
		  AliFMDParameters::kZeroSuppression|
		  AliFMDParameters::kAltroMap);

    if (!range) mask |= AliFMDParameters::kStripRange;
    if (!rate)  mask |= AliFMDParameters::kSampleRate;
    if (!peds)  mask |= AliFMDParameters::kPedestal;
    if (!gains) mask |= AliFMDParameters::kPulseGain;

    AliFMDParameters* pars = AliFMDParameters::Instance();
    pars->Init(kFALSE, mask);

    if (range)  pars->SetStripRange(range);
    if (rate)   pars->SetSampleRate(rate);
    if (peds)   pars->SetPedestal(peds);
    if (gains)  pars->SetGain(gains);
    
    Bool_t ret = AliFMDInput::Init();
    
    if (!fDets) return ret;

    fOut  = TFile::Open(Form("histo_%s", fRawFile.Data()), "RECREATE");

    Int_t m    = fELoss->GetXaxis()->GetNbins();
    Int_t mmin = fELoss->GetXaxis()->GetXmin();
    Int_t mmax = fELoss->GetXaxis()->GetXmax();
    for (Int_t d = 1; d <= 3; d++) {
      Int_t      nRng = (d == 1 ? 1 : 2);
      TObjArray* det  = 0;
      fDets->AddAt(det = new TObjArray(nRng), d-1);
      det->SetName(Form("FMD%d", d));
      TDirectory* detD = fOut->mkdir(det->GetName());
      for (Int_t q = 0; q < nRng; q++) { 
	Char_t r       = q == 0 ? 'I' : 'O';
	Int_t  nSec    = q == 0 ?  20 :  40;
	Int_t  nStr    = q == 0 ? 512 : 256;
	TObjArray* rng = 0;
	det->AddAt(rng = new TObjArray(nSec), q);
	rng->SetName(Form("FMD%d%c", d, r));
	TDirectory* rngD = detD->mkdir(rng->GetName());
	for (Int_t s = 0; s < nSec; s++) { 
	  TObjArray* sec = 0;
	  rng->AddAt(sec = new TObjArray(nStr), s);
	  sec->SetName(Form("FMD%d%c_%02d", d, r, s));
	  TDirectory* secD = rngD->mkdir(sec->GetName());
	  for (Int_t t = 0; t < nStr; t++) { 
	    secD->cd();
	    TH1* str = new TH1D(Form("FMD%d%c_%02d_%03d", d, r, s, t), 
				 Form("Scaled energy loss in FMD%d%c[%2d,%3d]",
				      d, r, s, t), m, mmin, mmax);
	    str->SetXTitle("#Delta E/#Delta E_{mip}");
	    // str->SetDirectory(secD);
	    sec->AddAt(str, t);
	  }
	}
      }
    }

    return ret;
  }

  //__________________________________________________________________
  Bool_t Begin(Int_t e)
  {
    fHasData = kFALSE;
    return AliFMDInput::Begin(e);
  }

  //__________________________________________________________________
  Bool_t ProcessRawCalibDigit(AliFMDDigit* digit)
  {
    if (!digit) return kTRUE;

    AliFMDParameters* parm = AliFMDParameters::Instance();
    UShort_t d             =  digit->Detector();
    Char_t   r             =  digit->Ring();
    UShort_t s             =  digit->Sector();
    UShort_t t             =  digit->Strip();
    Double_t gain          =  parm->GetPulseGain(d, r, s, t);
    Double_t ped           =  parm->GetPedestal(d, r, s, t);
    Double_t pedW          =  parm->GetPedestalWidth(d, r, s, t);
    Double_t adc           =  digit->Counts();
    Double_t threshold     =  pedW * fFactor;
    if (gain < 0.1 || gain > 10) return kTRUE;
    if (pedW > 10) { 
      Warning("ProcessRawCalibDigit", "FMD%d%c[%2d,%3d] is noisy: %f",
	      d, r, s, t, pedW);
      return kTRUE;
    }

    if (fFMDReader && fFMDReader->IsZeroSuppressed(d-1))
      adc += fFMDReader->NoiseFactor(d-1) * pedW;
    else 
      threshold            += ped;

    if (adc < threshold) return kTRUE;
    
    Double_t mult = (adc-ped) / (gain * parm->GetDACPerMIP());

    fHasData = kTRUE;
    fELoss->Fill(mult);

    // if (t >= 10) return kTRUE;
    TObjArray* det = static_cast<TObjArray*>(fDets->At(d-1));
    TObjArray* rng = static_cast<TObjArray*>(det->At(r == 'I' ? 0 : 1));
    TObjArray* sec = static_cast<TObjArray*>(rng->At(s));
    TH1*       str = static_cast<TH1*>(sec->At(t));
    str->Fill(mult);

    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t End()
  {
    if (fHasData) fGotNEvents++;
    return AliFMDInput::End();
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    std::cout << "A total of " << fGotNEvents << " with data" << std::endl;
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    fELoss->SetStats(kFALSE);
    fELoss->SetFillColor(kRed);
    fELoss->SetFillStyle(3001);
    fELoss->Scale(1. / fELoss->GetEntries());
    fELoss->DrawCopy("e1 bar");

    if (fDets && fOut) { 
      std::cout << "Flusing to disk ... " << std::flush;
      fOut->cd();
      fELoss->Write();
      fOut->Write();
      fOut->Close();
      std::cout << "done" << std::endl;
    }
    return kTRUE;
  }

  ClassDef(DrawCalibRaw,0);
};

//____________________________________________________________________
//
// EOF
//
