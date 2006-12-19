//____________________________________________________________________
//
// $Id$
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
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDUShortMap.h>
#include <AliFMDFloatMap.h>
#include <AliFMDRecPoint.h>
#include <AliESDFMD.h>
#include <AliLog.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>

/** @class DrawDigitsRecs
    @brief Draw digit ADC versus Rec point mult
    @code 
    Root> .L Compile.C
    Root> Compile("DrawDigitsRecs.C")
    Root> DrawDigitsRecs c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawDigitsRecs : public AliFMDInput
{
private:
  TH2D* fAdcVsSingle; // Histogram 
  AliFMDUShortMap fMap;
public:
  //__________________________________________________________________
  TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max) 
  {
    TArrayF bins(n+1);
    Float_t dp   = n / TMath::Log10(max / min);
    Float_t pmin = TMath::Log10(min);
    bins[0]      = min;
    for (Int_t i = 1; i < n+1; i++) {
      Float_t p = pmin + i / dp;
      bins[i]   = TMath::Power(10, p);
    }
    return bins;
  }
  //__________________________________________________________________
  DrawDigitsRecs(Int_t m=1100, Double_t amin=-0.5, Double_t amax=1099.5,
		 Int_t n=105, Double_t mmin=-0.5, Double_t mmax=20.5) 
  { 
    AddLoad(kDigits);
    AddLoad(kRecPoints);
    fAdcVsSingle = new TH2D("adcVsSingle", "ADC vs. Multiplicity (strip)", 
			    m, amin, amax, n, mmin, mmax);
    fAdcVsSingle->SetXTitle("ADC value");
    fAdcVsSingle->SetYTitle("Strip Multiplicity");
  }
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    fMap.Reset();
    return AliFMDInput::Begin(ev);
  }
  //__________________________________________________________________
  Bool_t ProcessDigit(AliFMDDigit* digit) 
  {
    // Cache the energy loss 
    if (!digit) return kFALSE;
    UShort_t det = digit->Detector();
    Char_t   rng = digit->Ring();
    UShort_t sec = digit->Sector();
    UShort_t str = digit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in digit", str));
      return kTRUE;
    }
    fMap(det, rng, sec, str) = digit->Counts();
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessRecPoint(AliFMDRecPoint* single)
  {
    if (!single) return kFALSE;
    UShort_t det = single->Detector();
    Char_t   rng = single->Ring();
    UShort_t sec = single->Sector();
    UShort_t str = single->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in single", str));
      return kFALSE;
    }
    fAdcVsSingle->Fill(fMap(det, rng, sec, str), single->Particles());
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

    fAdcVsSingle->SetStats(kFALSE);
    fAdcVsSingle->Draw("COLZ");
    return kTRUE;
  }

  ClassDef(DrawDigitsRecs,0);
  
};

//____________________________________________________________________
//
// EOF
//
