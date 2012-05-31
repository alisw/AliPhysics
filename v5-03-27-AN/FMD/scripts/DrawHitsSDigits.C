//____________________________________________________________________
//
// $Id: DrawHitsSDigits.C 20907 2007-09-25 08:44:03Z cholm $
//
// Script that contains a class to draw eloss from hits, versus ADC
// counts from sdigits, using the AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDSDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <AliLog.h>
#include <TMath.h>

/** @class DrawHitsSDigits
    @brief Draw hit energy loss versus sdigit ADC
    @code 
    Root> .L Compile.C
    Root> Compile("DrawHitsSDigits.C")
    Root> DrawHitsSDigits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawHitsSDigits : public AliFMDInput
{
private:
  TH2D* fElossVsAdc; // Histogram 
  AliFMDEdepMap fMap;
public:
  //__________________________________________________________________
  DrawHitsSDigits(Int_t n=900, Double_t emin=1e-3, Double_t emax=10, 
		 Int_t m=1100, Double_t amin=-0.5, Double_t amax=1099.5) 
  { 
    AddLoad(kSDigits);
    AddLoad(kHits);
    TArrayF eloss(MakeLogScale(n, emin, emax));
    TArrayF adcs(m+1);
    adcs[0] = amin;
    for (Int_t i = 1; i < m+1; i++) adcs[i] = adcs[i-1] + (amax-amin)/m;
    fElossVsAdc = new TH2D("bad", "#Delta E vs. ADC", 
			   eloss.fN-1, eloss.fArray, adcs.fN-1, adcs.fArray);
    fElossVsAdc->SetXTitle("#Delta E/#Delta x [MeV/cm]");
    fElossVsAdc->SetYTitle("ADC value");
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
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    Info("ProcessHit", "Processing hit %p", hit);
    // Cache the energy loss 
    if (!hit) return kFALSE;
    UShort_t det = hit->Detector();
    Char_t   rng = hit->Ring();
    UShort_t sec = hit->Sector();
    UShort_t str = hit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in hit", str));
      return kTRUE;
    }
    fMap(det, rng, sec, str).fEdep += hit->Edep();
    fMap(det, rng, sec, str).fN++;
    hit->Print();
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessSDigit(AliFMDSDigit* sdigit)
  {
    Info("ProcessSDigit", "Processing sdigit %p", sdigit);
    if (!sdigit) return kFALSE;
    UShort_t det = sdigit->Detector();
    Char_t   rng = sdigit->Ring();
    UShort_t sec = sdigit->Sector();
    UShort_t str = sdigit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in sdigit", str));
      return kFALSE;
    }
    fElossVsAdc->Fill(fMap(det, rng, sec, str).fEdep, sdigit->Counts());
    sdigit->Print();
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
    fElossVsAdc->SetStats(kFALSE);
    fElossVsAdc->Draw("COLZ");
    return kTRUE;
  }

  ClassDef(DrawHitsSDigits,0);
};

//____________________________________________________________________
//
// EOF
//
