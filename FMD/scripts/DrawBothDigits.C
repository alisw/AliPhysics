//____________________________________________________________________
//
// $Id: DrawBothDigits.C 20907 2007-09-25 08:44:03Z cholm $
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
#include <TH2F.h>
#include <AliFMDDigit.h>
#include <AliFMDSDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <AliLog.h>

/** @class DrawBothDigits
    @brief Draw hit energy loss versus digit ADC
    @code 
    Root> .L Compile.C
    Root> Compile("DrawBothDigits.C")
    Root> DrawBothDigits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawBothDigits : public AliFMDInput
{
private:
  TH2F* fTrackNos; // Histogram 
  AliFMDEdepMap fCache;
public:
  //__________________________________________________________________
  DrawBothDigits(Int_t max=300) 
    : AliFMDInput("galice.root")
  { 
    AddLoad(kDigits);
    AddLoad(kSDigits);
    fTrackNos = new TH2F("trackNos", "Track numbers", 
			 max+1, -1.5, max-.5, max+1, -1.5, max-.5);
    fTrackNos->SetXTitle("Digit track");
    fTrackNos->SetYTitle("SDigit track");
  }
  //__________________________________________________________________
  Bool_t Begin(Int_t evno)
  {
    fCache.Reset();
    return AliFMDInput::Begin(evno);
  }
  //__________________________________________________________________
  Bool_t ProcessSDigit(AliFMDSDigit* sdigit)
  {
    if (!sdigit) return kTRUE;
    AliFMDEdepHitPair& entry = fCache(sdigit->Detector(), 
				      sdigit->Ring(), 
				      sdigit->Sector(), 
				      sdigit->Strip());
    entry.fLabels.Set(sdigit->GetNTrack());
    Info("ProcessSDigit", "Got %d SDigit tracks", sdigit->GetNTrack());
    for (size_t i = 0; i < sdigit->GetNTrack(); i++) 
      entry.fLabels.fArray[i] = sdigit->GetTrack(i);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessDigit(AliFMDDigit* digit)
  {
    if (!digit) return kTRUE;
    AliFMDEdepHitPair& entry = fCache(digit->Detector(), 
				      digit->Ring(), 
				      digit->Sector(), 
				      digit->Strip());
    TArrayI& stracks = entry.fLabels;
    Info("ProcessDigit", "Got %d SDigit tracks, and %d Digit tracks", 
	 stracks.fN, digit->GetNTrack());
    for (Int_t i = 0; (i < stracks.fN || i < digit->GetNTrack()); i++) { 
      Int_t strack = (i < stracks.fN ? stracks.fArray[i] : -1);
      Int_t dtrack = (i < digit->GetNTrack() ? digit->GetTrack(i) : -1);
      fTrackNos->Fill(strack, dtrack);
    }
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
    fTrackNos->SetStats(kFALSE);
    // fTrackNos->Scale(1. / fAdc->GetEntries());
    fTrackNos->Draw("colz");
    return kTRUE;
  }

  ClassDef(DrawBothDigits,0);
};

//____________________________________________________________________
//
// EOF
//
