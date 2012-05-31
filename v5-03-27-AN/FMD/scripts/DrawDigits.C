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
#include <TH1D.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <AliLog.h>

/** @class DrawDigits
    @brief Draw hit energy loss versus digit ADC
    @code 
    Root> .L Compile.C
    Root> Compile("DrawDigits.C")
    Root> DrawDigits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawDigits : public AliFMDInput
{
private:
  TH1D* fAdc; // Histogram 
public:
  //__________________________________________________________________
  DrawDigits(Int_t m=1100, Double_t amin=-0.5, Double_t amax=1023.5) 
    : AliFMDInput("galice.root")
  { 
    AddLoad(kDigits);
    fAdc = new TH1D("adc", "ADC", m, amin, amax);
    fAdc->SetXTitle("ADC value");
    fAdc->Sumw2();
  }
  //__________________________________________________________________
  Bool_t ProcessDigit(AliFMDDigit* digit)
  {
    if (!digit) return kTRUE;
    fAdc->Fill(digit->Counts());
    digit->Print("l");
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
    fAdc->SetStats(kFALSE);
    fAdc->Scale(1. / fAdc->GetEntries());
    fAdc->Draw();
    return kTRUE;
  }

  ClassDef(DrawDigits,0);
};

//____________________________________________________________________
//
// EOF
//
