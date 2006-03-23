//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to compare the raw data written to the
// digits it's created from.
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <AliLog.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDUShortMap.h>
#include <AliFMDParameters.h>
#include <AliFMDGeometry.h>
#include <AliFMDRing.h>
#include <AliFMDDetector.h>
#include <iostream>
#include <TH3D.h>
#include <TStyle.h>
#include <TArrayF.h>
#include <TParticle.h>
#include <TCanvas.h>

class CheckAlign : public AliFMDInput
{
public:
  CheckAlign()
  {
    AddLoad(kHits);
    AddLoad(kDigits);
    AddLoad(kGeometry);
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    geom->Init();
    // geom->InitTransformations();
    Double_t iR  = geom->GetRing('I')->GetHighR()+5;
    Double_t oR  = geom->GetRing('O')->GetHighR()+5;
    Double_t z1l = geom->GetDetector(1)->GetInnerZ()-5;
    Double_t z1h = z1l+10;
    Double_t z2l = geom->GetDetector(2)->GetOuterZ()-5;
    Double_t z2h = geom->GetDetector(2)->GetInnerZ()+5;
    Double_t z3l = geom->GetDetector(3)->GetOuterZ()-5;
    Double_t z3h = geom->GetDetector(3)->GetInnerZ()+5;
    
    f1Hits   = new TH3D("hits1", "FMD1 hits", 
			300,-iR,iR,300,-iR,iR,100,z1l,z1h);
    f1Hits->SetMarkerColor(2);
    f1Hits->SetMarkerStyle(3);
      
    f2Hits   = new TH3D("hits2", "FMD2 hits", 
			300,-oR,oR,300,-oR,oR,100,z2l,z2h);
    f2Hits->SetMarkerColor(2);
    f2Hits->SetMarkerStyle(3);

    f3Hits   = new TH3D("hits3", "FMD3 hits", 
			300,-oR,oR,300,-oR,oR,100,z3l,z3h);
    f3Hits->SetMarkerColor(2);
    f3Hits->SetMarkerStyle(3);

    f1Digits = new TH3D("digits1", "FMD1 digits", 
			300,-iR,iR,300,-iR,iR,100,z1l,z1h); 
    f1Digits->SetMarkerColor(4);
    f1Digits->SetMarkerStyle(2);

    f2Digits = new TH3D("digits2", "FMD2 digits", 
			300,-oR,oR,300,-oR,oR,100,z2l,z2h);
    f2Digits->SetMarkerColor(4);
    f2Digits->SetMarkerStyle(2);

    f3Digits = new TH3D("digits3", "FMD3 hits", 
			300,-oR,oR,300,-oR,oR,100,z3l,z3h);
    f3Digits->SetMarkerColor(4);
    f3Digits->SetMarkerStyle(2);
  }
  Bool_t Init() 
  {
    Bool_t ret = AliFMDInput::Init();
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    geom->Init();
    geom->InitTransformations();
    AliFMDParameters* param = AliFMDParameters::Instance();
    param->Init();
    return ret;
  }
  
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*)
  {
    // Cache the energy loss 
    if (!hit) return kFALSE;
    UShort_t det = hit->Detector();
    UShort_t str = hit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in hit", str));
      return kTRUE;
    }
    switch (det) {
    case 1: f1Hits->Fill(hit->X(), hit->Y(), hit->Z()); break;
    case 2: f2Hits->Fill(hit->X(), hit->Y(), hit->Z()); break;
    case 3: f3Hits->Fill(hit->X(), hit->Y(), hit->Z()); break;
    }
    return kTRUE;
  }
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
    AliFMDParameters* param = AliFMDParameters::Instance();
    if (digit->Counts() < (param->GetPedestal(det, rng, sec, str) + 4 *
			   param->GetPedestalWidth(det, rng, sec, str)))
      return kTRUE;
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    Double_t x, y, z;
    geom->Detector2XYZ(det, rng, sec, str, x, y, z);
    switch (det) {
    case 1: f1Digits->Fill(x, y , z); break;
    case 2: f2Digits->Fill(x, y , z); break;
    case 3: f3Digits->Fill(x, y , z); break;
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

    TCanvas* c1 = new TCanvas("FMD1","FMD1");
    c1->cd();
    f1Hits->Draw();
    f1Digits->Draw("same");

    TCanvas* c2 = new TCanvas("FMD2","FMD2");
    c2->cd();
    f2Hits->Draw();
    f2Digits->Draw("same");

    TCanvas* c3 = new TCanvas("FMD3","FMD3");
    c3->cd();
    f3Hits->Draw();
    f3Digits->Draw("same");
    return kTRUE;
  }
protected:
  TH3D* f1Hits;
  TH3D* f2Hits;
  TH3D* f3Hits;
  TH3D* f1Digits;
  TH3D* f2Digits;
  TH3D* f3Digits;
};


//____________________________________________________________________
//
// EOF
//


  
  
